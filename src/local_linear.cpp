// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace {

void validate_inputs(const arma::mat& predictors,
                     const arma::vec& response,
                     const arma::umat& neighbors) {
    if (predictors.n_rows != response.n_elem) {
        Rcpp::stop("Predictor rows and response length differ.");
    }
    if (predictors.n_rows != neighbors.n_rows) {
        Rcpp::stop("One neighbor row is required per observation.");
    }
    if (neighbors.n_cols <= predictors.n_cols + 1) {
        Rcpp::stop("Neighborhood size must exceed d + 1.");
    }
    if (neighbors.max() >= predictors.n_rows) {
        Rcpp::stop("Neighbor index is outside the data range.");
    }
}

void validate_query_inputs(const arma::mat& predictors,
                           const arma::vec& response,
                           const arma::uvec& query_index,
                           const arma::umat& neighbors) {
    if (predictors.n_rows != response.n_elem) {
        Rcpp::stop("Predictor rows and response length differ.");
    }
    if (query_index.n_elem != neighbors.n_rows) {
        Rcpp::stop("One neighbor row is required per query index.");
    }
    if (query_index.n_elem == 0) {
        Rcpp::stop("At least one query index is required.");
    }
    if (neighbors.n_cols <= predictors.n_cols + 1) {
        Rcpp::stop("Neighborhood size must exceed d + 1.");
    }
    if (query_index.max() >= predictors.n_rows ||
        neighbors.max() >= predictors.n_rows) {
        Rcpp::stop("Query or neighbor index is outside the data range.");
    }
}

void validate_external_query_inputs(const arma::mat& predictors,
                                    const arma::vec& response,
                                    const arma::mat& query_predictors,
                                    const arma::umat& neighbors) {
    if (predictors.n_rows != response.n_elem) {
        Rcpp::stop("Predictor rows and response length differ.");
    }
    if (query_predictors.n_rows != neighbors.n_rows) {
        Rcpp::stop("One neighbor row is required per query point.");
    }
    if (query_predictors.n_rows == 0) {
        Rcpp::stop("At least one query point is required.");
    }
    if (query_predictors.n_cols != predictors.n_cols) {
        Rcpp::stop("Training and query dimensions differ.");
    }
    if (neighbors.n_cols <= predictors.n_cols + 1) {
        Rcpp::stop("Neighborhood size must exceed d + 1.");
    }
    if (neighbors.max() >= predictors.n_rows) {
        Rcpp::stop("Neighbor index is outside the data range.");
    }
}

arma::mat local_design(const arma::mat& predictors,
                       const arma::uvec& index,
                       const arma::rowvec& query) {
    const arma::uword k = index.n_elem;
    const arma::uword d = predictors.n_cols;
    arma::mat design(k, d + 1);
    design.col(0).ones();
    design.cols(1, d) = predictors.rows(index);
    design.cols(1, d).each_row() -= query;
    return design;
}

arma::vec ridge_coefficients(const arma::mat& design,
                             const arma::vec& response,
                             double lambda) {
    arma::mat cross = design.t() * design;
    arma::vec rhs = design.t() * response;
    for (arma::uword j = 1; j < cross.n_cols; ++j) {
        cross(j, j) += lambda;
    }

    arma::vec coefficient;
    const bool solved = arma::solve(
        coefficient,
        cross,
        rhs,
        arma::solve_opts::likely_sympd
    );
    if (!solved) {
        coefficient = arma::pinv(cross) * rhs;
    }
    return coefficient;
}

double ridge_sse(const arma::mat& design,
                 const arma::vec& response,
                 double lambda) {
    const arma::vec coefficient = ridge_coefficients(
        design, response, lambda
    );
    const arma::vec residual = response - design * coefficient;
    return std::max(0.0, arma::dot(residual, residual));
}

using distance_index = std::pair<double, arma::uword>;

struct forward_candidate_result {
    double total_sse;
    bool boundary_tie;
};

bool distance_index_less(const distance_index& left,
                         const distance_index& right) {
    if (left.first < right.first) return true;
    if (right.first < left.first) return false;
    return left.second < right.second;
}

std::vector<double> selected_squared_distances_block(
    const arma::mat& predictors,
    const arma::uvec& query_index,
    const arma::uvec& selected,
    arma::uword query_begin,
    arma::uword query_count
) {
    const arma::uword n = predictors.n_rows;
    const std::size_t rows = static_cast<std::size_t>(query_count);
    const std::size_t columns = static_cast<std::size_t>(n);
    if (rows > 0 && columns > std::vector<double>().max_size() / rows) {
        Rcpp::stop("Forward distance block is too large.");
    }
    std::vector<double> distances(rows * columns, 0.0);

    for (arma::uword query = 0; query < query_count; ++query) {
        const arma::uword query_row = query_index(query_begin + query);
        const std::size_t offset = static_cast<std::size_t>(query) * columns;
        for (arma::uword row = 0; row < n; ++row) {
            double squared_distance = 0.0;
            for (arma::uword column = 0; column < selected.n_elem; ++column) {
                const arma::uword variable = selected(column);
                const double difference =
                    predictors(query_row, variable) - predictors(row, variable);
                squared_distance += difference * difference;
            }
            distances[offset + row] = squared_distance;
        }
    }
    return distances;
}

forward_candidate_result forward_candidate_block(
    const arma::mat& predictors,
    const arma::vec& response,
    const arma::uvec& query_index,
    const arma::uvec& selected,
    arma::uword candidate,
    const std::vector<double>& base_distance,
    arma::uword query_begin,
    arma::uword query_count,
    arma::uword k,
    double lambda,
    double initial_sse
) {
    arma::uvec columns(selected.n_elem + 1);
    if (selected.n_elem > 0) {
        columns.head(selected.n_elem) = selected;
    }
    columns(selected.n_elem) = candidate;
    columns = arma::sort(columns);
    const arma::mat subset_predictors = predictors.cols(columns);
    const arma::uword n = predictors.n_rows;
    std::vector<distance_index> ranked(n);
    // Carry the running sum, rather than reducing block sums, to retain the
    // original per-candidate floating-point addition order over all queries.
    double total_sse = initial_sse;
    bool boundary_tie = false;

    for (arma::uword query = 0; query < query_count; ++query) {
        const arma::uword query_row = query_index(query_begin + query);
        const std::size_t offset = static_cast<std::size_t>(query) * n;
        const double query_value = predictors(query_row, candidate);
        for (arma::uword row = 0; row < n; ++row) {
            const double difference = query_value - predictors(row, candidate);
            double distance = base_distance[offset + row] +
                difference * difference;
            if (row == query_row) {
                distance = std::numeric_limits<double>::infinity();
            }
            ranked[row] = std::make_pair(distance, row);
        }

        std::nth_element(
            ranked.begin(), ranked.begin() + k, ranked.end(),
            distance_index_less
        );
        std::sort(
            ranked.begin(), ranked.begin() + k, distance_index_less
        );
        const double included_distance = ranked[k - 1].first;
        const double excluded_distance = ranked[k].first;
        const double tie_scale = std::max(
            1.0,
            std::max(
                std::abs(included_distance),
                std::abs(excluded_distance)
            )
        );
        if (std::abs(included_distance - excluded_distance) <=
            32.0 * std::numeric_limits<double>::epsilon() * tie_scale) {
            boundary_tie = true;
        }
        arma::uvec neighbors(k);
        for (arma::uword position = 0; position < k; ++position) {
            neighbors(position) = ranked[position].second;
        }
        const arma::mat design = local_design(
            subset_predictors,
            neighbors,
            subset_predictors.row(query_row)
        );
        total_sse += ridge_sse(
            design,
            response.elem(neighbors),
            lambda
        );
    }

    return forward_candidate_result{total_sse, boundary_tie};
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List local_variance_cpp(
    const arma::mat& predictors,
    const arma::vec& response,
    const arma::umat& neighbors,
    double lambda = 1e-6
) {
    validate_inputs(predictors, response, neighbors);
    const arma::uword n = predictors.n_rows;
    const arma::uword d = predictors.n_cols;
    const double residual_df = static_cast<double>(
        neighbors.n_cols - d - 1
    );
    double total_sse = 0.0;

    for (arma::uword i = 0; i < n; ++i) {
        const arma::uvec index = neighbors.row(i).t();
        const arma::mat design = local_design(
            predictors,
            index,
            predictors.row(i)
        );
        total_sse += ridge_sse(
            design,
            response.elem(index),
            lambda
        );
    }

    const double mean_sse = total_sse / static_cast<double>(n);
    return Rcpp::List::create(
        Rcpp::Named("risk") = mean_sse / residual_df,
        Rcpp::Named("mean_sse") = mean_sse,
        Rcpp::Named("residual_df") = residual_df
    );
}

// [[Rcpp::export]]
Rcpp::List query_variance_cpp(
    const arma::mat& predictors,
    const arma::vec& response,
    const arma::uvec& query_index,
    const arma::umat& neighbors,
    double lambda = 1e-6
) {
    validate_query_inputs(predictors, response, query_index, neighbors);
    const arma::uword m = query_index.n_elem;
    const arma::uword d = predictors.n_cols;
    const double residual_df = static_cast<double>(
        neighbors.n_cols - d - 1
    );
    double total_sse = 0.0;
    arma::vec contribution(m);

    for (arma::uword i = 0; i < m; ++i) {
        const arma::uvec index = neighbors.row(i).t();
        const arma::mat design = local_design(
            predictors,
            index,
            predictors.row(query_index(i))
        );
        const double local_sse = ridge_sse(
            design,
            response.elem(index),
            lambda
        );
        total_sse += local_sse;
        contribution(i) = local_sse / residual_df;
    }

    const double mean_sse = total_sse / static_cast<double>(m);
    return Rcpp::List::create(
        Rcpp::Named("risk") = mean_sse / residual_df,
        Rcpp::Named("mean_sse") = mean_sse,
        Rcpp::Named("residual_df") = residual_df,
        Rcpp::Named("contribution") = contribution
    );
}

// [[Rcpp::export]]
Rcpp::List forward_values_cpp(
    const arma::mat& predictors,
    const arma::vec& response,
    const arma::uvec& query_index,
    const arma::uvec& selected,
    const arma::uvec& candidates,
    int k,
    double lambda,
    double response_variance,
    int threads = 1
) {
    if (predictors.n_rows != response.n_elem) {
        Rcpp::stop("Predictor rows and response length differ.");
    }
    if (query_index.n_elem == 0 || query_index.max() >= predictors.n_rows) {
        Rcpp::stop("Query index is empty or outside the data range.");
    }
    if ((selected.n_elem > 0 && selected.max() >= predictors.n_cols) ||
        candidates.n_elem == 0 || candidates.max() >= predictors.n_cols) {
        Rcpp::stop("Selected or candidate index is outside the data range.");
    }
    if (k <= static_cast<int>(selected.n_elem) + 2 ||
        k >= static_cast<int>(predictors.n_rows)) {
        Rcpp::stop("Neighborhood size is invalid for forward evaluation.");
    }
    if (!std::isfinite(lambda) || lambda < 0 ||
        !std::isfinite(response_variance) || response_variance <= 0) {
        Rcpp::stop("Lambda or response variance is invalid.");
    }

    // A shared block preserves distance reuse across candidates without an
    // all-query N-by-N allocation. This is an execution constant, not a
    // statistical tuning parameter: every query and candidate is evaluated.
    const arma::uword query_block_size = 256;
    const arma::uword query_count = query_index.n_elem;
    arma::vec total_sse(candidates.n_elem, arma::fill::zeros);
    arma::vec values(candidates.n_elem, arma::fill::zeros);
    arma::uvec boundary_ties(candidates.n_elem, arma::fill::zeros);
    const arma::uword worker_count = std::max<arma::uword>(
        1,
        std::min<arma::uword>(
            candidates.n_elem,
            static_cast<arma::uword>(std::max(1, threads))
        )
    );

    for (arma::uword query_begin = 0; query_begin < query_count;) {
        Rcpp::checkUserInterrupt();
        const arma::uword block_count = std::min<arma::uword>(
            query_block_size, query_count - query_begin
        );
        const std::vector<double> base_distance = selected_squared_distances_block(
            predictors, query_index, selected, query_begin, block_count
        );
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&](arma::uword worker_index) {
            try {
                for (arma::uword position = worker_index;
                     position < candidates.n_elem;
                     position += worker_count) {
                    const forward_candidate_result result = forward_candidate_block(
                        predictors,
                        response,
                        query_index,
                        selected,
                        candidates(position),
                        base_distance,
                        query_begin,
                        block_count,
                        static_cast<arma::uword>(k),
                        lambda,
                        total_sse(position)
                    );
                    total_sse(position) = result.total_sse;
                    if (result.boundary_tie) boundary_ties(position) = 1U;
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
        };

        if (worker_count == 1) {
            worker(0);
        } else {
            std::vector<std::thread> workers;
            workers.reserve(worker_count);
            try {
                for (arma::uword index = 0; index < worker_count; ++index) {
                    workers.emplace_back(worker, index);
                }
            } catch (...) {
                // Join already-started workers even if thread creation fails.
                for (std::thread& thread : workers) {
                    if (thread.joinable()) thread.join();
                }
                throw;
            }
            for (std::thread& thread : workers) thread.join();
        }
        if (worker_error) std::rethrow_exception(worker_error);
        query_begin += block_count;
    }
    const double residual_df = static_cast<double>(
        k - static_cast<int>(selected.n_elem + 1) - 1
    );
    for (arma::uword position = 0; position < candidates.n_elem; ++position) {
        const double mean_sse = total_sse(position) / static_cast<double>(query_count);
        const double risk = mean_sse / residual_df;
        values(position) = (response_variance - risk) / response_variance;
    }
    return Rcpp::List::create(
        Rcpp::Named("values") = values,
        Rcpp::Named("boundary_ties") = boundary_ties
    );
}

// [[Rcpp::export]]
arma::mat query_slopes_cpp(
    const arma::mat& predictors,
    const arma::vec& response,
    const arma::mat& query_predictors,
    const arma::umat& neighbors,
    double lambda = 1e-6
) {
    validate_external_query_inputs(
        predictors, response, query_predictors, neighbors
    );
    const arma::uword m = query_predictors.n_rows;
    const arma::uword d = predictors.n_cols;
    arma::mat slopes(m, d, arma::fill::zeros);

    for (arma::uword i = 0; i < m; ++i) {
        const arma::uvec index = neighbors.row(i).t();
        const arma::mat design = local_design(
            predictors,
            index,
            query_predictors.row(i)
        );
        const arma::vec coefficient = ridge_coefficients(
            design,
            response.elem(index),
            lambda
        );
        slopes.row(i) = coefficient.subvec(1, d).t();
    }
    return slopes;
}
