// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "vendor/nanoflann/nanoflann.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cmath>
#include <exception>
#include <limits>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace {
using Neighbor = std::pair<double, arma::uword>;

// A bounded heap avoids the O(K) insertion cost of a sorted KNN result set.
// The row-index secondary key makes ties deterministic, including the kth
// boundary. The query itself is excluded by identity, not by zero distance.
class NeighborHeap {
    std::size_t capacity_;
    arma::uword self_;
  public:
    using DistanceType = double;
    using IndexType = arma::uword;
    std::vector<Neighbor> values;
    explicit NeighborHeap(std::size_t capacity) : capacity_(capacity), self_(0) {
        values.reserve(capacity);
    }
    void reset(arma::uword self) { values.clear(); self_ = self; }
    std::size_t size() const { return values.size(); }
    bool empty() const { return values.empty(); }
    bool full() const { return values.size() == capacity_; }
    double worstDist() const {
        return full() ? std::nextafter(values.front().first,
            std::numeric_limits<double>::infinity()) :
            std::numeric_limits<double>::infinity();
    }
    bool addPoint(double distance, arma::uword index) {
        if (index == self_) return true;
        const Neighbor point(distance, index);
        if (!full()) {
            values.push_back(point);
            std::push_heap(values.begin(), values.end());
        } else if (point < values.front()) {
            std::pop_heap(values.begin(), values.end());
            values.back() = point;
            std::push_heap(values.begin(), values.end());
        }
        return true;
    }
    void sort() { std::sort(values.begin(), values.end()); }
};

struct MatrixPoints {
    const arma::mat& X;
    explicit MatrixPoints(const arma::mat& x) : X(x) {}
    std::size_t kdtree_get_point_count() const { return X.n_rows; }
    double kdtree_get_pt(std::size_t row, std::size_t col) const {
        return X(row, col);
    }
    template<class BBox> bool kdtree_get_bbox(BBox&) const { return false; }
};
using Metric = nanoflann::L2_Simple_Adaptor<double, MatrixPoints, double, arma::uword>;
using Tree = nanoflann::KDTreeSingleIndexAdaptor<Metric, MatrixPoints, -1, arma::uword>;

float search_epsilon(double eps) {
    const double factor = (1.0 + eps) * (1.0 + eps);
    const double squared = factor - 1.0;
    if (!std::isfinite(eps) || eps < 0 || !std::isfinite(squared) ||
        squared > std::numeric_limits<float>::max()) {
        Rcpp::stop("eps must be a finite nonnegative distance tolerance.");
    }
    float bounded_factor = static_cast<float>(factor);
    if (static_cast<double>(bounded_factor) > factor)
        bounded_factor = std::nextafter(bounded_factor, 0.0f);
    // Round the factor before subtraction: stepping through eps values can
    // take millions of iterations when eps is much smaller than float(1)'s ULP.
    float value = bounded_factor - 1.0f;
    if (static_cast<double>(1.0f + value) > factor)
        value = std::nextafter(value, 0.0f);
    return value;
}

void validate_indexed(const arma::mat& X, const arma::uvec& query, int k) {
    if (X.n_rows < 2 || X.n_cols < 1 || !X.is_finite())
        Rcpp::stop("X must be a finite matrix with at least two rows and one column.");
    if (!query.n_elem || query.max() >= X.n_rows)
        Rcpp::stop("Query index is empty or outside the data range.");
    if (k < 1 || static_cast<arma::uword>(k) >= X.n_rows)
        Rcpp::stop("k must be between 1 and n - 1.");
}

bool boundary_tie(const std::vector<Neighbor>& values, arma::uword k) {
    if (values.size() <= k) return false;
    const double left = values[k - 1].first, right = values[k].first;
    const double scale = std::max(1.0, std::max(std::abs(left), std::abs(right)));
    return std::abs(left - right) <= 32 * std::numeric_limits<double>::epsilon() * scale;
}

double indexed_sse(const arma::mat& X, const arma::vec& y,
                   const std::vector<Neighbor>& nearest, arma::uword row,
                   arma::uword k, double lambda) {
    arma::uvec index(k);
    for (arma::uword j = 0; j < k; ++j) index(j) = nearest[j].second;
    arma::mat design(k, X.n_cols + 1);
    design.col(0).ones();
    design.cols(1, X.n_cols) = X.rows(index);
    design.cols(1, X.n_cols).each_row() -= X.row(row);
    const arma::vec local_y = y.elem(index);
    arma::mat cross = design.t() * design;
    arma::vec rhs = design.t() * local_y;
    for (arma::uword j = 1; j < cross.n_cols; ++j) cross(j, j) += lambda;
    arma::vec coefficient;
    const bool solved = arma::solve(coefficient, cross, rhs, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);
    if (!solved) coefficient = arma::pinv(cross) * rhs;
    const arma::vec residual = local_y - design * coefficient;
    return std::max(0.0, arma::dot(residual, residual));
}

struct RiskResult {
    double total_sse = 0;
    std::vector<double> contribution;
    std::vector<int> ties;
    bool any_tie = false;
};

RiskResult indexed_risk(const arma::mat& X, const arma::vec& y,
                        const arma::uvec& query, arma::uword k, double lambda,
                        float eps, bool keep, int block, bool main_thread,
                        const std::atomic<bool>* cancelled = nullptr) {
    MatrixPoints points(X);
    Tree tree(static_cast<int>(X.n_cols), points,
              nanoflann::KDTreeSingleIndexAdaptorParams(16));
    NeighborHeap nearest(std::min<arma::uword>(k + 1, X.n_rows - 1));
    const nanoflann::SearchParameters search(eps, false);
    RiskResult result;
    if (keep) result.contribution.resize(query.n_elem);
    if (main_thread) result.ties.resize(query.n_elem);
    const double df = static_cast<double>(k - X.n_cols - 1);
    arma::rowvec query_point(X.n_cols);
    for (arma::uword i = 0; i < query.n_elem; ++i) {
        if (i % static_cast<arma::uword>(block) == 0) {
            if (cancelled && cancelled->load()) return result;
            if (main_thread) Rcpp::checkUserInterrupt();
        }
        const arma::uword row = query(i);
        query_point = X.row(row);
        nearest.reset(row);
        tree.findNeighbors(nearest, query_point.memptr(), search);
        nearest.sort();
        if (nearest.size() < k) throw std::runtime_error("Indexed search returned too few neighbors.");
        const bool tied = boundary_tie(nearest.values, k);
        result.any_tie = result.any_tie || tied;
        if (main_thread) result.ties[i] = tied;
        const double sse = indexed_sse(X, y, nearest.values, row, k, lambda);
        result.total_sse += sse;
        if (keep) result.contribution[i] = sse / df;
    }
    return result;
}
} // namespace

// Diagnostic API only: the production risk path consumes one neighborhood
// immediately and never materializes this query-by-K output.
// [[Rcpp::export]]
Rcpp::List indexed_neighbors_cpp(const arma::mat& X, const arma::uvec& query_index,
                                 int k, double eps = 0) {
    validate_indexed(X, query_index, k);
    const float epsilon = search_epsilon(eps);
    MatrixPoints points(X);
    Tree tree(static_cast<int>(X.n_cols), points,
              nanoflann::KDTreeSingleIndexAdaptorParams(16));
    NeighborHeap nearest(std::min<arma::uword>(k + 1, X.n_rows - 1));
    Rcpp::IntegerMatrix indices(query_index.n_elem, k);
    Rcpp::NumericMatrix distances(query_index.n_elem, k);
    Rcpp::LogicalVector ties(query_index.n_elem);
    arma::rowvec point(X.n_cols);
    for (arma::uword i = 0; i < query_index.n_elem; ++i) {
        if (i % 256 == 0) Rcpp::checkUserInterrupt();
        point = X.row(query_index(i));
        nearest.reset(query_index(i));
        tree.findNeighbors(nearest, point.memptr(), nanoflann::SearchParameters(epsilon, false));
        nearest.sort();
        if (nearest.size() < static_cast<std::size_t>(k)) Rcpp::stop("Too few neighbors.");
        for (int j = 0; j < k; ++j) {
            indices(i, j) = nearest.values[j].second;
            distances(i, j) = nearest.values[j].first;
        }
        ties[i] = boundary_tie(nearest.values, k);
    }
    return Rcpp::List::create(Rcpp::Named("indices") = indices,
        Rcpp::Named("squared_distances") = distances, Rcpp::Named("boundary_ties") = ties);
}

// [[Rcpp::export]]
Rcpp::List indexed_risk_cpp(const arma::mat& X, const arma::vec& y,
    const arma::uvec& query_index, int k, double lambda = 1e-6, double eps = 0,
    bool keep_contribution = true, int query_block_size = 256) {
    validate_indexed(X, query_index, k);
    if (y.n_elem != X.n_rows || !y.is_finite()) Rcpp::stop("Response is invalid.");
    if (k <= static_cast<int>(X.n_cols) + 1) Rcpp::stop("Neighborhood must exceed d + 1.");
    if (!std::isfinite(lambda) || lambda < 0) Rcpp::stop("lambda must be finite and nonnegative.");
    if (query_block_size < 1) Rcpp::stop("query_block_size must be positive.");
    const RiskResult result = indexed_risk(X, y, query_index, k, lambda,
        search_epsilon(eps), keep_contribution, query_block_size, true);
    const double mean_sse = result.total_sse / static_cast<double>(query_index.n_elem);
    const double df = k - static_cast<double>(X.n_cols) - 1;
    return Rcpp::List::create(Rcpp::Named("risk") = mean_sse / df,
        Rcpp::Named("mean_sse") = mean_sse, Rcpp::Named("residual_df") = df,
        Rcpp::Named("contribution") = result.contribution,
        Rcpp::Named("boundary_ties") = result.ties);
}

// [[Rcpp::export]]
Rcpp::List indexed_forward_values_cpp(const arma::mat& X, const arma::vec& y,
    const arma::uvec& query_index, const arma::uvec& selected,
    const arma::uvec& candidates, int k, double lambda, double response_variance,
    int threads = 1, double eps = 0, int query_block_size = 256) {
    validate_indexed(X, query_index, k);
    if (y.n_elem != X.n_rows || !y.is_finite()) Rcpp::stop("Response is invalid.");
    if (!candidates.n_elem || candidates.max() >= X.n_cols ||
        (selected.n_elem && selected.max() >= X.n_cols)) Rcpp::stop("Invalid coordinates.");
    if (k <= static_cast<int>(selected.n_elem) + 2) Rcpp::stop("Neighborhood must exceed d + 1.");
    if (!std::isfinite(lambda) || lambda < 0 || !std::isfinite(response_variance) ||
        response_variance <= 0 || query_block_size < 1) Rcpp::stop("Invalid numerical controls.");
    const float epsilon = search_epsilon(eps);
    const arma::uword workers = std::max<arma::uword>(1,
        std::min<arma::uword>(candidates.n_elem, std::max(1, threads)));
    arma::vec scores(candidates.n_elem);
    arma::uvec ties(candidates.n_elem, arma::fill::zeros);
    std::exception_ptr error;
    std::mutex error_mutex;
    std::atomic<bool> cancelled(false);
    std::atomic<arma::uword> finished(0);
    std::condition_variable completion;
    std::mutex wait_mutex;
    auto worker = [&](arma::uword slot) {
        try {
            for (arma::uword j = slot; j < candidates.n_elem; j += workers) {
                if (cancelled.load()) break;
                arma::uvec columns(selected.n_elem + 1);
                if (selected.n_elem) columns.head(selected.n_elem) = selected;
                columns(selected.n_elem) = candidates(j);
                columns = arma::sort(columns);
                const arma::mat subset = X.cols(columns);
                const RiskResult result = indexed_risk(subset, y, query_index, k,
                    lambda, epsilon, false, query_block_size, workers == 1, &cancelled);
                const double risk = result.total_sse /
                    static_cast<double>(query_index.n_elem) /
                    static_cast<double>(k - columns.n_elem - 1);
                scores(j) = (response_variance - risk) / response_variance;
                ties(j) = result.any_tie;
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(error_mutex);
            if (!error) error = std::current_exception();
            cancelled.store(true);
        }
        ++finished;
        completion.notify_one();
    };
    Rcpp::checkUserInterrupt();
    if (workers == 1) worker(0);
    else {
        std::vector<std::thread> pool;
        try { for (arma::uword j = 0; j < workers; ++j) pool.emplace_back(worker, j); }
        catch (...) {
            cancelled.store(true);
            for (auto& thread : pool) if (thread.joinable()) thread.join();
            throw;
        }
        try {
            while (finished.load() < workers) {
                std::unique_lock<std::mutex> lock(wait_mutex);
                completion.wait_for(lock, std::chrono::milliseconds(50));
                Rcpp::checkUserInterrupt();
            }
        } catch (...) {
            cancelled.store(true);
            for (auto& thread : pool) thread.join();
            throw;
        }
        for (auto& thread : pool) thread.join();
    }
    if (error) std::rethrow_exception(error);
    Rcpp::checkUserInterrupt();
    return Rcpp::List::create(Rcpp::Named("values") = scores, Rcpp::Named("boundary_ties") = ties);
}
