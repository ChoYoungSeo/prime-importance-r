test_that("single-predictor ridge coefficients use MASS centering and RMS scale", {
  # A dimension-preserving copy of MASS is an independent numerical reference.
  reference <- MASS::lm.ridge
  original_body <- body(MASS::lm.ridge)
  body(reference) <- parse(text = gsub(
    "X[, -Inter]", "X[, -Inter, drop = FALSE]",
    paste(deparse(body(reference)), collapse = "\n"), fixed = TRUE
  ))[[1L]]
  X <- matrix(13 + 3 * seq(-1, 1, length.out = 60), ncol = 1L,
              dimnames = list(NULL, "x"))
  y <- 2 * X[, 1L] + sin(X[, 1L])
  for (rescale in c(TRUE, FALSE)) {
    data <- prime:::.prime_prepare(X, y, rescale, NULL, 1e-6)
    for (ridge_lambda in c(1e-6, 0.3)) {
      data$ridge_lambda <- ridge_lambda
      # Evaluate the actual initializer prefix through its fit assignment.
      # The remaining initializer remains the shared weight/BIC path.
      statements <- as.list(body(prime:::.prime_ridge_initial))[-1L]
      fit_position <- which(vapply(statements, function(x) {
        is.call(x) && identical(x[[1L]], quote(`<-`)) &&
          identical(x[[2L]], quote(fit))
      }, logical(1)))
      expect_length(fit_position, 1L)
      scope <- new.env(parent = environment(prime:::.prime_ridge_initial))
      scope$data <- data
      for (statement in statements[seq_len(fit_position)]) eval(statement, scope)
      scaled_y <- (data$y - mean(data$y)) / sqrt(data$response_variance)
      frame <- data.frame(y = scaled_y, data$X, check.names = FALSE)
      expected <- reference(y ~ ., data = frame, lambda = ridge_lambda)
      expect_equal(as.numeric(scope$fit$coef), as.numeric(expected$coef),
                   tolerance = 1e-12)
      initialized <- prime:::.prime_ridge_initial(data)
      models <- list(stats::lm(y ~ 1, data = frame), stats::lm(y ~ ., data = frame))
      expected_size <- which.min(vapply(models, stats::BIC, numeric(1))) - 1L
      expect_identical(initialized$chosen_size, expected_size)
    }
  }
  expect_identical(body(MASS::lm.ridge), original_body)
})

test_that("one predictor supports positive allocation and preserves empty-fit sentinel", {
  X <- matrix(seq(-1, 1, length.out = 100), ncol = 1L,
              dimnames = list(NULL, "x"))
  for (rescale in c(TRUE, FALSE)) {
    positive <- prime(X, 2 * X[, 1L] + X[, 1L]^2, rescale = rescale)
    expect_identical(positive$selected, 1L)
    expect_equal(unname(positive$effects), 1)
    expect_equal(positive$subset_values$projected, c(0, 1))
  }
  empty <- prime(X, sin(20 * pi * X[, 1L]))
  expect_length(empty$selected, 0L)
  expect_null(empty$subset_values)
  expect_equal(sum(empty$effects), 0)
  expect_equal(empty$signal, 0)
  expect_true(all(empty$screen$forward$trace$dimension == 1L))
  expect_true(all(!empty$screen$forward$trace$accepted))
})

test_that("ridge screening is invariant to response-colliding and non-syntactic names", {
  X <- cbind(alpha = seq(-1, 1, length.out = 100),
             beta = ((seq_len(100) * 37) %% 101) / 101)
  y <- 2 * X[, 1L] + 0.4 * X[, 2L]
  neutral <- prime(X, y)
  renamed <- X
  colnames(renamed) <- c("y", "non syntactic label")
  fitted <- prime(renamed, y)
  expect_identical(fitted$selected, neutral$selected)
  expect_identical(unname(fitted$effects), unname(neutral$effects))
  expect_identical(fitted$screen$forward$ridge_path,
                   neutral$screen$forward$ridge_path)
  expect_identical(names(fitted$effects), colnames(renamed))
  expect_identical(fitted$selected_names, colnames(renamed)[fitted$selected])

  # Internal neutral names preserve the complete numeric design for names
  # that the original formula interface already handled correctly.
  prepared <- prime:::.prime_prepare(X, y, TRUE, NULL, 1e-6)
  original_frame <- data.frame(y = prepared$y, prepared$X, check.names = FALSE)
  neutral_frame <- original_frame
  names(neutral_frame) <- c("y", "X1", "X2")
  expect_identical(as.numeric(stats::model.matrix(y ~ ., original_frame)),
                   as.numeric(stats::model.matrix(y ~ ., neutral_frame)))
  original_ridge <- MASS::lm.ridge(y ~ ., original_frame, lambda = 1e-6)
  neutral_ridge <- MASS::lm.ridge(y ~ ., neutral_frame, lambda = 1e-6)
  expect_identical(as.numeric(original_ridge$coef), as.numeric(neutral_ridge$coef))
})
