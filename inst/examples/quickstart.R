library(prime)

n <- 300
index <- seq_len(n)
frequencies <- sqrt(c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29))
base <- 2 * (outer(index, frequencies) %% 1) - 1
rho <- 0.9
X <- cbind(
  base[, 1],
  base[, 2],
  rho * base[, 1] + sqrt(1 - rho^2) * base[, 3],
  rho * base[, 2] + sqrt(1 - rho^2) * base[, 4],
  base[, 5],
  base[, 6],
  base[, 7],
  base[, 8]
)
colnames(X) <- paste0("X", seq_len(ncol(X)))
y <- X[, 1] - X[, 2] + X[, 5]^2 + 0.08 * sin(index * sqrt(31))

fit <- prime(X, y)
print(fit)
plot(fit)
