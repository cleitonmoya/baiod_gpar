# Generalized Poisson Distribution

library(gp)

set.seed(42)

# parâmetros da distribuição
N <- 1000
theta <- 5
lambda <- 0.5
mu <- theta/(1-lambda)

# gera uma amostra aleatória
x <- rgp(n=N, theta=theta, lambda=lambda, method="Branching")

# histograma e densidade estimada
hist(x, freq=FALSE, breaks=20)
lines(density(x), col="red")

# MLE
cat("mu:", mu)
cat("\nMLE:\n")
print(gp.mle(x))