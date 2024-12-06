# Generalized Poisson Distribution

library(gp)

set.seed(42)

# parâmetros da distribuição
N <- 1000
lambda <- 2
theta <- 0.5

mu <- lambda/(1-theta)
sigma2 <- lambda/(1-theta)^3

# gera uma amostra aleatória
# Note que, para biblioteca, GP(theta,lambda)
x <- rgp(n=N, lambda, theta, method="Branching")

# histograma e densidade estimada
hist(x, freq=FALSE, breaks=20)
lines(density(x), col="red")

# MLE
print(mu)
print(sigma2)
