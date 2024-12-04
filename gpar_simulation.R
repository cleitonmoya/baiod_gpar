# Generalized Poisson AR(1) - GPAR(1)
#
# X_t = S_t(X_{t-1}) + e_t, t=1, 2 , ...
#
# e_t    ~ GP(q*mu, xi)
# S_t(.) ~ QB(alpha, xi/mu, .)
# X_t    ~ GP(mu, xi)
#
# Constraints:
#   0 < alpha < 1
#   mu > 0
#   xi >= 0

library(gp)

graphics.off()  # fecha todos os gráficos
rm(list=ls())	# limpa o workspace
cat("\014") 	# limpa o console
set.seed(42)

# Parâmetros do processo
alpha <- 0.8
mu    <- 2
xi    <- 0.5
q <- 1-alpha
mu <- mu/(1-xi)
sigma2 <- mu/(1-xi)^3

# Simulação
T <- 100

# Quasi-binomial thinning pmf
dqb <- function(x, alpha, xi, n){
    q <- 1-alpha
    num <- alpha*q*choose(n,x) * (alpha + x*xi)^(x-1) * (q + (n-x)*xi)^(n-x-1)
    den <- (1 + n*xi)^(n-1)
    y <- num/den
    return(y)
}

# Quasi-binomial thinning cdf
pqb <- function(x, alpha, xi, n){
    x_ <- 0:x
    y <- sum(sapply(x_, dqb, alpha=alpha, xi=xi, n=n))
    return(y)
}

# Simula valores da quase-binomial
rqb <- function(size, alpha, xi, n){
    x <- 0:n
    w <- sapply(x, dqb, alpha=alpha, xi=xi, n=n)
    y <- sample(x, size=size, replace=TRUE, prob=w)
    return(y)
}

# Simula o processo GPAR(1)
x <- NULL
x[1] <- mu
for (t in 2:T){
    # Gera o erro aletório (função da biblioteca gp)
    e <- rgp(n=2, xi=xi, mu=q*mu, method="Branching")[1]
    
    # Operador quasi-binomial thinning (seq. de v.a. QB)
    s <- rqb(1, alpha, xi/mu, x[t-1])
    x[t] <- s + e
}

#####
# Gráfico
par(mar = c(2, 2, 0.2, 0.2))
plot(x, type='s')
abline(h=mu, col="gray", lty=2)

#####
# Exporta os dados
write.csv(x, file = "x.csv", row.names = FALSE)