# Generalized Poisson AR(1) - GPAR(1) Simulation
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

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

# (alpha, mu, sigma) <- (0.8, 2, 0.5)

# Parâmetros do processo
alpha <- 0.85
mu    <- 2
xi    <- 0.2 #0.2
q     <- 1-alpha

mean_x <- mu/(1-xi)
var_x <- mu/(1-xi)^3
sd_x <- sqrt(var_x)

# Simulação
T <- 120

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
    e <- rgp(n=2, q*mu, xi, method="Branching")[1]

    # Operador quasi-binomial thinning (seq. de v.a. QB)
    s <- rqb(1, alpha, xi/mu, x[t-1])
    x[t] <- s + e
}

#####
# Theoretical vs. sample mean and variance
printf("Theoretical mean: %.2f", mean_x)
printf("Sample mean: %.2f", mean(x))
printf("Theoretical var: %.2f", var_x)
printf("Sample var: %.2f", var(x))

#####
# Gráfico
par(mar = c(2, 2, 0.2, 0.2))
plot(x, type='s')
abline(h=mean_x, col="gray", lty=2)
abline(h=mean_x + 3*sd_x, col="gray", lty=3)
legend(x="topright", lty=c(2,3), col="gray", legend=c("E[X]", "E[X] + 3sigma[X]"))

#####
# Exporta os dados
write.csv(x, file = "data/x_100_mu0.csv", row.names = FALSE, col.names=FALSE)
