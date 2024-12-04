# Bayesian Inference of GPAR(1)
# GPAR(1): Generalized Poisson AR(1) 
#
# Model (Alzaid and Al-Osh, 1992):
#   X_t = S_t(X_{t-1}) + e_t, t=1, 2 , ...
#
#       e_t    ~ GP(q*mu, xi)    (Generalized Poisson)
#       S_t(.) ~ QB(alpha, xi/mu, .) (Quasi-Binomial)
#       X_t    ~ GP(mu, xi)      (Generalized Poisson)
#
#   Constraints:
#       0 < alpha < 1
#       q = 1-alpha
#       mu > 0
#       0 <= xi < 1
#
# Prior distriubtions:
#   alpha ~ Beta(a,b)
#   mu    ~ Gamma(c,d)
#   xi    ~ Beta(g,h)
#
# MCMC Method: 
#   Gibbs Sampling with Random Walk Metropolis
#
# Author: Cleiton Moya de Almeida
#####

graphics.off()  # close all plots
rm(list=ls())	# clear the (global) environment variables
cat("\014") 	# clear the console (Ctrl+L)
set.seed(42)

# Read the data
x_df <- read.csv("x.csv", header = TRUE)
x <- as.vector(x_df$x)

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

# Log-sum-exp auxiliary function
logsumexp <- function(x){
    c <- max(x)
    y <- c + log(sum((exp(x-c))))
    return(y)
}

# Auxiliary functions to compute full conditional posteriors
logf1 <- function(r, t, alpha, mu, xi, x){
    q <- 1-alpha
    y <- log(choose(x[t-1], r)) + (r-1)*log(alpha*mu + r*xi) + 
        (x[t-1]-r-1)*log(q*mu + (x[t-1]-r)*xi)
    return(y)
}

logf2 <- function(r, t, alpha, mu, xi, x){
    q <- 1-alpha
    y <- (x[t]-r-1)*log(q*mu + xi*(x[t]-r)) -xi*(x[t]-r) - 
        log(factorial(x[t]-r))
    return(y)
}

# Full conditional log-posterior of alpha
logpost_alpha <- function(alpha, mu, xi, x, a, b){
    q <- 1-alpha
    n <- length(x)
    logc1 <- (a+n-2)*log(alpha) + (b+2*n-3)*log(q) - q*mu*(n-1)
    
    soma1 <- 0
    for (t in 2:n){
        logc2 <- -(x[t-1]-1)*log(mu + x[t-1]*xi)
        Mt <- min(x[t], x[t-1])
        r_seq <- 0:Mt
        
        g <- function(r) logf1(r, t, alpha, mu, xi, x) + logf2(r, t, alpha, mu, xi, x)
        logsoma2 <- logsumexp(sapply(r_seq, g))
        soma1 <- soma1 + logc2 + logsoma2
    }
    logpost <- logc1 + soma1
    return(logpost)
}

# Full conditional log-posterior of mu
logpost_mu <- function(alpha, mu, xi, x, c, d){
    q <- 1-alpha
    n <- length(x)
    logc1 <- (c+2*n-3)*log(mu) - mu*(d + q*(n-1))
    
    soma1 <- 0
    for (t in 2:n){
        logc2 <- -(x[t-1]-1)*log(mu + x[t-1]*xi)
        Mt <- min(x[t], x[t-1])
        r_seq <- 0:Mt
        
        g <- function(r) logf1(r, t, alpha, mu, xi, x) + logf2(r, t, alpha, mu, xi, x)
        logsoma2 <- logsumexp(sapply(r_seq, g))
        soma1 <- soma1 + logc2 + logsoma2
    }
    logpost <- logc1 + soma1
    return(logpost)
}

# Full conditional log-posterior of xi
logpost_xi <- function(alpha, mu, xi, x, g, h){
    n <- length(x)
    logc1 <- (g-1)*log(xi) + (h-1)*log(1-xi)
    
    soma1 <- 0
    for (t in 2:n){
        logc2 <- -(x[t-1]-1)*log(mu + x[t-1]*xi)
        Mt <- min(x[t], x[t-1])
        r_seq <- 0:Mt
        
        g <- function(r) logf1(r, t, alpha, mu, xi, x) + logf2(r, t, alpha, mu, xi, x)
        logsoma2 <- logsumexp(sapply(r_seq, g))
        soma1 <- soma1 + logc2 + logsoma2
    }
    logpost <- logc1 + soma1
    return(logpost)
}

# Auxiliary functions for Metropolis
# Simulate a proposed value (random walk)
sim_xprop <- function(x, sigma){
    x_prop <- rnorm(n=1, mean=x, sd=sigma)
    return(x_prop)
}

# Acceptance probability of alpha_prop
logprob_accept_alpha <- function(alpha_prop, alpha, mu, xi, x, a, b){
    p1 <- logpost_alpha(alpha_prop, mu, xi, x, a, b)
    p2 <- logpost_alpha(alpha, mu, xi, x, a, b)
    if (p1 < p2){
        prob_accept <- min(1, exp(p1-p2))
    } else {
        prob_accept = 1
    }
    return(prob_accept)
}

# Acceptance probability of mu_prop
logprob_accept_mu <- function(mu_prop, mu, alpha, xi, x, c, d){
    p1 <- logpost_mu(alpha, mu_prop, xi, x, c, d)
    p2 <- logpost_mu(alpha, mu, xi, x, c, d)
    if (p1 < p2){
        prob_accept <- min(1, exp(p1-p2))
    } else {
        prob_accept = 1
    }
    return(prob_accept)
}

# Acceptance probability of xi_prop
logprob_accept_xi <- function(xi_prop, xi, alpha, mu, x, g, h){
    p1 <- logpost_xi(alpha, mu, xi_prop, x, g, h)
    p2 <- logpost_xi(alpha, mu, xi, x, g, h)
    if (p1 < p2){
        prob_accept <- min(1, exp(p1-p2))
    } else {
        prob_accept = 1
    }
    return(prob_accept)
}

# Samples alpha
ralpha <- function(p0, sigma, mu, xi, x, a, b){
    
    ac <- 0 # accepted flag
    alpha <- p0
    
    # simulate a proposed value for alpha
    alpha_prop <- sim_xprop(alpha, sigma)
    
    # verify the variable constrain
    if (alpha_prop > 0 & alpha_prop < 1){
        
        # acceptance probability
        prob_accept <- logprob_accept_alpha(alpha_prop, alpha, mu, xi, x, a, b)
        
        # acceptance criteria
        u <- runif(1)
        if  (u < prob_accept){
            alpha <- alpha_prop
            ac <- 1
        }
    }
    return(list(alpha=alpha, ac=ac))
}

# Samples mu
rmu <- function(lambda0, sigma, alpha, xi, x, c, d){
    
    ac <- 0 # accepted flag
    mu <- lambda0
    
    # simulate a proposed value for alpha
    mu_prop <- sim_xprop(mu, sigma)
    
    # verify the variable constrain
    if (mu_prop > 0){
        
        # acceptance probability
        prob_accept <- logprob_accept_mu(mu_prop, mu, alpha, xi, x, c, d)
        
        # acceptance criteria
        u <- runif(1)
        if  (u < prob_accept){
            mu <- mu_prop
            ac <- 1
        }
    }
    return(list(mu=mu, ac=ac))
}

# Samples xi
rxi <- function(theta0, sigma, alpha, mu, x, g, h){
    
    ac <- 0 # accepted flag
    xi <- theta0
    
    # simulate a proposed value for alpha
    xi_prop <- sim_xprop(xi, sigma)
    
    # verify the variable constrain
    if (xi_prop >= 0 & xi_prop < 1){
        
        # acceptance probability
        prob_accept <- logprob_accept_xi(xi_prop, xi, alpha, mu, x, g, h)
        
        # acceptance criteria
        u <- runif(1)
        if  (u < prob_accept){
            xi <- xi_prop
            ac <- 1
        }
    }
    return(list(xi=xi, ac=ac))
}

#####
# Gibbs Sampling (with random walking Metropopolis steps)

N <- 1200  # simulation steps
burnin <- 200

# Parameters initilization
alpha <- 0.1
mu    <- 10
xi    <- 0.1

# Model hyperparameters (priors)
# alpha      ~ Beta(a,b)
# mu ~ Gamma(c,d)
# xi  ~ Beta(g.h)
a <- 0.01
b <- 0.01
c <- 0.1
d <- 0.1
g <- 0.01
h <- 0.01

# Random walking Metropolis parameters
sigma_alpha <- 0.1
sigma_mu    <- 0.2
sigma_xi    <- 0.05

# Accepted steps counters
na_alpha <- 0
na_mu <- 0 
na_xi <- 0 

# Vectors with traces
Alpha <- numeric(N)
Mu <- numeric(N)
Xi <- numeric(N)

Alpha[1] <- alpha
Mu[1] <- mu
Xi[1] <- xi

start_time = proc.time() # para cálculo do tempo de execução
for (n in 2:N){
    
    if (n %% 100 == 0){
        printf("Step %d/%d: acceptance ratio: alpha: %.2f, mu: %.2f, xi: %.2f", n, N, na_alpha/n, na_mu/n, na_xi/n)
    }
    
    # Samples alpha
    res <- ralpha(alpha, sigma_alpha, mu, xi, x, a, b)
    alpha <- res$alpha 
    na_alpha <- na_alpha + res$ac
    Alpha[n] <- alpha
    
    # Samples mu
    res <- rmu(mu, sigma_mu, alpha, xi, x, c, d)
    mu <- res$mu
    na_mu <- na_mu +  res$ac
    Mu[n] <- mu
    
    # Samples xi
    res <- rxi(xi, sigma_xi, alpha, mu, x, g, h)
    xi <- res$xi
    na_xi <- na_xi +  res$ac
    Xi[n] <- xi
}

# Execution time
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.0f s\n", elapsed_time)

# Acceptance ratios
ac_p      <- na_alpha/N
ac_lambda <- na_mu/N
ac_theta  <- na_xi/N
printf("Acceptance rate - alpha: %.2f", ac_p)
printf("Acceptance rate - mu: %.2f", ac_lambda)
printf("Acceptance rate - xi: %.2f", ac_theta)

# Values after burn in
Alpha_final <- Alpha[(burnin+1):N]
Mu_final <- Mu[(burnin+1):N]
Xi_final <- Xi[(burnin+1):N]

# Mean values
alpha_mean <- mean(Alpha_final)
mu_mean <- mean(Mu_final)
xi_mean <- mean(Xi_final)
printf("alpha mean: %.2f", alpha_mean)
printf("mu mean: %.2f", mu_mean)
printf("xi mean: %.2f", xi_mean)

#####
# Trace plots
par(mfrow = c(3, 1), mar=c(4,3,2,2))
plot(Alpha, type="l", main="alpha", xlab="", ylab="", ylim=c(0,1))
abline(v=burnin, col="gray", lty=2)
plot(Mu, type="l", main="mu", xlab="", ylab="", ylim=c(0,10))
abline(v=burnin, col="gray", lty=2)
plot(Xi, type="l", main="xi", xlab="step", ylab="", ylim=c(0,1))
abline(v=burnin, col="gray", lty=2)

# Posterior densities
par(mfrow=c(1,1))
plot(density(Alpha_final), main="alpha")
plot(density(Mu_final), main="mu")
plot(density(Xi_final), main="xi")