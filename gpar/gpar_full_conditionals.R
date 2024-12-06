# GPAR(1) - Full conditional posteriors
# Author: Cleiton Moya de Almeida
#####

graphics.off()  # close all plots
rm(list=ls())	# clear the (global) environment variables
cat("\014") 	# clear the console (Ctrl+L)

# Read the data
x_df <- read.csv("x.csv", header = TRUE)
x <- as.vector(x_df$x)

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
    logc1 <- (g-1)*log(xi) - h*xi
    
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

# Model parameters
alpha <- 0.8
mu    <- 2
xi    <- 0.5

# Prior hyperparameters
a <- 0.01
b <- 0.01
c <- 0.1
d <- 0.1
g <- 0.1
h <- 0.1

# Plot the full conditional posteriors
alpha_seq      <- seq(0.01, 0.99, 0.01)
mu_seq <- seq(0.1, 10, 0.1)
xi_seq  <- seq(0.1, 0.9, 0.01)

h1 <- function(z) logpost_alpha(alpha=z, mu, xi, x, a, b)
h2 <- function(z) logpost_mu(alpha, mu=z, xi, x, c, d)
h3 <- function(z) logpost_xi(alpha, mu, xi=z, x, g, h)

y1 <- sapply(alpha_seq, h1)
y2 <- sapply(mu_seq, h2)
y3 <- sapply(xi_seq, h3)

plot(alpha_seq, y1, type="l", xlab="alpha", ylab="log-posterior density")
plot(mu_seq, y2, type="l", xlab="mu", ylab="log-posterior density")
plot(xi_seq, y3, type="l", xlab="xi", ylab="log-posterior density")