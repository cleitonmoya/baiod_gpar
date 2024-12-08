# GPAR(1) Maximum Likelihood Estimation
# Author: Cleiton Moya de Almeida
#####

graphics.off()  # close all plots
rm(list=ls())	# clear the (global) environment variables
cat("\014") 	# clear the console (Ctrl+L)

# Read the data
x_df <- read.csv("data/x_120_poinar.csv", header = FALSE)
x <- as.vector(x_df$V1)
#x_df <- read.csv("earthquakes.csv", col.names = c("Ano", "Intensidade"), header=FALSE)
#x <- as.vector(x_df$Intensidade)


# Log-sum-exp auxiliary function
logsumexp <- function(x){
    c <- max(x)
    y <- c + log(sum((exp(x-c))))
    return(y)
}

# Auxiliary function to computhe the likelihood
f1 <- function(r, t, alpha, mu, xi, x){
	q <- 1-alpha
	y <- choose(x[t-1], r) * (alpha*mu + r*xi)^(r-1) * 
		(q*mu + (x[t-1]-r)*xi)^(x[t-1]-r-1)
	return(y)
}

logf1 <- function(r, t, alpha, mu, xi, x){
    q <- 1-alpha
    y <- log(choose(x[t-1], r)) + (r-1)*log(alpha*mu + r*xi) + 
        (x[t-1]-r-1)*log(q*mu + (x[t-1]-r)*xi)
    return(y)
}

# Auxiliary function to computhe the likelihood
f2 <- function(r, t, alpha, mu, xi, x){
	q <- 1-alpha
	num <- (q*mu + xi*(x[t]-r))^(x[t]-r-1) * exp(-xi*(x[t]-r))
	den <- factorial(x[t]-r)
	y <- num/den
	return(y)
}

logf2 <- function(r, t, alpha, mu, xi, x){
    q <- 1-alpha
    y <- (x[t]-r-1)*log(q*mu + xi*(x[t]-r)) -xi*(x[t]-r) - 
        log(factorial(x[t]-r))
    return(y)
}

# Likelhood function (not in logscale, just for debug)
likelihood <- function(alpha, mu, xi, x){
	q <- 1-alpha
	n <- length(x)
	c1 <- alpha^(n-1) * (q*mu)^(2*(n-1)) * exp(-q*mu*(n-1))
	
	prod <- 1
	for (t in 2:n){
		c2 <- (mu+x[t-1]*xi)^(1-x[t-1])
		Mt <- min(x[t], x[t-1])
		soma <- 0
		for (r in 0:Mt){
			soma <- soma + f1(r,t,alpha,mu,xi,x)*f2(r,t,alpha,mu,xi,x)
		}
		prod <- prod*c2*soma
	}
	prod <- c1*prod
	
	return(prod)
}

# Loglikelihood function without logsumexp
loglikelihood_old <- function(alpha, mu, xi, x){
    q <- 1-alpha
    n <- length(x)
    logc1 <- (n-1)*log(alpha) + 2*(n-1)*log(q*mu) - q*mu*(n-1)
    
    soma1 <- 0
    for (t in 2:n){
        logc2 <- -(x[t-1]-1)*log(mu + x[t-1]*xi)
        soma2 <- 0
        Mt <- min(x[t], x[t-1])
        for (r in 0:Mt){
            soma2 <- soma2 + f1(r,t,alpha,mu,xi,x)*f2(r,t,alpha,mu,xi,x)
        }
        soma1 <- soma1 + logc2 + log(soma2)
    }
    logl <- logc1 + soma1
    return(logl)
}

# Loglikelihood function (with logsumexp)
loglikelihood <- function(alpha, mu, xi, x){
    q <- 1-alpha
    n <- length(x)
    logc1 <- (n-1)*log(alpha) + 2*(n-1)*log(q*mu) - q*mu*(n-1)
    
    soma1 <- 0
    for (t in 2:n){
        logc2 <- -(x[t-1]-1)*log(mu + x[t-1]*xi)
        Mt <- min(x[t], x[t-1])
        r_seq <- 0:Mt
        
        g <- function(r) logf1(r, t, alpha, mu, xi, x) + logf2(r, t, alpha, mu, xi, x)
        logsoma2 <- logsumexp(sapply(r_seq, g))
        soma1 <- soma1 + logc2 + logsoma2
    }
    logl <- logc1 + soma1
    return(logl)
}

#####
# Compute the the likelihood and loglikelihood for conference
# Model parameters
alpha <- 0.8
mu    <- 2
xi    <- 0.5

logl1 <- log(likelihood(alpha, mu, xi, x))
logl2 <- loglikelihood(alpha, mu, xi, x)
print(logl1)
print(logl2)

# Plot the loglikelihood function holding the two other parameters
# (conditional loglikelihoods)
p_seq      <- seq(0.01, 0.99, 0.01)
lambda_seq <- seq(0.1, 5, 0.01)
theta_seq  <- seq(0.1, 0.9, 0.01)

h1 <- function(z) loglikelihood(alpha=z, mu, xi, x)
h2 <- function(z) loglikelihood(alpha, mu=z, xi, x)
h3 <- function(z) loglikelihood(alpha, mu, xi=z, x)

y1 <- sapply(p_seq, h1)
y2 <- sapply(lambda_seq, h2)
y3 <- sapply(theta_seq, h3)

plot(p_seq, y1, type="l")
plot(lambda_seq, y2, type="l")
plot(theta_seq, y3, type="l")

#####
# Compute the Maximum Likelihood Estimation (MLE)
# We use the R unconstrained non-linear optmization function nlm.
# The natural parameters of the model (alpha, xi and mu) are constrained. 
# So, before applying the nlm function, we must convert the natural parameters 
# to unsconstrained working parameters.

# Natural to working parameters
npar2wpar <- function(mu, xi, alpha){
    lambda_w <- log(mu)     # constrain: mu > 0
    theta_w  <- qlogis(xi)   # constrain: 0 <= xi < 1
    p_w      <- qlogis(alpha)       # constrain: 0 < alpha < 1 (qlogis: logit function)
    wpar     <- c(lambda_w, theta_w, p_w)
    return(wpar)
}

# Working to natural parameters
wpar2npar <- function(wpar){
    lambda_w <- wpar[1]
    theta_w  <- wpar[2]
    p_w      <- wpar[3]
    mu   <- exp(lambda_w)
    xi    <- plogis(theta_w)
    alpha        <- plogis(p_w) # logistic function
    return(list(mu=mu, xi=xi, alpha=alpha))
}

# Negative loglikelihood
nll <- function(wpar, x){
    npar <- wpar2npar(wpar)
    nll <- -loglikelihood(npar$alpha, npar$mu, npar$xi, x)
    return(nll)
}

#####
wpar <- npar2wpar(mu, xi, alpha) # working parameters
res <- nlm(nll, wpar, x)            # optmization results
mle_est <- wpar2npar(res$estimate)  # mle estimates
print(mle_est)
