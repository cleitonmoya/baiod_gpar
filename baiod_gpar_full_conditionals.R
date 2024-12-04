# Bayesian Outlier Detection in GPAR(1)

# GPAR(1): Generalized Poisson AR(1) (Alzaid and Al-Osh, 1992):
#   X_t = S_t(X_{t-1}) + e_t, t=1, 2 , ...
#
#       e_t    ~ GP(q*mu, xi)        (Generalized Poisson)
#       S_t(.) ~ QB(alpha, xi/mu, .) (Quasi-Binomial)
#       X_t    ~ GP(mu, xi)          (Generalized Poisson)
#
#   Constraints:
#       0 < alpha < 1
#       q = 1-alpha
#       mu > 0
#       0 <= xi < 1
#
# Additive outliers modeling (Silva et al, 2018):
#	Y_t = X_t + eta_t*delta_t
#
#		eta_t|beta_t ~ Poisson(beta_t)
#		delta_t      ~ Bernoulli(p_t)
#
#	Constraints
#		beta_t > 0
#		0 < p_t < 1
#		eta_t, delta_t mutually independents
#
# Prior distributions:
#   alpha		 ~ Beta     (a, b)
#   mu  		 ~ Gamma    (c, d)
#   xi  		 ~ Beta     (g, h)
#   delta_t|p_t  ~ Bernoulli(p_t)
#	p_t			 ~ Beta     (l, m)
#	eta_t|beta_t ~ Poisson  (beta_t)
#	beta_t		 ~ Gamma    (v, w)
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

 
logf1 <- function(r, t, alpha, mu, xi, z){
	q <- 1-alpha
	logf1 <- log(choose(z[t-1], r)) + (r-1)*log(alpha*mu + r*xi) + 
		(z[t-1]-r-1)*log(q*mu + (z[t-1]-r)*xi)
	return(logf1)
}

log_f1 <- function(r, alpha, mu, xi, ys1){
    q <- 1-alpha
    logf1 <- log(choose(ys1, r)) + (r-1)*log(alpha*mu + r*xi) + 
        (ys1-r-1)*log(q*mu + (ys1-r)*xi)
    return(logf1)
}


logf2 <- function(r, t, alpha, mu, xi, z){
	q <- 1-alpha
	logf2 <- (z[t]-r-1)*log(q*mu + xi*(z[t]-r)) -xi*(z[t]-r) - 
		log(factorial(z[t]-r))
	return(logf2)
}

log_f2 <- function(r, alpha, mu, xi, ys){
    q <- 1-alpha
    logf2 <- (ys-r-1)*log(q*mu + xi*(ys-r)) -xi*(ys-r) - 
        log(factorial(ys-r))
    return(logf2)
}


# Full conditionallog-posterior of alpha
logpost_alpha <- function(alpha, mu, xi, delta, eta, y, a, b){
	z <- y - eta*delta
	q <- 1-alpha
	n <- length(y)
	logc1 <- (a+n-2)*log(alpha) + (b+2*n-3)*log(q) - mu*q*(n-1)
	
	soma1 <- 0
	for (t in 2:n){
		logc2 <- -(z[t-1]-1)*log(mu + z[t-1]*xi)
		Mt <- min(z[t], z[t-1])
		r_seq <- 0:Mt
		
		g <- function(r) logf1(r, t, alpha, mu, xi, z) + 
			logf2(r, t, alpha, mu, xi, z)
		logsoma2 <- logsumexp(sapply(r_seq, g))
		soma1 <- soma1 + logc2 + logsoma2
	}
	logpost <- logc1 + soma1
	return(logpost)
}

# Full conditional log-posterior of mu
logpost_mu <- function(alpha, mu, xi, delta, eta, y, c, d){
    z <- y - eta*delta
    q <- 1-alpha
    n <- length(y)
    logc1 <- (c+2*n-3)*log(mu) - mu*(d + q*(n-1))
    
    soma1 <- 0
    for (t in 2:n){
        logc2 <- -(z[t-1]-1)*log(mu + z[t-1]*xi)
        Mt <- min(z[t], z[t-1])
        r_seq <- 0:Mt
        
        g <- function(r) logf1(r, t, alpha, mu, xi, z) + 
            logf2(r, t, alpha, mu, xi, z)
        logsoma2 <- logsumexp(sapply(r_seq, g))
        soma1 <- soma1 + logc2 + logsoma2
    }
    logpost <- logc1 + soma1
    return(logpost)
}

# Full conditional log-posterior of xi
logpost_xi <- function(alpha, mu, xi, delta, eta, y, g, h){
    z <- y - eta*delta
    n <- length(y)
    logc1 <- (g-1)*log(xi) - h*xi
    
    soma1 <- 0
    for (t in 2:n){
        logc2 <- -(z[t-1]-1)*log(mu + z[t-1]*xi)
        Mt <- min(z[t], z[t-1])
        r_seq <- 0:Mt
        
        g <- function(r) logf1(r, t, alpha, mu, xi, z) + 
            logf2(r, t, alpha, mu, xi, z)
        logsoma2 <- logsumexp(sapply(r_seq, g))
        soma1 <- soma1 + logc2 + logsoma2
    }
    logpost <- logc1 + soma1
    return(logpost)
}

# logf, logfyj1, logfyj0: Used to compute logpost_delta_j
# ys:  y*[t]
# ys1: y*[t-1]
logf <- function(ys, ys1, alpha, mu, xi){
    q <- 1-alpha
    logc <- log(alpha) + 2*log(q*mu) - q*mu -(ys1-1)*log(mu + xi*ys1)
    
    Mt <- min(ys, ys1)
    if (Mt < 0){
    	Mt <- 0
    }
    r_seq <- 0:Mt
    g <- function(r) log_f1(r, alpha, mu, xi, ys1) + 
        log_f2(r, alpha, mu, xi, ys)
    logsoma <- logsumexp(sapply(r_seq, g))
    
    logf <- logc + logsoma
    return(logf)
}

# needs only delta[j-1], eta[j-1], eta[j] e eta[j+1]
logfyj1 <- function(j, alpha, mu, xi, delta, eta, p, y){
    
	#ys_3  <- y[j+1]-eta[j+1]
	#ys1_3 <- y[j]  -eta[j]
	#ys_4  <- y[j+1]
	#ys1_4 <- y[j]  -eta[j]
	
	#logf3 <- logf(ys_3, ys1_3, alpha, mu, xi)
	#logf4 <- logf(ys_4, ys1_4, alpha, mu, xi)
	
	#logf1 <- logf(y[j]  -eta[j]  , y[j-1]-eta[j-1], alpha, mu, xi)
    #logf2 <- logf(y[j]  -eta[j]  , y[j-1]         , alpha, mu, xi)
    logf3 <- logf(y[j+1]-eta[j+1], y[j]  -eta[j]  , alpha, mu, xi)
    logf4 <- logf(y[j+1]         , y[j]  -eta[j]  , alpha, mu, xi)
    
    g3 <- logf3 + log(p[j+1])
    g4 <- logf4 + log(1-p[j+1])
    if (delta[j-1] == 1){
    	logf1 <- logf(y[j]  -eta[j]  , y[j-1]-eta[j-1], alpha, mu, xi)
        out <- logf1 + logsumexp(c(g3,g4))
    } else {
    	logf2 <- logf(y[j]  -eta[j]  , y[j-1]         , alpha, mu, xi)
        out <- logf2 + logsumexp(c(g3,g4))
    }
    
    return(out)
}

# needs only delta[j-1], eta[j-1], eta[j+1]
logfyj0 <- function(j, alpha, mu, xi, delta, eta, p, y){
    
    logf1 <- logf(y[j]           , y[j-1]-eta[j-1], alpha, mu, xi)
    logf2 <- logf(y[j]           , y[j-1]         , alpha, mu, xi)
    logf3 <- logf(y[j+1]-eta[j+1], y[j]           , alpha, mu, xi)
    logf4 <- logf(y[j+1]         , y[j]           , alpha, mu, xi)
    
    g3 <- logf3 + log(p[j+1])
    g4 <- logf4 + log(1-p[j+1])
    if (delta[j-1] == 1){
        out <- logf1 + logsumexp(c(g3,g4))
    } else {
        out <- logf2 + logsumexp(c(g3,g4))
    }
    
    return(out)
}

# Full conditional log-posterior of delta_j
# needs only delta[j-1], eta[j-1], eta[j], eta[j+1]
# note it can be evaluated only for j=2:(n-1)
logpost_delta_j <- function(j, alpha, mu, xi, delta, eta, p, y){
    g1 <- log(p[j])   + logfyj1(j, alpha, mu, xi, delta, eta, p, y)
    g2 <- log(1-p[j]) + logfyj0(j, alpha, mu, xi, delta, eta, p, y)
    logpost <- g1 - logsumexp(c(g1, g2))
    return(logpost)
}

# Full conditional log-posterior of p_j
# needs only delta[j]
logpost_p_j <- function(j, delta, p_j, l, m){
	logpost <- (delta[j]+l-1)*log(p_j) + (m-delta[j])*log(1-p_j)
	return(logpost)
}

# Full conditional log-posterior of eta_j
logpost_eta_j <- function(j, alpha, mu, xi, delta, eta, eta_j, p, beta){
	logpi <- -beta[j] + eta_j*log(beta[j]) -log(factorial(eta_j))
	eta_ <- eta
	eta_[j] = eta_j
	logpost <- logpi + logfyj1(j, alpha, mu, xi, delta, eta_, p, y)
	return(logpost)
} 

# Full conditional log-posterior of beta_j
logpost_beta_j <- function(j, eta, beta_j, v, w){
	logpost <- -beta_j*(w+1) + (eta[j]+v-1)*log(beta_j)
	return(logpost)
}

# Read the data
x_df <- read.csv("x.csv", header = TRUE)
x <- as.vector(x_df$x)
n <- length(x)

# Insert the outliers
tau <- c(7, 26, 60, 90, 91)
eta_tau <- 9
y <- x
y[tau] <- y[tau] + eta_tau
plot(y, type="l")
points(tau, y[tau], col="red")

# Prepare delta and eta vectors
delta <- rep(0,n)
eta <- rep(0,n)
delta[tau] <- 1
eta[tau] <- eta_tau

# Prepare p and beta vectors
p    <- rep(0.05, n)
beta <- rep(10, n)

# Model parameters
alpha <- 0.8
mu    <- 2
xi    <- 0.5

# Prior hyperparameters
a <- 0.01 # alpha
b <- 0.01 # alpha
c <- 0.1  # mu
d <- 0.1  # mu	
g <- 0.01 # xi
h <- 0.01 # xi
l <- 5    # p_t 	
m <- 95   # p_t 
v <- 10   # eta_t
w <- 1    # eta_t

# Plot the full conditional posterior densities
alpha_seq <- seq(0.01, 0.99, 0.01)
mu_seq    <- seq(0.1, 10, 0.1)
xi_seq    <- seq(0.1, 0.9, 0.01)
j_seq     <- seq(2, n-1)
pj_seq	  <- seq(0.01, 0.99, 0.01)
eta7_seq  <- seq(0, y[7])
betaj_seq <- seq(0.1, 16, 0.1)

h1 <- function(z) logpost_alpha(alpha=z, mu, xi, delta, eta, y, a, b)
h2 <- function(z) logpost_mu(alpha, mu=z, xi, delta, eta, y, c, d)
h3 <- function(z) logpost_xi(alpha, mu, xi=z, delta, eta, y, g, h)

y1 <- sapply(alpha_seq, h1)
y2 <- sapply(mu_seq, h2)
y3 <- sapply(xi_seq, h3)

y4 <- rep(0, n)
for (j in j_seq){
    y4[j] <- exp(logpost_delta_j(j, alpha, mu, xi, delta, eta, p, y))
}

h5 <- function(z) logpost_p_j(7,  delta, p_j=z, l, m)
y5 <- sapply(pj_seq, h5)

h6 <- function(z) logpost_p_j(10, delta, p_j=z, l, m)
y6 <- sapply(pj_seq, h6)

h7 <- function(z) exp(logpost_eta_j(7, alpha, mu, xi, delta, eta, eta_j=z, p, beta))
y7 <- sapply(eta7_seq, h7)
y7 <- y7/sum(y7)

h8 <- function(z) logpost_beta_j(7, eta, beta_j=z, v, w)
y8 <- sapply(betaj_seq, h8)

h9 <- function(z) logpost_beta_j(10, eta, beta_j=z, v, w)
y9 <- sapply(betaj_seq, h9)

plot(alpha_seq, y1, type="l",xlab="alpha", ylab="log-posterior density")
plot(mu_seq, y2, type="l", xlab="mu", ylab="log-posterior density")
plot(xi_seq, y3, type="l", xlab="xi", ylab="log-posterior density")
plot(1:n, y4, type="h", xlab="t", ylab="posterior density", main="delta_t")
plot(pj_seq, y5, type="l", ylab="log-posterior density", main="p_7")
plot(pj_seq, y6, type="l", ylab="log-posterior density", main="p_10")
plot(eta7_seq, y7, type="b", xlab="eta_7", ylab="posterior pmf")
plot(betaj_seq, y8, type="l", xlab="beta_7", ylab="log-posterior density")
plot(betaj_seq, y9, type="l", xlab="beta_10", ylab="log-posterior density")