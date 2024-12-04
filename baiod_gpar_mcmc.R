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
	
	logf1 <- logf(y[j]  -eta[j]  , y[j-1]-eta[j-1], alpha, mu, xi)
    logf2 <- logf(y[j]  -eta[j]  , y[j-1]         , alpha, mu, xi)
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
logpost_eta_j <- function(j, alpha, mu, xi, delta, eta, eta_j, beta){
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
		printf("Step %d/%d: acceptance ratio: alpha: %.2f, mu: %.2f, xi: %.2f", 
			   n, N, na_alpha/n, na_mu/n, na_xi/n)
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

