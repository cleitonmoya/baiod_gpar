# Bayesian Outlier Detection in GPAR(1)
#
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
	T <- length(y)
	logc1 <- (a+T-2)*log(alpha) + (b+2*T-3)*log(q) - mu*q*(T-1)

	soma1 <- 0
	for (t in 2:T){
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
    T <- length(y)
    logc1 <- (c+2*T-3)*log(mu) - mu*(d + q*(T-1))

    soma1 <- 0
    for (t in 2:T){
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
    T <- length(y)
    logc1 <- (g-1)*log(xi) - h*xi

    soma1 <- 0
    for (t in 2:T){
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

    logf3 <- logf(y[j+1]-eta[j+1], y[j]           , alpha, mu, xi)
    logf4 <- logf(y[j+1]         , y[j]           , alpha, mu, xi)

    g3 <- logf3 + log(p[j+1])
    g4 <- logf4 + log(1-p[j+1])
    if (delta[j-1] == 1){
        logf1 <- logf(y[j]           , y[j-1]-eta[j-1], alpha, mu, xi)
        out <- logf1 + logsumexp(c(g3,g4))
    } else {
        logf2 <- logf(y[j]           , y[j-1]         , alpha, mu, xi)
        out <- logf2 + logsumexp(c(g3,g4))
    }

    return(out)
}

# Full conditional log-posterior of delta_j
# needs only delta[j-1], eta[j-1], eta[j], eta[j+1]
# note it can be evaluated only for j=2:(T-1)
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
logpost_eta_j <- function(j, alpha, mu, xi, delta, eta, eta_j, p, beta, y){
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
logprob_accept_alpha <- function(alpha_prop, alpha, mu, xi,
								 delta, eta, y, a, b){
	p1 <- logpost_alpha(alpha_prop, mu, xi, delta, eta, y, a, b)
	p2 <- logpost_alpha(alpha, mu, xi, delta, eta, y, a, b)
	if (p1 < p2){
		prob_accept <- min(1, exp(p1-p2))
	} else {
		prob_accept = 1
	}
	return(prob_accept)
}

# Acceptance probability of mu_prop
logprob_accept_mu <- function(mu_prop, mu, alpha, xi, delta, eta, y, c, d){
	p1 <- logpost_mu(alpha, mu_prop, xi, delta, eta, y, c, d)
	p2 <- logpost_mu(alpha, mu, xi, delta, eta, y, c, d)
	if (p1 < p2){
		prob_accept <- min(1, exp(p1-p2))
	} else {
		prob_accept = 1
	}
	return(prob_accept)
}

# Acceptance probability of xi_prop
logprob_accept_xi <- function(xi_prop, xi, alpha, mu, delta, eta, y, g, h){
	p1 <- logpost_xi(alpha, mu, xi_prop, delta, eta, y, g, h)
	p2 <- logpost_xi(alpha, mu, xi, delta, eta, y, g, h)
	if (p1 < p2){
		prob_accept <- min(1, exp(p1-p2))
	} else {
		prob_accept = 1
	}
	return(prob_accept)
}

# Acceptance probability of beta_j_prop
logprob_accept_beta_j <- function(j, beta_j_prop, beta_j, eta, v, w){
	p1 <- logpost_beta_j(j, eta, beta_j_prop, v, w)
	p2 <- logpost_beta_j(j, eta, beta_j, v, w)
	if (p1 < p2){
		prob_accept <- min(1, exp(p1-p2))
	} else {
		prob_accept = 1
	}
	return(prob_accept)
}

# Draw alpha (Random Walk Metropolis)
ralpha <- function(alpha0, sigma, mu, xi, delta, eta, y, a, b){

	ac <- 0 # accepted flag
	alpha <- alpha0

	# simulate a proposed value for alpha
	alpha_prop <- sim_xprop(alpha, sigma)

	# verify the variable constrain
	if (alpha_prop > 0 & alpha_prop < 1){

		# acceptance probability
		prob_accept <- logprob_accept_alpha(alpha_prop, alpha, mu, xi,
											delta, eta, y, a, b)

		# acceptance criteria
		u <- runif(1)
		if  (u < prob_accept){
			alpha <- alpha_prop
			ac <- 1
		}
	}
	return(list(alpha=alpha, ac=ac))
}

# Draw mu (Random Walk Metropolis)
rmu <- function(mu0, sigma, alpha, xi, delta, eta, y, c, d){

	ac <- 0 # accepted flag
	mu <- mu0

	# simulate a proposed value for alpha
	mu_prop <- sim_xprop(mu, sigma)

	# verify the variable constrain
	if (mu_prop > 0){

		# acceptance probability
		prob_accept <- logprob_accept_mu(mu_prop, mu, alpha, xi,
										 delta, eta, y, c, d)

		# acceptance criteria
		u <- runif(1)
		if  (u < prob_accept){
			mu <- mu_prop
			ac <- 1
		}
	}
	return(list(mu=mu, ac=ac))
}

# Draw xi (Random Walk Metropolis)
rxi <- function(xi0, sigma, alpha, mu, delta, eta, y, g, h){

	ac <- 0 # accepted flag
	xi <- xi0

	# simulate a proposed value for alpha
	xi_prop <- sim_xprop(xi, sigma)

	# verify the variable constrain
	if (xi_prop >= 0 & xi_prop < 1){

		# acceptance probability
		prob_accept <- logprob_accept_xi(xi_prop, xi, alpha, mu,
										 delta, eta, y, g, h)

		# acceptance criteria
		u <- runif(1)
		if  (u < prob_accept){
			xi <- xi_prop
			ac <- 1
		}
	}
	return(list(xi=xi, ac=ac))
}


# Draw delta (Bernoulli Distribution)
rdelta <- function(alpha, mu, xi, delta, eta, p, y){
	T <- length(y)
	j_seq <- seq(2, T-1)
	g <- function(z) exp(logpost_delta_j(j=z, alpha, mu, xi, delta, eta, p, y))
	prob_delta <- sapply(j_seq, g)
	prob_delta <- c(0, prob_delta, 0)
	delta[2:(T-1)] <- rbinom(T-2, 1, prob_delta)
	return(list(delta=delta, prob_delta=prob_delta))
}

# Draw eta_j (Discrete Distribution)
reta_j <- function (j, alpha, mu, xi, delta, eta, p, beta, y){
	etaj_seq <- seq(0, y[j])
	g <- function(z) exp(logpost_eta_j(j, alpha, mu, xi, delta, eta, eta_j=z, p, beta, y))
	eta_pmf <- sapply(etaj_seq, g)
	eta_pmf <- eta_pmf/sum(eta_pmf)
	eta_j <- sample(etaj_seq, 1, prob=eta_pmf)
	return(eta_j)
}

# Draw eta (only for delta_j = 1)
reta <- function(alpha, mu, xi, delta, eta, p, beta, y){
	j_seq <- which(delta==1) # sample eta lonly for j|delta_j=1
	if (length(j_seq) > 0) {
		g <- function(z) reta_j(j=z, alpha, mu, xi, delta, eta, p, beta, y)
		eta_1 <- sapply(j_seq, g)
		eta[j_seq] <- eta_1
	}
	return(eta)
}

# Draw p_j (Beta Distribution)
rp_j <- function(j, delta, l, m){
	pj <- rbeta(1, delta[j]+l, m-delta[j]+1)
	return(pj)
}

# Draw p
rp <- function(delta, l, m){
	T <- length(delta)
	j_seq = seq(1, T)
	g <- function(z) rp_j(j=z, delta, l, m)
	p <- sapply(j_seq, g)
	return(p)
}

# Draw beta_j (Random Walk Metropolis)
rbeta_j <- function(j, beta_j0, sigma, eta, v, w){
	ac <- 0 # accepted flag
	beta_j <- beta_j0

	# simulate a proposed value for beta_j
	beta_j_prop <- sim_xprop(beta_j, sigma)

	# verify the variable constrain
	if (beta_j_prop > 0){

		# acceptance probability
		prob_accept <- logprob_accept_beta_j(j, beta_j_prop, beta_j, eta, v, w)

		# acceptance criteria
		u <- runif(1)
		if  (u < prob_accept){
			beta_j <- beta_j_prop
			ac <- 1
		}
	}
	return(list(beta_j=beta_j, ac=ac))
}

# Draw beta
rbetavect <- function(x0, sigma, eta, v, w){
	beta <- x0
	T <- length(beta)
	j_seq <- seq(1, T)
	na <- 0

	for (j in j_seq){
		res <- rbeta_j(j, x0[j], sigma, eta, v, w)
		beta_j <- res$beta_j
		na_j <- res$ac
		na <- na + na_j
		beta[j] <- beta_j
	}
	ac <- na/T # taxa de aceitação média
	return(list(beta=beta, ac=ac))
}

# Read the data
x_df <- read.csv("x.csv", header = TRUE)
x <- as.vector(x_df$x)
T <- length(x)

# Insert the outliers
tau <- c(25, 75, 125, 175, 176)
eta_tau <- 10
y <- x
y[tau] <- y[tau] + eta_tau

# Simulation parameters
N      <- 1000 # number of steps
burnin <- 200

# Parameters intialization
alpha      <- 0.1
mu         <- 0.1
xi         <- 0.01
delta      <- rep(0, T)
prob_delta <- rep(0, T)
eta        <- rep(0, T)
p          <- rep(0.01, T)
beta       <- rep(5, T)

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
sigma_alpha  <- 0.1
sigma_mu     <- 0.5
sigma_xi     <- 0.02
sigma_beta   <- 3

# Accepted steps counters
na_alpha <- 0
na_mu <- 0
na_xi <- 0
ac_beta <- 0

# Vectors and matrices with traces
Alpha	   <- numeric(N)
Mu		   <- numeric(N)
Xi		   <- numeric(N)
Delta	   <- matrix(nrow=N, ncol=T)
Prob_delta <- matrix(nrow=N, ncol=T)
Eta 	   <- matrix(nrow=N, ncol=T)
P		   <- matrix(nrow=N, ncol=T)
Beta	   <- matrix(nrow=N, ncol=T)

Alpha[1] <- alpha
Mu[1] <- mu
Xi[1] <- xi
Delta[1, ] <- delta
Prob_delta[1, ] <- prob_delta
Eta[1, ] <- eta
P[1, ] <- p
Beta[1, ] <- beta

#####
start_time = proc.time() # para cálculo do tempo de execução
for (n in 2:N){

	if (n %% 100 == 0){
		printf("Step %d/%d: acceptance ratio: alpha: %.2f, mu: %.2f, xi: %.2f, beta: %.2f",
			   n, N, na_alpha/n, na_mu/n, na_xi/n, ac_beta)
	}

	# Draw alpha (Metropolis step)
	res <- ralpha(alpha, sigma_alpha, mu, xi, delta, eta, y, a, b)
	alpha <- res$alpha
	na_alpha <- na_alpha + res$ac
	Alpha[n] <- alpha

	# Draw mu (Metropolis step)
	res <- rmu(mu, sigma_mu, alpha, xi, delta, eta, y, c, d)
	mu <- res$mu
	na_mu <- na_mu +  res$ac
	Mu[n] <- mu

	# Draw xi (Metropolis step)
	res <- rxi(xi, sigma_xi, alpha, mu, delta, eta, y, g, h)
	xi <- res$xi
	na_xi <- na_xi +  res$ac
	Xi[n] <- xi

	# Draw delta (Bernoulli)
	res <- rdelta(alpha, mu, xi, delta, eta, p, y)
	delta <- res$delta
	prob_delta <-res$prob_delta
	Delta[n, ] <- delta
	Prob_delta[n, ] <- res$prob_delta

	# Draw eta (Discrete simulation)
	eta <- reta(alpha, mu, xi, delta, eta, p, beta, y)
	Eta[n, ] <- eta

	# Draw p (Beta)
	p <- rp(delta, l, m)
	P[n, ] <- p

	# Draw beta (Metropolis step)
	res <- rbetavect(beta, sigma_beta, eta, v, w)
	beta <- res$beta
	ac_beta <- (n-1)*(ac_beta)/n + res$ac/n
	Beta[n, ] <- beta
}
#####
# Execution time
end_time <- proc.time()
elapsed_time <- (end_time - start_time)[[3]]
printf("Execution time: %.0f s\n", elapsed_time)

# Acceptance ratios
# ac_p  <- na_alpha/N
# ac_mu <- na_mu/N
# ac_xi <- na_xi/N
# printf("Acceptance rate - alpha: %.2f", ac_p)
# printf("Acceptance rate - mu: %.2f", ac_mu)
# printf("Acceptance rate - xi: %.2f", ac_xi)
# printf("Acceptance rate - beta: %.2f", ac_beta)

# Values after burn in
Alpha_final <- Alpha[(burnin+1):N]
Mu_final <- Mu[(burnin+1):N]
Xi_final <- Xi[(burnin+1):N]
Delta_final <- Delta[(burnin+1):N, ]
Prob_delta_final <- Prob_delta[(burnin+1):N, ]
Eta_final <- Eta[(burnin+1):N, ]

# Mean values
alpha_mean <- mean(Alpha_final)
mu_mean <- mean(Mu_final)
xi_mean <- mean(Xi_final)
prob_delta_mean <- colMeans(Prob_delta_final)

# Estimatives using median
delta_est <- numeric(T)
eta_est <- numeric(T)

median_round_up <- function(x){
	q <- quantile(x, 0.5)
	med <- min(x[x >= q])
	return(med)
}

for (j in 1:T) {
    delta_est[j] <- median_round_up(Delta_final[, j])
    eta_est[j] <- median_round_up(Eta_final[, j])
}

# Estimated outliers using threshold probability
threshold <- 0.8
tau_est <- which(prob_delta_mean > threshold)

printf("alpha mean: %.2f", alpha_mean)
printf("mu mean: %.2f", mu_mean)
printf("xi mean: %.2f", xi_mean)
cat("tau_est: ", tau_est)
cat("\ntau: ", tau)
cat("\neta_est|tau_est: ", eta_est[tau_est])

#####
# Trace plots
par(mfrow = c(3, 1), mar=c(4,3,2,2))
plot(Alpha, type="l", main="alpha", xlab="", ylab="", ylim=c(0,1))
abline(v=burnin, col="gray", lty=2)
plot(Mu, type="l", main="mu", xlab="", ylab="")
abline(v=burnin, col="gray", lty=2)
plot(Xi, type="l", main="xi", xlab="step", ylab="", ylim=c(0,1))
abline(v=burnin, col="gray", lty=2)

par(mfrow = c(4, 1), mar=c(4,3,2,2))
plot(Prob_delta[, 25], type="l", main="delta_t", xlab="", ylab="")
lines(Prob_delta[, 75], type="l", main="delta_t", xlab="", ylab="", col=2)
lines(Prob_delta[, 125], type="l", main="delta_t", xlab="", ylab="", col=3)
abline(v=burnin, col="gray", lty=2)
plot(P[, 25], type="l", main="p_t", xlab="", ylab="")
abline(v=burnin, col="gray", lty=2)
plot(Eta[, 25], type="s", main="eta_t", xlab="", ylab="")
abline(v=burnin, col="gray", lty=2)
plot(Beta[, 25], type="l", main="beta_t", xlab="", ylab="")
abline(v=burnin, col="gray", lty=2)


#####
#Outlier detection
par(mfrow = c(3, 1), mar=c(4,3,2,2))
plot(y, type="s", main="Time series and outliers", xlab="", ylab="", font.main=1)
points(tau, y[tau], col="blue", pch=0)
points(tau_est, y[tau_est], col="red", pch=4)
legend(x="topright", pch=c(0, 4), col=c("blue", "red"), legend=c("label", "detected"))
plot(prob_delta_mean, type="h", main=expression(delta[t]~posterior~probability~(mean)), xlab="", ylab="")
legend(x="topright", col="red", lty=2, legend="threshold")
abline(h=threshold, col="red", lty=2)
y_ <- rep(0, T)
y_[tau_est] <- eta_est[tau_est]
plot(1:T, y_, type="h", main=expression(eta[t] * " | " * (delta[t]==1) ~ (median)), xlab="t", ylab="")

#####
# Posterior densities
par(mfrow=c(1,1))
plot(density(Alpha_final), main="alpha")
plot(density(Mu_final), main="mu")
plot(density(Xi_final), main="xi")

