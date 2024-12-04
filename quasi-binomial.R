# Quasi-binomial Distribution

set.seed(42)

# Quasi-binomial thinning pmf
dqb <- function(x, p, theta, n){
    q <- 1-p
    num <- p*q*choose(n,x) * (p + x*theta)^(x-1) * (q + (n-x)*theta)^(n-x-1)
    den <- (1 + n*theta)^(n-1)
    y <- num/den
    return(y)
}

# Quasi-binomial thinning cdf
pqb <- function(x, p, theta, n){
    x_ <- 0:x
    y <- sum(sapply(x_, dqb, p=p, theta=theta, n=n))
    return(y)
}

# Simula valores da quase-binomial
rqb <- function(size, p, theta, n){
    x <- 0:n
    w <- sapply(x, dqb, p=p, theta=theta, n=n)
    y <- sample(x, size=size, replace=TRUE, prob=w)
    return(y)
}
    
n <- 100
p <- 0.5
theta <- 0.2
lambda <- 0.5

x <- 0:n

# Funções PMF e CDF
y_pmf <- sapply(x, dqb, p=p, theta=theta/lambda, n=n)
y_cdf <- sapply(x, pqb, p=p, theta=theta/lambda, n=n)
plot(x, y_pmf)
plot(x, y_cdf)

#####
# Amostra valores
y <- rqb(100, p, theta/lambda, n)
hist(y)
