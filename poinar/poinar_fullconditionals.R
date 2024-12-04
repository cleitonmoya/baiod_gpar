# Modelo INAR(1) Poisson
# Inferência Baysiana
# Gibbs Sampling com passos de Metropolis Hastings

# Lê o arquivo CSV
setwd("C:/Users/cleit/OneDrive/Documentos/Projetos/R/sim_est/prova2")
tab <- read.csv("Y1.txt", header = FALSE)
Y <- tab$V1


# Função auxiliar para cálculo das condicionais completas
funT <- function(i, t, alpha, lamb, Y){
    p1 <- lamb^(Y[t]-i)/factorial(Y[t]-i)
    p2 <-choose(Y[t-1],i)
    p3 <- alpha^i*(1-alpha)^(Y[t-1]-i)
    p <- p1*p2*p3
    return(p)
}


# Condicional completa para alpha
logpost_alpha <- function(alpha, lamb, Y, a, b){
    n <- length(Y)
    logp1 <- (a-1)*log(alpha) + (b-1)*log(1-alpha)
    s1 <- 0
    for (t in 2:n){
        Mt <- min(Y[t],Y[t-1])
        if (Mt > 0){
            i_ <- 0:Mt
            y_ <- sapply(i_, function(i) funT(i, t, alpha, lamb, Y))
            s2 <- log(sum(y_))
        } else {
            s2 <- 0
        }
        s1 <- s1 + s2
    }
    p <- s1 + logp1
    return(p)
}


# Condicional completa de lambda
logpost_lamb <- function(lamb, alpha, Y, c, d){
    n <- length(Y)
    logp1 <- -(d+n)*lamb + (c-1)*log(lamb)
    
    s1 = 0
    for (t in 2:n){
        Mt = min(Y[t],Y[t-1])
        if (Mt > 0){
            i_ <- 0:Mt
            y_ <- sapply(i_, function(i) funT(i, t, alpha, lamb, Y))
            s2 <- log(sum(y_))
        } else {
            s2 <- 0
        }
        s1 = s1 + s2
    }
    p <- s1 + logp1
    return(p)
}

# hiperparâmetros das distribuições a priori
a <- 0.001
b <- 0.001
c <- 0.001
d <- 0.001

#####
# Densidade de alpha
lamb <- 1
alpha_ <- seq(0.10, 1.00, by=0.01)
logpdf_alpha <- sapply(alpha_, function(x) logpost_alpha(x,lamb,Y,a,b))
plot(alpha_, logpdf_alpha, type="l", xlab="alpha", ylab="Density")
print(alpha_[which.max(logpdf_alpha)])

#####
#Densidade de lambda
alpha <- 0.85
lamb_ <- seq(0.10, 5.00, by=0.01)
logpdf_lamb <- sapply(lamb_, function(x) logpost_lamb(x,alpha,Y,a,b))
plot(lamb_, logpdf_lamb, type="l", xlab="lambda", ylab="Density")
print(lamb_[which.max(logpdf_lamb)])
