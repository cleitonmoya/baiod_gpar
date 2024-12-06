# Modelo INAR(1) Poisson
# Inferência Baysiana
# Gibbs Sampling com passos de Metropolis Hastings


rm(list=ls()) # limpa o workspace
set.seed(42)

# Lê o arquivo CSV
setwd("C:/Users/cleit/OneDrive/Documentos/Projetos/R/sim_est/prova2")
tab <- read.csv("Y2.txt", header = FALSE)
Y <- tab$V1

# Função auxiliar para impressão de strings com variáveis numéricas
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}


#####
# Condicionais completas 

# Função auxiliar para cálculo das condicionais completas
funT <- function(i, t, alpha, lamb, Y){
    p1 <- lamb^(Y[t]-i)/factorial(Y[t]-i)
    p2 <-choose(Y[t-1],i)
    p3 <- alpha^i*(1-alpha)^(Y[t-1]-i)
    p <- p1*p2*p3
    return(p)
}

# Condicional completa para alpha (log-pdf não normalizada)
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

# Condicional completa de lambda (log-pdf não normalizada)
logpost_lamb <- function(lamb, alpha, Y, c, d){
    n <- length(Y)
    logp1 <- -(d+n)*lamb + (c-1)*log(lamb)
    
    s1 = 0
    for (t in 2:n){
        Mt = min(Y[t],Y[t-1])
        #if (Mt > 0){
            i_ <- 0:Mt
            y_ <- sapply(i_, function(i) funT(i, t, alpha, lamb, Y))
            s2 <- log(sum(y_))
        #} else {
        #    s2 <- 0
        #}
        s1 = s1 + s2
    }
    p <- s1 + logp1
    return(p)
}


#####
# Funções auxiliares - Metropolis Hastings

# Simula um valor da distribuição proposta
simula_xprop <- function(mu, sigma){
    x <- rnorm(n=1, mean=mu, sd=sigma)
    return(x)
}

# Calcula a probabiliadde de aceitação para alpha
log_prob_aceit_alpha <- function(x, x_prop, sigma, lamb, Y, a, b){
    p1 <- logpost_alpha(x_prop, lamb, Y, a, b)
    p2 <- logpost_alpha(x, lamb, Y, a, b)
    if (p1 < p2){
        prob_aceit <- min(1, exp(p1-p2))
    } else {
    prob_aceit = 1
    }
    return(prob_aceit)
}

# Calcula a probabiliadde de aceitação para lambda
log_prob_aceit_lamb<- function(x, x_prop, sigma, alpha, Y, a, b){
    p1 <- logpost_lamb(x_prop, alpha, Y, a, b)
    p2 <- logpost_lamb(x, alpha, Y, a, b)
    if (p1 < p2){
        prob_aceit <- min(1, exp(p1-p2))
    } else {
        prob_aceit = 1
    }
    return(prob_aceit)
}

# Gera uma amostra de alpha pelo método de Metropolis-Hastings
amostra_alpha <- function(sigma, x0, lamb, Y, a, b){
    
    aceito <- 0
    x <- x0
    
    # simula valor proposto
    x_prop <- simula_xprop(x, sigma)
    if (x_prop > 0 & x_prop < 1){
    
        # probabilidade de aceitação
        prob_aceit <- log_prob_aceit_alpha(x, x_prop, sigma, lamb, Y, a, b)
        
        # critério de aceição
        u = runif(1)
        if  (u < prob_aceit){
            x <- x_prop
            aceito <- 1
        }
    }
    return(list(x=x,aceito=aceito))
}

# Gera uma amostra de lambda pelo método de Metropolis-Hastings
amostra_lamb<- function(sigma, x0, alpha, Y, a, b){
    
    aceito <- 0
    x <- x0
    
    # simula valor proposto
    x_prop <- simula_xprop(x, sigma)
    if (x_prop > 0){
        
        # probabilidade de aceitação
        prob_aceit <- log_prob_aceit_lamb(x, x_prop, sigma, alpha, Y, a, b)
        
        # critério de aceição
        u = runif(1)
        if  (u < prob_aceit){
            x <- x_prop
            aceito <- 1
        }
    }
    return(list(x=x,aceito=aceito))
}


#####
# Algoritmo Gibbs Sampling com passos de Metropolis
N <- 1000 # número de passos
burnin <- 200

# hiperparâmetros das distribuições a priori
a <- 0.001
b <- 0.001
c <- 0.001
d <- 0.001

# parâmetros do passeio aleatório
sigma_alpha <- 0.02
sigma_lamb <- 0.2 #0.5

# inicialização
alpha <- 0.5
lamb <- 0.5

# número de aceitações de alpha
na_alpha <- 0 
na_lamb <- 0 

# inicialização dos vetores
Alpha <- numeric(N)
Lamb <- numeric(N)
Alpha[1] = alpha
Lamb[1] = lamb

start_time = proc.time() # para cálculo do tempo de execução
for (n in 2:N){
    
    if (n %% 200 == 0){
        printf("Passo %d/%d: aceit_alpha=%.2f, aceit_lamb=%.2f", n, N, na_alpha/n, na_lamb/n)
    }

    # gera uma amostra de alpha
    x <- amostra_alpha(sigma_alpha, alpha, lamb, Y, a, b)
    alpha <- x$x 
    aceit_alpha <- x$aceito
    na_alpha <- na_alpha + aceit_alpha
    Alpha[n] <- alpha
    
    # gera uma amostra de lambda
    x <- amostra_lamb(sigma_lamb, lamb, alpha, Y, c, d)
    lamb <- x$x
    aceit_lamb <- x$aceito
    na_lamb <- na_lamb + aceit_lamb
    Lamb[n] <- lamb
}

# Tempo de execução
end_time = proc.time()
elapsed_time = (end_time - start_time)[[3]]
printf("Tempo de execução: %.0f s\n", elapsed_time)

# Taxas de aceitação
taxa_aceit_alpha = na_alpha/N
taxa_aceit_lamb = na_lamb/N
printf("Taxa de aceitação de alpha: %.2f", taxa_aceit_alpha)
printf("Taxa de aceitação de lambda: %.2f", taxa_aceit_lamb)

# Valores médios
alpha_mean <- mean(Alpha[burnin+1:T])
lambda_mean <- mean(Lamb[burnin+1:T])
printf("Valor médio de alpha: %.2f", alpha_mean)
printf("Valor médio de lambda: %.2f", lambda_mean)

#####
# Gráficos
par(mfrow = c(2, 1), mar=c(4,3,2,2))
plot(Alpha, type="l", main="Alpha", xlab="", ylab="")
plot(Lamb, type="l", main="Lambda", xlab="Iteração", ylab="")