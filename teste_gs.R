# Teste gs

# parâmetros na forma do artigo
alpha <- 0.8
mu    <- 2
qsi   <- 0.5

# fórmula do artigo
gs_artigo <- function(s, xt1, alpha, mu, qsi){
    alpha * choose(xt1, s) * (alpha + s*(qsi/mu))^(s-1) * 
        (1 - alpha -s*(qsi/mu))^(xt1-s)
}

gs_cleiton <- function(s, xt1, alpha, mu, qsi){
    mu*alpha*(1-alpha) * choose(xt1, s) * (mu*alpha + s*qsi)^(s-1) *
        ((1-alpha)*mu + (xt1-s)*qsi)^(xt1-s-1) /
        (mu+xt1*qsi)^(xt1-1)
}

gs_cleiton2 <- function(s, xt1, alpha, mu, qsi){
    alpha*(1-alpha) * choose(xt1, s) * (alpha + s*qsi/mu)^(s-1) *
        (1-alpha + (xt1-s)*qsi/mu)^(xt1-s-1) /
        (1+xt1*qsi/mu)^(xt1-1)
}

xt1 <- 10
s <- 5

g1 <- gs_artigo(s, xt1, alpha, mu, qsi)
g2 <- gs_cleiton(s, xt1, alpha, mu, qsi)
g3 <- gs_cleiton2(s, xt1, alpha, mu, qsi)
print(g1)
print(g2)
print(g3)