library(coda)

rm(list=ls()) # limpa o workspace
set.seed(10)

# Lê o arquivo CSV
tab <- read.csv("D:/OneDrive/Documentos/Projetos/Python/sim_est/prova2/baiod/results/Y1_5300_1/Alpha.txt", header = FALSE)
x1 <- mcmc(tab$V1)

tab <- read.csv("D:/OneDrive/Documentos/Projetos/Python/sim_est/prova2/baiod/results/Y1_5300_2/Alpha.txt", header = FALSE)
x2 <- mcmc(tab$V1)

tab <- read.csv("D:/OneDrive/Documentos/Projetos/Python/sim_est/prova2/baiod/results/Y1_5300_3/Alpha.txt", header = FALSE)
x3 <- mcmc(tab$V1)

lista = list(x1, x2, x3)
listaMcmc = mcmc.list(lista)
res1 = gelman.diag(listaMcmc)$psrf
res = res1[1]