#!/usr/bin/env Rscript

################ Read Command line inputs ###########
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("At least 2 arguments must be supplied (input file).n", call.=FALSE)
}
############### Inputs #########################
N <- as.integer(args[1])
Ns <- as.integer(args[2])
###############################################
Lp <- 0.0210870
Ln <- 0.0169807
kT <- 2.495   # kT at 300K

Veff <- read.table('./Veff.txt')$V1
Mrat <- read.table('./mratio.txt')$V1
TI <- read.table('./TI.txt')$V1/kT

dG <- -log(Mrat*Veff/( ((N+1)^2) * (Lp*Ln)^3)) + TI
write.table(x = dG,row.names = FALSE,col.names = FALSE,file = 'dG.txt')
P <- exp(-dG)/sum(exp(-na.omit(dG)))
P[is.na(P)] <- 0
S <- sample(x = c(0: (length(P)-1) ),prob = P,replace = TRUE,size = Ns)
write.table(x = S,row.names = FALSE,col.names = FALSE,file = 'boltzmann.txt')

