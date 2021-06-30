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
Gm <- as.numeric(args[3])
###############################################
Lp <- 0.0210870
Ln <- 0.0169807
kT <- 2.479   # kT at 298K
print(Gm)
#Gm <- 0.1   # Sample at temperature = T/Gm 
	    #(Eg: For Gm = 0.1 ->  T/Gm = 298/0.1 = 2980 K)

Veff <- read.table('./Veff.txt')$V1
Mrat <- read.table('./mratio.txt')$V1
TI <- read.table('./TI.txt')$V1/kT

Pr <- read.table('./Pr.txt')$V1

dG <- -log(Mrat*Veff/( ((N+1)^2) * (Lp*Ln)^3)) + TI -log(Pr)
write.table(x = dG,row.names = FALSE,col.names = FALSE,file = 'dG.txt')
dG <- dG -min(na.omit(dG))
P <- exp(-Gm*dG)/sum(exp(-na.omit(Gm*dG)))
P[is.na(P)] <- 0
S <- sample(x = c(0: (length(P)-1) ),prob = P,replace = TRUE,size = Ns)


Pr <- ( exp(-dG)/sum(exp(-na.omit(dG))) ) / (P)
Pr[is.na(Pr)] <- 0
Pr[is.nan(Pr)] <- 0

write.table(x = S,row.names = FALSE,col.names = FALSE,file = 'boltzmann.txt')
write.table(x = Pr,row.names = FALSE,col.names = FALSE,file = 'Pr.txt')
write.table(x = P,row.names = FALSE,col.names = FALSE,file = 'P.txt')

