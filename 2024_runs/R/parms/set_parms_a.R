library(rjags)
library(R2jags)
# library(ggplot2)
library(lattice)
library(tidyverse)

library(here)

# Source files
source("IPM_dat.R")
makeData()
saveRDS(jags.data.m, '2024_runs/input_data.rds')
source(here('zuur/MCMCSupportHighstatV2.R'))
source(here('zuur/HighstatLibV7.R'))


# JAGS settings ----
parms = c("tau.proc",  "tau.obs", "sigma2.proc",
          "N", 'Np',
          'rec', 'recp', 'sigma.rec',
          "eps", "osa", "osa_mean", "osa_sd",
          "posa", "nll",
          "m",
          'tau.co',
          'sigma.Z', 'Z0',
          'Minf', 'sigma.laa',
          "beta", "alpha", 
          "gamma", "delta",
          "epsilon",
          "I2.rep", "I3.rep", "I4.rep", "I.rep")


# MCMC settings
ni <- 20000; nt <- 6; nb <- 5000; nc <- 3

jags.data.m$Ni <- 3 # ages - N is abundance and i is the index for age
# jags.data.m$M <- 3 # maturity in matrix matM - this may not be needed

jags.data.m$N2end <- 18
jags.data.m$N2start <- 19
jags.data.m$N3end <- 11
jags.data.m$N3start <- 12

source('2024_runs/R/models/model_a.R')
tC = model
tC.txt = "model"
# run model----
jagsout <- jags(jags.data.m, 
                parameters = parms,
                n.iter=ni, n.burnin = nb, 
                n.chains=nc, n.thin=nt, 
                model.file = textConnection(tC))
out <- jagsout$BUGSoutput
# warnings()

# the below are just various output that I comment on or off for convenience
out$sims.list$N2 <- out$sims.list$N[,,1]
out$sims.list$N3 <- out$sims.list$N[,,2]
out$sims.list$N4 <- out$sims.list$N[,,3]

## extract raw values from chains
raw <- ls_out(out)
ls_all <- ls_med(raw)
med <- ls_all$ls_med
cri <- ls_all$ls_cri
pri <- ls_all$ls_pri

source(here('2024_runs/R/get_stats.R'))
all.list <- list(raw=raw, 
                 med=med, cri=cri, pri=pri, 
                 stats=get_stats(jagsout) #This sometimes doesn't work ??
                 )
saveRDS(all.list, '2024_runs/data/models/model_a.rds')

