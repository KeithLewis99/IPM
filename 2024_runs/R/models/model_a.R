model <- '
 model {
#PRIORS
###### Need to check that these are reasonable
## Prior for sd of process - N2-N3[t] uninformative
# a is process variance for age group, tp for time period pre/post collapse

   sigma.proc ~ dunif(0.01, 20)
   sigma2.proc <- pow(sigma.proc, 2)
   tau.proc <- pow(sigma.proc, -2)

## Prior for sd of observation - I2-I4[t] - uninformative
   sigma.obs ~ dunif(0.01, 20)
   sigma2.obs <- pow(sigma.obs, 2)
   tau.obs <- pow(sigma.obs, -2)

## Prior for sd of maturity - mat2-mat3[t] - uninformative
   sigma.mat ~ dunif(0.01, 1)
   sigma2.mat <- pow(sigma.mat, 2)
   tau.m <- pow(sigma.mat, -2)
   
# Prior for sd of rec - N[t, 0] - uninformative
   sigma.rec ~ dunif(0.01, 20)
   sigma2.rec <- pow(sigma.rec, 2)
   tau.rec <- pow(sigma.rec, -2)
   
## Prior for catchability
   # sigma.q ~ dunif(.01, 1)
   # sigma2.q <- pow(sigma.q, 2)
   # tau.q <- pow(sigma.q, -2)

### Priors for Initial values for N2-N4[t] informative - based on actual values
   rec[1] ~ dnorm(15, 1/9)
   N[1,1] ~ dnorm(12.8, 1/9)
   N[1,2] ~ dnorm(11.3, 1/9)
   N[1,3] ~ dnorm(8.2, 1/9)

## Priors for maturity at age
   # theta0 ~ dnorm(2., .5)
   # theta1 ~ dnorm(5., 2.)
   for (a in 1:Ni){
      for(t in 1:(n.occasions)){
         # m[t, a] <- logit(theta0 + theta1*a)
         m[t,a] ~ dbeta(100*matMp[t,a], 100*(1-matMp[t,a]))
      }
   }
   
### From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
### priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   alpha ~ dnorm(0, 100^-2)      # int
   beta ~ dnorm(0, 100^-2)       # larval abund
   for(a in 1:Ni){
      # alpha[a] ~ dnorm(0, 100^-2)      # int
      # beta[a] ~ dnorm(0, 100^-2)       # larval abund
      gamma[a] ~ dunif(0.01, 1)         # limit to 0 to 1
      delta[a] ~ dunif(1.29, 2.56)     # to make delta > max tice after scaling
      epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
   }
   
## Catchability
   # for(a in 1:Ni){
   #    logq[a] ~ dnorm(0, tau.q)
   #    q[a] <- exp(logq[a])
   # }
   
## Mortality
   sigma.Z ~ dunif(0.01, 1)
   sigma2.Z <- pow(sigma.Z, 2)
   tau.Z <- pow(sigma.Z, -2)
   for(a in 1:Ni){
      Z0[a] ~ dlnorm(0, tau.Z)
      # Z0[a] <- exp(logZ0[a])
   }
   
## 2010 offset for age 2s
   # N2offset ~ dlnorm(0, 100^-2)
   # N2offset <- exp(logN2offset)

# LIKELIHOODS
## State process
   
# Survivability/Mortality

   for(a in 1:Ni){
      for(t in 1:10){
         # s[t,a] <- exp(-Z[t, a])
         Z[t,a] <- gamma[a]*4*TI[t]/delta[a]*(1-TI[t]/delta[a]) +
                   Z0[a]
      }
      for(t in 11:(n.occasions-3)){
         # s[t,a] <- exp(-Z[t, a])
         Z[t,a] <- gamma[a]*4*TI[t]/delta[a]*(1-TI[t]/delta[a]) +
                   epsilon[a]*CO[t] + #Added condition back in
                   Z0[a]
      }
      for(t in (n.occasions-2):n.occasions){
         # s[t,a] <- exp(-Z[t, a])
         Z[t,a] <- gamma[a]*4*TI[t]/delta[a]*(1-TI[t]/delta[a]) +
                   Z0[a]
      }
   }
   
# Expected Abundance

   # for(t in 2:25){
   #    N[t,1] <- rec[t-1] - Z[t,1] #* s[t-1, 1]
   #    for(a in 2:Ni){
   #       N[t, a] <- N[t-1, a-1] - Z[t,a] #* s[t-1, a]
   #    }
   # }
   # N[26, 1] <- rec[25] - Z[25,1] #* s[25, 1] + N2offset #checking the 2010 problem
   # for(a in 2:Ni){
   #    N[26, a] <- N[25, a-1] - Z[25,a] #* s[25, a]
   # }
   # for(t in 27:n.occasions){
   #    N[t,1] <- rec[t-1] - Z[t,1] #* s[t-1, 1]
   #    for(a in 2:Ni){
   #       N[t, a] <- N[t-1, a-1] - Z[t,a] #* s[t-1, a]
   #    }
   # }
   
   for(t in 2:17){
      N[t,1] <- rec[t-1] - Z[t,1] # these rec are estimated directly
      for(a in 2:Ni){
         N[t, a] <- Np[t-1, a-1] - Z[t,a]
      }
   }
   for(t in 18:n.occasions){
      N[t,1] <- recp[t-1] - Z[t,1] # these rec are functional
      for(a in 2:Ni){
         N[t, a] <- Np[t-1, a-1] - Z[t,a]
      }
   }
   
   for(t in 2:17){ #prior for recruitment before LD
      rec[t] ~ dlnorm(log(rec[t-1]), tau.rec)
   }
   for(t in 18:(n.occasions-1)){ #abundance based on LD once available
      rec[t] <- alpha + beta*LD[t-1] # age 2 based on Larval Density for these years
   }
   
   # for(t in 2:17){ #prior for recruitment before LD
   #    # N[t,1] ~ dnorm(LDmu, LDsd)
   # }
   # LDmu <- mean(LD[17:n.occasions])
   # LDsd <- sd(LD[17:n.occasions])


## Process Error

   for(t in 1:(n.occasions-1)){
      recp[t] ~ dnorm(rec[t], tau.proc)
      for(a in 1:Ni){
         Np[t,a] ~ dnorm(N[t,a], tau.proc)
      }
   }
   for(a in 1:Ni){
      Np[n.occasions,a] ~ dnorm(N[n.occasions,a], tau.proc)
   }
   # for(t in 18:n.occasions-1){
   #    recp[t] ~ dnorm(rec[t], tau.proc)
   # }


# ######################################

# raw resids
   for (a in 1:Ni){
      for (t in 1:n.occasions){
         eps[t, a] <- matI[t,a] - N[t,a]
      }
   }

   # one step ahead resids - could set a p
   for (a in 1:Ni){
     posa[1, a] ~ dnorm(0, 1/10)
     # posa[2, a] ~ dnorm(0, 1/10)
     osa[1, a] ~ dnorm(0, 1/10)
     # osa[2, a] ~ dnorm(0, 1/10)
     osa_mean[1, a] ~ dnorm(0, 1/10)
     # osa_mean[2, a] ~ dnorm(0, 1/10)
     # osa_sd[1, a] ~ dnorm(0, 1/10)
     # osa_sd[2, a] ~ dnorm(0, 1/10)
     osa_sd[1, a] ~ dunif(0.1, 10)
     # osa_sd[2, a] ~ dunif(0.1, 10)
     # pe[1, a] ~ dnorm(0, 1/10)
   }

for (a in 1:Ni){
   for (t in 3:n.occasions){
        osa_mean[t-1,a] <- mean(N[1:(t-1), a])
        osa[t-1,a] <- N[t,a] - osa_mean[t-1,a]
        osa_sd[t-1,a] <- sd(osa[1:(t-1), a]) # this is the SD for t-1
        posa[t-1,a] <- osa[t-2,a]/osa_sd[t-2,a]
   }
}

## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


   for (a in 1:Ni){  
      for (t in 1:n.occasions) {
         matI[t,a] ~ dnorm(Np[t,a], tau.obs)  # sampled observation
         nll[t,a] <- (-log(dnorm(matI[t,a], Np[t,a], tau.obs)))
      }
   }

   for (t in 1:n.occasions) {
      I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
   }





# Assessing the fit of the state-space model
   ## 1. Compute fit statistics for observed data.
   ### 1.1 Discrepancy meansure: mean absolute error

   for (t in 1:n.occasions) {
      I.exp[t] <- log(exp(N[t,1]) + exp(N[t,2]) + exp(N[t,3]))
      Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
   }
   Dmape.obs <- sum(Dssm.obs)

   # ## 1.2 Test statistic: number of turns or switches - jaggedness
   for (t in 1:(n.occasions-2)){
      Tt1.obs[t] <- step(I[t+2] - I[t+1])
      Tt2.obs[t] <- step(I[t+1] - I[t])
      # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
      # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
      # Tt22.obs[t] <- step(I2[t+1] - I2[t])
      # Tt23.obs[t] <- step(I3[t+1] - I3[t])
      # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
      # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
      Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
   }
   Tturn.obs <- sum(Tt3.obs)


   ## 2.1 Simulated data
for (t in 1:n.occasions){
      I2.rep[t] ~ dnorm(N[t,1], tau.obs)
      I3.rep[t] ~ dnorm(N[t,2], tau.obs)
      I4.rep[t] ~ dnorm(N[t,3], tau.obs)
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}


   Dmape.rep <- sum(Dssm.rep)


   ##Test statistic: number of turns or switches - jaggedness
   for (t in 1:(n.occasions-2)){
      Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
      Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
      Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
   }
   Tturn.rep <- sum(Tt3.rep)

}'
