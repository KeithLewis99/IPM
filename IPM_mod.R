# Model ln scale: N2-N4 + forecast for each age----

cap.v7 = '
 model {
#PRIORS
      ###### Need to check that these are reasonable
 ## Prior for sd of process - N2[t] uninformative
   sigma.proc2 ~ dunif(0.01, 20)     
   sigma2.proc2 <- pow(sigma.proc2, 2) 
   tau.proc2 <- pow(sigma.proc2, -2) 

 ## Prior for sd of process - N3[t] uninformative
   sigma.proc3 ~ dunif(0.01, 20)     
   sigma2.proc3 <- pow(sigma.proc3, 2) 
   tau.proc3 <- pow(sigma.proc3, -2) 

 ## Prior for sd of process - N4[t] uninformative
    sigma.proc4 ~ dunif(0.01, 20)     
   sigma2.proc4 <- pow(sigma.proc4, 2) 
   tau.proc4 <- pow(sigma.proc4, -2) 

 ## Prior for sd of observation - I2-I4[t] - uninformative
   sigma.obs ~ dunif(0.01, 20)   
   sigma2.obs <- pow(sigma.obs, 2) 
   tau.obs <- pow(sigma.obs, -2) 
 
   ### Priors for Initial values for N2-N4[t] informative - based on actual values
   N2[1] ~ dnorm(8.5, 1/9)    
   N3[1] ~ dnorm(8.9, 1/9)    
    N4[1] ~ dnorm(6, 1/9)    
   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
   for (t in 2:4){
   N2[t] ~ dnorm(8.5, 1/9)
   }

# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   alpha2 ~ dnorm(0, 100^-2)      # int
   beta2 ~ dnorm(0, 100^-2)       # larval abund
   gamma2 ~ dunif(0, 100)         # tices-max rate of increase
   delta2 ~ dgamma(11.5, 5.7)     # tice-width
   #gamma2 ~ dnorm(0, 100^-2)     # condition # for CO
   #epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO

   # various model formulations - keep for eventual DIC 
   for (t in 5:n.occasions) {
      #mu2[t] <- alpha2 + beta2*LD[t-2]
      #mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*CO[t-1]
      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2)
      #mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2) + epsilon2*CO[t-1]

      N2[t] ~ dnorm(mu2[t], tau.proc2)
      
      # Divya suggested this formulation
      #N2[t] <- alpha + beta*LD[t-2] + gamma*TI[t]*(1-TI[t]/delta) + pe2[t]
      #pe2[t] ~ dnorm(0, tau.proc2)
   }     


   ### N3
   # Priors for killing N3
   #alpha3 ~ dnorm(0, 100^-2)          # int
   alpha3 ~ dunif(0, 1)                # int
   gamma3 ~ dunif(0, 100)              #tices-max rate of increase
   # gamma3 ~ dunif(0, 3.65)           # was uniform
   delta3 ~ dgamma(11.5, 5.7)          #tice-width
   epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
   mu3[1] ~ dunif(0,1)                 # need this for a mu3[1]
   
   # calculate a survival for N2 -> N3; mu3[2:24]
   for(t in 2:n.occasions){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
    }

   # calculate the N3
   for(t in 2:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))+mu3[t], tau.proc3)
   # Divyas formulation
      #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
      #pe3[t+1] ~ dnorm(0, tau.proc3)
      }

   ###N4
   # Priors for killing N4
   alpha4 ~ dunif(0, 1)             # int
   gamma4 ~ dunif(0, 100)           # tices-max rate of increase
   delta4 ~ dgamma(11.5, 5.7)       # tice-width
   epsilon4 ~ dnorm(0, 100^-2)      # condition # for CO
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:n.occasions){
      mu4[t] <- alpha4 + gamma4*TI[t]*(1-TI[t]/delta4) + epsilon4*CO[t-1]
          }
      
   # calculate the N4
   for(t in 2:n.occasions-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
            }

 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index 
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.

   for (t in 1:n.occasions) {          
      I2[t] ~ dnorm(N2[t], tau.obs)       # sampled observation
      I3[t] ~ dnorm(N3[t], tau.obs)
      I4[t] ~ dnorm(N4[t], tau.obs)
      I[t] ~ dnorm(log(exp(I2[t]) + exp(I3[t]) + exp(I4[t])), tau.obs)
   }

# Assessing the fit of the state-space model
 ## 1. Compute fit statistics for observed data.
  ### 1.1 Discrepancy meansure: mean absolute error
   for (t in 1:n.occasions) {
   I.exp[t] <- log(exp(N2[t]) + exp(N3[t]) + exp(N4[t]))
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
   #    y2.rep[t] ~ dnorm(N2[t], tau.obs) 
   #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
   I2.rep[t] ~ dnorm(N2[t], tau.obs)
   I3.rep[t] ~ dnorm(N3[t], tau.obs)
   I4.rep[t] ~ dnorm(N4[t], tau.obs)
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
