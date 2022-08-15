# Model ln scale: N2-N4 + forecast for each age; 1999-2000----

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
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
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
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t,3], tau.proc[3])
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




# Model ln scale: N2-N4 + forecast for each age----
# this has common parameters for N3:N4 for the forecast model but not N2

cap.v8 = '
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
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   # Divyas formulation
      #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
      #pe3[t+1] ~ dnorm(0, tau.proc3)
      }

   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:n.occasions){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
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




# Model ln scale: N2-N4 + forecast for each age----
# This is a demographic model.  The forecast is applied only to N2 but only demographic values are applied to N3.

cap.v9 = '
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

   # calculate the N3
   for(t in 2:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t])), tau.proc3)
   # Divyas formulation
      #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
      #pe3[t+1] ~ dnorm(0, tau.proc3)
      }

   ###N4


   # calculate the N4
   for(t in 2:n.occasions-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95)), tau.proc4)
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

#######################################################
# Below are a series of attempts to use time series and smoothing techniques, i.e., an attempt to get rid of the pattern that has the early years below median and later years above median. Use a rdm walk on N2.  But this makes no sense biologically - at least it runs.
#######################################################

# Model ln scale: N2-N4 but with a rdm walk----
## no forecast on N2, only rdm walk.  Forecast on N3 but not N4.

cap.v10 = '
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

# LIKELIHOODS
 ## State process
   ### N2
   # random walk
   u ~ dnorm(0, 1/16) # this is based on a SD of 2 so Var = 4- then, I doubled it for no good reason other than to see what would happen

   # various model formulations - keep for eventual DIC 
   for (t in 2:n.occasions) {
      N2[t] ~ dnorm(N2[t-1] + u, tau.proc2)
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
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +  epsilon3*CO[t-1]
    }

   # calculate the N3
   for(t in 2:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   # Divyas formulation
      #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
      #pe3[t+1] ~ dnorm(0, tau.proc3)
      }

   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:n.occasions){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
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


 

# # Model ln scale: N2-N4 + forecast for each age ----
# # and with a rdm walk (previous one does not have forecast)

# cap.v10 = '
#  model {
# #PRIORS
#       ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2[t] uninformative
#    sigma.proc2 ~ dunif(0.01, 20)
#    sigma2.proc2 <- pow(sigma.proc2, 2)
#    tau.proc2 <- pow(sigma.proc2, -2)
# 
#  ## Prior for sd of process - N3[t] uninformative
#    sigma.proc3 ~ dunif(0.01, 20)
#    sigma2.proc3 <- pow(sigma.proc3, 2)
#    tau.proc3 <- pow(sigma.proc3, -2)
# 
#  ## Prior for sd of process - N4[t] uninformative
#     sigma.proc4 ~ dunif(0.01, 20)
#    sigma2.proc4 <- pow(sigma.proc4, 2)
#    tau.proc4 <- pow(sigma.proc4, -2)
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)
#    sigma2.obs <- pow(sigma.obs, 2)
#    tau.obs <- pow(sigma.obs, -2)
# 
#    ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N2[1] ~ dnorm(8.5, 1/9)
#    N3[1] ~ dnorm(8.9, 1/9)
#     N4[1] ~ dnorm(6, 1/9)
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
#    for (t in 2:4){
#    N2[t] ~ dnorm(8.5, 1/9)
#    }
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    alpha2 ~ dnorm(0, 100^-2)      # int
#    beta2 ~ dnorm(0, 100^-2)       # larval abund
#    gamma2 ~ dunif(0, 100)         # tices-max rate of increase
#    delta2 ~ dgamma(11.5, 5.7)     # tice-width
#    #epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO
# u ~ dnorm(0, 1/16) # this is based on a SD of 2 so Var = 4- then, I doubled it for no good reason other than to see what would happen
# 
#    # various model formulations - keep for eventual DIC
#    for (t in 5:n.occasions) {
#      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2) + u
#      N2[t] ~ dnorm(mu2[t], tau.proc2)
#    }
# 
# 
#    ### N3
#    # Priors for killing N3
#    #alpha3 ~ dnorm(0, 100^-2)          # int
#    alpha3 ~ dunif(0, 1)                # int
#    gamma3 ~ dunif(0, 100)              #tices-max rate of increase
#    # gamma3 ~ dunif(0, 3.65)           # was uniform
#    delta3 ~ dgamma(11.5, 5.7)          #tice-width
#    epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
#    mu3[1] ~ dunif(0,1)                 # need this for a mu3[1]
# 
#    # calculate a survival for N2 -> N3; mu3[2:24]
#    for(t in 2:n.occasions){
#       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
#     }
# 
#    # calculate the N3
#    for(t in 2:n.occasions-1){
#       N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
#    # Divyas formulation
#       #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
#       #pe3[t+1] ~ dnorm(0, tau.proc3)
#       }
# 
#    ###N4
#    mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
# 
#    # calculate a survival for N3 -> N4; mu4[2:24]
#    for(t in 2:n.occasions){
#       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
#           }
# 
#    # calculate the N4
#    for(t in 2:n.occasions-1){
#       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
#             }
# 
#  ## Observation
#    ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
#    #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
#    #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
# 
#    for (t in 1:n.occasions) {
#       I2[t] ~ dnorm(N2[t], tau.obs)       # sampled observation
#       I3[t] ~ dnorm(N3[t], tau.obs)
#       I4[t] ~ dnorm(N4[t], tau.obs)
#       I[t] ~ dnorm(log(exp(I2[t]) + exp(I3[t]) + exp(I4[t])), tau.obs)
#    }
# 
# # Assessing the fit of the state-space model
#  ## 1. Compute fit statistics for observed data.
#   ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#    I.exp[t] <- log(exp(N2[t]) + exp(N3[t]) + exp(N4[t]))
#    Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
# 
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
# 
# 
#  ## 2.1 Simulated data
#     for (t in 1:n.occasions){
#    #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
#    #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I2.rep[t] ~ dnorm(N2[t], tau.obs)
#    I3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I4.rep[t] ~ dnorm(N4[t], tau.obs)
#    I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
#    Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
#     }
#    Dmape.rep <- sum(Dssm.rep)
# 
# 
# ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# 
# }'
# 
# parms2 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs",
#             "N2",  "N3", "N4",
#             "mu2", "alpha2", "beta2",  "gamma2", "delta2",
#             "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
#             "mu4",
#             "Dssm.obs", "Dmape.obs",  "Tturn.obs",
#             "Dssm.rep", "Dmape.rep",  "Tturn.rep",
#             "I2.rep", "I3.rep", "I4.rep", "I.rep"
# )
# ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
# 
# # run model
# #source("IPM_mod.R")
# ssm26 <- jags(jags.data, parameters=parms2, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(cap.v10))
# ssm26




# # Model ln scale: N2-N4 + forecast for each age----
# # Autoregression 

# cap.v10 = '
#  model {
# #PRIORS
#       ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2[t] uninformative
#    sigma.proc2 ~ dunif(0.01, 20)     
#    sigma2.proc2 <- pow(sigma.proc2, 2) 
#    tau.proc2 <- pow(sigma.proc2, -2) 
# 
#  ## Prior for sd of process - N3[t] uninformative
#    sigma.proc3 ~ dunif(0.01, 20)     
#    sigma2.proc3 <- pow(sigma.proc3, 2) 
#    tau.proc3 <- pow(sigma.proc3, -2) 
# 
#  ## Prior for sd of process - N4[t] uninformative
#     sigma.proc4 ~ dunif(0.01, 20)     
#    sigma2.proc4 <- pow(sigma.proc4, 2) 
#    tau.proc4 <- pow(sigma.proc4, -2) 
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)   
#    sigma2.obs <- pow(sigma.obs, 2) 
#    tau.obs <- pow(sigma.obs, -2) 
#  
#    ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N2[1] ~ dnorm(8.5, 1/9)    
#    N3[1] ~ dnorm(8.9, 1/9)    
#     N4[1] ~ dnorm(6, 1/9)    
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
#    for (t in 2:4){
#    N2[t] ~ dnorm(8.5, 1/9)
#    }
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    alpha2 ~ dnorm(0, 100^-2)      # int
#    beta2 ~ dnorm(0, 100^-2)       # larval abund
#    gamma2 ~ dunif(0, 100)         # tices-max rate of increase
#    delta2 ~ dgamma(11.5, 5.7)     # tice-width
#    #epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO
# 
#
# for (i in 1:p) {
# zeta[i] ~ dnorm(0, 100^-2)
# }
# 
#    # various model formulations - keep for eventual DIC 
#    for (t in 5:n.occasions) {
#      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2) + inprod(zeta, I2[(t-p):(t-1)])
#      N2[t] ~ dnorm(mu2[t], tau.proc2)
#    }     
# 
# 
#    ### N3
#    # Priors for killing N3
#    #alpha3 ~ dnorm(0, 100^-2)          # int
#    alpha3 ~ dunif(0, 1)                # int
#    gamma3 ~ dunif(0, 100)              #tices-max rate of increase
#    # gamma3 ~ dunif(0, 3.65)           # was uniform
#    delta3 ~ dgamma(11.5, 5.7)          #tice-width
#    epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
#    mu3[1] ~ dunif(0,1)                 # need this for a mu3[1]
#    
#    # calculate a survival for N2 -> N3; mu3[2:24]
#    for(t in 2:n.occasions){
#       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
#     }
# 
#    # calculate the N3
#    for(t in 2:n.occasions-1){
#       N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
#    # Divyas formulation
#       #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
#       #pe3[t+1] ~ dnorm(0, tau.proc3)
#       }
# 
#    ###N4
#    mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
# 
#    # calculate a survival for N3 -> N4; mu4[2:24]
#    for(t in 2:n.occasions){
#       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
#           }
#       
#    # calculate the N4
#    for(t in 2:n.occasions-1){
#       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
#             }
# 
#  ## Observation
#    ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
#    #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index 
#    #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
# 
#    for (t in 1:n.occasions) {          
#       I2[t] ~ dnorm(N2[t], tau.obs)       # sampled observation
#       I3[t] ~ dnorm(N3[t], tau.obs)
#       I4[t] ~ dnorm(N4[t], tau.obs)
#       I[t] ~ dnorm(log(exp(I2[t]) + exp(I3[t]) + exp(I4[t])), tau.obs)
#    }
# 
# # Assessing the fit of the state-space model
#  ## 1. Compute fit statistics for observed data.
#   ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#    I.exp[t] <- log(exp(N2[t]) + exp(N3[t]) + exp(N4[t]))
#    Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
#    
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
# 
# 
#  ## 2.1 Simulated data
#     for (t in 1:n.occasions){
#    #    y2.rep[t] ~ dnorm(N2[t], tau.obs) 
#    #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I2.rep[t] ~ dnorm(N2[t], tau.obs)
#    I3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I4.rep[t] ~ dnorm(N4[t], tau.obs)
#    I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
#    Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
#     }
#    Dmape.rep <- sum(Dssm.rep)
# 
#      
# ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# 
# }'
# 
# jags.data$p <- 1
# 
# parms2 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
#             "N2",  "N3", "N4",
#             "mu2", "alpha2", "beta2",  "gamma2", "delta2",
#             "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
#             "mu4", 
#             "Dssm.obs", "Dmape.obs",  "Tturn.obs", 
#             "Dssm.rep", "Dmape.rep",  "Tturn.rep",
#             "I2.rep", "I3.rep", "I4.rep", "I.rep"
# ) 
# ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
# 
# # run model
# #source("IPM_mod.R")
# ssm26 <- jags(jags.data, parameters=parms2, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(cap.v10))
# ssm26
# 
# 
# plot(1999:2023, jags.data$I2)
# lines(1999:2022, apply(ssm26$BUGSoutput$sims.list$mu, 2, 'median'))





# # Model ln scale: N2-N4 + forecast for each age----
## As above but with model averaging (MA)

# cap.v11 = '
#  model {
# #PRIORS
#       ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2[t] uninformative
#    sigma.proc2 ~ dunif(0.01, 20)
#    sigma2.proc2 <- pow(sigma.proc2, 2)
#    tau.proc2 <- pow(sigma.proc2, -2)
# 
#  ## Prior for sd of process - N3[t] uninformative
#    sigma.proc3 ~ dunif(0.01, 20)
#    sigma2.proc3 <- pow(sigma.proc3, 2)
#    tau.proc3 <- pow(sigma.proc3, -2)
# 
#  ## Prior for sd of process - N4[t] uninformative
#     sigma.proc4 ~ dunif(0.01, 20)
#    sigma2.proc4 <- pow(sigma.proc4, 2)
#    tau.proc4 <- pow(sigma.proc4, -2)
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)
#    sigma2.obs <- pow(sigma.obs, 2)
#    tau.obs <- pow(sigma.obs, -2)
# 
#    ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N2[1] ~ dnorm(8.5, 1/9)
#    N3[1] ~ dnorm(8.9, 1/9)
#     N4[1] ~ dnorm(6, 1/9)
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
#    for (t in 2:4){
#    N2[t] ~ dnorm(8.5, 1/9)
#    }
# 
# # set up residuals
# for(t in 1:4){
# eps[t] <- alpha2 - N2[t] # this is probably wrong should be N2 - alpha2
# }
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    alpha2 ~ dnorm(0, 100^-2)      # int
#    beta2 ~ dnorm(0, 100^-2)       # larval abund
#    gamma2 ~ dunif(0, 100)         # tices-max rate of increase
#    delta2 ~ dgamma(11.5, 5.7)     # tice-width
#    #epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO
# 
# for (i in 1:q) {
# eta[i] ~ dnorm(0, 100^-2)
# }
# 
# 
# # various model formulations - keep for eventual DIC
#    for (t in 5:n.occasions) {
#      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2) + ma_mean[t]
#      ma_mean[t] <- inprod(eta, eps[(t-q):(t-1)])
#      N2[t] ~ dnorm(mu2[t], tau.proc2)
#      eps[t] <- N2[t] -  mu2[t]
#    }
# 
# 
#    ### N3
#    # Priors for killing N3
#    #alpha3 ~ dnorm(0, 100^-2)          # int
#    alpha3 ~ dunif(0, 1)                # int
#    gamma3 ~ dunif(0, 100)              #tices-max rate of increase
#    # gamma3 ~ dunif(0, 3.65)           # was uniform
#    delta3 ~ dgamma(11.5, 5.7)          #tice-width
#    epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
#    mu3[1] ~ dunif(0,1)                 # need this for a mu3[1]
# 
#    # calculate a survival for N2 -> N3; mu3[2:24]
#    for(t in 2:n.occasions){
#       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
#     }
# 
#    # calculate the N3
#    for(t in 2:n.occasions-1){
#       N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
#    # Divyas formulation
#       #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
#       #pe3[t+1] ~ dnorm(0, tau.proc3)
#       }
# 
#    ###N4
#    mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
# 
#    # calculate a survival for N3 -> N4; mu4[2:24]
#    for(t in 2:n.occasions){
#       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
#           }
# 
#    # calculate the N4
#    for(t in 2:n.occasions-1){
#       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
#             }
# 
#  ## Observation
#    ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
#    #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
#    #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
# 
#    for (t in 1:n.occasions) {
#       I2[t] ~ dnorm(N2[t], tau.obs)       # sampled observation
#       I3[t] ~ dnorm(N3[t], tau.obs)
#       I4[t] ~ dnorm(N4[t], tau.obs)
#       I[t] ~ dnorm(log(exp(I2[t]) + exp(I3[t]) + exp(I4[t])), tau.obs)
#    }
# 
# # Assessing the fit of the state-space model
#  ## 1. Compute fit statistics for observed data.
#   ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#    I.exp[t] <- log(exp(N2[t]) + exp(N3[t]) + exp(N4[t]))
#    Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
# 
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
# 
# 
#  ## 2.1 Simulated data
#     for (t in 1:n.occasions){
#    #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
#    #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I2.rep[t] ~ dnorm(N2[t], tau.obs)
#    I3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I4.rep[t] ~ dnorm(N4[t], tau.obs)
#    I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
#    Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
#     }
#    Dmape.rep <- sum(Dssm.rep)
# 
# 
# ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# }'
# #
# jags.data$q <- 1
# 
# parms2 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs",
#             "N2",  "N3", "N4",
#             "mu2", "alpha2", "beta2",  "gamma2", "delta2",
#             "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
#             "mu4",
#             "Dssm.obs", "Dmape.obs",  "Tturn.obs",
#             "Dssm.rep", "Dmape.rep",  "Tturn.rep",
#             "I2.rep", "I3.rep", "I4.rep", "I.rep"
# )
# ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
# 
# # run model
# #source("IPM_mod.R")
# ssm26 <- jags(jags.data, parameters=parms2, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(cap.v11))
# ssm26
# 
# 
# plot(1999:2023, jags.data$I2)
# lines(1999:2022, apply(ssm26$BUGSoutput$sims.list$mu, 2, 'median'))





# # Model ln scale: N2-N4 + forecast for each age----
## Autoregression- model averaging (MA)

# cap.v12 = '
#  model {
# #PRIORS
#       ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2[t] uninformative
#    sigma.proc2 ~ dunif(0.01, 20)
#    sigma2.proc2 <- pow(sigma.proc2, 2)
#    tau.proc2 <- pow(sigma.proc2, -2)
# 
#  ## Prior for sd of process - N3[t] uninformative
#    sigma.proc3 ~ dunif(0.01, 20)
#    sigma2.proc3 <- pow(sigma.proc3, 2)
#    tau.proc3 <- pow(sigma.proc3, -2)
# 
#  ## Prior for sd of process - N4[t] uninformative
#     sigma.proc4 ~ dunif(0.01, 20)
#    sigma2.proc4 <- pow(sigma.proc4, 2)
#    tau.proc4 <- pow(sigma.proc4, -2)
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)
#    sigma2.obs <- pow(sigma.obs, 2)
#    tau.obs <- pow(sigma.obs, -2)
# 
#    ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N2[1] ~ dnorm(8.5, 1/9)
#    N3[1] ~ dnorm(8.9, 1/9)
#     N4[1] ~ dnorm(6, 1/9)
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
#    for (t in 2:4){
#    N2[t] ~ dnorm(8.5, 1/9)
#    }
# 
# # set up residuals
# for(t in 1:4){
# eps[t] <- alpha2 - N2[t] # probably should be N2-alpha2
# }
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    alpha2 ~ dnorm(0, 100^-2)      # int
#    beta2 ~ dnorm(0, 100^-2)       # larval abund
#    gamma2 ~ dunif(0, 100)         # tices-max rate of increase
#    delta2 ~ dgamma(11.5, 5.7)     # tice-width
#    #epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO
# 
# 
# for (i in 1:p) {
#    eta[i] ~ dnorm(0, 100^-2)
# }
# 
# for (i in 1:q) {
#    zeta[i] ~ dnorm(0, 100^-2)
# }
# 
# # various model formulations - keep for eventual DIC
#    for (t in 5:n.occasions) {
#      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2) + inprod(eta, N2[(t-p):(t-1)]) + inprod(zeta, eps[(t-q):(t-1)])
#      N2[t] ~ dnorm(mu2[t], tau.proc2)
#      eps[t] <- N2[t] - alpha2 - mu2[t]
#    }
# 
# 
#    ### N3
#    # Priors for killing N3
#    #alpha3 ~ dnorm(0, 100^-2)          # int
#    alpha3 ~ dunif(0, 1)                # int
#    gamma3 ~ dunif(0, 100)              #tices-max rate of increase
#    # gamma3 ~ dunif(0, 3.65)           # was uniform
#    delta3 ~ dgamma(11.5, 5.7)          #tice-width
#    epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
#    mu3[1] ~ dunif(0,1)                 # need this for a mu3[1]
# 
#    # calculate a survival for N2 -> N3; mu3[2:24]
#    for(t in 2:n.occasions){
#       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
#     }
# 
#    # calculate the N3
#    for(t in 2:n.occasions-1){
#       N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
#    # Divyas formulation
#       #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
#       #pe3[t+1] ~ dnorm(0, tau.proc3)
#       }
# 
#    ###N4
#    mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
# 
#    # calculate a survival for N3 -> N4; mu4[2:24]
#    for(t in 2:n.occasions){
#       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
#           }
# 
#    # calculate the N4
#    for(t in 2:n.occasions-1){
#       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
#             }
# 
#  ## Observation
#    ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
#    #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
#    #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
# 
#    for (t in 1:n.occasions) {
#       I2[t] ~ dnorm(N2[t], tau.obs)       # sampled observation
#       I3[t] ~ dnorm(N3[t], tau.obs)
#       I4[t] ~ dnorm(N4[t], tau.obs)
#       I[t] ~ dnorm(log(exp(I2[t]) + exp(I3[t]) + exp(I4[t])), tau.obs)
#    }
# 
# # Assessing the fit of the state-space model
#  ## 1. Compute fit statistics for observed data.
#   ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#    I.exp[t] <- log(exp(N2[t]) + exp(N3[t]) + exp(N4[t]))
#    Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
# 
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
# 
# 
#  ## 2.1 Simulated data
#     for (t in 1:n.occasions){
#    #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
#    #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I2.rep[t] ~ dnorm(N2[t], tau.obs)
#    I3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I4.rep[t] ~ dnorm(N4[t], tau.obs)
#    I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
#    Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
#     }
#    Dmape.rep <- sum(Dssm.rep)
# 
# 
# ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# 
# }'
# 
# jags.data$p <- 1
# jags.data$q <- 1
# 
# parms2 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs",
#             "N2",  "N3", "N4",
#             "mu2", "alpha2", "beta2",  "gamma2", "delta2",
#             "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
#             "mu4",
#             "Dssm.obs", "Dmape.obs",  "Tturn.obs",
#             "Dssm.rep", "Dmape.rep",  "Tturn.rep",
#             "I2.rep", "I3.rep", "I4.rep", "I.rep"
# )
# ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
# 
# # run model
# #source("IPM_mod.R")
# ssm26 <- jags(jags.data, parameters=parms2, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(cap.v12))
# ssm26
# 
# 
# plot(1999:2023, jags.data$I2)
# lines(1999:2022, apply(ssm26$BUGSoutput$sims.list$mu, 2, 'median'))




 
# Model ln scale: N2-N4 + forecast for each age - NO ALPHA for age 2----
# as above but with NO ALPHA variable for age 2.  I removed this because I think that the alpha is forcing everything to deviate from a grand mean and this means that lower values (early 2000s) and higher values (2010s) are all below or above the line respectively.

# cap.v13 = '
#  model {
# #PRIORS
#       ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2[t] uninformative
#    sigma.proc2 ~ dunif(0.01, 20)
#    sigma2.proc2 <- pow(sigma.proc2, 2)
#    tau.proc2 <- pow(sigma.proc2, -2)
# 
#  ## Prior for sd of process - N3[t] uninformative
#    sigma.proc3 ~ dunif(0.01, 20)
#    sigma2.proc3 <- pow(sigma.proc3, 2)
#    tau.proc3 <- pow(sigma.proc3, -2)
# 
#  ## Prior for sd of process - N4[t] uninformative
#     sigma.proc4 ~ dunif(0.01, 20)
#    sigma2.proc4 <- pow(sigma.proc4, 2)
#    tau.proc4 <- pow(sigma.proc4, -2)
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)
#    sigma2.obs <- pow(sigma.obs, 2)
#    tau.obs <- pow(sigma.obs, -2)
# 
#    ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N2[1] ~ dnorm(8.5, 1/9)
#    N3[1] ~ dnorm(8.9, 1/9)
#     N4[1] ~ dnorm(6, 1/9)
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
#    for (t in 2:4){
#    N2[t] ~ dnorm(8.5, 1/9)
#    }
# 
# # set up residuals
# for(t in 1:4){
# eps[t] <- alpha2 - N2[t] #probably wrong N2 - alpha
# }
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    alpha2 ~ dnorm(0, 100^-2)      # int
#    beta2 ~ dnorm(0, 100^-2)       # larval abund
#    gamma2 ~ dunif(0, 100)         # tices-max rate of increase
#    delta2 ~ dgamma(11.5, 5.7)     # tice-width
#    #epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO
# 
# for (i in 1:p) {
#    eta[i] ~ dnorm(0, 100^-2)
# }
# for (i in 1:q) {
#    zeta[i] ~ dnorm(0, 100^-2)
# }
# 
# # various model formulations - keep for eventual DIC
#    for (t in 5:n.occasions) {
#      mu2[t] <- beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2) + ar_mean[t] + ma_mean[t]
#      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
#      ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
#      N2[t] ~ dnorm(mu2[t], tau.proc2)
#      eps[t] <- N2[t] - alpha2 - mu2[t]
#    }
# 
# 
#    ### N3
#    # Priors for killing N3
#    #alpha3 ~ dnorm(0, 100^-2)          # int
#    alpha3 ~ dunif(0, 1)                # int
#    gamma3 ~ dunif(0, 100)              #tices-max rate of increase
#    # gamma3 ~ dunif(0, 3.65)           # was uniform
#    delta3 ~ dgamma(11.5, 5.7)          #tice-width
#    epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
#    mu3[1] ~ dunif(0,1)                 # need this for a mu3[1]
# 
#    # calculate a survival for N2 -> N3; mu3[2:24]
#    for(t in 2:n.occasions){
#       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
#     }
# 
#    # calculate the N3
#    for(t in 2:n.occasions-1){
#       N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
#    # Divyas formulation
#       #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
#       #pe3[t+1] ~ dnorm(0, tau.proc3)
#    }
# 
#    ###N4
#    mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
# 
#    # calculate a survival for N3 -> N4; mu4[2:24]
#    for(t in 2:n.occasions){
#       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
#           }
# 
#    # calculate the N4
#    for(t in 2:n.occasions-1){
#       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
#             }
# 
#  ## Observation
#    ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
#    #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
#    #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
# 
#    for (t in 1:n.occasions) {
#       I2[t] ~ dnorm(N2[t], tau.obs)       # sampled observation
#       I3[t] ~ dnorm(N3[t], tau.obs)
#       I4[t] ~ dnorm(N4[t], tau.obs)
#       I[t] ~ dnorm(log(exp(I2[t]) + exp(I3[t]) + exp(I4[t])), tau.obs)
#    }
# 
# # Assessing the fit of the state-space model
#  ## 1. Compute fit statistics for observed data.
#   ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#    I.exp[t] <- log(exp(N2[t]) + exp(N3[t]) + exp(N4[t]))
#    Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
# 
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
# 
# 
#  ## 2.1 Simulated data
#     for (t in 1:n.occasions){
#    #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
#    #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I2.rep[t] ~ dnorm(N2[t], tau.obs)
#    I3.rep[t] ~ dnorm(N3[t], tau.obs)
#    I4.rep[t] ~ dnorm(N4[t], tau.obs)
#    I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
#    Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
#     }
#    Dmape.rep <- sum(Dssm.rep)
# 
# 
# ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# 
# }'
# 
# jags.data$p <- 3
# jags.data$q <- 3
# 
# parms2 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs",
#             "N2",  "N3", "N4",
#             "mu2", "alpha2", "beta2",  "gamma2", "delta2",
#             "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
#             "mu4",
#             "Dssm.obs", "Dmape.obs",  "Tturn.obs",
#             "Dssm.rep", "Dmape.rep",  "Tturn.rep",
#             "I2.rep", "I3.rep", "I4.rep", "I.rep"
# )
# 
# vars_vAR <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs")
# vars_Nyear <- c("N2[10]", "N3[10]", "N4[10]")
# vars_N2 <- c("mu2[10]","alpha2", "beta2",  "gamma2", "delta2")
# vars_N3 <- c("mu3[10]", "alpha3", "gamma3", "delta3", "epsilon3", "mu4[10]")
# vars_N4 <- c(NA)
# 
# ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
# 
# # run model
# #source("IPM_mod.R")
# ssm26 <- jags(jags.data, parameters=parms2, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(cap.v13))
# ssm26
# 
# 
# plot(1999:2023, jags.data$I2)
# lines(2003:2023, apply(ssm26$BUGSoutput$sims.list$mu2, 2, 'median'))





# Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998

cap.v20 = '
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
   # for (t in 2:4){
   # N2[t] ~ dnorm(8.5, 1/9)
   # }

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


for (t in 1:N2end) { #18
      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2)
   }     

for (t in 2:N2end) { #18
      N2[t] ~ dnorm(mu2[t], tau.proc2)
   }     

for (t in N2start:n.occasions) { #19
      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2)
      N2[t] ~ dnorm(mu2[t], tau.proc2)
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
   for(t in 2:N3end){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) 
    }

   for(t in N3start:n.occasions){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
    }

   # calculate the N3
   for(t in 2:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
      }
   # for(t in N3start:n.occasions-1){
   #    N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   #    }


   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:N3end){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
          }

   for(t in N3start:n.occasions){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
          }
      
   # calculate the N4
   for(t in 2:N3end-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
            }
   for(t in N3start:n.occasions-1){
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


## Extend model time + ARMA (age 2 only)----
### ln scale: N2-N4 + forecast for N2 only
#### This model has the ARMA smoothers but only for 1985-2002 and it is ONLY for the age 2 capelin. 

# cap.v22 = '
#  model {
# #PRIORS
# ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2[t] uninformative
#    sigma.proc2 ~ dunif(0.01, 20)
#    sigma2.proc2 <- pow(sigma.proc2, 2)
#    tau.proc2 <- pow(sigma.proc2, -2)
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)
#    sigma2.obs <- pow(sigma.obs, 2)
#    tau.obs <- pow(sigma.obs, -2)
# 
#    ### Priors for Initial values for N2-N4[t] informative - based on actual values
# 
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# 
# # # set up residuals
#  for(t in 1:max(p:q)){
#  eps[t] <- alpha2 - N2[t] # probably wrong N2[t] - alpha2
#  mu2[t] ~ dunif(0,10)
#  N2[t] ~ dnorm(8.5, 1/9)
#  }
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    alpha2 ~ dnorm(0, 100^-2)      # int
#    beta2 ~ dnorm(0, 100^-2)       # larval abund
#    gamma2 ~ dunif(0, 100)         # tices-max rate of increase
#    delta2 ~ dgamma(11.5, 5.7)     # tice-width
# 
#  for (i in 1:q) {
#     eta[i] ~ dnorm(0, 100^-2)
#  }
# 
#  for (i in 1:p) {
#     zeta[i] ~ dnorm(0, 100^-2)
#  }
# 
# for (t in (q+1):18) {
#      ar_mean[t] <- inprod(zeta, N2[(t-p):(t-1)])
#      ma_mean[t] <- inprod(eta, eps[(t-q):(t-1)])
# }
# 
# 
# for (t in (q+1):18) {
#      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + ar_mean[t] + ma_mean[t]
#      N2[t] ~ dnorm(mu2[t], tau.proc2)
#      eps[t] <- N2[t] - mu2[t]
# }
# 
#    for (t in 1:18) {
#       I2[t] ~ dnorm(N2[t], tau.obs)       # sampled observation
#     }
# }'
# 




#22: Extend model time ----
### ln scale: N2-N4 + forecast for each age
#### This model combines the 1985-1998 data but with ARMA smoothers for N2 ONLY.
#### Also, N2 are dealt with from 1985-2002:2003-2023 but N3 & N4 are dealt with from 1985-1995: 1996-2023.

cap.v22 = '
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
   N3[1] ~ dnorm(8.9, 1/9)
   N4[1] ~ dnorm(6, 1/9)

# # set up residuals, and initial values of N2 and mu2
 for(t in 1:max(p, q)){
   eps[t] <- alpha2 - N2[t] # probably wrong N2[t] - alpha2
   N2[t] ~ dnorm(8.5, 1/9)
   mu2[t] ~ dunif(0,10)           
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

  for (i in 1:p) {
    eta[i] ~ dnorm(0, 100^-2)
 }

   for (i in 1:q) {
    zeta[i] ~ dnorm(0, 100^-2)
   }
 
for (t in (max(p:q)+1):18) {
      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + ar_mean[t] + ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
   }

for (t in 19:n.occasions) {
      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + ar_mean[t] + ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
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
   for(t in 2:11){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
    }

   for(t in 12:n.occasions){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
    }

   # calculate the N3
   for(t in 2:11-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }
   for(t in 12:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }

   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:11){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
    }

   for(t in 12:n.occasions){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
    }

   # calculate the N4
   for(t in 2:11-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
   }

   for(t in 12:n.occasions-1){
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
   
## 1.2 Test statistic: number of turns or switches - jaggedness
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


# 23: Extend model time ----
### ln scale: N2-N4 + forecast for each age but including LD and CO
#### This model combines the 1985-1998 data but with ARMA smoothers for N2.

cap.v23 = '
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
   N3[1] ~ dnorm(8.9, 1/9)
   N4[1] ~ dnorm(6, 1/9)

# # set up residuals, and initial values of N2 and mu2
 for(t in 1:max(p, q)){
   eps[t] <- alpha2 - N2[t] # probably wrong N2[t] - alpha2
   N2[t] ~ dnorm(8.5, 1/9)
   mu2[t] ~ dunif(0,10)           
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
   epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO

  for (i in 1:p) {
    eta[i] ~ dnorm(0, 100^-2)
 }

   for (i in 1:q) {
    zeta[i] ~ dnorm(0, 100^-2)
   }
 
for (t in (max(p:q)+1):18) {
      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + ar_mean[t] #+ ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      #ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
   }

for (t in 19:n.occasions) {
      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2)  + ar_mean[t] #+ ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      #ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
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
   for(t in 2:11){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
    }

   for(t in 12:n.occasions){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
    }

   # calculate the N3
   for(t in 2:11-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }
   for(t in 12:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }

   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:11){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
    }

   for(t in 12:n.occasions){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
    }

   # calculate the N4
   for(t in 2:11-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
   }

   for(t in 12:n.occasions-1){
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
   
## 1.2 Test statistic: number of turns or switches - jaggedness
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





# 24: Extend model time ----
### ln scale: N2-N4: DEMOGRAPHIC MODEL ONLY

cap.v24 = '
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
   # for (t in 2:4){
   # N2[t] ~ dnorm(8.5, 1/9)
   # }
   
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
   
   
   for (t in 2:N2end) { # t = 18
      mu2[t] <- alpha2
      N2[t] ~ dnorm(mu2[t], tau.proc2)
   }     
   
   for (t in N2start:n.occasions) { #t = 19
      mu2[t] <- alpha2 + beta2*LD[t-2]
      N2[t] ~ dnorm(mu2[t], tau.proc2)
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
   for(t in 2:N3end){
      mu3[t] <- alpha3 
   }
   
   for(t in N3start:n.occasions){
      mu3[t] <- alpha3
   }
   
   # calculate the N3
   for(t in 2:n.occasions-1){
      #N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t])*mu3[t]), tau.proc3)
   }

   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
   
   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:N3end){
      mu4[t] <- alpha3
   }
   
   for(t in N3start:n.occasions){
      mu4[t] <- alpha3
   }
   
   # calculate the N4
   for(t in 2:N3end-1){
      #N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95)*mu4[t]), tau.proc4)
   }
   for(t in N3start:n.occasions-1){
#      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95)*mu4[t]), tau.proc4)
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





# 25: Extend model time ----
### ln scale: N2-N4 
#### SMOOTHER MODEL: This model combines the 1985-1998 data but with ARMA smoothers for age 2.

cap.v25 = '
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
   N3[1] ~ dnorm(8.9, 1/9)
   N4[1] ~ dnorm(6, 1/9)

# # set up residuals, and initial values of N2 and mu2
 for(t in 1:max(p, q)){
   eps[t] <- alpha2 - N2[t] # probably wrong N2[t] - alpha2
   N2[t] ~ dnorm(8.5, 1/9)
   mu2[t] ~ dunif(0,10)           
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

# AR
  for (i in 1:p) {
    eta[i] ~ dnorm(0, 100^-2)
 }

# MA
   for (i in 1:q) {
    zeta[i] ~ dnorm(0, 100^-2)
   }
 
for (t in (max(p,q)+1):18) {
      mu2[t] <- alpha2 + ar_mean[t] + ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
   }

for (t in 19:n.occasions) {
      mu2[t] <- alpha2 + ar_mean[t] + ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
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
   for(t in 2:11){
      mu3[t] <- alpha3
    }

   for(t in 12:n.occasions){
      mu3[t] <- alpha3
    }

   # calculate the N3
   for(t in 2:11-1){
      #N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t])*mu3[t]), tau.proc3)
   }
   for(t in 12:n.occasions-1){
      #N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t])*mu3[t]), tau.proc3)
   }

   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:11){
      mu4[t] <- alpha3
    }

   for(t in 12:n.occasions){
      mu4[t] <- alpha3
    }

   # calculate the N4
   for(t in 2:11-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
   }

   for(t in 12:n.occasions-1){
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
   
## 1.2 Test statistic: number of turns or switches - jaggedness
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





# 26: Extend model time ----
### ln scale: N2-N4 
#### SMOOTHER MODEL: This model combines the 1985-1998 data but with ARMA smoothers for age 2 and 3. b==6

cap.v26 = '
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
   #N3[1] ~ dnorm(8.9, 1/9)
   N4[1] ~ dnorm(6, 1/9)

# # set up residuals, and initial values of N2 and mu2
 for(t in 1:max(p, q)){
   eps[t] <- alpha2 - N2[t] # probably wrong N2[t] - alpha2
   N2[t] ~ dnorm(8.5, 1/9)
   mu2[t] ~ dunif(0,10)           
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

  for (i in 1:p) {
    eta[i] ~ dnorm(0, 100^-2)
 }

   for (i in 1:q) {
    zeta[i] ~ dnorm(0, 100^-2)
   }
 
for (t in (max(p:q)+1):18) {
      mu2[t] <- alpha2 + ar_mean[t] + ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
   }

for (t in 19:n.occasions) {
      mu2[t] <- alpha2 + ar_mean[t] + ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
   }

### N3
   # Priors for killing N3
   #alpha3 ~ dnorm(0, 100^-2)          # int
   alpha3 ~ dunif(0, 1)                # int
   gamma3 ~ dunif(0, 100)              #tices-max rate of increase
   # gamma3 ~ dunif(0, 3.65)           # was uniform
   delta3 ~ dgamma(11.5, 5.7)          #tice-width
   epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
   

# # set up residuals, and initial values of N2 and mu2
 for(t in 1:max(p, q)){
   eps3[t] <- alpha3 - N3[t]
   N3[t] ~ dnorm(8.5, 1/9)
   mu3[t] ~ dunif(0,10)           
 }
 
  for (i in 1:p) {
    eta3[i] ~ dnorm(0, 100^-2)
 }

   for (i in 1:q) {
    zeta3[i] ~ dnorm(0, 100^-2)
   }
   
   
   # calculate a survival for N2 -> N3; mu3[2:24]
   for(t in (max(p:q)+1):11){
      mu3[t] <- alpha3 + ar_mean3[t] #+ ma_mean3[t]
      ar_mean3[t] <- inprod(eta3, N3[(t-p):(t-1)])
      #ma_mean3[t] <- inprod(zeta3, eps3[(t-q):(t-1)])
      eps3[t] <- N3[t] - mu3[t]
    }

   for(t in 12:n.occasions){
      mu3[t] <- alpha3 + ar_mean3[t] #+ ma_mean3[t]
      ar_mean3[t] <- inprod(eta3, N3[(t-p):(t-1)])
      #ma_mean3[t] <- inprod(zeta3, eps3[(t-q):(t-1)])
      eps3[t] <- N3[t] - mu3[t]
    }

   # calculate the N3
   for(t in (max(p:q)+1):11-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }
   for(t in 12:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }

   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:11){
      mu4[t] <- alpha3
    }

   for(t in 12:n.occasions){
      mu4[t] <- alpha3
    }

   # calculate the N4
   for(t in 2:11-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
   }

   for(t in 12:n.occasions-1){
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
   
## 1.2 Test statistic: number of turns or switches - jaggedness
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




# 27: Extend model time ----
### ln scale: N2-N4 + forecast for each age but including LD and CO
#### This model combines the 1985-1998 data but with ARMA smoothers.b==6

cap.v27 = '
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
   N3[1] ~ dnorm(8.9, 1/9)
   N4[1] ~ dnorm(6, 1/9)

# # set up residuals, and initial values of N2 and mu2
 for(t in 1:max(p, q)){
   eps[t] <- alpha2 - N2[t] # probably wrong N2[t] - alpha2
   N2[t] ~ dnorm(8.5, 1/9)
   mu2[t] ~ dunif(0,10)           
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
   epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO

  for (i in 1:p) {
    eta[i] ~ dnorm(0, 100^-2)
 }

   for (i in 1:q) {
    zeta[i] ~ dnorm(0, 100^-2)
   }
 
for (t in (max(p:q)+1):11) {
      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + ar_mean[t] #+ ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      #ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
   }

for (t in 12:18) { # has to be 12 and not 11 bc no CO value in 1994.
      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + epsilon2*CO[t-1] + ar_mean[t] #+ ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      #ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
   }


for (t in 19:n.occasions) {
      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2) + epsilon2*CO[t-1] + ar_mean[t] #+ ma_mean[t]
      ar_mean[t] <- inprod(eta, N2[(t-p):(t-1)])
      #ma_mean[t] <- inprod(zeta, eps[(t-q):(t-1)])
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
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
   for(t in 2:11){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
    }

   for(t in 12:n.occasions){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
    }

   # calculate the N3
   for(t in 2:11-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }
   for(t in 12:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }

   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:11){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
    }

   for(t in 12:n.occasions){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
    }

   # calculate the N4
   for(t in 2:11-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
   }

   for(t in 12:n.occasions-1){
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

## 1.2 Test statistic: number of turns or switches - jaggedness
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



# 28: Extend model time ----
### ln scale: N2-N4 + forecast for each age but including LD and CO
#### This model combines the 1985-1998 data but with a Random Walk: b==5

cap.v28 = '
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
   N3[1] ~ dnorm(8.9, 1/9)
   N4[1] ~ dnorm(6, 1/9)

# # set up residuals, and initial values of N2 and mu2
 #for(t in 1:max(p, q)){
   eps[1] <- alpha2 - N2[1] # probably wrong N2[t] - alpha2
   N2[1] ~ dnorm(8.5, 1/9)
   mu2[1] ~ dunif(0,10)           
 #}

# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   alpha2 ~ dnorm(0, 100^-2)      # int
   beta2 ~ dnorm(0, 100^-2)       # larval abund
   gamma2 ~ dunif(0, 100)         # tices-max rate of increase
   delta2 ~ dgamma(11.5, 5.7)     # tice-width
   epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO
    u ~ dnorm(0, 1/16) # this is based on a SD of 2 so Var = 4- then, I doubled it for no good reason other than to see what would happen

#for (t in (max(p:q)+1):11) {
for (t in 2:11) {
#      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + u
#      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + N2[t-1] + u
      mu2[t] <-  N2[t-1] + u
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
   }

for (t in 12:18) { # has to be 12 and not 11 bc no CO value in 1994.
     # mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + epsilon2*CO[t-1] + u
#      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + epsilon2*CO[t-1] + N2[t-1] + u
      mu2[t] <-  N2[t-1] + u
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
   }


for (t in 19:n.occasions) {
     #mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2) + epsilon2*CO[t-1] + u
#      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2) + epsilon2*CO[t-1] + N2[t-1] + u
      mu2[t] <- N2[t-1] + u
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      eps[t] <- N2[t] - mu2[t]
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
   for(t in 2:11){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + u
    }

   for(t in 12:n.occasions){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +  epsilon3*CO[t-1] + u
    }

   # calculate the N3
   for(t in 2:11-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }
   for(t in 12:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }

   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]

   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:11){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
    }

   for(t in 12:n.occasions){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
    }

   # calculate the N4
   for(t in 2:11-1){
 #  N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
      N4[t+1] ~ dnorm(ifelse(log(exp(N3[t])) > 0, log(exp(N3[t])*(1-0.95))*mu4[t], 1), tau.proc4)
   }

   for(t in 12:n.occasions-1){
      #N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
      N4[t+1] ~ dnorm(ifelse(log(exp(N3[t])) > 0, log(exp(N3[t])*(1-0.95))*mu4[t], 1), tau.proc4)
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

## 1.2 Test statistic: number of turns or switches - jaggedness
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





# 29: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998
#### Change Point regression b==2

cap.v29 = '
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
   # for (t in 2:4){
   # N2[t] ~ dnorm(8.5, 1/9)
   # }
   
   # LIKELIHOODS
   ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   alpha2[1] ~ dnorm(0, 100^-2)      # int
   alpha2[2] ~ dnorm(0, 100^-2)      # int
   beta2[1] ~ dnorm(0, 100^-2)       # larval abund
   beta2[2] ~ dnorm(0, 100^-2)       # larval abund
   gamma2 ~ dunif(0, 100)         # tices-max rate of increase
   delta2 ~ dgamma(11.5, 5.7)     # tice-width
   #gamma2 ~ dnorm(0, 100^-2)     # condition # for CO
   #epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO
   t_1 ~ dunif(1985, 2021)
   mu2[1] ~ dunif(0,10)
   
   
   for (t in 2:18) { #t = 18
      mu2[t] <- alpha2[J[t]]+ beta2[J[t]]*(t[t]-t_1)
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      J[t] <- 1 + step(t[t] - t_1)
   }     
   
   for (t in 19:n.occasions) { #t = 19
      mu2[t] <- alpha2[J[t]] + beta2[J[t]]*(t[t]-t_1)*LD[t-2]
      N2[t] ~ dnorm(mu2[t], tau.proc2)
      J[t] <- 1 + step(t[t] - t_1)
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
   for(t in 2:N3end){
      mu3[t] <- alpha3 
   }
   
   for(t in N3start:n.occasions){
      mu3[t] <- alpha3
   }
   
   # calculate the N3
   for(t in 2:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }

   
   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
   
   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:N3end){
      mu4[t] <- alpha3
   }
   
   for(t in N3start:n.occasions){
      mu4[t] <- alpha3
   }
   
   # calculate the N4
   for(t in 2:N3end-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
   }
   for(t in N3start:n.occasions-1){
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






# 30: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998
#### And Soviet data - this may be a NULL

cap.v30 = '
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
   # for (t in 2:4){
   # N2[t] ~ dnorm(8.5, 1/9)
   # }
   
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
   
   
   for (t in 1:N2end) { #18
      mu2[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2)
   }     
   
   for (t in 2:N2end) { #18
      N2[t] ~ dnorm(mu2[t], tau.proc2)
   }     
   
   for (t in N2start:n.occasions) { #19
      mu2[t] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2)
      N2[t] ~ dnorm(mu2[t], tau.proc2)
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
   for(t in 2:N3end){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) 
   }
   
   for(t in N3start:n.occasions){
      mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
   }
   
   # calculate the N3
   for(t in 2:n.occasions-1){
      N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   }
   # for(t in N3start:n.occasions-1){
   #    N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
   #    }
   
   
   ###N4
   mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
   
   # calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:N3end){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
   }
   
   for(t in N3start:n.occasions){
      mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
   }
   
   # calculate the N4
   for(t in 2:N3end-1){
      N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
   }
   for(t in N3start:n.occasions-1){
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



# 31: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form

cap.v31 = '
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
   N[1,1] ~ dnorm(8.5, 1/9)    
   N[1,2] ~ dnorm(8.9, 1/9)    
   N[1,3] ~ dnorm(6, 1/9)    
   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
   # for (t in 2:4){
   # N2[t] ~ dnorm(8.5, 1/9)
   # }

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

for (a in 1:Ni){
   for (t in 1:18) { #18
      mu[t,a] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2)
   }     
}

for (a in 1:Ni){
   for (t in 2:18) { #18
      N[t,a] ~ dnorm(mu[t,a], tau.proc2)
   }     
}

for (a in 1:Ni){
   for (t in 19:n.occasions) { #19
      mu[t,a] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2)
      N[t,a] ~ dnorm(mu[t,a], tau.proc2)
   }     
}


# ### N3
#    # Priors for killing N3
#    #alpha3 ~ dnorm(0, 100^-2)          # int
#    alpha3 ~ dunif(0, 1)                # int
#    gamma3 ~ dunif(0, 100)              #tices-max rate of increase
#    # gamma3 ~ dunif(0, 3.65)           # was uniform
#    delta3 ~ dgamma(11.5, 5.7)          #tice-width
#    epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
#    mu3[1] ~ dunif(0,1)                 # need this for a mu3[1]
# 
# for(a in )   
#    # calculate a survival for N2 -> N3; mu3[2:24]
#    for(t in 2:18){
#       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) 
#     }
# 
#    for(t in 19:n.occasions){
#       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
#     }
# 
#    # calculate the N3
#    for(t in 2:n.occasions-1){
#       N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
#       }
# 
# 
#    ###N4
#    mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
# 
#    # calculate a survival for N3 -> N4; mu4[2:24]
#    for(t in 2:18){
#       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
#           }
# 
#    for(t in 19:n.occasions){
#       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
#           }
#       
#    # calculate the N4
#    for(t in 2:18-1){
#       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
#             }
#    for(t in 19:n.occasions-1){
#       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
#             }

 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index 
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.

for (a in 1:Ni){
   for (t in 1:n.occasions) {          
      matI[t,a] ~ dnorm(N[t,a], tau.obs)       # sampled observation
      # matI[2,t] ~ dnorm(N3[t], tau.obs)
      # matI[3,t] ~ dnorm(N4[t], tau.obs)
#      I[t] ~ dnorm(log(exp(I2[t]) + exp(I3[t]) + exp(I4[t])), tau.obs)
   }
}


}'



# 32: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics.

# cap.v32 = '
#  model {
# #PRIORS
# ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2[t] uninformative
#    sigma.proc2 ~ dunif(0.01, 20)     
#    sigma2.proc2 <- pow(sigma.proc2, 2) 
#    tau.proc2 <- pow(sigma.proc2, -2) 
# 
#  ## Prior for sd of process - N3[t] uninformative
#    sigma.proc3 ~ dunif(0.01, 20)     
#    sigma2.proc3 <- pow(sigma.proc3, 2) 
#    tau.proc3 <- pow(sigma.proc3, -2) 
# 
#  ## Prior for sd of process - N4[t] uninformative
#     sigma.proc4 ~ dunif(0.01, 20)     
#    sigma2.proc4 <- pow(sigma.proc4, 2) 
#    tau.proc4 <- pow(sigma.proc4, -2) 
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)   
#    sigma2.obs <- pow(sigma.obs, 2) 
#    tau.obs <- pow(sigma.obs, -2) 
#  
#    ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N[1,1] ~ dnorm(8.5, 1/9)    
#    N[1,2] ~ dnorm(8.9, 1/9)    
#    N[1,3] ~ dnorm(6, 1/9)    
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
#    # for (t in 2:4){
#    # N2[t] ~ dnorm(8.5, 1/9)
#    # }
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    alpha2 ~ dnorm(0, 100^-2)      # int
#    beta2 ~ dnorm(0, 100^-2)       # larval abund
#    gamma2 ~ dunif(0, 100)         # tices-max rate of increase
#    delta2 ~ dgamma(11.5, 5.7)     # tice-width
#    #gamma2 ~ dnorm(0, 100^-2)     # condition # for CO
#    #epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO
# 
# for (a in 1:Ni){
#    for (t in 1:18) { #18
#       mu[t,a] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2)
#    }     
# }
# 
# for (a in 1:Ni){
#    for (t in 2:18) { #18
#       N[t,a] ~ dnorm(mu[t,a], tau.proc2)
#    }     
# }
# 
# for (a in 1:Ni){
#    for (t in 19:n.occasions) { #19
#       mu[t,a] <- alpha2 + beta2*LD[t-2] + gamma2*TI[t]*(1-TI[t]/delta2)
#       N[t,a] ~ dnorm(mu[t,a], tau.proc2)
#    }     
# }
# 
# 
# # ### N3
# #    # Priors for killing N3
# #    #alpha3 ~ dnorm(0, 100^-2)          # int
# #    alpha3 ~ dunif(0, 1)                # int
# #    gamma3 ~ dunif(0, 100)              #tices-max rate of increase
# #    # gamma3 ~ dunif(0, 3.65)           # was uniform
# #    delta3 ~ dgamma(11.5, 5.7)          #tice-width
# #    epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
# #    mu3[1] ~ dunif(0,1)                 # need this for a mu3[1]
# # 
# # for(a in )   
# #    # calculate a survival for N2 -> N3; mu3[2:24]
# #    for(t in 2:18){
# #       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) 
# #     }
# # 
# #    for(t in 19:n.occasions){
# #       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
# #     }
# # 
# #    # calculate the N3
# #    for(t in 2:n.occasions-1){
# #       N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
# #       }
# # 
# # 
# #    ###N4
# #    mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
# # 
# #    # calculate a survival for N3 -> N4; mu4[2:24]
# #    for(t in 2:18){
# #       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
# #           }
# # 
# #    for(t in 19:n.occasions){
# #       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
# #           }
# #       
# #    # calculate the N4
# #    for(t in 2:18-1){
# #       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
# #             }
# #    for(t in 19:n.occasions-1){
# #       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
# #             }
# 
#  ## Observation
#    ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
#    #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index 
#    #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
# 
# for (a in 1:Ni){
#    for (t in 1:n.occasions) {          
#       matI[t,a] ~ dnorm(N[t,a], tau.obs)       # sampled observation
#       # matI[2,t] ~ dnorm(N3[t], tau.obs)
#       # matI[3,t] ~ dnorm(N4[t], tau.obs)
#       
#    }
# }
# 
# for (t in 1:n.occasions) { 
# I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs)
# }
# 
#    # Assessing the fit of the state-space model
#    ## 1. Compute fit statistics for observed data.
#    ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#       I.exp[t] <- log(exp(N[t,1]) + exp(N[t,2]) + exp(N[t,3]))
#       Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
#    
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
#    
#    
#    ## 2.1 Simulated data
#    for (t in 1:n.occasions){
#       #    y2.rep[t] ~ dnorm(N2[t], tau.obs) 
#       #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#       I2.rep[t] ~ dnorm(N[t,1], tau.obs)
#       I3.rep[t] ~ dnorm(N[t,2], tau.obs)
#       I4.rep[t] ~ dnorm(N[t,3], tau.obs)
#       I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
#       Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
#    }
#    Dmape.rep <- sum(Dssm.rep)
#    
#    
#    ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# 
# 
# }'



# # 33: Extend model time----
# ## ln scale: N2-N4 + forecast for each age
# ### Added in the capelin data from 1985-1998 but in matrix form
# #### AS 31 but with I and all the diagnostics - this adds parameters to the loops
# 
# cap.v33 = '
#  model {
# #PRIORS
# ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2[t] uninformative
#    # sigma.proc2 ~ dunif(0.01, 20)
#    # sigma2.proc2 <- pow(sigma.proc2, 2)
#    # tau.proc2 <- pow(sigma.proc2, -2)
# 
#  # ## Prior for sd of process - N3[t] uninformative
# #    sigma.proc3 ~ dunif(0.01, 20)
#  #   sigma2.proc3 <- pow(sigma.proc3, 2)
#   #  tau.proc3 <- pow(sigma.proc3, -2)
#  #
#  # ## Prior for sd of process - N4[t] uninformative
#  #    sigma.proc4 ~ dunif(0.01, 20)
#  #   sigma2.proc4 <- pow(sigma.proc4, 2)
#  #   tau.proc4 <- pow(sigma.proc4, -2)
# 
# for(a in 1:2){
#     sigma.proc[a] ~ dunif(0.01, 20)
#    sigma2.proc[a] <- pow(sigma.proc[a], 2)
#    tau.proc[a] <- pow(sigma.proc[a], -2)
# }
# 
# #for(a in 2:Ni){
#    sigma.proc3 ~ dunif(0.01, 20)
#    sigma2.proc3 <- pow(sigma.proc3, 2)
#    tau.proc3 <- pow(sigma.proc3, -2)
# #}
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)
#    sigma2.obs <- pow(sigma.obs, 2)
#    tau.obs <- pow(sigma.obs, -2)
# 
#    ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N[1,1] ~ dnorm(8.5, 1/9)
#    N[1,2] ~ dnorm(8.9, 1/9)
#    N[1,3] ~ dnorm(6, 1/9)
# 
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# # for (t in 2:4){
#  #   N2[t] ~ dnorm(8.5, 1/9)
#   #  }
# 
# for(a in 1:Ni){
# mu[1, a] ~ dunif(0,10)
# }
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    for(a in 1:2){
#    alpha[a] ~ dnorm(0, 100^-2)      # int
#    beta[a] ~ dnorm(0, 100^-2)       # larval abund
#    gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
#    delta[a] ~ dgamma(11.5, 5.7)     # tice-width
#    epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
# }
# 
#    # gamma2 ~ dnorm(0, 100^-2)     # condition # for CO
#     epsilon2 ~ dnorm(0, 100^-2)   # condition # for CO
# 
# # mu:N2; 1985:2002
#    # for (t in 1:18) { #18
#    #    mu[t,1] <- alpha[1] + gamma[1]*TI[t]*(1-TI[t]/delta[1])
#    # }
# 
# # mu:N3-4; 1985:2002
# for (a in 1:2){
#    for (t in 2:18) { #18
#       mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
#    }
# }
# 
# 
# # N: N2-3; 1985:2002
#    for (t in 2:18) { #18
#       N[t,1] ~ dnorm(mu[t,1], tau.proc[1])
#    }
# 
#    for (t in 2:18-1) { #18
#       N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2])
#    }
# 
#  for (t in 2:18-1) { #18
#       N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,2], tau.proc[2])
#  }
# 
# aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
# 
# # mu and N2: 2003-present
# for (t in 19:n.occasions) { #19
#       mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
#       N[t,1] ~ dnorm(mu[t,1], tau.proc[1])
#    }
# 
# # mu: N3 and N4; 2003-present
# ##
# for (a in 2:Ni){
#   for (t in 19:n.occasions) { #19
#     mu[t,a] <- alpha[2] + gamma[2]*TI[t]*(1-TI[t]/delta[2]) + epsilon[2]*CO[t-1]
#    }
# }
# 
# 
# #for (a in 2:Ni){
#  #  for (t in 19:n.occasions-1) { #19
#   #    N[t+1,a] ~ dnorm(log(exp(N[t,a-1])*(1-m[t]))*mu[t,a], tau.proc[a])
#    #}
# #}
# #N3; 2003-present
#    for (t in 19:n.occasions-1) { #19
#       N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2])
#    }
# #N4; 2003-present
#    for (t in 19:n.occasions-1) { #19
#       N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,2], tau.proc[2])
#    }


# ### N3
#    # Priors for killing N3
#    #alpha3 ~ dnorm(0, 100^-2)          # int
#    alpha3 ~ dunif(0, 1)                # int
#    gamma3 ~ dunif(0, 100)              #tices-max rate of increase
#    # gamma3 ~ dunif(0, 3.65)           # was uniform
#    delta3 ~ dgamma(11.5, 5.7)          #tice-width
#    epsilon3 ~ dnorm(0, 100^-2)         # condition # for CO
#    mu3[1] ~ dunif(0,1)                 # need this for a mu3[1]
#
# for(a in )
#    # calculate a survival for N2 -> N3; mu3[2:24]
#    for(t in 2:18){
#       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
#     }
#
#    for(t in 19:n.occasions){
#       mu3[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) +             epsilon3*CO[t-1]
#     }
#
#    # calculate the N3
#    for(t in 2:n.occasions-1){
#       N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))*mu3[t], tau.proc3)
#       }
#
#
#    ###N4
#    mu4[1] ~ dunif(0,1)              # need this for a mu3[1]
#
#    # calculate a survival for N3 -> N4; mu4[2:24]
#    for(t in 2:18){
#       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3)
#           }
#
#    for(t in 19:n.occasions){
#       mu4[t] <- alpha3 + gamma3*TI[t]*(1-TI[t]/delta3) + epsilon3*CO[t-1]
#           }
#
#    # calculate the N4
#    for(t in 2:18-1){
#       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
#             }
#    for(t in 19:n.occasions-1){
#       N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
#             }

 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.

# for (a in 1:Ni){
#    for (t in 1:n.occasions) {
#       matI[t,a] ~ dnorm(N[t,a], tau.obs)       # sampled observation
#       # matI[2,t] ~ dnorm(N3[t], tau.obs)
#       # matI[3,t] ~ dnorm(N4[t], tau.obs)
# 
#    }
# }
# 
# for (t in 1:n.occasions) {
# I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs)
# }
# 
#    # Assessing the fit of the state-space model
#    ## 1. Compute fit statistics for observed data.
#    ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#       I.exp[t] <- log(exp(N[t,1]) + exp(N[t,2]) + exp(N[t,3]))
#       Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
# 
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
# 
# 
#    ## 2.1 Simulated data
#    for (t in 1:n.occasions){
#       #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
#       #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#       I2.rep[t] ~ dnorm(N[t,1], tau.obs)
#       I3.rep[t] ~ dnorm(N[t,2], tau.obs)
#       I4.rep[t] ~ dnorm(N[t,3], tau.obs)
#       I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
#       Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
#    }
#    Dmape.rep <- sum(Dssm.rep)
# 
# 
#    ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# 
# }'

# 33: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code
##### no eps

cap.v33 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
for(a in 1:2){
    sigma.proc[a] ~ dunif(0.01, 20)
   sigma2.proc[a] <- pow(sigma.proc[a], 2)
   tau.proc[a] <- pow(sigma.proc[a], -2)
}

#for(a in 2:Ni){
   # sigma.proc3 ~ dunif(0.01, 20)
   # sigma2.proc3 <- pow(sigma.proc3, 2)
   # tau.proc3 <- pow(sigma.proc3, -2)
#}

 ## Prior for sd of observation - I2-I4[t] - uninformative
   sigma.obs ~ dunif(0.01, 20)
   sigma2.obs <- pow(sigma.obs, 2)
   tau.obs <- pow(sigma.obs, -2)

### Priors for Initial values for N2-N4[t] informative - based on actual values
   N[1,1] ~ dnorm(8.5, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)

   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }

# Priors for initial value of mu - may not need this
for(a in 1:Ni){
   mu[1, a] ~ dunif(0,10)
}

# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   for(a in 1:2){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}


# mu:N2-3; 1985:2002
for (a in 1:2){
   for (t in 2:18) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
   }
}


# N: N2-3; 1985:2002
   for (t in 2:18) { #18
      N[t,1] ~ dnorm(mu[t,1], tau.proc[1])
   }

   for (t in 2:18-1) { #18
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2])
   }

 for (t in 2:18-1) { #18
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,2], tau.proc[2])
 }

######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
      mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
      N[t,1] ~ dnorm(mu[t,1], tau.proc[1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    mu[t,a] <- alpha[2] + gamma[2]*TI[t]*(1-TI[t]/delta[2]) + epsilon[2]*CO[t-1]
   }
}


#N3; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2])
   }

#N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,2], tau.proc[2])
   }


 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.

for (a in 1:Ni){
   for (t in 1:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs)       # sampled observation
   }
}

for (t in 1:n.occasions) {
I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs)
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
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
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





# 34: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code
##### eps added for N2-N4

# cap.v34 = '
#  model {
# #PRIORS
# ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2-N3[t] uninformative
# for(a in 1:Ni){
#     sigma.proc[a] ~ dunif(0.01, 20)
#    sigma2.proc[a] <- pow(sigma.proc[a], 2)
#    tau.proc[a] <- pow(sigma.proc[a], -2)
# }
# 
# #for(a in 2:Ni){
#    # sigma.proc3 ~ dunif(0.01, 20)
#    # sigma2.proc3 <- pow(sigma.proc3, 2)
#    # tau.proc3 <- pow(sigma.proc3, -2)
# #}
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)
#    sigma2.obs <- pow(sigma.obs, 2)
#    tau.obs <- pow(sigma.obs, -2)
# 
# ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N[1,1] ~ dnorm(8.5, 1/9)
#    N[1,2] ~ dnorm(8.9, 1/9)
#    N[1,3] ~ dnorm(6, 1/9)
# 
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# # for (t in 2:4){
#  #   N2[t] ~ dnorm(8.5, 1/9)
#   #  }
# 
# # Priors for initial value of mu - may not need this
# for(a in 1:Ni){
#    mu[1, a] ~ dunif(0,10)
#    #eps[1,a] <- mu[1,a] - N[1,a] # probably wrong N[t] - mu
# }
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    for(a in 1:Ni){
#    alpha[a] ~ dnorm(0, 100^-2)      # int
#    beta[a] ~ dnorm(0, 100^-2)       # larval abund
#    gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
#    delta[a] ~ dgamma(11.5, 5.7)     # tice-width
#    epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
# }
# 
# 
# # mu:N2-3; 1985:2002
# for (a in 1:Ni){
#    for (t in 2:18) { #18
#       mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
#    }
# }
# 
# 
# # N: N2-3; 1985:2002
#    for (t in 2:18) { #18
#       N[t,1] ~ dnorm(mu[t,1], tau.proc[1])
#    }
# 
#    for (t in 2:18-1) { #18
#       N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2])
#    }
# 
#  for (t in 2:18-1) { #18
#       N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[3])
#  }
# 
# ######################################
# 
# # mu and N2: 2003-present
# for (t in 19:n.occasions) { #19
#       mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
#       N[t,1] ~ dnorm(mu[t,1], tau.proc[1])
#    }
# 
# # mu: N3 and N4; 2003-present
# ##
# for (a in 2:Ni){
#   for (t in 19:n.occasions) { #19
#     mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
#    }
# }
# 
# 
# #N3; 2003-present
#    for (t in 19:n.occasions-1) { #19
#       N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2])
#    }
# 
# #N4; 2003-present
#    for (t in 19:n.occasions-1) { #19
#       N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[3])
#    }
# 
# 
# for (a in 1:Ni){
#      for (t in 1:n.occasions){
#      eps[t,a] <- mu[t,a] - N[t,a]
#      }
# }
# 
# 
#  ## Observation
#    ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
#    #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
#    #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
# 
# for (a in 1:Ni){
#    for (t in 1:n.occasions) {
#       matI[t,a] ~ dnorm(N[t,a], tau.obs)       # sampled observation
#    }
# }
# 
# for (t in 1:n.occasions) {
#      I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs)
# }
# 
#    # Assessing the fit of the state-space model
#    ## 1. Compute fit statistics for observed data.
#    ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#       I.exp[t] <- log(exp(N[t,1]) + exp(N[t,2]) + exp(N[t,3]))
#       Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
# 
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
# 
# 
#    ## 2.1 Simulated data
#    for (t in 1:n.occasions){
#       #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
#       #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#       I2.rep[t] ~ dnorm(N[t,1], tau.obs)
#       I3.rep[t] ~ dnorm(N[t,2], tau.obs)
#       I4.rep[t] ~ dnorm(N[t,3], tau.obs)
#       I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
#       Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
#    }
#    Dmape.rep <- sum(Dssm.rep)
# 
# 
#    ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# 
# }'
# 
# 
# 
# # 35: Extend model time----
# ## ln scale: N2-N4 + forecast for each age
# ### Added in the capelin data from 1985-1998 but in matrix form
# #### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code
# ##### eps added for N2-N4 and osa
# #####  minimize code and reduce tau.proc to 2 values
# 
# cap.v34 = '
#  model {
# #PRIORS
# ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2-N3[t] uninformative
# for(a in 1:2){
#     sigma.proc[a] ~ dunif(0.01, 20)
#    sigma2.proc[a] <- pow(sigma.proc[a], 2)
#    tau.proc[a] <- pow(sigma.proc[a], -2)
# }
# 
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#    sigma.obs ~ dunif(0.01, 20)
#    sigma2.obs <- pow(sigma.obs, 2)
#    tau.obs <- pow(sigma.obs, -2)
# 
# ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N[1,1] ~ dnorm(8.5, 1/9)
#    N[1,2] ~ dnorm(8.9, 1/9)
#    N[1,3] ~ dnorm(6, 1/9)
# 
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# # for (t in 2:4){
#  #   N2[t] ~ dnorm(8.5, 1/9)
#   #  }
# 
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    for(a in 1:Ni){
#    alpha[a] ~ dnorm(0, 100^-2)      # int
#    beta[a] ~ dnorm(0, 100^-2)       # larval abund
#    gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
#    delta[a] ~ dgamma(11.5, 5.7)     # tice-width
#    epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
# }
# 
# 
# # mu:N2-3; 1985:2002
# for (a in 1:Ni){
#    for (t in 1:18) { #18
#       mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
#    }
# }
# 
# 
# # N: N2-3; 1985:2002
#    for (t in 2:18) { #18
#       N[t,1] ~ dnorm(mu[t,1], tau.proc[1])
#    }
# 
#    for (t in 2:18-1) { #18
#       N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2]) #N3
#       N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2])
#    } #N4
# 
# 
# ######################################
# 
# # mu and N2: 2003-present
# for (t in 19:n.occasions) { #19
#       mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
#       N[t,1] ~ dnorm(mu[t,1], tau.proc[1])
#    }
# 
# # mu: N3 and N4; 2003-present
# ##
# for (a in 2:Ni){
#   for (t in 19:n.occasions) { #19
#     mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
#    }
# }
# 
# 
# #N3 & N4; 2003-present
#    for (t in 19:n.occasions-1) { #19
#       N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2]) #N3
#       N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2]) #N4
#    }
# 
# 
# # process error
# for (a in 1:Ni){
#      for (t in 1:n.occasions){
#      # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
#      #eps[t,a] <- N[t,a] - matI[t,a]
#      #eps[t,a] <- N[t,a] - mean(N[,a])
#      eps[t,a] <- N[t,a] - mu[t,a]
#      }
# }
# 
# # one step ahead resids - could set a p
# for (a in 1:Ni){
#      osa[1, a] ~ dnorm(0, 1/10)
# }
# 
# for (a in 1:Ni){
#      for (t in 2:n.occasions){
#      osa_mean[t,a] <- mean(N[1:(t-1), a])
#      osa[t,a] <- N[t,a] - osa_mean[t,a]
#      }
# }
# 
#  ## Observation
#    ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
#    #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
#    #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
# 
# for (a in 1:Ni){
#    for (t in 1:n.occasions) {
#       matI[t,a] ~ dnorm(N[t,a], tau.obs)       # sampled observation
#    }
# }
# 
# for (t in 1:n.occasions) {
#      # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs)
#      I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
# }
# 
#    # Assessing the fit of the state-space model
#    ## 1. Compute fit statistics for observed data.
#    ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#       I.exp[t] <- log(exp(N[t,1]) + exp(N[t,2]) + exp(N[t,3]))
#       Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
# 
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
# 
# 
#    ## 2.1 Simulated data
#    for (t in 1:n.occasions){
#       #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
#       #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#       I2.rep[t] ~ dnorm(N[t,1], tau.obs)
#       I3.rep[t] ~ dnorm(N[t,2], tau.obs)
#       I4.rep[t] ~ dnorm(N[t,3], tau.obs)
#       I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs)
#       Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
#    }
#    Dmape.rep <- sum(Dssm.rep)
# 
# 
#    ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# 
# 
# 
# }'


# 36: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) - these correspond to the data series and NOT the pre/post collapse
#### added osa and posa resids

cap.v34 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, 
for(a in 1:2){
     for(tp in 1:2){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:2){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
   N[1,1] ~ dnorm(8.5, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)

   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   for(a in 1:Ni){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}


# mu:N2-3; 1985:2002
for (a in 1:Ni){
   for (t in 1:18) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
   }
}


# N: N2-3; 1985:2002
   for (t in 2:18) { #18
      N[t,1] ~ dnorm(mu[t,1], tau.proc[1,1])
   }

   for (t in 2:18-1) { #18
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[1, 2]) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[1, 2])
   } #N4


######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
      mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
      N[t,1] ~ dnorm(mu[t,1], tau.proc[2, 1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
   }
}


#N3 & N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2, 2]) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2, 2]) #N4
   }


# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     #eps[t,a] <- N[t,a] - matI[t,a]
     #eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids - could set a p
for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     posa[1, a] ~ dnorm(0, 1/10)
}

for (a in 1:Ni){
     for (t in 2:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     posa[t,a] <- osa[t,a]/sd(osa[,a])
     }
}
 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}


for (t in 8:n.occasions){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
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


# 37: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids

# cap.v35 = '
#  model {
# #PRIORS
# ###### Need to check that these are reasonable
#  ## Prior for sd of process - N2-N3[t] uninformative
#  # a is process variance for age group, tp for time period pre/post collapse
# for(a in 1:2){
#      for(tp in 1:3){
#         sigma.proc[tp, a] ~ dunif(0.01, 20)
#         sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
#         tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
#      }
# }
# 
# 
#  ## Prior for sd of observation - I2-I4[t] - uninformative
#  for(tp in 1:3){
#    sigma.obs[tp] ~ dunif(0.01, 20)
#    sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
#    tau.obs[tp] <- pow(sigma.obs[tp], -2)
#  }
# 
# 
# ### Priors for Initial values for N2-N4[t] informative - based on actual values
#    N[1,1] ~ dnorm(8.5, 1/9)
#    N[1,2] ~ dnorm(8.9, 1/9)
#    N[1,3] ~ dnorm(6, 1/9)
# 
#    ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# # for (t in 2:4){
#  #   N2[t] ~ dnorm(8.5, 1/9)
#   #  }
# 
# 
# # LIKELIHOODS
#  ## State process
#    ### N2
#    # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
#    # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
#    for(a in 1:Ni){
#    alpha[a] ~ dnorm(0, 100^-2)      # int
#    beta[a] ~ dnorm(0, 100^-2)       # larval abund
#    gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
#    delta[a] ~ dgamma(11.5, 5.7)     # tice-width
#    epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
# }
# 
# 
# # mu:N2-3; 1985:1990
# for (a in 1:Ni){
#    for (t in 1:6) { #18
#       mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
#    }
# }
# 
# 
# # N: N2-3; 1985:1990
#    for (t in 2:6) { #18
#       N[t,1] ~ dnorm(mu[t,1], tau.proc[1,1])
#    }
# 
#    for (t in 2:6-1) { #18
#       N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[1, 2]) #N3
#       N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[1, 2])
#    } #N4
# 
# 
# for (a in 1:Ni){
#    for (t in 7:18) { #18
#       mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
#    }
# }
# 
# # mu:N2-3; 1991:2002
# # N: N2-3; 1985:1990
#    for (t in 7:18) { #18
#       N[t,1] ~ dnorm(mu[t,1], tau.proc[1,1])
#    }
# 
#    for (t in 7:18-1) { #18
#       N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2, 2]) #N3
#       N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2, 2])
#    } #N4
# 
# ######################################
# 
# # mu and N2: 2003-present
# for (t in 19:n.occasions) { #19
#       mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
#       N[t,1] ~ dnorm(mu[t,1], tau.proc[2, 1])
#    }
# 
# # mu: N3 and N4; 2003-present
# ##
# for (a in 2:Ni){
#   for (t in 19:n.occasions) { #19
#     mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
#    }
# }
# 
# 
# #N3 & N4; 2003-present
#    for (t in 19:n.occasions-1) { #19
#       N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2, 2]) #N3
#       N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2, 2]) #N4
#    }
# 
# 
# # process error
# for (a in 1:Ni){
#      for (t in 1:n.occasions){
#      # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
#      #eps[t,a] <- N[t,a] - matI[t,a]
#      #eps[t,a] <- N[t,a] - mean(N[,a])
#      eps[t,a] <- N[t,a] - mu[t,a]
#      }
# }
# 
# # one step ahead resids - could set a p
# for (a in 1:Ni){
#      osa[1, a] ~ dnorm(0, 1/10)
#      posa[1, a] ~ dnorm(0, 1/10)
# }
# 
# for (a in 1:Ni){
#      for (t in 2:n.occasions){
#      osa_mean[t,a] <- mean(N[1:(t-1), a])
#      osa[t,a] <- N[t,a] - osa_mean[t,a]
#      posa[t,a] <- osa[t,a]/sd(osa[,a])
#      }
# }
#  ## Observation
#    ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
#    #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
#    #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
# 
# 
# for (a in 1:Ni){
#    for (t in 1:7) {
#       matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
#    }
# }
# 
# 
# for (a in 1:Ni){
#    for (t in 8:n.occasions) {
#       matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
#    }
# }
# 
# 
# for (t in 1:7) {
#      # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
#      I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
# }
# 
# 
# for (t in 8:n.occasions) {
#      # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
#      I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
# }
# 
# 
#    # Assessing the fit of the state-space model
#    ## 1. Compute fit statistics for observed data.
#    ### 1.1 Discrepancy meansure: mean absolute error
#    for (t in 1:n.occasions) {
#       I.exp[t] <- log(exp(N[t,1]) + exp(N[t,2]) + exp(N[t,3]))
#       Dssm.obs[t] <- abs((I[t] - I.exp[t])/I[t])
#    }
#    Dmape.obs <- sum(Dssm.obs)
# 
#    # ## 1.2 Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.obs[t] <- step(I[t+2] - I[t+1])
#       Tt2.obs[t] <- step(I[t+1] - I[t])
#       # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
#       # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
#       # Tt22.obs[t] <- step(I2[t+1] - I2[t])
#       # Tt23.obs[t] <- step(I3[t+1] - I3[t])
#       # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
#       # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
#       Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
#    }
#    Tturn.obs <- sum(Tt3.obs)
# 
# 
#    ## 2.1 Simulated data
# for (t in 1:7){
#       #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
#       #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#       I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
#       I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
#       I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
#       I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
#       Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
# }
# 
# 
# for (t in 8:n.occasions){
#       #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
#       #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
#       I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
#       I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
#       I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
#       I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
#       Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
# }
# 
#    Dmape.rep <- sum(Dssm.rep)
# 
# 
#    ##Test statistic: number of turns or switches - jaggedness
#    for (t in 1:(n.occasions-2)){
#       Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
#       Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
#       Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
#    }
#    Tturn.rep <- sum(Tt3.rep)
# 
# }'

# 37: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids

cap.v35 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, tp for time period pre/post collapse
for(a in 1:2){
     for(tp in 1:3){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:3){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
   N[1,1] ~ dnorm(8.5, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)

   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   for(a in 1:Ni){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}


# mu:N2-3; 1985:1990

for (a in 1:Ni){
   for (t in 1:18) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
   }
}

# N: N2-3; 1985:1990
   for (t in 2:18) { #18
      N[t,1] ~ dnorm(mu[t,1], ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
   }

   for (t in 2:18-1) { #18
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2])) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2]))
   } #N4



######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
      mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
      N[t,1] ~ dnorm(mu[t,1], tau.proc[2, 1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
   }
}


#N3 & N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2, 2]) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2, 2]) #N4
   }


# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     #eps[t,a] <- N[t,a] - matI[t,a]
     #eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids - could set a p
for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     posa[1, a] ~ dnorm(0, 1/10)
}

for (a in 1:Ni){
     for (t in 2:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     posa[t,a] <- osa[t,a]/sd(osa[,a])
     }
}
 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}


for (t in 8:n.occasions){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
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

# 38: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids
#### add TI, CO, and LD in time appropriate periods

cap.v36 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, tp for time period pre/post collapse
for(a in 1:2){
     for(tp in 1:3){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:3){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
   N[1,1] ~ dnorm(8.5, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)

   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   for(a in 1:Ni){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}


#mu:N2-3; 1985:1990 - would prefer to use this code but subset error 
# for (a in 1:Ni){
#    for (t in 1:18) { #18
#       mu[t,a] <- ifelse(t<12,
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]),
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1])
#       # this is OK for now but need CO by age and this may not be appro for age 2
#    }
# }

for (a in 1:Ni){
   for (t in 1:11) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
   }
}

for (a in 1:Ni){
   for (t in 12:18) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
   }
}


# N2-3; 1985:1990
   for (t in 2:18) { #18
      N[t,1] ~ dnorm(mu[t,1], 
      ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
   }

   for (t in 2:18-1) { #18
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2])) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2]))
   } #N4



# ######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
      mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
      N[t,1] ~ dnorm(mu[t,1], tau.proc[2, 1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
   }
}


#N3 & N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2, 2]) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2, 2]) #N4
   }


# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     #eps[t,a] <- N[t,a] - matI[t,a]
     #eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids - could set a p
for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     posa[1, a] ~ dnorm(0, 1/10)
}

for (a in 1:Ni){
     for (t in 2:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     posa[t,a] <- osa[t,a]/sd(osa[,a])
     }
}
 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}


for (t in 8:n.occasions){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
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

# 39: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids
#### add rdm walk

cap.v37 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, tp for time period pre/post collapse
for(a in 1:2){
     for(tp in 1:3){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:3){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
   N[1,1] ~ dnorm(8.5, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)
   # mu[1,1] ~ dnorm(0, 10)
   # mu[1,2] ~ dnorm(0, 10)
   # mu[1,3] ~ dnorm(0, 10)
   mu[1,1] ~ dunif(0,10)
   mu[1,2] ~ dunif(0,10)
   mu[1,3] ~ dunif(0,10)


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   for(a in 1:Ni){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}
   u ~ dnorm(0, 1/16) # this is based on a SD of 2 so Var = 4- then, I doubled it for no good reason other than to see what would happen



   for (t in 2:n.occasions) { #18
#     mu[t,1] <- alpha[1] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + N[t-1, 1] 
     # mu[t,1] <- ifelse(N[t-1, 1] > 0, alpha[1] + N[t-1, 1] + u, alpha[1])
      mu[t,1] <- N[t-1, 1] + u
   }

for (a in 2:Ni){
   for (t in 2:n.occasions) { #18
     #mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) 
     #mu[t,a] <- ifelse(N[t-1, a] > 0, alpha[a] + N[t-1, a] + u, alpha[a] + u)
#     mu[t,a] <- ifelse(N[t-1, a] > 0, N[t-1, a] + u, u)
      mu[t,a] <- N[t-1, a] + u
   }
}


# N2-3; 1985:1990
   for (t in 2:n.occasions) { #18
      N[t,1] ~ dnorm(mu[t,1], 
      ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
   }

   for (t in 2:n.occasions) { #18
    #  N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2])) #N3
#      N[t+1,3] ~ dnorm(ifelse(N[t,2] > 0, log(exp(N[t,2])*(1-0.95))*mu[t,3], 2), ifelse(t<=6, tau.proc[1, 2], tau.proc[2, 2]))
     N[t,2] ~ dnorm(mu[t,2], ifelse(t<=6, tau.proc[1,2], tau.proc[2,2])) 
     N[t,3] ~ dnorm(mu[t,3], ifelse(t<=6, tau.proc[1,2], tau.proc[2,2])) 
   } #N4

######----------

# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     #eps[t,a] <- N[t,a] - matI[t,a]
     #eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids - could set a p
for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     posa[1, a] ~ dnorm(0, 1/10)
}

for (a in 1:Ni){
     for (t in 2:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     #posa[t,a] <- (osa[t,a]+0.001)/sd(osa[,a])
     }
}


 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}

for (t in 8:n.occasions){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
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


# 40: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids
#### add TI, CO, and LD in time appropriate periods
#### replace rdm walk with AR

cap.v38 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, tp for time period pre/post collapse
for(a in 1:2){
     for(tp in 1:3){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:3){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
   N[1,1] ~ dnorm(8.5, 1/9)
   N[2,1] ~ dnorm(8.9, 1/9)
   #N[3,1] ~ dnorm(6, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[2,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)
   N[2,3] ~ dnorm(6, 1/9)
   mu[1,1] ~ dunif(0,10)
   mu[2,1] ~ dunif(0,10)
   mu[1,2] ~ dunif(0,10)
   mu[2,2] ~ dunif(0,10)
   mu[1,3] ~ dunif(0,10)
   mu[2,3] ~ dunif(0,10)

   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   for(a in 1:Ni){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}
   
# AR
#for (a in 1:Ni){
  for (i in 1:p) {
  #  eta[i,a] ~ dnorm(0, 100^-2)
   eta[i] ~ dnorm(0, 100^-2)
 }
#}


#mu:N2-3; 1985:1990 - would prefer to use this code but subset error 
# for (a in 1:Ni){
#    for (t in 1:18) { #18
#       mu[t,a] <- ifelse(t<12,
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]),
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1])
#       # this is OK for now but need CO by age and this may not be appro for age 2
#    }
# }

for (a in 1:Ni){
   for (t in (p+1):11) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + ar_mean[t,a]
      ar_mean[t,a] <- inprod(eta, N[(t-p):(t-1), a])
   }
}

for (a in 1:Ni){
   for (t in 12:18) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1] 
   }
}


# N2-3; 1985:1990
   for (t in (p+1):18) { #18
      N[t,1] ~ dnorm(mu[t,1], 
      ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
   }

   for (t in (p+1):18-1) { #18
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2])) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2]))
   } #N4



# ######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
      mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1] 
      N[t,1] ~ dnorm(mu[t,1], tau.proc[2, 1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
  }
}


#N3 & N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2, 2]) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2, 2]) #N4
   }


# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     #eps[t,a] <- N[t,a] - matI[t,a]
     #eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids - could set a p
for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     posa[1, a] ~ dnorm(0, 1/10)
}

for (a in 1:Ni){
     for (t in 2:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     posa[t,a] <- osa[t,a]/sd(osa[,a])
     }
}
 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}


for (t in 8:n.occasions){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
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

# 41: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids
#### add TI, CO, and LD in time appropriate periods
#### replace rdm walk with MA

cap.v39 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, tp for time period pre/post collapse
for(a in 1:2){
     for(tp in 1:3){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:3){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
   N[1,1] ~ dnorm(8.5, 1/9)
   N[2,1] ~ dnorm(8.9, 1/9)
   #N[3,1] ~ dnorm(6, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[2,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)
   N[2,3] ~ dnorm(6, 1/9)
   mu[1,1] ~ dunif(0,10)
   mu[2,1] ~ dunif(0,10)
   mu[1,2] ~ dunif(0,10)
   mu[2,2] ~ dunif(0,10)
   mu[1,3] ~ dunif(0,10)
   mu[2,3] ~ dunif(0,10)
   # eps[1,1] ~ dunif(0,10)
   # eps[2,1] ~ dunif(0,10)
   # eps[1,2] ~ dunif(0,10)
   # eps[2,2] ~ dunif(0,10)
   # eps[1,3] ~ dunif(0,10)
   # eps[2,3] ~ dunif(0,10)

   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
  for(a in 1:Ni){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}
   
# AR
#for (a in 1:Ni){
#  for (i in 1:p) {
  #  eta[i,a] ~ dnorm(0, 100^-2)
 #  eta[i] ~ dnorm(0, 100^-2)
# }
#}

# MA
   for (i in 1:q) {
    zeta[i] ~ dnorm(0, 100^-2)
   }


#mu:N2-3; 1985:1990 - would prefer to use this code but subset error 
# for (a in 1:Ni){
#    for (t in 1:18) { #18
#       mu[t,a] <- ifelse(t<12,
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]),
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1])
#       # this is OK for now but need CO by age and this may not be appro for age 2
#    }
# }

for (a in 1:Ni){
   for (t in (max(p,q)+1):11) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])  + ma_mean[t,a]
      ma_mean[t,a] <- inprod(zeta, eps[(t-q):(t-1), a])
   }
}

for (a in 1:Ni){
   for (t in 12:18) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]  + ma_mean[t,a]
      ma_mean[t,a] <- inprod(zeta, eps[(t-q):(t-1), a])
   }
}


# N2-3; 1985:1990
   for (t in (max(p,q)+1):18) { #18
      N[t,1] ~ dnorm(mu[t,1], 
      ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
   }

   for (t in (max(p,q)+1):18-1) { #18
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2])) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2]))
   } #N4



# ######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
      mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1] + ma_mean[t,1]
      ma_mean[t,1] <- inprod(zeta, eps[(t-q):(t-1), 1])
      N[t,1] ~ dnorm(mu[t,1], tau.proc[2, 1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1] + ma_mean[t,a]
      ma_mean[t,a] <- inprod(zeta, eps[(t-q):(t-1), a])
  }
}


#N3 & N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2, 2]) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2, 2]) #N4
   }


# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     #eps[t,a] <- N[t,a] - matI[t,a]
     #eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids - could set a p
for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     posa[1, a] ~ dnorm(0, 1/10)
}

for (a in 1:Ni){
     for (t in 2:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     posa[t,a] <- osa[t,a]/sd(osa[,a])
     }
}
 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}


for (t in 8:n.occasions){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
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

# 42: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids
#### add TI, CO, and LD in time appropriate periods
#### replace rdm walk with ARMA

cap.v40 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, tp for time period pre/post collapse
for(a in 1:2){
     for(tp in 1:3){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:3){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
for(a in 1:Ni){
     for(t in 1:max(p,q)){
     N[t,a] ~ dnorm(8, 1/9)
     mu[t,a] ~ dunif(0,10)
   }
}


### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
  for(a in 1:Ni){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
  }

for(t in 1:18){
   u[t] ~ dnorm(0, 1/16) # this is based on a SD of 2 so Var = 4- then, I doubled it for no good reason other than to see what would happen
}
   
# AR
#for (a in 1:Ni){
  for (i in 1:p) {
  #  eta[i,a] ~ dnorm(0, 100^-2)
  eta[i] ~ dnorm(0, 100^-2)
 }
#}

# MA
   for (i in 1:q) {
    zeta[i] ~ dnorm(0, 100^-2)
   }


#mu:N2-3; 1985:1990 - would prefer to use this code but subset error 
# for (a in 1:Ni){
#    for (t in 1:18) { #18
#       mu[t,a] <- ifelse(t<12,
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]),
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1])
#       # this is OK for now but need CO by age and this may not be appro for age 2
#    }
# }

for (a in 1:Ni){
   for (t in (max(p,q)+1):11) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + ma_mean[t,a]         #+ ar_mean[t,a]
      ar_mean[t,a] <- inprod(eta, N[(t-p):(t-1), a])
       ma_mean[t,a] <- inprod(zeta, eps[(t-q):(t-1), a])
   }
}

for (a in 1:Ni){
   for (t in 12:18) { #18
      mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1] + ma_mean[t,a]
      #+ ar_mean[t,a]
      ar_mean[t,a] <- inprod(eta, N[(t-p):(t-1), a])
      ma_mean[t,a] <- inprod(zeta, eps[(t-q):(t-1), a])
   }
}


# N2-3; 1985:1990
   for (t in (max(p,q)+1):18) { #18
      N[t,1] ~ dnorm(mu[t,1], 
      ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
   }

   for (t in (max(p,q)+1):18-1) { #18
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2])) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2]))
   } #N4



# ######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
      mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1] + ma_mean[t,1]
      N[t,1] ~ dnorm(mu[t,1], tau.proc[2, 1])
      ma_mean[t,1] <- inprod(zeta, eps[(t-q):(t-1), 1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
  }
}


#N3 & N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2], tau.proc[2, 2]) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3], tau.proc[2, 2]) #N4
   }


# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     #eps[t,a] <- N[t,a] - matI[t,a]
     #eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids - could set a p
for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     posa[1, a] ~ dnorm(0, 1/10)
}

for (a in 1:Ni){
     for (t in 2:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     posa[t,a] <- osa[t,a]/sd(osa[,a])
     }
}
 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}

     Dmape.rep <- sum(Dssm.rep)
}'




# 43: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids
#### add TI, CO, and LD in time appropriate periods
#### replace rdm walk with AR but no alpha and for all N2 - tried N3 and N4 but get lm.fit problems

cap.v41 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, tp for time period pre/post collapse
for(a in 1:2){
     for(tp in 1:3){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:3){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
   N[1,1] ~ dnorm(8.5, 1/9)
   N[2,1] ~ dnorm(8.9, 1/9)
   #N[3,1] ~ dnorm(6, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[2,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)
   N[2,3] ~ dnorm(6, 1/9)
   mu[1,1] ~ dunif(0,10)
   mu[2,1] ~ dunif(0,10)
   mu[1,2] ~ dunif(0,10)
   mu[2,2] ~ dunif(0,10)
   mu[1,3] ~ dunif(0,10)
   mu[2,3] ~ dunif(0,10)
   alpha[1] ~ dnorm(0, 100^-2)  
   alpha[2] ~ dnorm(0, 100^-2)  
   alpha[3] ~ dnorm(0, 100^-2)  

   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   for(a in 1:Ni){
#   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}
   
# AR
#for (a in 1:Ni){
  for (i in 1:p) {
  #  eta[i,a] ~ dnorm(0, 100^-2)
   eta[i] ~ dnorm(0, 100^-2)
 }
#}


#mu:N2-3; 1985:1990 - would prefer to use this code but subset error 
# for (a in 1:Ni){
#    for (t in 1:18) { #18
#       mu[t,a] <- ifelse(t<12,
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]),
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1])
#       # this is OK for now but need CO by age and this may not be appro for age 2
#    }
# }


   for (t in (p+1):11) { #18
      #mu[t,1] <- alpha[1] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + ar_mean[t,1]
      mu[t,1] <- gamma[1]*TI[t]*(1-TI[t]/delta[1]) + ar_mean[t,1]
       ar_mean[t,1] <- inprod(eta, N[(t-p):(t-1), 1])
   }


for (a in 2:Ni){
   for (t in (p+1):11) { #18
      #mu[t,a] <-  alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + ar_mean[t,a]
       mu[t,a] <-  alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) 
      #mu[t,a] <-  gamma[a]*TI[t]*(1-TI[t]/delta[a]) + ar_mean[t,a]
#     #  ar_mean[t,a] <- inprod(eta, N[(t-p):(t-1), a])
   }
}


   for (t in 12:18) { #18
#      mu[t,1] <- alpha[1] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1] + ar_mean[t,1]
            mu[t,1] <- gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1] + ar_mean[t,1]
       ar_mean[t,1] <- inprod(eta, N[(t-p):(t-1), 1])
   }


for (a in 2:Ni){
   for (t in 12:18) { #18
      mu[t,a] <- alpha[a] +gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1] 
      #+ ar_mean[t,a]
       #ar_mean[t,a] <- inprod(eta, N[(t-p):(t-1), a])
   }
}

# N2-3; 1985:1990
   for (t in (p+1):18) { #18
      N[t,1] ~ dnorm(mu[t,1]+0.01, 
      ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
   }

   for (t in (p+1):18-1) { #18
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2]+0.01, ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2])) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3]+0.01, ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2]))
   } #N4



# ######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
#      mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1] + ar_mean[t,1]
    
       mu[t,1] <- beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1] + ar_mean[t,1]

       ar_mean[t,1] <- inprod(eta, N[(t-p):(t-1), 1])
      N[t,1] ~ dnorm(mu[t,1]+0.01, tau.proc[2, 1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1] 
    #+ ar_mean[t,a]
     #ar_mean[t,a] <- inprod(eta, N[(t-p):(t-1), a])
  }
}


#N3 & N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(log(exp(N[t,1])*(1-m[t]))*mu[t,2]+0.01, tau.proc[2, 2]) #N3
      N[t+1,3] ~ dnorm(log(exp(N[t,2])*(1-0.95))*mu[t,3]+0.01, tau.proc[2, 2]) #N4
   }


# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     #eps[t,a] <- N[t,a] - matI[t,a]
     #eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids - could set a p
for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     posa[1, a] ~ dnorm(0, 1/10)
}

for (a in 1:Ni){
     for (t in 2:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     posa[t,a] <- osa[t,a]/sd(osa[,a])
     }
}
 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}


for (t in 8:n.occasions){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
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


# 44: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids
#### add rdm walk in time appropriate periods
##### NOTE THAT POSA NOT WORKING

cap.v42 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, tp for time period pre/post collapse
for(a in 1:2){
     for(tp in 1:3){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:3){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
  N[1,1] ~ dnorm(8.5, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)
   # mu[1,1] ~ dnorm(0, 10)
   # mu[1,2] ~ dnorm(0, 10)
   # mu[1,3] ~ dnorm(0, 10)
   mu[1,1] ~ dunif(0,10)
   mu[1,2] ~ dunif(0,10)
   mu[1,3] ~ dunif(0,10)
   
   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   for(a in 1:Ni){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}

  u ~ dnorm(0, 1/16) # this is based on a SD of 2 so Var = 4- then, I doubled it for no good reason other than to see what would happen


#mu:N2-3; 1985:1990 - would prefer to use this code but subset error 
# for (a in 1:Ni){
#    for (t in 1:18) { #18
#       mu[t,a] <- ifelse(t<12,
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]),
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1])
#       # this is OK for now but need CO by age and this may not be appro for age 2
#    }
# }


for (a in 1:Ni){
   for (t in 2:11) { #18
      mu[t,a] <- N[t-1, a] + u 
      #+ gamma[a]*TI[t]*(1-TI[t]/delta[a])
   }
}

for (a in 1:Ni){
   for (t in 12:18) { #18
      mu[t,a] <- N[t-1, a] + u
      #+ u + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
   }
}


# N2-3; 1985:1990
   for (t in 2:18) { #18
      N[t,1] ~ dnorm(mu[t,1], 
      ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
   }

   for (t in 2:18-1) { #18
      N[t+1,2] ~ dnorm(mu[t,2], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2])) #N3
      N[t+1,3] ~ dnorm(mu[t,3], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2]))
   } #N4


# ######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
      #mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
      mu[t,1] <- N[t-1, 1] + u 
      N[t,1] ~ dnorm(mu[t,1], tau.proc[3, 1])
     # ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
     # N[t,1] ~ dnorm(mu[t,1], tau.proc[2, 1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    #mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
    mu[t,a] <- N[t-1, a] + u
   }
}


#N3 & N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(mu[t,2], tau.proc[3, 2]) #N3
      N[t+1,3] ~ dnorm(mu[t,3], tau.proc[3, 2]) #N4
   }


# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     # eps[t,a] <- N[t,a] - matI[t,a]
     # eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids - could set a p
for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     posa[1, a] ~ dnorm(0, 1/10)
}

for (a in 1:Ni){
     for (t in 2:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     posa[t,a] <- osa[t,a]
     #/sd(osa[,a])
     }
}
 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}

#    # Assessing the fit of the state-space model
#    ## 1. Compute fit statistics for observed data.
#    ### 1.1 Discrepancy meansure: mean absolute error
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}


for (t in 8:n.occasions){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
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


# 45: Extend model time----
## ln scale: N2-N4 + forecast for each age
### Added in the capelin data from 1985-1998 but in matrix form
#### AS 31 but with I and all the diagnostics - this adds parameters to the loops and cleaned up much excess code Ni = N age[i]
##### eps added for N2-N4 and osa
#####  minimize code and reduce tau.proc to 2 values
##### extend tau.proc and tau.obs to two time periods (tp) corresponding to the pre/post collapse
#### added osa and posa resids
#### add rdm walk in time appropriate periods
#### add TI, CO, and LD in time appropriate periods

cap.v43 = '
 model {
#PRIORS
###### Need to check that these are reasonable
 ## Prior for sd of process - N2-N3[t] uninformative
 # a is process variance for age group, tp for time period pre/post collapse
for(a in 1:2){
     for(tp in 1:3){
        sigma.proc[tp, a] ~ dunif(0.01, 20)
        sigma2.proc[tp,a] <- pow(sigma.proc[tp, a], 2)
        tau.proc[tp, a] <- pow(sigma.proc[tp, a], -2)
     }
}


 ## Prior for sd of observation - I2-I4[t] - uninformative
 for(tp in 1:3){
   sigma.obs[tp] ~ dunif(0.01, 20)
   sigma2.obs[tp] <- pow(sigma.obs[tp], 2)
   tau.obs[tp] <- pow(sigma.obs[tp], -2)
 }


### Priors for Initial values for N2-N4[t] informative - based on actual values
  N[1,1] ~ dnorm(8.5, 1/9)
   N[1,2] ~ dnorm(8.9, 1/9)
   N[1,3] ~ dnorm(6, 1/9)
   # mu[1,1] ~ dnorm(0, 10)
   # mu[1,2] ~ dnorm(0, 10)
   # mu[1,3] ~ dnorm(0, 10)
   mu[1,1] ~ dunif(0,10)
   mu[1,2] ~ dunif(0,10)
   mu[1,3] ~ dunif(0,10)
   
   ### Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
# for (t in 2:4){
 #   N2[t] ~ dnorm(8.5, 1/9)
  #  }


# LIKELIHOODS
 ## State process
   ### N2
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019 - mostly uninformative but see TI-width - not sure here
   for(a in 1:Ni){
   alpha[a] ~ dnorm(0, 100^-2)      # int
   beta[a] ~ dnorm(0, 100^-2)       # larval abund
   gamma[a] ~ dunif(0, 100)         # tices-max rate of increase
   delta[a] ~ dgamma(11.5, 5.7)     # tice-width
   epsilon[a] ~ dnorm(0, 100^-2)   # condition # for CO
}

  u ~ dnorm(0, 1/16) # this is based on a SD of 2 so Var = 4- then, I doubled it for no good reason other than to see what would happen


#mu:N2-3; 1985:1990 - would prefer to use this code but subset error 
# for (a in 1:Ni){
#    for (t in 1:18) { #18
#       mu[t,a] <- ifelse(t<12,
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]),
#       alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1])
#       # this is OK for now but need CO by age and this may not be appro for age 2
#    }
# }


for (a in 1:Ni){
   for (t in 2:11) { #18
      mu[t,a] <- N[t-1, a] + gamma[a]*TI[t]*(1-TI[t]/delta[a])
   }
}

for (a in 1:Ni){
   for (t in 12:18) { #18
      mu[t,a] <- N[t-1, a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
   }
}


# N2-3; 1985:1990
   for (t in 2:18) { #18
      N[t,1] ~ dnorm(mu[t,1], 
      ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
   }

   for (t in 2:18-1) { #18
      N[t+1,2] ~ dnorm(mu[t,2], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2])) #N3
      N[t+1,3] ~ dnorm(mu[t,3], ifelse( t<=6, tau.proc[1, 2], tau.proc[2, 2]))
   } #N4


# ######################################

# mu and N2: 2003-present
for (t in 19:n.occasions) { #19
      #mu[t,1] <- alpha[1] + beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
      mu[t,1] <- N[t-1, 1] +  beta[1]*LD[t-2] + gamma[1]*TI[t]*(1-TI[t]/delta[1]) + epsilon[1]*CO[t-1]
      N[t,1] ~ dnorm(mu[t,1], tau.proc[3, 1])
     # ifelse(t<=6, tau.proc[1,1], tau.proc[2,1]))
     # N[t,1] ~ dnorm(mu[t,1], tau.proc[2, 1])
   }

# mu: N3 and N4; 2003-present
##
for (a in 2:Ni){
  for (t in 19:n.occasions) { #19
    #mu[t,a] <- alpha[a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
    mu[t,a] <- N[t-1, a] + gamma[a]*TI[t]*(1-TI[t]/delta[a]) + epsilon[a]*CO[t-1]
   }
}


#N3 & N4; 2003-present
   for (t in 19:n.occasions-1) { #19
      N[t+1,2] ~ dnorm(mu[t,2], tau.proc[3, 2]) #N3
      N[t+1,3] ~ dnorm(mu[t,3], tau.proc[3, 2]) #N4
   }


# process error
for (a in 1:Ni){
     for (t in 1:n.occasions){
     # eps[t,a] <- mu[t,a] - N[t,a] # probably wrong N[t] - mu
     # eps[t,a] <- N[t,a] - matI[t,a]
     # eps[t,a] <- N[t,a] - mean(N[,a])
     eps[t,a] <- N[t,a] - mu[t,a]
     }
}

# one step ahead resids (osa) - set priors because you cannot calculate these for the first couple of values (I think because the if t = 1, then t-1 is zero which is nonsensical and if t = 2, then there is only one value for the SD which at least in R, results in a NA)

for (a in 1:Ni){
     osa[1, a] ~ dnorm(0, 1/10)
     osa[2, a] ~ dnorm(0, 1/10)
     osa_sd[1, a] ~ dnorm(0, 1/10)
     osa_sd[2, a] ~ dnorm(0, 1/10)
    # posa[1, a] ~ dnorm(0, 1/10)
    # posa[2, a] ~ dnorm(0, 1/10)
}

# calculate the osa and the Pearson osa (posa).  I cannot get this to work in JAGS but it seems easy enough in R.  So proceed 
for (a in 1:Ni){
     for (t in 3:n.occasions){
     osa_mean[t,a] <- mean(N[1:(t-1), a])
     osa[t,a] <- N[t,a] - osa_mean[t,a]
     osa_sd[t,a] <- sd(osa[1:(t-1), a])
     #osa_sd[t,a] <- sd(1/osa[1:(t-1), a])
     # posa[t,a] <- osa[t,a]/sd(osa[,a]) # this works - rest do not
      #posa[t,a] <- osa[t,a] %*%(osa_sd[t,a])
     }
}

# OK - so you cannot divide a matrix by a matrix but you can multiple a matrix by its inverse. However, I got an error (not positive definite and in JAGS, inverse only works on symetrical matrices).  So multiplying the below works but this is not what I want.  
# for (a in 1:Ni){
 #    for (t in 1:n.occasions){
  #    posa[t,1] <- osa[t,1] %*% pow(osa_sd[t,1], 1) # works but not right answer
      #posa[t,1] <- osa[t,1] %*% pow(osa_sd[t,1], -1) # results in invalid parent values
     #posa[t,1] <- osa[t,1] %*% pow(osa_sd[t,1], 1) # as above but with osa_sd as 1/osa_sd - again, invalid parent values
     #posa[t,1] <- pow(osa_sd[t,1], -1) %*% osa[t,1]   # results in invalid parent values

   #  }
# }


 ## Observation
   ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts - eliminateed this for now
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.


for (a in 1:Ni){
   for (t in 1:7) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[1])       # sampled observation
   }
}


for (a in 1:Ni){
   for (t in 8:n.occasions) {
      matI[t,a] ~ dnorm(N[t,a], tau.obs[2])       # sampled observation
   }
}


for (t in 1:7) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[1])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}


for (t in 8:n.occasions) {
     # I[t] ~ dnorm(log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3])), tau.obs[2])
     I[t] <- log(exp(matI[t,1]) + exp(matI[t,2]) + exp(matI[t,3]))
}

#    # Assessing the fit of the state-space model
#    ## 1. Compute fit statistics for observed data.
#    ### 1.1 Discrepancy meansure: mean absolute error
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
for (t in 1:7){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[1])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[1])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[1])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[1])
      Dssm.rep[t] <- abs((I.rep[t] - I.exp[t])/I.rep[t])
}


for (t in 8:n.occasions){
      #    y2.rep[t] ~ dnorm(N2[t], tau.obs)
      #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
      I2.rep[t] ~ dnorm(N[t,1], tau.obs[2])
      I3.rep[t] ~ dnorm(N[t,2], tau.obs[2])
      I4.rep[t] ~ dnorm(N[t,3], tau.obs[2])
      I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.obs[2])
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
