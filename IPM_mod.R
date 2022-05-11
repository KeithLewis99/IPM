# Model ln scale: tice N3 mortality and INdex SE and split N2----
# try to fix the priors ito variance

cap.v7 = '
 model {
 #PRIORS
 # Prior for sd of process - N2[t] based
    sigma.proc2 ~ dunif(0.01, 20)     
   sigma2.proc2 <- pow(sigma.proc2, 2) 
   tau.proc2 <- pow(sigma.proc2, -2) 

 # Prior for sd of process - N3[t] based
    sigma.proc3 ~ dunif(0.01, 20)     
   sigma2.proc3 <- pow(sigma.proc3, 2) 
   tau.proc3 <- pow(sigma.proc3, -2) 

 # Prior for sd of process - N4[t] based
    sigma.proc4 ~ dunif(0.01, 20)     
   sigma2.proc4 <- pow(sigma.proc4, 2) 
   tau.proc4 <- pow(sigma.proc4, -2) 

 # Prior for sd of sampled index - I[t]
   # sigma.ind ~ dunif(0.01, 1000)   
   # sigma2.ind <- pow(sigma.ind, 2) 
   # tau.ind <- pow(sigma.ind, -2) 
 
 # Prior for sd of observation - I[t]
   sigma.obs ~ dunif(0.01, 20)   
   sigma2.obs <- pow(sigma.obs, 2) 
   tau.obs <- pow(sigma.obs, -2) 
 

   
 # Prior for mean larval density - based on range of values plus extra
   #sigma.LD ~ dunif(0.01, 10)       
   # sigma.LD ~ dunif(0.01, 10)       
   # sigma2.LD <- pow(sigma.LD, 2) 
   # tau.LD <- pow(sigma.LD, -2) 
   #LD_lat ~ dnorm (LD, LD.tau) # latent larval density based on observed larval density and precision of the LD
   
 ### Initial values for N2 and N3
   N2[1] ~ dnorm(8.5, 1/9)    # Prior for initial population size - based on N2[1] above  
   N3[1] ~ dnorm(8.9, 1/9)    # Prior for initial population size - based on N1 above 
    N4[1] ~ dnorm(6, 1/9)    # Prior for initial population size - based on N1 above 
 # Values for N2[2-4]: required because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
   for (t in 2:4){
   N2[t] ~ dnorm(8.5, 1/9)
   }

 # LIKELIHOODS
 ## State process 
# make N2 from below and make a fraction immature
   # for (t in 1:(n.occasions)){ 
   # Nt2[t] ~ dnorm(N2[t], tau.proc)  # so this produces total N2 from below but with process error
   # Ni2[t] <- Nt2[t]*(1-m[t])        # immature N2
   # }

# Priors for killing N3
   #alpha1 ~ dnorm(0, 100^-2)                 # int
   alpha1 ~ dunif(0, 1)                 # int
   gamma1 ~ dunif(0, 100)                    #tices-max rate of increase
   # gamma1 ~ dunif(0, 3.65)                 # was uniform
   delta1 ~ dgamma(11.5, 5.7)               #tice-width
   epsilon1 ~ dnorm(0, 100^-2)              # condition # for CO
   mu3[1] ~ dunif(0,1)                    # need this for a mu3[1]
   
   
# calculate a survival for N2 -> N3; mu3[2:24]
   for(t in 2:n.occasions){
      mu3[t] <- alpha1 + gamma1*TI[t]*(1-TI[t]/delta1) + epsilon1*CO[t-1]
      #mu3[t] <- alpha1
   }

# calculate the N3
   for(t in 2:n.occasions-1){
         #N3[t+1] ~ dnorm(Ni2[t]*S[t], tau.proc)
         #N3[t+1] ~ dnorm(N2[t]*(1-m[t])*S[t], tau.proc)
         N3[t+1]  ~ dnorm(log(exp(N2[t])*(1-m[t]))+mu3[t], tau.proc3)
         #N3[t+1] <- N2[t]*(1-m[t])*exp(mu3[t]) + pe3[t+1]
         #pe3[t+1] ~ dnorm(0, tau.proc3)
      }

# Priors for killing N4
   alpha2 ~ dunif(0, 1)                 # int
   gamma2 ~ dunif(0, 100)                    #tices-max rate of increase
   delta2 ~ dgamma(11.5, 5.7)               #tice-width
   epsilon2 ~ dnorm(0, 100^-2)              # condition # for CO
   mu4[1] ~ dunif(0,1)                    # need this for a mu3[1]
   
   
# calculate a survival for N3 -> N4; mu4[2:24]
   for(t in 2:n.occasions){
      mu4[t] <- alpha2 + gamma2*TI[t]*(1-TI[t]/delta2) + epsilon2*CO[t-1]
    }

# calculate the N4
   for(t in 2:n.occasions-1){
         N4[t+1]  ~ dnorm(log(exp(N3[t])*(1-0.95))*mu4[t], tau.proc4)
      }


 ## Observation
 ### see Schaub and Kerry pg 263 - this is for estimated indices instead of counts
   #### y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   #### N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.

   for (t in 1:n.occasions) {          
      #y2[t] ~ dnorm(N2[t], tau.proc) # the true index
      #y3[t] ~ dnorm(N3[t], tau.proc)

      I2[t] ~ dnorm(N2[t], tau.obs)       # sampled observation
      I3[t] ~ dnorm(N3[t], tau.obs)
      I4[t] ~ dnorm(N4[t], tau.obs)
      
      #I2[t] ~ dnorm(y2[t], SE)       # sampled observation
      #I3[t] ~ dnorm(y3[t], SE)
      I[t] <- log(exp(I2[t]) + exp(I3[t]) + exp(I4[t]))
   }


 ## larval density
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2019
   # this is a reasonable approach given the correlation between I2 and I
  alpha ~ dnorm(0, 100^-2) # int
  beta ~ dnorm(0, 100^-2) # larval abund
  gamma ~ dunif(0, 100)  #tices-max rate of increase
 delta ~ dgamma(11.5, 5.7) #tice-width
  #gamma ~ dnorm(0, 100^-2) # condition # for CO
  #epsilon ~ dnorm(0, 100^-2) # condition # for CO
  #sigma ~ dunif(0, 100) 
 
  
   for (t in 5:n.occasions) {
      #mut[t] <- alpha # want this so that I can have a null model but it wont run?????
      #mu[t] <- alpha + beta*LD[t-2]
      #mu[t] <- alpha + beta*LD[t-2] + gamma*CO[t-1]
      mu2[t] <- alpha + beta*LD[t-2] + gamma*TI[t]*(1-TI[t]/delta)
      #mu[t] <- alpha + beta*LD[t-2] + gamma*TI[t]*(1-TI[t]/delta) + epsilon*CO[t-1]
      #N2[t] ~ dnorm(mu[t], tau.proc) # 

      #N2[t] <- alpha + beta*LD[t-2] + gamma*TI[t]*(1-TI[t]/delta) + pe2[t]
      N2[t] ~ dnorm(mu2[t], tau.proc2)

      #pe2[t] ~ dnorm(0, tau.proc2)

      #N2[t] <- mu[t] # no process error
      
   }     
     
# Assessing the fit of the state-space model
   ## 1. Compute fit statistics for observed data.
   ## 1.1 Discrepancy meansure: mean absolute error
   # for (t in 1:n.occasions) {
   # I.exp[t] <- log(exp(N2[t]) + exp(N3[t]))
   # }
   # 
   # ## 1.2 Test statistic: number of turns or switches - jaggedness
   # for (t in 1:(n.occasions-2)){
   #    Tt1.obs[t] <- step(I[t+2] - I[t+1])
   #    Tt2.obs[t] <- step(I[t+1] - I[t])
   #    # Tt12.obs[t] <- step(I2[t+2] - I2[t+1])
   #    # Tt13.obs[t] <- step(I3[t+2] - I3[t+1])
   #    # Tt22.obs[t] <- step(I2[t+1] - I2[t])
   #    # Tt23.obs[t] <- step(I3[t+1] - I3[t])
   #    # Tt1.obs[t] <- log(exp(Tt12.obs[t]) + exp(Tt13.obs[t]))
   #    # Tt2.obs[t] <- log(exp(Tt22.obs[t]) + exp(Tt23.obs[t]))
   #    Tt3.obs[t] <- equals(Tt1.obs[t] + Tt2.obs[t], 1)
   # }
   # Tturn.obs <- sum(Tt3.obs)
   # 
   # 
   # ## 2.1 Simulated data
    for (t in 1:n.occasions){
   #    y2.rep[t] ~ dnorm(N2[t], tau.obs) 
   #    y3.rep[t] ~ dnorm(N3[t], tau.obs)
       I2.rep[t] ~ dnorm(N2[t], tau.obs)
       I3.rep[t] ~ dnorm(N3[t], tau.obs)
       I4.rep[t] ~ dnorm(N4[t], tau.obs)
   #   I.rep[t] ~ dnorm(log(exp(I2.rep[t]) + exp(I3.rep[t]) + exp(I4.rep[t])), tau.ind)
    }
   # 
   # 
   # ##Test statistic: number of turns or switches - jaggedness
   # for (t in 1:(n.occasions-2)){
   #    Tt1.rep[t] <- step(I.rep[t+2] - I.rep[t+1])
   #    Tt2.rep[t] <- step(I.rep[t+1] - I.rep[t])
   #    Tt3.rep[t] <- equals(Tt1.rep[t] + Tt2.rep[t], 1)
   # }
   # Tturn.rep <- sum(Tt3.rep)

}'
