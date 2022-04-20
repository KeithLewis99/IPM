



# Model ln scale: tice N3 mortality and INdex SE and split N2----
# try to fix the priors ito variance

cap.v7 = '
 model {
    sigma.proc ~ dunif(0.01, 20)        # Prior for sd of state process - based on max sd of 10.8
   sigma2.proc <- pow(sigma.proc, 2) 
   tau.proc <- pow(sigma.proc, -2) 
 
   sigma.obs ~ dunif(0.01, 20)       # Prior for sd of observation process 
   sigma2.obs <- pow(sigma.obs, 2) 
   tau.obs <- pow(sigma.obs, -2) 
 
   Sld ~ dunif(0, 1)    # Prior for mean larval density - based on range of values plus extra
   sigma.LD ~ dunif(0.01, 10)        # Prior for sd of larval density
   sigma2.LD <- pow(sigma.LD, 2) 
   tau.LD <- pow(sigma.LD, -2) 
 
  # NEEDED?? # S2 ~ dunif(0, 1)    #Prior for mean larval density - based on range of values plus extra - 
 
 
 # Likelihood 
 ## State process 
 ### Initial values for N2 and N3
   N2[1] ~ dnorm(8.5, 1/9)    # Prior for initial population size - based on N2[1] above   
   N3[1] ~ dnorm(8.9, 1/9)    # Prior for initial population size - based on N1 above 
 
 # these values are because 2000-2003 are NA for LD and therefore, N2 cant be calculated in this model formulation. VAlues are from the data, variance is made up
   for (t in 2:4){
   N2[t] ~ dnorm(8.5, 1/9)
   }
   

# trying to find a way to produce N3 from N2 (based on immature N2 and survival to age 3)
   for (t in 1:(n.occasions-1)){ 
   Nt2[t] ~ dnorm(N2[t], tau.proc)  # so this produces total N2 from below but with process error
   Ni2[t] <- Nt2[t]*(1-m[t])        # immature N2
   Nm2[t] <- Nt2[t]*(m[t])          # mature N2

# trying to kill off immature N2 so that not all become N3
      sur[t] ~ dunif(0, 1)  # prior for survival rate
      n3[t] <- Ni2[t]*sur[t] # the survival rate of N2
    N3[t+1] ~ dnorm(n3[t], tau.proc) # the number of N3 with process error

   # this produces the pool of N2 that could become N3
   #  N3[t+1] ~ dnorm(X2[t], tau.proc) 
    # N2[t] <- N2[t]*(1-m[t])
     # X2[t] <- alpha + gamma*TI[t]*(1-TI[t]/delta) # how to cause mortality in N3???
#    # N3[t+1] ~ dnorm(X2[t], tau.proc) 
     #N3[t+1] ~ dlnorm(N2[t]*(1-m[t]), tau.proc) 
     #N3[t+1] ~ dlnorm(N2[t] + (1-m[t]) + Sld, tau.proc) 
     
     #sigma.sur ~ dunif(0, 10)       # Prior for sd of survival
      #sigma2.sur <- pow(sigma.sur, 2) 
      #tau.sur <- pow(sigma.sur, -2) 
   } 

 ## Observation process - aim to make this age disaggregated in next phase
 # see Schaub and Kerry pg 263 - this is for estimated indices instead of counts

   # y[t] is the "true" index that is sampled by I[t] - tau.obs is the sampling error of the index
   # N[t] is the "true" population (process) where the tau is the additional residual error - i may have tehse confused.
   sigma.ind ~ dunif(0.01, 1000)       # Prior for sd of index
   sigma2.ind <- pow(sigma.ind, 2) 
   tau.ind <- pow(sigma.ind, -2) 

   for (t in 1:n.occasions) { 
   # the estiamted index
   I2[t] ~ dnorm(y2[t], tau.ind) 
     I3[t] ~ dnorm(y3[t], tau.ind)

   y2[t] ~ dnorm(N2[t], tau.obs) 
     y3[t] ~ dnorm(N3[t], tau.obs)
   }

 ## larval density
   # From Murphy the equation relating R = LD*S is R = 0.40x + 2.80
   # priors from Lewis et al. 2020
   # this is a reasonable approach given the correlation between I2 and I
  alpha ~ dnorm(0, 100^-2) # int
  beta ~ dnorm(0, 100^-2) # larval abund
  gamma ~ dunif(0, 100)  #tices-max rate of increase
  delta ~ dgamma(11.5, 5.7) #tice-width
  #gamma ~ dnorm(0, 100^-2) # condition # for CO
  #epsilon ~ dnorm(0, 100^-2) # condition # for CO
  sigma ~ dunif(0, 100) 
   for (t in 5:n.occasions) {
      #mut[t] <- alpha # want this so that I can have a null model but it wont run?????
      #mu[t] <- alpha + beta*LD[t-2]
      #mu[t] <- alpha + beta*LD[t-2] + gamma*CO[t-1]
      mu[t] <- alpha + beta*LD[t-2] + gamma*TI[t]*(1-TI[t]/delta)
      #mu[t] <- alpha + beta*LD[t-2] + gamma*TI[t]*(1-TI[t]/delta) + epsilon*CO[t-1]
      N2[t] ~ dnorm(mu[t], sigma^-2)
   }
}'