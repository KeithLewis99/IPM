# the purpose of this file is simply to have a repository for the large number of combinations of model inputs (model name and variables that will be used in outputs) and parameter values that may be of interest.


# model----

if (b==1){ # model with separate parms for each age
     parms = parms1
     tC = cap.v7
     tC.txt = "cap.v7"
     vars_vAR <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs")
     vars_Nyear <- c("N2[10]", "N3[10]", "N4[10]")
     vars_N2 <- c("mu2[10]","alpha2", "beta2",  "gamma2", "delta2")
     vars_N3 <- c("mu3[10]", "alpha3", "gamma3", "delta3", "epsilon3")
     vars_N4 <- c("mu4[10]", "alpha4", "gamma4", "delta4", "epsilon4")
     
} else if (b==2) { # model separate parms for N2 and N3:N4  - EXT
     parms = parms2
     # tC = cap.v8
     # tC.txt = "cap.v8"
     tC = cap.v22
     tC.txt = "cap.v22"
     vars_vAR <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs")
     vars_Nyear <- c("N2[10]", "N3[10]", "N4[10]")
     vars_N2 <- c("mu2[10]","alpha2", "beta2",  "gamma2", "delta2")
     vars_N3 <- c("mu3[10]", "alpha3", "gamma3", "delta3", "epsilon3", "mu4[10]")
     vars_N4 <- c(NA)
     
} else if (b==3) { # demographic model - no forecast for N4
     parms = parms3
     tC = cap.v9
     tC.txt = "cap.v9"
     vars_vAR <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs")
     vars_Nyear <- c("N2[10]", "N3[10]", "N4[10]")
     vars_N2 <- c("mu2[10]","alpha2", "beta2",  "gamma2", "delta2")
     vars_N3 <- c(NA)
     vars_N4 <- c(NA)
     
} else if (b==4) { # model separate parms for N2 and N3:N4
     parms = parms4
     tC = cap.v10
     tC.txt = "cap.v10"
     vars_vAR <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs")
     vars_Nyear <- c("N2[10]", "N3[10]", "N4[10]")
     vars_N2 <- c("u")
     vars_N3 <- c("mu3[10]", "alpha3", "gamma3", "delta3", "epsilon3", "mu4[10]")
     vars_N4 <- c(NA)
}


# Parameters -----
## cap.v7: Model ln scale: N2-N4 + forecast for each age; 1999-2000
### b == 1
parms1 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
            "N2",  "N3", "N4",
            "mu2", "alpha2", "beta2",  "gamma2", "delta2",
            "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
            "mu4", "alpha4", "gamma4", "delta4", "epsilon4",
            "Dssm.obs", "Dmape.obs",  "Tturn.obs", 
            "Dssm.rep", "Dmape.rep",  "Tturn.rep",
            "I2.rep", "I3.rep", "I4.rep", "I.rep"
) 

# cap.v8: Model ln scale: N2-N4 + forecast for each age
## this has common parameters for N3:N4 for the forecast model but not N2
### b == 2
parms2 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
            "N2",  "N3", "N4", "eta", "zeta", "eps", "ar_mean", "ma_mean",
            "mu2", "alpha2", "beta2",  "gamma2", "delta2",
            "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
            "mu4",
            "Dssm.obs", "Dmape.obs",  "Tturn.obs", 
            "Dssm.rep", "Dmape.rep",  "Tturn.rep",
            "I2.rep", "I3.rep", "I4.rep", "I.rep"
) 


# cap.v9: This model may have been deleted but the idea here is to have no forecast model except for age 2 and just have demographic values generate age 3 and age 4.
### b == 3 demographic model - no forecast for N4
parms3 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
            "N2",  "N3", "N4",
            "mu2", "alpha2", "beta2",  "gamma2", "delta2",
            "Dssm.obs", "Dmape.obs",  "Tturn.obs", 
            "Dssm.rep", "Dmape.rep",  "Tturn.rep",
            "I2.rep", "I3.rep", "I4.rep", "I.rep"
) 

#  cap.v10: Model ln scale: N2-N4 but with a rdm walk----
## no forecast on N2, only rdm walk.  Forecast on N3 but not N4.
### b == 4
parms4 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
            "N2",  "N3", "N4", "u",
            "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
            "mu4", 
            "Dssm.obs", "Dmape.obs",  "Tturn.obs", 
            "Dssm.rep", "Dmape.rep",  "Tturn.rep",
            "I2.rep", "I3.rep", "I4.rep", "I.rep"
) 


# # cap.v10: Model ln scale: N2-N4 + forecast for each age ----
# # and with a rdm walk (previous one does not have forecast)
### b == 5
parms5 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
            "N2",  "N3", "N4",
            "mu2", "alpha2", "beta2",  "gamma2", "delta2", "u",
            "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
            "mu4", 
            "Dssm.obs", "Dmape.obs",  "Tturn.obs", 
            "Dssm.rep", "Dmape.rep",  "Tturn.rep",
            "I2.rep", "I3.rep", "I4.rep", "I.rep"
) 


#  , "pe3", "pe2",
# "I.exp", "I2.rep", "I3.rep", "I4.rep", "I.rep","I2", "I3", "I4", "I",
# "Tt1.obs", "Tt2.obs", "Tt3.obs", "Tt1.rep", "Tt2.rep", "Tt3.rep",