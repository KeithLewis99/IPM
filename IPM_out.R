
# required packages
library(rjags)
library(R2jags)
library(ggplot2)



# Start----
rm(list=ls())

# Source files
source("IPM_dat.R")
source("IPM_fun.R")
source("IPM_mod.R")



# Model ln scale: tice N3 mortality and INdex SE and split N2----
# try to fix the priors ito variance

# JAGS settings
parms <- c("Sld", "N2",  "N3", "tau.proc", "tau.obs", "tau.LD", "tau.ind", "I2", "I3", "y2", "y3", "alpha", "beta", "gamma")

# MCMC settings
#ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
ni <- 200000; nt <- 30; nb <- 30000; nc <- 3

# run model
ssm26 <- jags(jags.data, parameters=parms, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(cap.v7))
ssm26

# create ouput
out <- ssm26$BUGSoutput 

ssm26_dic <- out$DIC

## extract raw values from chains
raw <- ls_out(out)
str(raw)

#extract medians, credible intervals, and prediction intervals
calc <- ls_med(raw)
df_calc <- do.call(rbind, calc) # this doesn't work
write(df_calc, "out2.csv")
str(calc)
#cbind(N2_med, N3_med, N2_med+N3_med)

## figures
# N2: observation median v process median
plot(calc$I2_med, calc$N2_med)
# N3: observation median v process median
plot(calc$I3_med, calc$N3_med)
# Observation median over time
plot(seq(1999:2023), calc$I_med)

# observation median v real data - relation is perfect - is this OK?
plot(calc$I2_med, jd$I2)
plot(calc$I3_med, jd$I3)
plot(calc$N2_med, jd$I2)
plot(calc$N3_med, jd$I3)


## IPM plot----
year <- 1999:2021
ly <- length(year)
forecast <- 2022:2023
lf <- length(forecast)

tmp_plot <- ipm_plot(calc)
tmp_plot
ggsave("tmp_plot1.pdf")

     # Other plots
     tp = out$sims.list$tau.proc
     tp_med = apply(tp,2,'median') # median values of y_pred
     tp_ci = apply(tp,2,'quantile', c(0.1, 0.9)) # median values of y_pred

# Diagnostics----
#Source
source('C:/Users/lewiske/Documents/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('C:/Users/lewiske/Documents/R/zuur_rcode/HighstatLibV7.R')
library(lattice)

filepath_gen <- "biomass_cond_ag1_2_DIC_R3" 
filepath <- paste0(filepath_gen, "/recruitment_1")


print(out, intervals=c(0.025, 0.975), digits = 3)
out$mean


     #overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
     mean(out$sims.list$FitNew > out$sims.list$Fit)
     # mean = 0.546333

## Mixing ----     
# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
# vars for forecast model     
     source('C:/Users/lewiske/Documents/R/zuur_rcode/MCMCSupportHighstatV2.R')
var1 <- c('alpha', 'beta')
#vars1 <- c('alpha', 'beta', 'gamma', "delta", "epsilon", "sigma")

MyBUGSChains(out, var1)
mix1 <- MyBUGSChains(out, var1)

#ggsave(MyBUGSChains(out, vars1), filename = paste0("Bayesian/", filepath, "/chains-forecast.pdf"), width=10, height=8, units="in")

# vars for state space and demographic vars
vars2 <- c("Sld")
MyBUGSChains(out, vars2)
mix2 <- MyBUGSChains(out, vars2)
#ggsave(MyBUGSChains(out, vars2), filename = paste0("Bayesian/", filepath, "/chains-demographic.pdf"), width=10, height=8, units="in")


# vars for variances
vars3 <- c("tau.proc", "tau.obs", "tau.LD", "tau.ind")
MyBUGSChains(out, vars3)
mix3 <- MyBUGSChains(out, vars3)
#ggsave(MyBUGSChains(out, vars3), filename = paste0("Bayesian/", filepath, "/chains-variance.pdf"), width=10, height=8, units="in")


##autocorrelation----
MyBUGSACF(out, var1)
autocorr1 <- MyBUGSACF(out, var1)
#ggsave(MyBUGSACF(out, vars1), filename = paste0("Bayesian/", filepath, "/auto_corr-forecast.pdf"), width=10, height=8, units="in")

MyBUGSACF(out, vars2)
autocorr2 <- MyBUGSACF(out, vars2)
#ggsave(MyBUGSACF(out, vars2), filename = paste0("Bayesian/", filepath, "/auto_corr-demographic.pdf"), width=10, height=8, units="in")

MyBUGSACF(out, vars3)
autocorr3 <- MyBUGSACF(out, vars3)
#ggsave(MyBUGSACF(out, vars3), filename = paste0("Bayesian/", filepath, "/auto_corr-autocorrelation.pdf"), width=10, height=8, units="in")


## Model Validation ----
#(see Zuuer et al. 2013 for options for calculating Pearson residuals)
# # Residual diagnostics
# E1 <- out$mean$PRes # Pearson resids
# F1 <- out$mean$expY # Expected values
# N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
# D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD
# 
# pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
# par(mfrow = c(2,2), mar = c(5,5,2,2))
# plot(x=F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals")
# abline(h = 0, lty = 2)
# # see notes in Mortality model
# plot(y = N2, x = F1, xlab = "Fitted values", ylab = "Observed data") # should follow the line
# abline(coef = c(0,1), lty = 2)
# par(mfrow = c(1,1))
# dev.off()
# 
# # Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
# pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
# par(mfrow = c(2,2), mar = c(5,5,2,2))
# MyVar <- c("tice.std", "meandCond_lag.std")
# df_diag <- as.data.frame(model_data)
# df_diag <- cbind(df_diag, E1)
# plot(df_diag$ST, df_diag$EI, xlab = "Surface Tows", ylab = "Pearson resids")
# plot(df_diag$PS, df_diag$EI, xlab = "Pseudo-calanus", ylab = "Pearson resids")
# par(mfrow = c(1,1))
# dev.off()
