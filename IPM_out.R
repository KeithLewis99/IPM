
# required packages
library(rjags)
library(R2jags)
library(ggplot2)
library(lattice)

# Start----
#rm(list=ls())


# Source files
source("IPM_dat.R")
source("IPM_fun.R")
source("IPM_mod.R")
source('C:/Users/lewiske/Documents/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('C:/Users/lewiske/Documents/R/zuur_rcode/HighstatLibV7.R')

# JAGS settings ----

parms1 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
           "N2",  "N3", "N4",
           "mu2", "alpha2", "beta2",  "gamma2", "delta2",
           "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
           "mu4", "alpha4", "gamma4", "delta4", "epsilon4",
           "Dssm.obs", "Dmape.obs",  "Tturn.obs", 
           "Dssm.rep", "Dmape.rep",  "Tturn.rep",
           "I2.rep", "I3.rep", "I4.rep", "I.rep"
           ) 

parms2 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
            "N2",  "N3", "N4",
            "mu2", "alpha2", "beta2",  "gamma2", "delta2",
            "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
            "mu4", 
            "Dssm.obs", "Dmape.obs",  "Tturn.obs", 
            "Dssm.rep", "Dmape.rep",  "Tturn.rep",
            "I2.rep", "I3.rep", "I4.rep", "I.rep"
) 

parms3 <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
            "N2",  "N3", "N4",
            "mu2", "alpha2", "beta2",  "gamma2", "delta2",
            "Dssm.obs", "Dmape.obs",  "Tturn.obs", 
            "Dssm.rep", "Dmape.rep",  "Tturn.rep",
            "I2.rep", "I3.rep", "I4.rep", "I.rep"
) 

#  , "pe3", "pe2",
# "I.exp", "I2.rep", "I3.rep", "I4.rep", "I.rep","I2", "I3", "I4", "I",
# "Tt1.obs", "Tt2.obs", "Tt3.obs", "Tt1.rep", "Tt2.rep", "Tt3.rep",

# model----
b <- 1
if (b==1){ # model with separate parms for each age
    parms = parms1
    tC = cap.v7
    tC.txt = "cap.v7"
    vars_vAR <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs")
    vars_Nyear <- c("N2[10]", "N3[10]", "N4[10]")
    vars_N2 <- c("mu2[10]","alpha2", "beta2",  "gamma2", "delta2")
    vars_N3 <- c("mu3[10]", "alpha3", "gamma3", "delta3", "epsilon3")
    vars_N4 <- c("mu4[10]", "alpha4", "gamma4", "delta4", "epsilon4")

} else if (b==2) { # model separate parms for N2 and N3:N4
    parms = parms2
    tC = cap.v8
    tC.txt = "cap.v8"
    vars_vAR <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs")
    vars_Nyear <- c("N2[10]", "N3[10]", "N4[10]")
    vars_N2 <- c("mu2[10]","alpha2", "beta2",  "gamma2", "delta2")
    vars_N3 <- c("mu3[10]", "alpha3", "gamma3", "delta3", "epsilon3", "mu4[10]")
    vars_N4 <- c(NA)
} else if (b==3) { # demographic model - no forecast for N3:N4
    parms = parms3
    tC = cap.v9
    tC.txt = "cap.v9"
    vars_vAR <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs")
    vars_Nyear <- c("N2[10]", "N3[10]", "N4[10]")
    vars_N2 <- c("mu2[10]","alpha2", "beta2",  "gamma2", "delta2")
    vars_N3 <- c(NA)
    vars_N4 <- c(NA)

}

# MCMC settings
ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
#ni <- 200000; nt <- 30; nb <- 30000; nc <- 3
#ni <- 2000000; nt <- 150; nb <- 300000; nc <- 3

# run model
#source("IPM_mod.R")
ssm26 <- jags(jags.data, parameters=parms, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(tC))
ssm26

# JAGS ouput ----
out <- ssm26$BUGSoutput 
str(out$sims.list)

# DIC to dashboard
ssm26_dic <- out$DIC


## extract raw values from chains
raw <- ls_out(out)
str(raw)

#extract medians, credible intervals, and prediction intervals
ls_all <- ls_med(raw)
calc <- ls_all$ls_med
cri <- ls_all$ls_cri
pri <- ls_all$ls_pri
str(ls_all,1)
str(calc, 1)
str(cri, 1)
str(pri, 1)

cbind(1999:2023, calc$Nt2, jd$I2, calc$N3, jd$I3, calc$mu3)


# figures ----
# N2: observation median v process median
 plot(jd$I2, calc$N2)
# N3: observation median v process median
 plot(jd$I3, calc$N3)
# N4: observation median v process median
 plot(jd$I4, calc$N4)
# Observation median over time
 plot(seq(1999:2023), ls_all$N_med)


## IPM plot----
# variables for IPM plots
year <- 1999:2021
ly <- length(year)
forecast <- 2022:2023
lf <- length(forecast)

#source("IPM_fun.R")

# combined N2-N4[t]
tmp_plot <- ipm_plot(df_med = ls_all$N_med, df_cri = ls_all$N_ci, df_pri = ls_all$Pr_ci, df_dat = df_cap[15:39,]) # ignore warnings - all legit NAs although df_cap needs to be updated.
tmp_plot <- tmp_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = log(exp(I2) + exp(I3)), x = year),
                                      shape = 16, size = 2)
tmp_plot

# N2[t] - create plot, then add the capelin data
tmpN2_plot <- ipm_plot(df_med = calc$N2, df_cri = cri$N2_cri, df_pri = pri$I2.rep_pri, df_dat = df_cap[15:39,]) # ignore warnings - all legit NAs although df_cap needs to be updated.
tmpN2_plot <- tmpN2_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = I2, x = year),
                                      shape = 16, size = 2)
tmpN2_plot


# N3[t]
tmpN3_plot <- ipm_plot(df_med = calc$N3, df_cri = cri$N3_cri, df_pri = pri$I3.rep_pri, df_dat = df_cap[15:39,]) # ignore warnings - all legit NAs although df_cap needs to be updated.
#tmpN3_plot <- tmpN3_plot + 
    
tmpN3_plot <- tmpN3_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = I3, x = year),
                                      shape = 16, size = 2)
tmpN3_plot 



tmpN4_plot <- ipm_plot(df_med = calc$N4, df_cri = cri$N4_cri, df_pri = pri$I4.rep_pri, df_dat = df_cap[15:39,]) # ignore warnings - all legit NAs although df_cap needs to be updated.
#tmpN3_plot <- tmpN3_plot + 

tmpN4_plot <- tmpN4_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = I4, x = year),
                                      shape = 16, size = 2)
tmpN4_plot 

# ggsave("N4_plot.pdf")



# Diagnostics----
# these are just for when figures need to be saved to folders
filepath_gen <- "biomass_cond_ag1_2_DIC_R3"
filepath <- paste0(filepath_gen, "/recruitment_1")

# check convergence
out$summary[rownames(out$summary), c("Rhat")]
out$summary[, 8:9] # these all look really good suggesting good convergence.  This line is easier to read than the one above but I wanted to know how to extract that info just in case.

# calculations for effective sample size - n.eff should be > # of chains *100 (check on this). It seems fine with 2M runs.  But tab_neff shows the samples less than 300.
N_samples <- nc*(ni-nb)/nt
neff <- nc*100  #n.eff should be >nc*100
tab_neff <- out$summary[rownames(out$summary), c("n.eff")][out$summary[rownames(out$summary), c("n.eff")]< 300]
tab_neffa <- out$summary[rownames(out$summary), c("n.eff")][out$summary[rownames(out$summary), c("n.eff")]> 300]

tab_neff


## Mixing ----
# Asess mixing of chains to see if one MCMC goes badly Zuur et al. 2013, pg 83
# vars for forecast model
    # source('C:/Users/lewiske/Documents/R/zuur_rcode/MCMCSupportHighstatV2.R')

# vars for variances
#MyBUGSChains(out, vars_vAR)
mix_var <- MyBUGSChains(out, vars_vAR)
#ggsave(MyBUGSChains(out, vars3), filename = paste0("Bayesian/", filepath, "/chains-variance.pdf"), width=10, height=8, units="in")

# vars for state space and demographic vars
#MyBUGSChains(out, vars_Nyear)
mix_vars_Nyear <- MyBUGSChains(out, vars_Nyear)
#ggsave(MyBUGSChains(out, vars2), filename = paste0("Bayesian/", filepath, "/chains-demographic.pdf"), width=10, height=8, units="in")


### N2 ----
#MyBUGSChains(out, vars_N2)
mix_N2 <- MyBUGSChains(out, vars_N2)
#ggsave(MyBUGSChains(out, vars1), filename = paste0("Bayesian/", filepath, "/chains-forecast.pdf"), width=10, height=8, units="in")

### N3 ----
#MyBUGSChains(out, vars_N3)
if(b ==1 | b==2){
    mix_N3 <- MyBUGSChains(out, vars_N3)    
}

#ggsave(MyBUGSChains(out, vars1), filename = paste0("Bayesian/", filepath, "/chains-forecast.pdf"), width=10, height=8, units="in")

### N4 ----
#MyBUGSChains(out, vars_N4)
if(b == 1){
    mix_N4 <- MyBUGSChains(out, vars_N4)    
}

#ggsave(MyBUGSChains(out, vars1), filename = paste0("Bayesian/", filepath, "/chains-forecast.pdf"), width=10, height=8, units="in")

##autocorrelation ----
#MyBUGSACF(out, vars_vAR)
autocorr_vars_vAR <- MyBUGSACF(out, vars_vAR)
#ggsave(MyBUGSACF(out, vars1), filename = paste0("Bayesian/", filepath, "/auto_corr-forecast.pdf"), width=10, height=8, units="in")

#MyBUGSACF(out, vars_Nyear)
autocorr_vars_Nyear <- MyBUGSACF(out, vars_Nyear)
#ggsave(MyBUGSACF(out, vars2), filename = paste0("Bayesian/", filepath, "/auto_corr-demographic.pdf"), width=10, height=8, units="in")

### N2 ----
#MyBUGSACF(out, vars_N2)
autocorr_N2 <- MyBUGSACF(out, vars_N2)
#ggsave(MyBUGSACF(out, vars3), filename = paste0("Bayesian/", filepath, "/auto_corr-autocorrelation.pdf"), width=10, height=8, units="in")

### N3 ----
#MyBUGSACF(out, vars_N3)
if(b ==1 | b ==2){
    autocorr_N3 <- MyBUGSACF(out, vars_N3)    
}

#ggsave(MyBUGSACF(out, vars3), filename = paste0("Bayesian/", filepath, "/auto_corr-autocorrelation.pdf"), width=10, height=8, units="in")

### N4 ----
#MyBUGSACF(out, vars_N4)
if (b==1){
    autocorr_N4 <- MyBUGSACF(out, vars_N4)    
}

#ggsave(MyBUGSACF(out, vars3), filename = paste0("Bayesian/", filepath, "/auto_corr-autocorrelation.pdf"), width=10, height=8, units="in")


## Model Validation ----
# #(see Schuab and Kery 2022, pg 272-282 - note I tried to do a lot of this in R but abandonded it as I wasn't getting proper results - not sure how to do the algebra in R

### Posterior Predictive Checks (PPC) ----
#### Mean absolute percentage error see formula
Dmape.obs <- out$sims.list$Dmape.obs # this is a GOF but i'm 
Dmape.rep <- out$sims.list$Dmape.rep
# Dmape.GOF <- 100/7000*sum(Dmape.obs) # no idea if this is "right"

#### Bayesian p-value
pB <- mean(Dmape.rep > Dmape.obs)
pB # this shoudl be ~ 0.5 but not near 0 or 1

##### plot of PPC
p <- ggplot()
p <- p + geom_point(aes(x = Dmape.obs, y = Dmape.rep))
p <- p + geom_abline(intercept = 0, slope = 1)
p <- p + xlab("Discrepancy observed data") + ylab("Discrepancy replicate data")
p <- p + xlim(0, max(Dmape.obs))
p <- p + theme_bw()
p <- p + annotate(geom = "text", x = 4, y = 0.5, label = bquote(p[B]), colour = "red")
p <- p + annotate(geom = "text", x = 4.6, y = 0.5, label = paste("=", round(pB, 2)), colour = "black")
p
obs_v_rep <- p

# number of switches, i.e., jaggedness pg 274 & 279 in S&K
hist(out$sims.list$Tturn.rep, xlim = c(7, 22), xlab = "Number of switches \n (replicated data)")
abline(v= out$mean$Tturn.obs, col = "red")

# need simulated data----

### - Residuals for Covariates----
#### (see Zuur et al. 2013 for options for calculating Pearson residuals) - note that I am opting to do a lot of this outside of JAGS due to run time issues.
#### Residual diagnostics

##### raw residuals - I reason that the N2 is the process which is what the linear model is predicting and the mu is the fitted value
##### this is mu for N2 which missed the first 4 years.
if (b ==1 ){
    resN2 <- raw$N2[5:25] - raw$mu2
    resN3 <- raw$N3 - raw$mu3
    resN4 <- raw$N4 - raw$mu4
} else if (b == 2){
    resN2 <- raw$N2[5:25] - raw$mu2
    resN3 <- raw$N3 - raw$mu3
} else if (b == 3){
    resN2 <- raw$N2[5:25] - raw$mu2
}

#### Pearson residuals
presN2 <- sweep(resN2, 1, raw$tau.obs, "/")  # I do not understand why the "1" works as this indicates rowwise division but it seems to work based on the work below, i.e. change values of z and w to manually get same results as presN2 
if (b ==1|b==2){
    presN3 <- sweep(resN3, 1, raw$tau.obs, "/")    
}

if (b==1){
    presN4 <- sweep(resN4, 1, raw$tau.obs, "/")    
}



str(resN2)
str(raw$tau.obs)
str(as.vector(raw$tau.obs))
str(presN2)
# z <- 7500 # row
# w<- 21 # column
#     
# x <- resN2[z,w]
# y <- raw$tau.obs[z]
# x/y
# presN2[z,w]


# get median values
# Cooks' D - Zuur pg 58 is a leave-one-observation-out measure of influence
if(b ==1){
    # raw resids
    resN2_mean <- apply(resN2, 2, 'mean')
    resN3_mean <- apply(resN3, 2, 'mean')
    resN4_mean <- apply(resN4, 2, 'mean')
    
    # Pearson resids but I don't think we need these
    presN2_mean <- apply(presN2, 2, 'mean')
    presN3_mean <- apply(presN3, 2, 'mean')
    presN4_mean <- apply(presN4, 2, 'mean')
    
    dN2 <- presN2_mean^2
    dN3 <- presN3_mean^2
    dN4 <- presN4_mean^2
    
} else if (b ==2 ){
    # raw resids
    resN2_mean <- apply(resN2, 2, 'mean')
    resN3_mean <- apply(resN3, 2, 'mean')

    # Pearson resids but I don't think we need these
    presN2_mean <- apply(presN2, 2, 'mean')
    presN3_mean <- apply(presN3, 2, 'mean')
    
    dN2 <- presN2_mean^2
    dN3 <- presN3_mean^2

} else if (b ==3){
    # raw resids
    resN2_mean <- apply(resN2, 2, 'mean')
    # Pearson resids but I don't think we need these
    presN2_mean <- apply(presN2, 2, 'mean')
    dN2 <- presN2_mean^2
}



# x<- 21
# dN2
# presN2_mean[x]^2

#pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=calc$N2[5:25], y = resN2_mean, xlab = "Fitted values", ylab = "Raw residuals")
abline(h = 0, lty = 2)
# # see notes in Mortality model
plot(y = jd$I2, x = calc$N2, xlab = "Fitted values", ylab = "Observed data") # should follow the line
abline(coef = c(0,1), lty = 2)
#normality
histogram(resN2_mean)
# Cook's D
plot(y = dN2, x = 1:21, xlab = "Observation", ylab = "Cook's D") # should follow the line

par(mfrow = c(1,1))
dev.off()

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
#pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(jags.data$LD[3:23], resN2_mean, xlab = "Larval Density", ylab = "Pearson resids")
plot(jags.data$TI[5:25], resN2_mean, xlab = "Ice retreat", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()


# Posteriors & Priors ----
priormean <- 0
priorsd <- 100
alpha <- out$sims.list$alpha2
prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
limits <- c(min(alpha)-0.3, max(alpha) + 0.3)
x_label <- "Intercept"
bin_1 <- mean(alpha)/100

df_quant <- quantile(alpha, c(0.025, 0.975))
df_cred <- subset(alpha, alpha > df_quant[1] & alpha < df_quant[2])


p1 <- postPriors(df = alpha, df2 = prior, df3 =df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)

p1


# End----

### archived code-----
### overdispersion
# I don't think I need this code as there are no distributions that can be overdispersed in the current model (only Guassian and gamma)
# - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
# I think that this is only needed for distributions that can be overdispersed like Poisson.  Gamma has a dispersion parameter and normal assumes constant variance.
# But confirm with PDR.

# # But, this may not be a problem for nomral and uniform distirbutions - seems to be mostly a Poisson and perhaps binomial thing.
# mean(out$sims.list$FitNew > out$sims.list$Fit)
# # mean = 0.546333
# 
# #squared resids
# presN2_sq <- presN2^2
# presN2_fit <- rowSums(presN2_sq)
# 
# # pluggin in different rows - results seem to be the same for the rowSums approach as for when these are calculated with subscripts.
# sum(presN2_sq[3000,])
# presN2_fit[3000]
# 
# # need to simulate data for the New values
# mean(out$sims.list$FitNew > presN2_fit)


## variance/mean
# I don't think I need this either and S&K just do it for the Poisson model 
# 


#####My attempt to do mean absolute percentage error in R BUT i CAN'T GET IT TO WORK SO DOING IT IN JAGS FOR NOW

# this is mu for N2 which missed the first 4 years.

# I.exp <- out$sims.list$I.exp  # expected value of I based on sum of Ns
# I.rep <- out$sims.list$I.rep  # simulated value of I based on sum of simulated I values
# 
# 
# ss.exp <- I.exp[, 1:c(ly)] # don't want predicted values
# lcap <- length(df_cap$abundance_med) - 2 # get length of capelin series but subtract the NAs that were added for the prediction
# ss.obs <- log(df_cap$abundance_med[15:lcap]*1000) # put appropriate years on log scale and appropriate scale
# 
# ### Mean absolute percentage error see formula in S&K pg. 274.
#  n = jd$n.occasions
#  Dmape.obs <- 100/jd$n.occasions[1] *sum(abs((ss.obs-ss.exp)/ss.obs), na.rm = T)
#  Dmape.obs.t <- 100/jd$n.occasions[1] *rowSums(abs((ss.obs-ss.exp)/ss.obs), na.rm = T)
# # sum(Dmape.obs)
# # 

# ss.rep <- I.rep[1:c(ly)]
# Dmape.rep <- 100/jd$n.occasions[1] *sum(abs((ss.rep-ss.exp)/ss.rep), na.rm = T)
# Dmape.rep.t <- 100/jd$n.occasions[1] *colSums(abs((ss.rep-ss.exp)/ss.rep), na.rm = T)
# 
