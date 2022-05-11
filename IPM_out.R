
# required packages
library(rjags)
library(R2jags)
library(ggplot2)
library(lattice)

# Start----
rm(list=ls())

# Source files
source("IPM_dat.R")
source("IPM_fun.R")
source("IPM_mod.R")
source('C:/Users/lewiske/Documents/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('C:/Users/lewiske/Documents/R/zuur_rcode/HighstatLibV7.R')

# JAGS settings ----
parms <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs", 
           "N2",  "N3", "N4",
           "I2.rep", "I3.rep", "I4.rep",
           "mu2", "alpha2", "beta2",  "gamma2", "delta2",
           "mu3", "alpha3", "gamma3", "delta3", "epsilon3",
           "mu4", "alpha4", "gamma4", "delta4", "epsilon4"
                  #diagnostics
           ) 

# "Ni2","Nt2","gamma", "delta",  "I.exp", "I.rep",
# "Tturn.obs", "Tturn.rep"      
#  , "pe3", "pe2",

# "sigma", "I2", "I3", "I",
#parms <- c("N2",  "N3", "mu", "tau.proc", "tau.obs", "tau.LD", "tau.ind", "I2", "I3", "I", "I2.rep", "I3.rep", "I.exp", "I.rep", "alpha", "beta", "gamma",  "Tturn.obs", "Tturn.rep") #"sigma",


# MCMC settings
ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
#ni <- 200000; nt <- 30; nb <- 30000; nc <- 3
#ni <- 2000000; nt <- 150; nb <- 300000; nc <- 3

# run model
#source("IPM_mod.R")
ssm26 <- jags(jags.data, parameters=parms, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(cap.v7))
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

# combined N2-N4[t]
tmp_plot <- ipm_plot(df_med = ls_all$N_med, df_cri = ls_all$N_ci, df_pri = ls_all$Pr_ci, df_dat = df_cap[15:39,]) # ignore warnings - all legit NAs although df_cap needs to be updated.
tmp_plot <- tmp_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = log(exp(I2) + exp(I3)), x = year),
                                      shape = 16, size = 1.5)
tmp_plot

# N2[t] - create plot, then add the capelin data
tmpN2_plot <- ipm_plot(df_med = calc$N2, df_cri = cri$N2_cri, df_pri = pri$I2.rep_pri, df_dat = df_cap[15:39,]) # ignore warnings - all legit NAs although df_cap needs to be updated.
tmpN2_plot <- tmpN2_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = I2, x = year),
                                      shape = 16, size = 1.5)
tmpN2_plot


# N3[t]
tmpN3_plot <- ipm_plot(df_med = calc$N3, df_cri = cri$N3_cri, df_pri = pri$I3.rep_pri, df_dat = df_cap[15:39,]) # ignore warnings - all legit NAs although df_cap needs to be updated.
#tmpN3_plot <- tmpN3_plot + 
    
tmpN3_plot <- tmpN3_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = I3, x = year),
                                      shape = 16, size = 1.5)
tmpN3_plot 



tmpN4_plot <- ipm_plot(df_med = calc$N4, df_cri = cri$N4_cri, df_pri = pri$I4.rep_pri, df_dat = df_cap[15:39,]) # ignore warnings - all legit NAs although df_cap needs to be updated.
#tmpN3_plot <- tmpN3_plot + 

tmpN4_plot <- tmpN4_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = I4, x = year),
                                      shape = 16, size = 1.5)
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
vars_vAR <- c("tau.proc2", "tau.proc3", "tau.proc4", "tau.obs")
MyBUGSChains(out, vars_vAR)
mix_var <- MyBUGSChains(out, vars_vAR)
#ggsave(MyBUGSChains(out, vars3), filename = paste0("Bayesian/", filepath, "/chains-variance.pdf"), width=10, height=8, units="in")

# vars for state space and demographic vars
vars_Nyear <- c("N2[10]", "N3[10]", "N4[10]")
MyBUGSChains(out, vars_Nyear)
mix_vars_Nyear <- MyBUGSChains(out, vars_Nyear)
#ggsave(MyBUGSChains(out, vars2), filename = paste0("Bayesian/", filepath, "/chains-demographic.pdf"), width=10, height=8, units="in")


### N2 ----
vars_N2 <- c("mu2[10]","alpha2", "beta2",  "gamma2", "delta2")
MyBUGSChains(out, vars_N2)
mix_N2 <- MyBUGSChains(out, vars_N2)
#ggsave(MyBUGSChains(out, vars1), filename = paste0("Bayesian/", filepath, "/chains-forecast.pdf"), width=10, height=8, units="in")

### N3 ----
vars_N3 <- c("mu3[10]", "alpha3", "gamma3", "delta3", "epsilon3")
MyBUGSChains(out, vars_N3)
mix_N3 <- MyBUGSChains(out, vars_N3)
#ggsave(MyBUGSChains(out, vars1), filename = paste0("Bayesian/", filepath, "/chains-forecast.pdf"), width=10, height=8, units="in")

### N2 ----
vars_N4 <- c("mu4[10]", "alpha4", "gamma4", "delta4", "epsilon4")
MyBUGSChains(out, vars_N4)
mix_N4 <- MyBUGSChains(out, vars_N4)
#ggsave(MyBUGSChains(out, vars1), filename = paste0("Bayesian/", filepath, "/chains-forecast.pdf"), width=10, height=8, units="in")

##autocorrelation ----
MyBUGSACF(out, vars_vAR)
autocorr_vars_vAR <- MyBUGSACF(out, vars_vAR)
#ggsave(MyBUGSACF(out, vars1), filename = paste0("Bayesian/", filepath, "/auto_corr-forecast.pdf"), width=10, height=8, units="in")

MyBUGSACF(out, vars_Nyear)
autocorr_vars_Nyear <- MyBUGSACF(out, vars_Nyear)
#ggsave(MyBUGSACF(out, vars2), filename = paste0("Bayesian/", filepath, "/auto_corr-demographic.pdf"), width=10, height=8, units="in")

### N2 ----
MyBUGSACF(out, vars_N2)
autocorr_N2 <- MyBUGSACF(out, vars_N2)
#ggsave(MyBUGSACF(out, vars3), filename = paste0("Bayesian/", filepath, "/auto_corr-autocorrelation.pdf"), width=10, height=8, units="in")

### N3 ----
MyBUGSACF(out, vars_N3)
autocorr_N3 <- MyBUGSACF(out, vars_N3)
#ggsave(MyBUGSACF(out, vars3), filename = paste0("Bayesian/", filepath, "/auto_corr-autocorrelation.pdf"), width=10, height=8, units="in")

### N4 ----
MyBUGSACF(out, vars_N4)
autocorr_N4 <- MyBUGSACF(out, vars_N4)
#ggsave(MyBUGSACF(out, vars3), filename = paste0("Bayesian/", filepath, "/auto_corr-autocorrelation.pdf"), width=10, height=8, units="in")


## Model Validation ----
#(see Zuur et al. 2013 for options for calculating Pearson residuals) - note that I am opting to do a lot of this outside of JAGS due to run time issues.
# # Residual diagnostics

# this is mu for N2 which missed the first 4 years.
# resN2 <- raw$N2 - raw$mu2
# sigmaJ <- raw$sigmaJ
# presN2 <- resN2/as.vector(sigmaJ) # I do not see why this needs as.vector but it seems to work

#  pluggin in different columns and rows - results seem to be the same for the apply approach as for when these are calculated with subscripts.
# resN2[,21][7500]/sigmaJ[7500]
# presN2[,21][7500]
# presN2_med = apply(presN2,2,'median')

# plot(calc$N2_med[5:25], presN2_med)



# E1 <- out$mean$PRes # Pearson resids
# F1 <- out$mean$expY # Expected values
# N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
# D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD
#
# #pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
# # par(mfrow = c(2,2), mar = c(5,5,2,2))
# # plot(x=calc$N2_med[5:25], y = presN2_med, xlab = "Fitted values", ylab = "Pearson residuals")
# # abline(h = 0, lty = 2)
# # # see notes in Mortality model
# # plot(y = calc$I2_med[5:25], x = calc$N2_med[5:25], xlab = "Fitted values", ylab = "Observed data") # should follow the line
# # abline(coef = c(0,1), lty = 2)
# # par(mfrow = c(1,1))
# # dev.off()
# # # 
# # # # Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
# # # pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
# # par(mfrow = c(2,2), mar = c(5,5,2,2))
# # #MyVar <- c("tice.std", "meandCond_lag.std")
# # #df_diag <- as.data.frame(model_data)
# # #df_diag <- cbind(df_diag, E1)
# # plot(jags.data$LD[3:23], presN2_med, xlab = "Larval Density", ylab = "Pearson resids")
# # plot(jags.data$TI[5:25], presN2_med, xlab = "Ice retreat", ylab = "Pearson resids")
# # par(mfrow = c(1,1))
# # # dev.off()
# 

# ### overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
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


## variance/mean ----
# 
# 
## State-space ----
# 
# 
## Model Validation ----
# #(see Schuab and Kery 2022, pg 272-282 - note that I am opting to do a lot of this outside of JAGS due to run time issues.  
# 
# # this is mu for N2 which missed the first 4 years.
# # I <- out$sims.list$I
# # I.exp <- out$sims.list$I.exp
# # I.rep <- out$sims.list$I.rep
# # 
# # I.exp_median <- out$median$I
# # 
# # #ss.exp <- I.exp[, 1:c(ly)] 
# # #lcap <- length(df_cap$abundance_med) - 2
# # #ss.obs <- log(df_cap$abundance_med[15:lcap]*1000)
# # 
# # 
# # # Mean absolute percentage error        
# #         # n = jd$n.occasions
# # #Dmape.obs <- 100/jd$n.occasions[1] *sum(abs((ss.obs-ss.exp)/ss.obs), na.rm = T)
# # Dmape.obs <- 100/jd$n.occasions[1] *rowSums(abs((I-I.exp)/I), na.rm = T)
# # 
# # #ss.rep <- calc$I.rep[1:c(ly)] 
# # 
# # #Dmape.rep <- 100/jd$n.occasions[1] *sum(abs((ss.rep-ss.exp)/ss.rep), na.rm = T)
# # Dmape.rep <- 100/jd$n.occasions[1] *rowSums(abs((I.rep-I.exp)/I.rep), na.rm = T)
# # 
# # # Bayesian p-value
# # pB <- mean(Dmape.rep > Dmape.obs)
# # 
# # p <- ggplot()
# # p <- p + geom_point(aes(x = Dmape.obs, y = Dmape.rep))
# # p <- p + geom_abline(intercept = 0, slope = 1)
# # p <- p + xlab("Discrepancy observed data") + ylab("Discrepancy replicate data")
# # p <- p + xlim(0, 60)
# # p <- p + theme_bw()
# # p <- p + annotate(geom = "text", x = 50, y = 2, label = bquote(p[B]), colour = "red")
# # p <- p + annotate(geom = "text", x = 55, y = 2, label = paste("=", round(pB, 2)), colour = "black")
# # p
# # 
# # 
# # # need simulated data        
# # 
# # # number of switches
# # out$sims.list$Tturn.obs
# # par(mfrow = c(1,2), mar = c(5,5,2,2))
# # hist(out$sims.list$Tturn.obs)
# # hist(out$sims.list$Tturn.rep)
# 
# 
# 
# #install.packages('IPMbook')
# 
# 
# 
# # plotGOF <- function(jagsout, obs, rep, main=NA, showP=TRUE,
# #                     ylab="Discrepancy replicate data", xlab="Discrepancy observed data",
# #                     pch=16, cex = 0.8, col=1){
# #         OBS <- jagsout$sims.list[[obs]]
# #         REP <- jagsout$sims.list[[rep]]
# #         lim <- quantile(c(OBS, REP), c(0.0001, 0.999))
# #         plot(OBS, REP, pch=pch, cex=cex, ylim=lim, xlim=lim,
# #              ylab=ylab, xlab=xlab, main=main, axes=FALSE, col=col)
# #         axis(1); axis(2)
# #         segments(lim[1], lim[1], lim[2], lim[2], lty=3)
# #         bp <- round(mean(REP > OBS),2)
# #         if(showP){
# #                 loc <- ifelse(bp < 0.5, "topleft", "bottomright")
# #                 legend(loc, legend=bquote(p[B]==.(bp)), bty="n")
# #         }
# #         return(invisible(bp))
# # }
# # 
# # plotGOF(ssm26, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co, 0.3))
# 
# 

# End----