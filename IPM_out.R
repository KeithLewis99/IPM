
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
# model - the value of b will determine what model and parameters from IPM_JAGS-settings.R
b <- 9
smoother <- "no"
matrix <- "yes"

source("IPM_JAGS-settings.R")


# MCMC settings
ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
#ni <- 200000; nt <- 60; nb <- 30000; nc <- 3
# ni <- 2000000; nt <- 600; nb <- 300000; nc <- 3
# ni <- 5000000; nt <- 1000; nb <- 300000; nc <- 3 # this produces really nice ACFs!!!!

# these are just preliminary values for p (auto-regression: AR) and q (MA: moving average)
if(smoother == "yes"){
    jags.data$p <- 2
    jags.data$q <- 2
}

jags.data.m$Ni <- 3

# these are values to make the JAGS code more generalized, i.e., that the indices are not hard coded.  Currently applies only to cap.v20.

if(disaggregated == "1985-present"){
    jags.data$N2end <- 18
    jags.data$N2start <- 19
    jags.data$N3end <- 11
    jags.data$N3start <- 12
} else {
    jags.data$N2end <- 4
    jags.data$N2start <- 5
    jags.data$N3end <- jags.data$n.occasions
    jags.data$N3start <- 11
}

# run model----
#source("IPM_mod.R")
if(matrix == "no") {
    ssm26 <- jags(jags.data, parameters=parms, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(tC))
    ssm26
    out <- ssm26$BUGSoutput 

    
} else if(matrix == "yes"){
    ssm27 <- jags(jags.data.m, parameters=parms, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(tC))
    ssm27
    out <- ssm27$BUGSoutput 
    out$sims.list$N2 <- out$sims.list$N[,,1]
    out$sims.list$N3 <- out$sims.list$N[,,2]
    out$sims.list$N4 <- out$sims.list$N[,,3]
    out$sims.list$N <- NULL
    out$sims.list$mu2 <- out$sims.list$mu[,,1]
    out$sims.list$mu3 <- out$sims.list$mu[,,2]
    out$sims.list$mu4 <- out$sims.list$mu[,,3]
    out$sims.list$eps2 <- out$sims.list$eps[,,1]
    out$sims.list$eps3 <- out$sims.list$eps[,,2]
    out$sims.list$eps4 <- out$sims.list$eps[,,3]
    out$sims.list$posa2 <- out$sims.list$posa[,,1]
    out$sims.list$posa3 <- out$sims.list$posa[,,2]
    out$sims.list$posa4 <- out$sims.list$posa[,,3]
}

str(out$sims.list)

# JAGS output ----
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

if(disaggregated == "1985-present") {
    cbind(1985:2023, calc$Nt2, jd$I2, calc$N3, jd$I3, calc$mu3)
} else {
    cbind(1999:2023, calc$Nt2, jd$I2, calc$N3, jd$I3, calc$mu3)
} 


# DIC to dashboard
ssm26_dic <- out$DIC

# Bayesian R-squared - see scratch pad



# figures ----

if(disaggregated == "1985-present") {
    year <- 1985:2021
 } else {
     year <- 1999:2021
 } 

ly <- length(year)
forecast <- 2022:2023
lf <- length(forecast)

# N2: observation median v process median
plot(jd$I2, calc$N2)
# N3: observation median v process median
plot(jd$I3, calc$N3)
# N4: observation median v process median
plot(jd$I4, calc$N4)
# Observation median over time
plot(c(year,forecast), ls_all$N_med)

## IPM plot----
# variables for IPM plots

# set capelin years
if(disaggregated == "1985-present") {
    cap <- df_cap 
} else {
    cap <- df_cap[15:39,]
}

# combined N2-N4[t]
tmp_plot <- ipm_plot(df_med = ls_all$N_med, df_cri = ls_all$N_ci, df_pri = ls_all$Pr_ci, df_dat = cap) # ignore warnings - all legit NAs although df_cap needs to be updated. 
tmp_plot <- tmp_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = log(exp(I2) + exp(I3)), x = year),
                                      shape = 16, size = 2)
tmp_plot



# N2[t] - create plot, then add the capelin data
tmpN2_plot <- ipm_plot(df_med = calc$N2, df_cri = cri$N2_cri, df_pri = pri$I2.rep_pri, df_dat = cap) # ignore warnings - all legit NAs although df_cap needs to be updated.
tmpN2_plot <- tmpN2_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = I2, x = year),
                                      shape = 16, size = 2)
tmpN2_plot


# N3[t]
#source("IPM_fun.R")
tmpN3_plot <- ipm_plot(df_med = calc$N3, df_cri = cri$N3_cri, df_pri = pri$I3.rep_pri, df_dat = cap) # ignore warnings - all legit NAs although df_cap needs to be updated.
#tmpN3_plot <- tmpN3_plot + 
    
tmpN3_plot <- tmpN3_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = I3, x = year),
                                      shape = 16, size = 2)
tmpN3_plot 



tmpN4_plot <- ipm_plot(df_med = calc$N4, df_cri = cri$N4_cri, df_pri = pri$I4.rep_pri, df_dat = cap) # ignore warnings - all legit NAs although df_cap needs to be updated.
#tmpN3_plot <- tmpN3_plot + 

tmpN4_plot <- tmpN4_plot + geom_point(data = df_dis_tabLog,
                                      aes(y = I4, x = year),
                                      shape = 16, size = 2)
tmpN4_plot 


# Process error----
## One-step ahead resids----
## see AugerMethe 2021
### couldn't get this to run in JAGS
plot(density(out$sims.list$osa[,,1]))
str(out$sims.list$osa[,,1])

# this give the standard deviation for each t - I don't think this is right bc its sd[t] which is not in the formula
# sd_osa <-  apply(out$sims.list$osa[,,1],2,'sd')
# str(sd_osa)
# 
# # this gives the sd for the whole matrix 
# # then, the Pearson resids for the matrix/a single sd and then the median - this also doesn't seem right because its the sd of all the various iterations
# sd_osa <- sd(out$sims.list$osa[,,1])
# posa <- sweep(out$sims.list$osa[,,1], 2, sd_osa, FUN="/")
# str(posa)
# posa_med <-  apply(out$sims.list$osa[,,1],2,'median')
# 
# # this is the median of the osa divided by a single sd - also doesn't seem right bc its a median value divided by the single sd
# tmp <- apply(out$sims.list$osa[,,1],2,'median')
# posa <- tmp/sd_osa

# sd for each row - calculate teh sd across the row, then divide each osa by that sd - I think that this is the right one and the intent of Auger-Methe.  I can also do this in JAGS now but this helped me to figure out how to write the JAGS code

osa_sd <-  apply(out$sims.list$osa[,,1],1,'sd')
str(t(osa_sd))
posa <- sweep(out$sims.list$osa[,,1], 1, osa_sd, FUN="/")
posa_med <- apply(posa, 2, 'median')
str(out$sims.list$osa[,,1])
str(posa)
str(posa_med)
plot(density(posa[,1]))

## Bubble Plot - osa/posa resids
# this shows that the JAGS approach and the R approach are equivalent
plot(posa_med, calc$posa2)

posa_df <- as.data.frame(cbind(calc$posa2, calc$posa3, calc$posa4))
posa_df <- cbind(posa_df, 1985:2023)
posa_df <- posa_df %>% rename(posa2 = V1, posa3 = V2, posa4 = V3, year = '1985:2023')
posa_wide <- pivot_longer(posa_df, cols = c("posa2", "posa3", "posa4"), names_to = "age", values_to = "pres")
posa_wide$sign <- NA

#eps_wide <- pivot_longer(eps, cols = c("eps2", "eps3", "eps4"), names_to = "age", values_to = "pres")
posa_wide$sign <- NA

## create positive and negative colours
for(i in seq_along(posa_wide$pres)){
     if(posa_wide$pres[i] > 0){
          posa_wide$sign[i] <- "pos"
     } else {
          posa_wide$sign[i] <- "neg"
     }
}
head(posa_wide)

## Bubble Plot - posa resids
p <- ggplot(data = posa_wide, aes(x = year, y = age, size = pres, colour = sign))
p <- p + geom_point()
p


# qqplot
qqnorm(calc$N2)
qqline(calc$N2)

plot(1985:2023, calc$posa2)
qqnorm(calc$posa2)
qqline(calc$posa2)

acf(calc$posa2)

# Bubble Plot - raw resids
## bind median residuals and then pivot them to make a long data set
eps <- as.data.frame(cbind(calc$eps2, calc$eps3, calc$eps4))
eps <- cbind(eps, 1985:2023)
eps <- eps %>% rename(eps2 = V1, eps3 = V2, eps4 = V3, year = '1985:2023')
eps_wide <- pivot_longer(eps, cols = c("eps2", "eps3", "eps4"), names_to = "age", values_to = "pres")
eps_wide$sign <- NA

## create positive and negative colours
for(i in seq_along(eps_wide$pres)){
     if(eps_wide$pres[i] > 0){
          eps_wide$sign[i] <- "pos"
     } else {
          eps_wide$sign[i] <- "neg"
     }
}
head(eps_wide)

## Bubble Plot - raw resids
p <- ggplot(data = eps_wide, aes(x = year, y = age, size = pres, colour = sign))
p <- p + geom_point()
p



# trends in process error
## create a data frame with the years and residuals for each year
calc$year <- rep(NA, 39)
calc$year <- 1985:2023
eps_trend <- as.data.frame(cbind(calc$year, calc$eps2, calc$eps3, calc$eps4))
eps_trend <- eps_trend %>% rename(year = V1, eps2 = V2, eps3 = V3, eps4 = V4)
str(calc)
str(eps_trend)

# create a dataframe with the years and credible intervals for each year
eps2_cri <- as.data.frame(cbind(year = 1985:2023, min = cri$eps2_cri[1,], max = cri$eps2_cri[2,]))
eps3_cri <- as.data.frame(cbind(year = 1985:2023, min = cri$eps3_cri[1,], max = cri$eps3_cri[2,]))
eps4_cri <- as.data.frame(cbind(year = 1985:2023, min = cri$eps4_cri[1,], max = cri$eps4_cri[2,]))


# eps2: create a plot 
p <- ggplot()
p <- p + geom_point(data = eps_trend, aes(x = year, y = eps2))
p <- p + geom_ribbon(data=eps2_cri, aes(x = year,
                    ymax = max, 
                    ymin = min),
                alpha = 0.5, fill = "grey")
p

# eps3: create a plot 
p <- ggplot()
p <- p + geom_point(data = eps_trend, aes(x = year, y = eps3))
p <- p + geom_ribbon(data=eps3_cri, aes(x = year,
                                        ymax = max, 
                                        ymin = min),
                     alpha = 0.5, fill = "grey")
p

# eps4: create a plot 
p <- ggplot()
p <- p + geom_point(data = eps_trend, aes(x = year, y = eps4))
p <- p + geom_ribbon(data=eps4_cri, aes(x = year,
                                        ymax = max, 
                                        ymin = min),
                     alpha = 0.5, fill = "grey")
p


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

low <- min(out$sims.list$Tturn.rep)
high <- max(out$sims.list$Tturn.rep)
# number of switches, i.e., jaggedness pg 274 & 279 in S&K
hist(out$sims.list$Tturn.rep, xlim = c(low, high), xlab = "Number of switches \n (replicated data)")
abline(v= out$mean$Tturn.obs, col = "red")

# need simulated data----

### - Residuals for Covariates----
#### (see Zuur et al. 2013 for options for calculating Pearson residuals) - note that I am opting to do a lot of this outside of JAGS due to run time issues.
#### Residual diagnostics

if(disaggregated == "1985-present") {
    i <- 1:39
} else {
    i <- 1:25
} 

##### raw residuals - I reason that the N2 is the process which is what the linear model is predicting and the mu is the fitted value
##### this is mu for N2 which missed the first 4 years.
if (b == 1|b == 6|b == 5){
    resN2 <- raw$N2[,i] - raw$mu2
    resN3 <- raw$N3[,i] - raw$mu3
    resN4 <- raw$N4[,i] - raw$mu4
} else if (b == 2){
    resN2 <- raw$N2[,i] - raw$mu2
    resN3 <- raw$N3[,i] - raw$mu3
} else if (b == 3){
    resN2 <- raw$N2[i] - raw$mu2
}

#### Pearson residuals
presN2 <- sweep(resN2, 1, raw$tau.obs, "/")  # I do not understand why the "1" works as this indicates rowwise division but it seems to work based on the work below, i.e. change values of z and w to manually get same results as presN2 
if (b ==1|b==2|b==6|b==5){
    presN3 <- sweep(resN3, 1, raw$tau.obs, "/")    
}

if (b==1|b==6|b==5){
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
if(b ==1|b==6|b==5){
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
plot(x=calc$N2[i], y = resN2_mean, xlab = "Fitted values", ylab = "Raw residuals")
abline(h = 0, lty = 2)
# # see notes in Mortality model
plot(y = jd$I2, x = calc$N2, xlab = "Fitted values", ylab = "Observed data") # should follow the line
abline(coef = c(0,1), lty = 2)
#normality
histogram(resN2_mean)

if(disaggregated == "1985-present"){
    N2xaxis <- 1:39
    N3xaxis <- 1:39
    N4xaxis <- 1:39
} else {
    N2xaxis <- 1:25
    N3xaxis <- 1:25
    N4xaxis <- 1:25
}
# Cook's D
plot(y = dN2, x = N2xaxis, xlab = "Observation", ylab = "Cook's D") # should follow the line

par(mfrow = c(1,1))
dev.off()


if(disaggregated == "1985-present"){ # set resN2_mean
    LDxaxis <- 17:37
    TIxaxis <- 1:39
} else {
    LDxaxis <- 3:23
    TIxaxis <- 5:25
}

# Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
#pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(jags.data$LD[LDxaxis], resN2_mean[LDxaxis], xlab = "Larval Density", ylab = "Pearson resids")
plot(jags.data$TI[TIxaxis], resN2_mean[TIxaxis], xlab = "Ice retreat", ylab = "Pearson resids")
par(mfrow = c(1,1))
dev.off()


# Posteriors & Priors ----
#alpha

a2 <- post_param(param = "alpha", priormean = 0, priorsd = 100, jags = out$sims.list$alpha2, x_label = "Intercept - alpha") 

# pa2 <- postPriors(df = a2$jags, df2 = a2$prior, df3 = a2$df_cred, limits=a2$limits, x_label=a2$x_label, priormean=a2$priormean, priorsd=a2$priorsd, by_bin = a2$bin_1)
# 
# pa2

#beta
b2 <- post_param(param = "beta", priormean = 0, priorsd = 100, jags = out$sims.list$beta2, x_label = "Larval Density - beta") 

# pb2 <- postPriors(df = b2$jags, df2 = b2$prior, df3 = b2$df_cred, limits = b2$limits, x_label=b2$x_label, priormean=b2$priormean, priorsd=b2$priorsd, by_bin = a2$bin_1)

#gamma - MRI
#source("IPM_fun.R")
g2 <- post_param(param = "gamma", jags = out$sims.list$gamma2) # max rate of increase

 pg2 <- postPriors(df = g2$jags, df2 = g2$prior, df3 = g2$df_cred, limits=g2$limits, x_label=g2$x_label, by_bin = g2$bin_1)
# 
 pg2

#delta - width
d2 <- post_param(param = "delta", priormean = 11.5, priorsd = 5.7, jags = out$sims.list$delta2) 

# pd2 <- postPriors(df = d2$jags, df2 = d2$prior, df3 = d2$df_cred, limits=d2$limits, x_label=d2$x_label, priormean=d2$priormean, priorsd=d2$priorsd, by_bin = d2$bin_1)
# 
# pd2


# process error-----
# this is for a process error that is constant within time periods.
# this cuts the array on the "a" which is the tau.proc for different ages. Then the apply cuts it by the z which is the time series.
tp1 <- apply(out$sims.list$tau.proc[,1,],2,'median')
tp2 <- apply(out$sims.list$tau.proc[,2,],2,'median')
cri_tp1 <- apply(out$sims.list$tau.proc[,1,],2,'quantile', c(0.025, 0.975))
cri_tp2 <- apply(out$sims.list$tau.proc[,2,],2,'quantile', c(0.025, 0.975))

plot(1:4, cri$tau.proc_cri)

# these are the estiamtes for a fully time varying process error
tp1 <- apply(out$sims.list$tau.proc[,,1],2,'median')
tp2 <- apply(out$sims.list$tau.proc[,,2],2,'median')
cri_tp1 <- apply(out$sims.list$tau.proc[,,1],2,'quantile', c(0.025, 0.975))
cri_tp2 <- apply(out$sims.list$tau.proc[,,2],2,'quantile', c(0.025, 0.975))


# remember that tau is the precision - t
plot(1985:2023, 1/tp1, type= 'n', ylim = c(0,500))
lines(1985:2023, 1/cri_tp1[2,])
lines(1985:2023, 1/cri_tp1[1,])

plot(1985:2023, 1/tp2, type= 'n', ylim = c(0,500))
lines(1985:2023, 1/cri_tp2[2,])
lines(1985:2023, 1/cri_tp2[1,])

# trying the same with the variance
s2 <- apply(out$sims.list$sigma2.proc[,,1],2,'median')
s3 <- apply(out$sims.list$sigma2.proc[,,2],2,'median')
cri_s2 <- apply(out$sims.list$sigma2.proc[,,1],2,'quantile', c(0.025, 0.975))
cri_s3 <- apply(out$sims.list$sigma2.proc[,,2],2,'quantile', c(0.025, 0.975))

plot(1985:2023, s2, type= 'n', ylim = c(0,300))
lines(1985:2023, 1/cri_s2[2,])
lines(1985:2023, 1/cri_s2[1,])

plot(1985:2023, s3, type= 'n', ylim = c(0,300))
lines(1985:2023, 1/cri_s3[2,])
lines(1985:2023, 1/cri_s3[1,])

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
