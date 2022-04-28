
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

# Model ln scale: tice N3 mortality and INdex SE and split N2----
# try to fix the priors ito variance

# JAGS settings
parms <- c("Sld", "N2",  "N3", "mu", "tau.proc", "tau.obs", "tau.LD", "tau.ind", "I2", "I3", "I", "I2.rep", "I3.rep", "I.exp", "I.rep", "y2", "y3", "alpha", "beta", "gamma", "sigma", "Tturn.obs", "Tturn.rep")

# MCMC settings
ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
#ni <- 200000; nt <- 30; nb <- 30000; nc <- 3

# run model
source("IPM_mod.R")
ssm26 <- jags(jags.data, parameters=parms, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(cap.v7))
ssm26

# create ouput
out <- ssm26$BUGSoutput 

# Send DIC to dashboard
ssm26_dic <- out$DIC

## extract raw values from chains
raw <- ls_out(out)
str(raw)

#extract medians, credible intervals, and prediction intervals
calc <- ls_med(raw)
str(calc)
#df_calc <- do.call(rbind, calc) # this doesn't work
#write(df_calc, "out2.csv")
#cbind(N2_med, N3_med, N2_med+N3_med)


# calculations for effective sample size - n.eff should be > # of chains *100
Neff <- nc*(ni-nb)/nt
neff <- nc*100  #n.eff should be >nc*100

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


#source("IPM_fun.R")
tmp_plot <- ipm_plot(df1 = calc, df2 = df_cap[15:39,]) # ignore warnings - all legit NAs although df_cap needs to be updated.
tmp_plot
ggsave("tmp_plot1.pdf")

     # Other plots
     tp = out$sims.list$tau.proc
     tp_med = apply(tp,2,'median') # median values of y_pred
     tp_ci = apply(tp,2,'quantile', c(0.1, 0.9)) # median values of y_pred

# Diagnostics----
# these are just for when figures need to be saved to folders
filepath_gen <- "biomass_cond_ag1_2_DIC_R3" 
filepath <- paste0(filepath_gen, "/recruitment_1")

# print
print(out, intervals=c(0.025, 0.975), digits = 3)
out$mean


## N2 ----
### Mixing ----     
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


###autocorrelation ----
MyBUGSACF(out, var1)
autocorr1 <- MyBUGSACF(out, var1)
#ggsave(MyBUGSACF(out, vars1), filename = paste0("Bayesian/", filepath, "/auto_corr-forecast.pdf"), width=10, height=8, units="in")

MyBUGSACF(out, vars2)
autocorr2 <- MyBUGSACF(out, vars2)
#ggsave(MyBUGSACF(out, vars2), filename = paste0("Bayesian/", filepath, "/auto_corr-demographic.pdf"), width=10, height=8, units="in")

MyBUGSACF(out, vars3)
autocorr3 <- MyBUGSACF(out, vars3)
#ggsave(MyBUGSACF(out, vars3), filename = paste0("Bayesian/", filepath, "/auto_corr-autocorrelation.pdf"), width=10, height=8, units="in")


### Model Validation ----
#(see Zuuer et al. 2013 for options for calculating Pearson residuals) - note that I am opting to do a lot of this outside of JAGS due to run time issues.  
# # Residual diagnostics

# this is mu for N2 which missed the first 4 years.
resN2 <- raw$N2[, 5:25] - raw$mu
sigmaJ <- raw$sigmaJ
presN2 <- resN2/as.vector(sigmaJ) # I do not see why this needs as.vector but it seems to work 

#  pluggin in different columns and rows - results seem to be the same for the apply approach as for when these are calculated with subscripts. 
resN2[,21][7500]/sigmaJ[7500]
presN2[,21][7500]
presN2_med = apply(presN2,2,'median')

plot(calc$N2_med[5:25], presN2_med)



# E1 <- out$mean$PRes # Pearson resids
# F1 <- out$mean$expY # Expected values
# N2 <- out$mean$N2   # N2 - observed values? Why do these fill in the NAs of df$ln_biomass_med???
# D <- out$mean$D     # this is SSQ - but i'm looking for Cook'sD
# 
#pdf(paste0("Bayesian/", filepath, "/fit_obs.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x=calc$N2_med[5:25], y = presN2_med, xlab = "Fitted values", ylab = "Pearson residuals")
abline(h = 0, lty = 2)
# see notes in Mortality model
plot(y = calc$I2_med[5:25], x = calc$N2_med[5:25], xlab = "Fitted values", ylab = "Observed data") # should follow the line
abline(coef = c(0,1), lty = 2)
par(mfrow = c(1,1))
dev.off()
# 
# # Residuals v covariates Zuur et al. 2013, pg 59-60: look for no patterns; patterns may indicated non-linear
# pdf(paste0("Bayesian/", filepath, "/resid_covar.pdf"))
par(mfrow = c(2,2), mar = c(5,5,2,2))
#MyVar <- c("tice.std", "meandCond_lag.std")
#df_diag <- as.data.frame(model_data)
#df_diag <- cbind(df_diag, E1)
plot(jags.data$LD[3:23], presN2_med, xlab = "Larval Density", ylab = "Pearson resids")
plot(jags.data$TI[5:25], presN2_med, xlab = "Ice retreat", ylab = "Pearson resids")
par(mfrow = c(1,1))
# dev.off()



#overdispersion - values close to 0.5 indicate a good fit pg. 77 Zuur et al. 2013 Beginners guide to GLM and GLMM
# But, this may not be a problem for nomral and uniform distirbutions - seems to be mostly a Poisson and perhaps binomial thing.  
mean(out$sims.list$FitNew > out$sims.list$Fit)
# mean = 0.546333

#squared resids
presN2_sq <- presN2^2
presN2_fit <- rowSums(presN2_sq)

# pluggin in different rows - results seem to be the same for the rowSums approach as for when these are calculated with subscripts. 
sum(presN2_sq[3000,])
presN2_fit[3000]

# need to simulate data for the New values
mean(out$sims.list$FitNew > presN2_fit)


# variance/mean


## State-space ----


### Model Validation ----
#(see Schuab and Kery 2022, pg 272-282 - note that I am opting to do a lot of this outside of JAGS due to run time issues.  

# this is mu for N2 which missed the first 4 years.
I <- out$sims.list$I
I.exp <- out$sims.list$I.exp
I.rep <- out$sims.list$I.rep

I.exp_median <- out$median$I

#ss.exp <- I.exp[, 1:c(ly)] 
#lcap <- length(df_cap$abundance_med) - 2
#ss.obs <- log(df_cap$abundance_med[15:lcap]*1000)


# Mean absolute percentage error        
        # n = jd$n.occasions
#Dmape.obs <- 100/jd$n.occasions[1] *sum(abs((ss.obs-ss.exp)/ss.obs), na.rm = T)
Dmape.obs <- 100/jd$n.occasions[1] *rowSums(abs((I-I.exp)/I), na.rm = T)

#ss.rep <- calc$I.rep[1:c(ly)] 

#Dmape.rep <- 100/jd$n.occasions[1] *sum(abs((ss.rep-ss.exp)/ss.rep), na.rm = T)
Dmape.rep <- 100/jd$n.occasions[1] *rowSums(abs((I.rep-I.exp)/I.rep), na.rm = T)

# Bayesian p-value
pB <- mean(Dmape.rep > Dmape.obs)

p <- ggplot()
p <- p + geom_point(aes(x = Dmape.obs, y = Dmape.rep))
p <- p + geom_abline(intercept = 0, slope = 1)
p <- p + xlab("Discrepancy observed data") + ylab("Discrepancy replicate data")
p <- p + xlim(0, 60)
p <- p + theme_bw()
p <- p + annotate(geom = "text", x = 50, y = 2, label = bquote(p[B]), colour = "red")
p <- p + annotate(geom = "text", x = 55, y = 2, label = paste("=", round(pB, 2)), colour = "black")
p


# need simulated data        

# number of switches
out$sims.list$Tturn.obs
par(mfrow = c(1,2), mar = c(5,5,2,2))
hist(out$sims.list$Tturn.obs)
hist(out$sims.list$Tturn.rep)



#install.packages('IPMbook')



# plotGOF <- function(jagsout, obs, rep, main=NA, showP=TRUE,
#                     ylab="Discrepancy replicate data", xlab="Discrepancy observed data",
#                     pch=16, cex = 0.8, col=1){
#         OBS <- jagsout$sims.list[[obs]]
#         REP <- jagsout$sims.list[[rep]]
#         lim <- quantile(c(OBS, REP), c(0.0001, 0.999))
#         plot(OBS, REP, pch=pch, cex=cex, ylim=lim, xlim=lim,
#              ylab=ylab, xlab=xlab, main=main, axes=FALSE, col=col)
#         axis(1); axis(2)
#         segments(lim[1], lim[1], lim[2], lim[2], lty=3)
#         bp <- round(mean(REP > OBS),2)
#         if(showP){
#                 loc <- ifelse(bp < 0.5, "topleft", "bottomright")
#                 legend(loc, legend=bquote(p[B]==.(bp)), bty="n")
#         }
#         return(invisible(bp))
# }
# 
# plotGOF(ssm26, "Dmape.obs", "Dmape.rep", main="State-space model", col=alpha(co, 0.3))


# End----