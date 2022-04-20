
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
parms <- c("Sld", "N2",  "N3", "tau.proc", "tau.obs", "tau.LD", "tau.ind", "I2", "I3", "y2", "y3")

# MCMC settings
ni <- 20000; nt <- 6; nb <- 5000; nc <- 3

# run model
ssm26 <- jags(jags.data, parameters=parms, n.iter=ni, n.burnin = nb, n.chains=nc, n.thin=nt, model.file = textConnection(cap.v7))
ssm26

# create ouput
out <- ssm26$BUGSoutput 

## extract raw values from chains
raw <- ls_out(out)
str(raw)

#extract medians, credible intervals, and prediction intervals
calc <- ls_med(raw)
df_calc <- do.call(rbind, calc) # this doesn't work
write(df_calc, "out2.csv")
str(calc)
cbind(N2_med, N3_med, N2_med+N3_med)

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

