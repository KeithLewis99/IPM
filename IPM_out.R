
# required packages

library(rjags)
library(R2jags)

rm(list=ls())


# Source files
source("IPM_fun.R")
source("IPM_dat.R")
source("IPM_mod.R")
#save <- "no"

parms <- c("Sld", "N2",  "N3", "tau.proc", "tau.obs", "tau.LD", "tau.ind", "I2", "I3", "y2", "y3")

# MCMC settings
ni <- 20000; nt <- 6; nb <- 5000; nc <- 3

jags.data <- ls_jag("yes", "yes")
jd <- as.data.frame(jags.data)

ssm26 <- jags(jags.data, parameters=parms, model.file = textConnection(cap.v7))
ssm26
out <- ssm26$BUGSoutput 

## extract raw values from chains
raw <- ls_out(out)
str(raw)

#extract medians, credible intervals, and prediction intervals

calc <- ls_med(raw)
str(calc)
cbind(N2_med, N3_med, N2_med+N3_med)

## figures ----
plot(calc$I2_med, calc$N2_med)
plot(calc$I3_med, calc$N3_med)
plot(seq(1999:2023), calc$I_med)

plot(calc$I2_med, jd$I2)
plot(calc$I3_med, jd$I3)


# IPM plot
year <- 1999:2021
ly <- length(year)
forecast <- 2022:2023
lf <- length(forecast)

source("IPM_fun.R")
ipm_plot(calc)


# Other plots
ggplot(jd, aes(x = lag(LD,2) , y = I2)) + geom_point()

# 2018 doesn't match but condition the best ever that year, 2021 doesn't match but survey incomplete, 2010 doesn't match but survey known to be low.keith
p <- ggplot(jd, aes(x = lag(log10(LD),2) , y = log10(I2), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p)


p <- ggplot(jd, aes(x = log10(LD) , y = lead(log10(I3), 1), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p)


p <- ggplot(jd, aes(x = log10(I2) , y = lead(log10(I3), 1), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p)


tp = out$sims.list$tau.proc
tp_med = apply(tp,2,'median') # median values of y_pred
tp_ci = apply(tp,2,'quantile', c(0.1, 0.9)) # median values of y_pred

