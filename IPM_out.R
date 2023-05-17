# The purpose of this file is to take the bundled data from IPM_dat.R and the models from IPM_mod.R to make an age-structured, state-space model for capelin.  Model outputs are fed to IPM_out_diag.Rmd which is run through the IPM_master.R file which renders the output, i.e., makes all the figures.

## IPM_JAGS_settings is a helper file to indicate what model and what paramaters JAGS should estimate.

# Start----
# required packages
library(rjags)
library(R2jags)
library(ggplot2)
library(lattice)


# Source files
source("IPM_dat.R")
source("IPM_fun.R")
source("IPM_mod.R")
source('C:/Users/lewiske/Documents/R/zuur_rcode/MCMCSupportHighstatV2.R')
source('C:/Users/lewiske/Documents/R/zuur_rcode/HighstatLibV7.R')

# JAGS settings ----
# model - the value of b will determine what model and parameters from IPM_JAGS-settings.R
b <- 10
smoother <- "no"
matrix <- "yes"

source("IPM_JAGS-settings.R")


# MCMC settings
#ni <- 1000; nt <- 6; nb <- 50; nc <- 3
 ni <- 20000; nt <- 6; nb <- 5000; nc <- 3
# ni <- 200000; nt <- 60; nb <- 30000; nc <- 3
# ni <- 2000000; nt <- 600; nb <- 300000; nc <- 3
# ni <- 5000000; nt <- 1000; nb <- 300000; nc <- 3 # this produces really nice ACFs!!!!

# these are just preliminary values for p (auto-regression: AR) and q (MA: moving average)
if(smoother == "yes"){
    jags.data.m$p <- 2
    jags.data.m$q <- 2
}


# set index for loops
jags.data.m$Ni <- 3 # ages - N is abundance and i is the index for age
jags.data.m$M <- 3 # maturity in matrix matM - this may not be needed

# these are values to make the JAGS code more generalized, i.e., that the indices are not hard coded.  Currently applies only to cap.v20.
# jags.data.m$n.occasions <- 6
# jags.data.m$matI <- jags.data.m$matI[1:6, 1:3]
# jags.data.m$m <- jags.data.m$m[1:6]
# jags.data.m$LD <- jags.data.m$LD[1:6]
# jags.data.m$TI <- jags.data.m$TI[1:6]
# jags.data.m$CO <- jags.data.m$CO[1:6]
# year <- 1985:1995

# This was an attempt to make the start and end dates generalizable - not sure how effective it was,,,,
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
    #out$sims.list$N <- NULL
    out$sims.list$mu2 <- out$sims.list$mu[,,1]
    #out$sims.list$mu3 <- out$sims.list$mu[,,2]
    #out$sims.list$mu4 <- out$sims.list$mu[,,3]
    # out$sims.list$eps2 <- out$sims.list$eps[,,1]
    # out$sims.list$eps3 <- out$sims.list$eps[,,2]
    # out$sims.list$eps4 <- out$sims.list$eps[,,3]
    # out$sims.list$osa2 <- out$sims.list$osa[,,1]
    # out$sims.list$osa3 <- out$sims.list$osa[,,2]
    # out$sims.list$osa4 <- out$sims.list$osa[,,3]
    # out$sims.list$osa_sd2 <- out$sims.list$osa_sd[,,1]
    # out$sims.list$osa_sd3 <- out$sims.list$osa_sd[,,2]
    # out$sims.list$osa_sd4 <- out$sims.list$osa_sd[,,3]
}

# all of this is (L100 - L419) is to find out why i'm getting negative values which we then, can't take the log of - starting to wonder if log is a good idea

#Shaekel only
# looking for values < 0 in JAGS output - abundance * maturity
exp(out$sims.list$N[,,1])*(1-out$sims.list$m[,,1])
str(out$sims.list$N[,,1]*(1-out$sims.list$m[,,1]))
tmp <- exp(out$sims.list$N[,,1])*(1-out$sims.list$m[,,1])
tmp[tmp < 0]

# looking for values < 0 in JAGS output - abundance * maturity - catch (not in all models)
exp(out$sims.list$N[,,1])*out$sims.list$m[,,1] - out$sims.list$C[,,1]
str(out$sims.list$N[,,1]*(out$sims.list$m[,,1]))
tmp <- exp(out$sims.list$N[,,1])*out$sims.list$m[,,1] - out$sims.list$C[,,1]
tmp <- log(exp(out$sims.list$N[,,1])*out$sims.list$m[,,1]- out$sims.list$C[,,1])
str(tmp)
rownames(tmp) <- 1:length(tmp[,1])

tmp[rowSums(is.nan(tmp[,1:39]))>1,]

# search for problem runs
x <-1480
y <- 25
z <- 1

exp(out$sims.list$N[x,y,z])*out$sims.list$m[x,y,z] - out$sims.list$C[x,y,z]
exp(out$sims.list$N[x,y,z])*out$sims.list$m[x,y,z]
exp(out$sims.list$N[x,y,z])
out$sims.list$m[x,y,z]
log(exp(out$sims.list$N[x,y,z])*out$sims.list$m[x,y,z])
log(exp(out$sims.list$N[x,y,z]))
log(out$sims.list$C[x,y,z])
plot(density(log(exp(out$sims.list$N[,y,z])*out$sims.list$m[,y,z])))

# doing the same but with the actual data
exp(jags.data.m$matI[25,1])*jags.data.m$matM[25,1]
log(exp(jags.data.m$matI[25,1])*jags.data.m$matM[25,1])
jags.data.m$matCAA[25,1]
log(jags.data.m$matCAA[25,1])
jags.data.m$matI-log(jags.data.m$matCAA)


plot(density(out$sims.list$si[,1]))
plot(density(out$sims.list$sm[,1]))
plot(density(out$sims.list$si[,2]))
plot(density(out$sims.list$sm[,2]))

apply(out$sims.list$si, 2, 'median')
apply(out$sims.list$sm, 2, 'median')

# calculations but with assummed mortalities
exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + (exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1])*.5    

exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 
(exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1])*.5
log((exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1])*.5)


exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + (exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1])*.5

exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]*.5


log(exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + (exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1])*.5)

log(exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]*.5)

(exp(jags.data.m$matI[,2])*jags.data.m$matM[,2]-jags.data.m$matCAA[,2])*.5
log((exp(jags.data.m$matI[,2])*jags.data.m$matM[,2]-jags.data.m$matCAA[,2])*.5)


# look at survival values - all 
str(out$sims.list)
#plot(density(head(out$sims.list$s[,,1])))
plot(density(out$sims.list$s[,,1]))
head(out$sims.list$s[,,1])
apply(head(out$sims.list$s[,,1]), 2, 'median')
apply(out$sims.list$s[,,1], 2, 'median')
apply(out$sims.list$s[,,1], 2, 'mean')
plot(density(head(out$sims.list$logit_s[,,1])))
plot(density(out$sims.list$logit_s[,,1]))

apply(head(out$sims.list$m[,,1]), 2, 'median')
apply(head(out$sims.list$m[,,2]), 2, 'median')
plot(density(out$sims.list$m[,,1]))
plot(density(out$sims.list$m[,,2]))

#plot(density(head(out$sims.list$s[,,2])))
plot(density(out$sims.list$s[,,2]))
apply(head(out$sims.list$s[,,2]), 2, 'median')
apply(head(out$sims.list$s[,,2]), 2, 'mean')
plot(density(head(out$sims.list$logit_s[,,2])))

apply(out$sims.list$alpha, 2, 'median')


age2 = apply(out$sims.list$N[,,1], 2, 'median')
age3 = apply(out$sims.list$N[,,2], 2, 'median')
age4 = apply(out$sims.list$N[,,3], 2, 'median')

# cohort graphs ----
tmp2 <- as.data.frame(cbind(year = 1985:2023, N2=age2, N3 = lead(age3), N4=lead(age4, 2)))
tmp2$immN <- log(exp(tmp2$N2)*(1-jags.data.m$matM[,1]))
tmp2$matN <- log(exp(tmp2$N2)*jags.data.m$matM[,1])

tmp2_long <- pivot_longer(tmp2, cols = c("N2", "matN", "immN", "N3", "N4"))
level_order <- c("N2", "matN", "immN", "N3", "N4")
tmp2_long$name <- factor(tmp2_long$name, levels=level_order)

# whole time series
p <- ggplot(data = tmp2_long, aes (x = year, y = value, fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "darkgreen", "red"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p


# post collapse
p <- ggplot(data = tmp2_long[31:195,], aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "darkgreen", "red"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p

# post collapse; N2-N4 only
p <- ggplot(data = tmp2_long[31:195,] %>% filter(name == "N2" | name == "N3" | name == "N4"), aes (x = year, y = value, fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("lightgoldenrod2", "darkgreen", "red"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p

# Mariano exercise
tmp <- as.data.frame(cbind(year = 1985:2023, I2 = jags.data.m$matI[,1], I3 = lead(jags.data.m$matI[,2], 1), I4 = lead(jags.data.m$matI[,3], 2)))

str(tmp)

tmp3_long <- tmp %>% 
     pivot_longer(cols = c("I2", "I3", "I4"))
str(tmp3_long)

tmp3_long$exp <- rep(NA, 117)
M <- 0.6
N_ <- NA
for (i in seq_along(tmp3_long$exp)){
     if(tmp3_long$name[i] == "I2"){
          N_[i] <- tmp3_long$value[i]
     } else if (tmp3_long$name[i] == "I3") {
          N_[i] <- tmp3_long$value[i-1]*exp(-M)
     } else{
          N_[i] <- tmp3_long$value[i-2]*exp(-M*2)
     }
     tmp3_long$exp[i] <- N_[i]
}

str(tmp3_long)
head(tmp3_long, 20)

p <- ggplot()
p <- p + geom_bar(data = tmp3_long, aes (x = name, y = value, fill = factor(name), group=1), stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("lightgoldenrod2", "darkgreen", "red"))
#p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
#p <- p + geom_line(aes(x = name, y = exp)) 
p <- p + facet_wrap(~ year)
p

# trying for the above but with immature sepearted from amture
tmp_mat <- as.data.frame(cbind(
   year = 1985:2023, 
   I2 = jags.data.m$matI[,1]*jags.data.m$matM[,1], 
   I3 = lead(jags.data.m$matI[,2],1)*lead(jags.data.m$matM[,1], 1), 
   I4 = lead(jags.data.m$matI[,3],1)*lead(jags.data.m$matM[,1], 2)))

tmp_mat_long <- tmp_imm %>% 
   pivot_longer(cols = c("I2", "I3", "I4"))
str(tmp_mat_long)

tmp3_long$mat <- "all"
tmp_mat_long$mat <- "mat"
tmp_mat_long$exp <- NA
tmp_test <- rbind(tmp3_long, tmp_mat_long)
tmp_test$mat <- as.factor(tmp_test$mat)
str(tmp_test)

p <- ggplot(data = tmp_test, aes (x = name, y = value, fill=mat, colour = mat, alpha=mat))
p <- p + geom_col(position=position_identity())
p <- p + scale_colour_manual(values=c("lightblue", "pink"))
p <- p + scale_fill_manual(values=c("lightblue","pink"))
p <- p + theme_bw()
p <- p + geom_line(data = tmp_test[!is.na(tmp_test$exp),], aes(x = name, y = exp)) 
p <- p + facet_wrap(~ year)
p

p <- ggplot(data = tmp_test %>%
          filter(year == 1985), aes (x = name, y = value, fill=mat, colour = mat, alpha=mat))
p <- p + geom_col(position=position_identity())
p <- p + scale_colour_manual(values=c("lightblue", "pink"))
p <- p + scale_fill_manual(values=c("lightblue","pink"))
p <- p + theme_bw()
p <- p + geom_line(data = tmp_test[!is.na(tmp_test$exp),], aes(x = name, y = exp)) 
#p <- p + facet_wrap(~ year)
p

p <- ggplot(data = tmp_test[!is.na(tmp_test$exp),], aes(x = name, y = exp))
p <- p + geom_line()
p


p <- ggplot()
p <- p + geom_col(data = tmp_test[tmp_test$mat == "all",] %>%
                     filter(year == 1985), aes (x = name, y = value, fill=mat, colour = mat, alpha=mat), position=position_identity())
p <- p + geom_col(data = tmp_test[tmp_test$mat == "mat",] %>%
                     filter(year == 1985), aes (x = name, y = value, fill=mat, colour = mat, alpha=mat), position=position_identity())
p <- p + scale_colour_manual(values=c("lightblue", "pink"))
p <- p + scale_fill_manual(values=c("lightblue","pink"))
p <- p + theme_bw()
p <- p + geom_line(data = tmp_test[!is.na(tmp_test$exp),], aes(x = name, y = exp)) 
p <- p + facet_wrap(~ year)
p

p <- ggplot()
p <- p + geom_col(data = tmp_test[tmp_test$mat == "all",], aes (x = name, y = value, fill=mat, colour = mat, alpha=mat), position=position_identity())
p <- p + geom_col(data = tmp_test[tmp_test$mat == "mat",], aes (x = name, y = value, fill=mat, colour = mat, alpha=mat), position=position_identity())
p <- p + scale_colour_manual(values=c("lightblue", "pink"))
p <- p + scale_fill_manual(values=c("lightblue","pink"))
p <- p + theme_bw()
p <- p + geom_line(data = tmp_test[!is.na(tmp_test$exp),], aes(x = name, y = exp)) 
p <- p + facet_wrap(~ year)
p


# compare I and N----

tmp2 <- as.data.frame(cbind(year = 1985:2023, N2=age2, N3 = lead(age3), N4=lead(age4, 2)))
tmp2$immN <- log(exp(tmp2$N2)*(1-jags.data.m$matM[,1]))
tmp2$matN <- log(exp(tmp2$N2)*jags.data.m$matM[,1])

names_matI <- c()
tmp4 <- as.data.frame(jags.data.m$matI) |> rename(I2 = V1, I3 = V2, I4 = V3)
tmp4$immI <- log(exp(tmp4$I2)*(1-jags.data.m$matM[,1]))
tmp4$matI <- log(exp(tmp4$I2)*jags.data.m$matM[,1])

tmp4
tmp5 <- as.data.frame(cbind(tmp2, tmp4))

tmp5_long <- pivot_longer(tmp5, cols = c("N2", "matN", "immN", "N3", "N4", "I2", "I3", "I4", "immI", "matI"))
#level_order <- c("N2", "matN", "immN", "N3", "N4", "I2", "I3", "I4", "immI", "matI")
level_order <- c("I2", "immI", "N2", "immN" , "I3", "N3")

tmp6_long <- tmp5_long %>%
   filter(name %in% level_order) 

tmp6_long$name <- factor(tmp6_long$name, levels=level_order)

p <- ggplot(data = tmp6_long, aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "red", "darkgreen", "gray"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p


p <- ggplot(data = tmp6_long[1:36,], aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "red", "darkgreen", "gray"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p


p <- ggplot(data = tmp6_long[37:234,], aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "red", "darkgreen", "gray"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p


# JAGS output ----
## extract raw values from chains
raw <- ls_out(out)
str(raw)

plot(density(raw$logit_s[,,1]/raw$gamma[1]*(raw$delta[1]/2)*(1-(raw$delta[1]/2)/raw$delta[1])))

tmp <- (exp(raw$N[,29:32,1])*(1-raw$m[,29:32,1])+ raw$Ntb[,15:18,1])*0.5
tmp1 <- log((exp(raw$N[,30:31,1])*(1-raw$m[,30:31,1])+ raw$Ntb[,16:17,1])*.5)
tmp2 <- log((exp(raw$N[,30:31,1])*(1-raw$m[,30:31,1]))*.5)
tmp <- exp(raw$N[,13,1])*(1-raw$m[,13,1])-raw$C[,13,1]
tmp <- cbind(exp(raw$N[,14,1]), raw$m[,14,1], raw$C[,14,1])
subset(tmp, tmp[,1]*tmp[,2]<tmp[,3])
str(tmp)
test <- ifelse(tmp[,1]*tmp[,2]-tmp[,3]<=1, 1, tmp[,1]*tmp[,2]-tmp[,3])
tt <- cbind(tmp, test)
tt[tt[,4] < 1]
tt[tt[,4] < 1,]
tt[,5] <- rep(NA, length(tt))
str(tt)
log(tt[,4])

plot(density(tmp[tmp<1000]))
x <- "[,3]"
raw$m[raw$m <= 0]
raw[[x]][raw[[x]] <= 0]
exp(raw[[x]][raw[[x]] <= 0])
raw$C[,13,1][raw$C[,13,1] <=0]
raw[[x]][is.na(raw[[x]])]
raw[[x]][is.nan(raw[[x]])]

tmp[tmp<0]
log(tmp)
mean(log(tmp))
var(log(tmp))
plot(density)
# this worked by taking the mu1 out of the N2 equation.  It was still calculated.  

raw$mu[raw$mu < 0]

#extract medians, credible intervals, and prediction intervals
#source("IPM_fun.R")
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


# Residuals ----
## One-step ahead resids----
## see AugerMethe 2021
### couldn't get this to run in JAGS
plot(density(out$sims.list$osa[,,1]))
str(out$sims.list$osa[,,1])

# these are equivalent.  I can calculate them in R and in JAGS but I can't get POSA
raw$osa_sd[1,5,1]
sd(raw$osa[1, 1:4,1])

# these are equivalent.  I can calculate them in R and JAGS and I can get POSA
sd(raw$osa_sd[1,1:39,1])
sd(raw$osa[1,,1])

# 
raw$posa[1,5,1]
raw$osa[1,5,1]/sd(raw$osa[1, 1:4,1]) # this is equivalent to next line but not posa
raw$osa[1,5,1]/raw$osa_sd[1,5,1]
raw$osa[1,5,1]/sd(raw$osa[1, 1:39,1]) # this produces the same as posa 


raw$osa[1,3,1]/raw$osa_sd[1,3,1]


# I think this is on the right track but I can't get it to work in JAGS
#tmpM <- matrix(NA, c(7500, 39)) # create an emptly matrix
# out$sims.list$posa <- array(NA, c(7500, 39, 3)) # create an empty array for the Pearson OSA resids
# str(out$sims.list$posa)
# str(out$sims.list$osa)
# str(out$sims.list$osa_sd)

# create POSA and fill the array
# for(i in 1:3){
#      out$sims.list$posa[,,i] <- 
#           out$sims.list$osa[,,i]/out$sims.list$osa_sd[,,i]
# }

str(out$sims.list$posa)

# calculate teh median value for the POSA
posa_med <- matrix(NA, nrow=37, ncol=3)
for(i in 1:3){
     posa_med[,i] <- apply(out$sims.list$posa[,,i], 2, 'median')     
}
str(posa_med)
posa_df <- as.data.frame(posa_med)
posa_df <- cbind(posa_df, 1987:2023)
posa_df <- posa_df %>% rename(posa2 = V1, posa3 = V2, posa4 = V3, year = '1987:2023')
str(posa_df)
head(posa_df)

# pivot the dataframe to longer so that its easier to graph
## this does the same as the lines above but much more efficiently
posa_long <- pivot_longer(posa_df, cols = c("posa2", "posa3", "posa4"), names_to = c("age"), values_to = "pres")
head(posa_long)
posa_long$sign <- NA

for(i in seq_along(posa_long$age)){
     if(posa_long$pres[i] > 0){
          posa_long$sign[i] <- "pos"
     } else {
          posa_long$sign[i] <- "neg"
     }          
}     

# p <- ggplot(data = posa_long, aes(x = year, y = age, size = pres, colour = sign))
# p <- p + geom_point()
# p
# posa <- p

# bind process value and the posa: these are for the POSA v Nx plots
posaN2 <- as.data.frame(cbind(N2 = jags.data.m$matI[3:length(yearF),1], posa = posa_med[,1]))
posaN3 <- as.data.frame(cbind(N3 = jags.data.m$matI[3:length(yearF),2], posa = posa_med[,2]))
posaN4 <- as.data.frame(cbind(N4 = jags.data.m$matI[3:length(yearF),3], posa = posa_med[,3]))


# p <- ggplot(data = posaN2, aes(x = N2, y = posa))
# p <- p + geom_point()
# p
# plot(posa_df$year, posa_df$posa2)
# qqnorm(posa_df$posa2)
# qqline(posa_df$posa2)

# qqplot
# qqnorm(calc$N2)
# qqline(calc$N2)
# 
# plot(1985:2023, calc$posa2)
# qqnorm(calc$posa2)
# qqline(calc$posa2)
# 
# acf(calc$posa2)

## raw resids ----
## Observation error
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
# p <- ggplot(data = eps_wide, aes(x = year, y = age, size = pres, colour = sign))
# p <- p + geom_point()
# p
# raw_resid <- p

# trends in observation error
## create a data frame with the years and residuals for each year
calc$year <- rep(NA, 39)
calc$year <- 1985:2023
eps_trend <- as.data.frame(cbind(calc$year, calc$eps2, calc$eps3, calc$eps4))
eps_trend <- eps_trend %>% rename(year = V1, eps2 = V2, eps3 = V3, eps4 = V4)
str(calc)
str(eps_trend)

# create a dataframe with the years and credible intervals for each year but don't match the graphs in the Trends page on the dashboards
eps2_cri <- as.data.frame(cbind(year = 1985:2023, min = cri$eps2_cri[1,], max = cri$eps2_cri[2,]))
eps3_cri <- as.data.frame(cbind(year = 1985:2023, min = cri$eps3_cri[1,], max = cri$eps3_cri[2,]))
eps4_cri <- as.data.frame(cbind(year = 1985:2023, min = cri$eps4_cri[1,], max = cri$eps4_cri[2,]))


# Test the above.  Reproduce the JAGS results in R - make sure its doing what you think it is.
# get the mean of N2 by row - creates 7500 means
N2_mean = apply(out$sims.list$N2,1,'mean')
plot(density(N2_mean))


# create a matrix full of the acoustic survey results
N2_mat <- matrix(jags.data.m$matI[,1], nrow = 7500, ncol = 39, byrow = T)
str(N2_mat)
head(N2_mat)

#subtract each N2 value from the mean
N2_diff <- N2_mat - N2_mean

# get the median value
N2_med <- apply(N2_diff, 2, 'median', (na.rm=T))
head(N2_med)

# results are equal to JAGS 
plot(N2_med, calc$eps2)

# and it doesn't equate to this crude approach but this seems to match the graphs in the Trends page on the dashboards
plot(jags.data.m$matI[,1]-calc$N2, N2_med)

plot(1:39, jags.data.m$matI[,1]-calc$N2)

## process error----
## see AugerMethe 2021 pg 29; 
### note that process error (pe) is always going to be the difference between the process output and the process equation.  In the case of Model 47, this is simple because mu2 is also N[t] which generates N[t+1].  But its probably best to generate this in JAGS but I have also done this in R to determine if there is any equivalency
#### But these are normally distributed around 0 with a variance tau.proc.


#source("IPM_fun.R")
#### First, do this in R - subtract the N2s from the mu's
pe2 <- out$sims.list$N2[,2:39] - out$sims.list$mu2[,1:38]
str(pe2)


# N2 - get the median and Credibile intervals and plot
pe2_med <- apply(pe2, 2, 'median')
pe2_cri <- apply(pe2, 2, 'quantile', c(0.025, 0.975))
pe2_plot <- pe_plot(df_med = pe2_med, df_cri = pe2_cri, df_pri = pe2_cri)
pe2_plot

# N3 - get the median and Credibile intervals and plot
# pe3 <- out$sims.list$N3[,2:39] - out$sims.list$mu3[,1:38]
# str(pe3)
# pe3_med <- apply(pe3, 2, 'median')
# pe3_cri <- apply(pe3, 2, 'quantile', c(0.025, 0.975))
# pe3_plot <- pe_plot(df_med = pe3_med, df_cri = pe3_cri, df_pri = pe3_cri)
# pe3_plot
# 
# # N3 - get the median and Credibile intervals and plot
# pe4 <- out$sims.list$N4[,2:39] - out$sims.list$mu4[,1:38]
# str(pe4)
# pe4_med <- apply(pe4, 2, 'median')
# pe4_cri <- apply(pe4, 2, 'quantile', c(0.025, 0.975))
# pe4_plot <- pe_plot(df_med = pe4_med, df_cri = pe4_cri, df_pri = pe4_cri)
# pe4_plot


# Bubble plot for R method
# year_pe <- 1986:2023
# out$sims.list$N4[,2:39]
# pe_trend <- as.data.frame(cbind(year_pe, pe2_med, pe3_med, pe4_med))
# str(pe_trend)
# pe_trend_long <- pivot_longer(pe_trend, cols = c("pe2_med", "pe3_med", "pe4_med"), names_to = "age", values_to = "pe")
# str(pe_trend_long)
# 
# ## create positive and negative colours
# pe_trend_long$sign <- NA
# for(i in seq_along(pe_trend_long$pe)){
#      if(pe_trend_long$pe[i] > 0){
#           pe_trend_long$sign[i] <- "pos"
#      } else {
#           pe_trend_long$sign[i] <- "neg"
#      }
# }
# 
# p <- ggplot(data = pe_trend_long, aes(x = year_pe, y = age, size = pe, colour = sign))
# p <- p + geom_point()
# p

## Process error----
###JAGS approach: Should be equivalent to the above but easier to get - changes can be made in JAGS models and not in R.

str(calc$pe)
year_pe1 <- 1985:2023
pe_trend1 <- as.data.frame(cbind(year_pe1, calc$pe[,1], calc$pe[,2], calc$pe[,3]))
pe_trend1 <- pe_trend1 %>% rename(year = year_pe1, pe2 = V2, pe3 = V3, pe4 = V4)
str(pe_trend1)
pe_trend_long1 <- pivot_longer(pe_trend1, cols = c("pe2", "pe3", "pe4"), names_to = "age", values_to = "pe")
str(pe_trend_long1)

## create positive and negative colours
pe_trend_long1$sign <- NA
for(i in seq_along(pe_trend_long1$pe)){
     if(pe_trend_long1$pe[i] > 0){
          pe_trend_long1$sign[i] <- "pos"
     } else {
          pe_trend_long1$sign[i] <- "neg"
     }
}

# this is the same as the R equivalent.  
p <- ggplot(data = pe_trend_long1, aes(x = year, y = age, size = pe, colour = sign))
p <- p + geom_point()
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
#dev.off()
mix_vars_Nyear <- MyBUGSChains(out, vars_Nyear)
#ggsave(MyBUGSChains(out, vars2), filename = paste0("Bayesian/", filepath, "/chains-demographic.pdf"), width=10, height=8, units="in")


### N2 ----
#MyBUGSChains(out, vars_N2)
mix_N2 <- MyBUGSChains(out, vars_N2)
#ggsave(MyBUGSChains(out, vars1), filename = paste0("Bayesian/", filepath, "/chains-forecast.pdf"), width=10, height=8, units="in")

### N3 ----
#MyBUGSChains(out, vars_N3)
if(b ==1 | b==2 | b==10){
    mix_N3 <- MyBUGSChains(out, vars_N3)    
}

#ggsave(MyBUGSChains(out, vars1), filename = paste0("Bayesian/", filepath, "/chains-forecast.pdf"), width=10, height=8, units="in")

### N4 ----
#MyBUGSChains(out, vars_N4)
if(b == 1| b==10){
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
if(b ==1 | b ==2 | b==10){
    autocorr_N3 <- MyBUGSACF(out, vars_N3)    
}

#ggsave(MyBUGSACF(out, vars3), filename = paste0("Bayesian/", filepath, "/auto_corr-autocorrelation.pdf"), width=10, height=8, units="in")

### N4 ----
#MyBUGSACF(out, vars_N4)
if (b==1| b==10){
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
#hist(out$sims.list$Tturn.rep, xlim = c(low, high), xlab = "Number of switches \n (replicated data)")
#abline(v= out$mean$Tturn.obs, col = "red")

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
#### May not need this anymore because its done in JAGS
# if (b == 1|b == 6|b == 5){
#     resN2 <- raw$N2[,i] - raw$mu2
#     resN3 <- raw$N3[,i] - raw$mu3
#     resN4 <- raw$N4[,i] - raw$mu4
# } else if (b == 2){
#     resN2 <- raw$N2[,i] - raw$mu2
#     resN3 <- raw$N3[,i] - raw$mu3
# } else if (b == 3){
#     resN2 <- raw$N2[i] - raw$mu2
# } 
# 
# #### Pearson residuals
# presN2 <- sweep(resN2, 1, raw$tau.obs, "/")  # I do not understand why the "1" works as this indicates rowwise division but it seems to work based on the work below, i.e. change values of z and w to manually get same results as presN2 
# if (b ==1|b==2|b==6|b==5){
#     presN3 <- sweep(resN3, 1, raw$tau.obs, "/")    
# }
# 
# if (b==1|b==6|b==5){
#     presN4 <- sweep(resN4, 1, raw$tau.obs, "/")    
# }
# 
# 
# 
# str(resN2)
# str(raw$tau.obs)
# str(as.vector(raw$tau.obs))
# str(presN2)
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
#plot(x=calc$N2[i], y = resN2_mean, xlab = "Fitted values", ylab = "Raw residuals")
plot(x=calc$N2[i], y = eps$eps2, xlab = "Fitted values", ylab = "Raw residuals")

abline(h = 0, lty = 2)
# # see notes in Mortality model
plot(y = jd$I2, x = calc$N2, xlab = "Fitted values", ylab = "Observed data") # should follow the line
abline(coef = c(0,1), lty = 2)
#normality
#histogram(resN2_mean)
histogram(eps$eps2)

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
plot(y = posa_df$posa2^2, x = N2xaxis[3:length(yearF)], xlab = "Observation", ylab = "Cook's D") # should follow the line

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
#plot(jags.data$LD[LDxaxis], resN2_mean[LDxaxis], xlab = "Larval Density", ylab = "Pearson resids")
plot(jags.data$LD[LDxaxis], posa_df$posa2[LDxaxis], xlab = "Larval Density", ylab = "Pearson resids")

#plot(jags.data$TI[TIxaxis], resN2_mean[TIxaxis], xlab = "Ice retreat", ylab = "Pearson resids")
plot(jags.data$TI[TIxaxis], posa_df$posa2[TIxaxis], xlab = "Ice retreat", ylab = "Pearson resids")

par(mfrow = c(1,1))
dev.off()


# Posteriors & Priors ----
#alpha

#a2 <- post_param(param = "alpha", priormean = 0, priorsd = 100, jags = out$sims.list$alpha2, x_label = "Intercept - alpha") 

a2 <- post_param(param = "alpha[1]", priormean = 0, priorsd = 100, jags = out$sims.list$alpha[,1], x_label = "Intercept - alpha") 

# pa2 <- postPriors(df = a2$jags, df2 = a2$prior, df3 = a2$df_cred, limits=a2$limits, x_label=a2$x_label, priormean=a2$priormean, priorsd=a2$priorsd, by_bin = a2$bin_1)
# 
# pa2

#beta

b2 <- post_param(param = "beta[1]", priormean = 0, priorsd = 100, jags = out$sims.list$beta[,1], x_label = "Larval Density - beta") 

# pb2 <- postPriors(df = b2$jags, df2 = b2$prior, df3 = b2$df_cred, limits = b2$limits, x_label=b2$x_label, priormean=b2$priormean, priorsd=b2$priorsd, by_bin = a2$bin_1)

#gamma - MRI
#source("IPM_fun.R")
g2 <- post_param(param = "gamma[1]", jags = out$sims.list$gamma[,1]) # max rate of increase

pg2 <- postPriors(df = g2$jags, df2 = g2$prior, df3 =           g2$df_cred, limits=g2$limits, x_label=g2$x_label, by_bin =    g2$bin_1)
# 
 pg2

#delta - width
d2 <- post_param(param = "delta[1]", priormean = 11.5, priorsd = 5.7, jags = out$sims.list$delta[,1]) 

# pd2 <- postPriors(df = d2$jags, df2 = d2$prior, df3 = d2$df_cred, limits=d2$limits, x_label=d2$x_label, priormean=d2$priormean, priorsd=d2$priorsd, by_bin = d2$bin_1)
# 
# pd2


# process error-----
# this is for a process error that is constant within time periods.
# this cuts the array on the "a" which is the tau.proc for different ages. Then the apply cuts it by the z which is the time series.
# tp1 <- apply(out$sims.list$tau.proc[,1,],2,'median')
# tp2 <- apply(out$sims.list$tau.proc[,2,],2,'median')
# cri_tp1 <- apply(out$sims.list$tau.proc[,1,],2,'quantile', c(0.025, 0.975))
# cri_tp2 <- apply(out$sims.list$tau.proc[,2,],2,'quantile', c(0.025, 0.975))

# plot(1:4, cri$tau.proc_cri)
# 
# # these are the estiamtes for a fully time varying process error
# tp1 <- apply(out$sims.list$tau.proc[,,1],2,'median')
# tp2 <- apply(out$sims.list$tau.proc[,,2],2,'median')
# cri_tp1 <- apply(out$sims.list$tau.proc[,,1],2,'quantile', c(0.025, 0.975))
# cri_tp2 <- apply(out$sims.list$tau.proc[,,2],2,'quantile', c(0.025, 0.975))
# 
# 
# # remember that tau is the precision - t
# plot(1985:2023, 1/tp1, type= 'n', ylim = c(0,500))
# lines(1985:2023, 1/cri_tp1[2,])
# lines(1985:2023, 1/cri_tp1[1,])
# 
# plot(1985:2023, 1/tp2, type= 'n', ylim = c(0,500))
# lines(1985:2023, 1/cri_tp2[2,])
# lines(1985:2023, 1/cri_tp2[1,])
# 
# # trying the same with the variance
# s2 <- apply(out$sims.list$sigma2.proc[,,1],2,'median')
# s3 <- apply(out$sims.list$sigma2.proc[,,2],2,'median')
# cri_s2 <- apply(out$sims.list$sigma2.proc[,,1],2,'quantile', c(0.025, 0.975))
# cri_s3 <- apply(out$sims.list$sigma2.proc[,,2],2,'quantile', c(0.025, 0.975))
# 
# plot(1985:2023, s2, type= 'n', ylim = c(0,300))
# lines(1985:2023, 1/cri_s2[2,])
# lines(1985:2023, 1/cri_s2[1,])
# 
# plot(1985:2023, s3, type= 'n', ylim = c(0,300))
# lines(1985:2023, 1/cri_s3[2,])
# lines(1985:2023, 1/cri_s3[1,])


#2010 problem----
## plot of the age disaggregaged data
tmp <- df_dis_summ %>%
     filter(age == 2| age == 3| age == 4)

p <- ggplot(data = tmp, aes(x = year, y = log(abun), colour = age))
p <- p + geom_point()
p <- p + scale_colour_manual(values = c("black", "red", "pink"))
p
ggsave("2010_problem.png", width = 7, height = 5)

# as above but from 1985
tmp <- pivot_longer(df_dis_tabLog[,1:4], cols = c("I2", "I3", "I4"), names_to = "age", values_to = "Ln_abund")
p <- ggplot(data = tmp, aes(x = year, y = Ln_abund, colour = age))
p <- p + geom_point()
p <- p + scale_colour_manual(values = c("black", "red", "pink"))
p


# plot the equivalent age disaggregated values from teh process equation
# 1999-pres
tmp <- as.data.frame(cbind(year = 1985:2023, N2 = calc$N2, N3 = calc$N3, N4 = calc$N4))
tmp <- pivot_longer(tmp, cols = c("N2", "N3", "N4"), names_to = c("age"), values_to = "ln_abun")

p <- ggplot(data = tmp, aes(x = year, y = ln_abun, colour = age))
p <- p + geom_point()
p <- p + scale_colour_manual(values = c("black", "red", "pink"))
p
ggsave("2010_problem_SSM.png", width = 7, height = 5)

# plot the equivalent age disaggregated values from teh process equation - 1999-present
tmp <- subset(tmp, year >=1999)

p <- ggplot(data = tmp, aes(x = year, y = ln_abun, colour = age))
p <- p + geom_point()
p <- p + scale_colour_manual(values = c("black", "red", "pink"))
p
ggsave("2010_problem_SSM_1999.png", width = 7, height = 5)


# Divya table
str(df_dis_tabLog)
tmp <- as.data.frame(cbind(year = 1985:2023, jags.data.m$matI, jags.data.m$m))
tmp <- rename(tmp, I2 = V2, I3 = V3, I4 = V4, m = V5)
tmp$s2 <- lead(exp(tmp$I3), 1)/(exp(tmp$I2) - exp(tmp$I2)*tmp$m)
str(tmp)

plot(tmp$year, tmp$s2)
abline(h=1)

# S ----
calc$s

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
