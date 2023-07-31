# # all of this is (L100 - L419) is to find out why i'm getting negative values which we then, can't take the log of - starting to wonder if log is a good idea
# 
# #Shaekel only ----
# # looking for values < 0 in JAGS output - abundance * maturity
# exp(out$sims.list$N[,,1])*(1-out$sims.list$m[,,1])
# str(out$sims.list$N[,,1]*(1-out$sims.list$m[,,1]))
# tmp <- exp(out$sims.list$N[,,1])*(1-out$sims.list$m[,,1])
# tmp[tmp < 0]
# 
# rownames(tmp) <- 1:length(tmp[,1])
# 
# tmp[rowSums(is.nan(tmp[,1:39]))>1,]

# problem runs
## search for problem runs
x <-1480
y <- 25
z <- 1

# the below is to try and determine why i'm getting negative values 
# tmp <- exp(out$sims.list$N[,,1])*(1-out$sims.list$m[,,1])
# tmp[1480, 25] # a negative value
# exp(out$sims.list$N[x,y,z])*(1-out$sims.list$m[x,y,z]) # N*(1-m)
# exp(out$sims.list$N[x,y,z]) # N
# 1-out$sims.list$m[x,y,z] # 1-m
# w <- exp(out$sims.list$N[x,y,z]) # N
# u <- 1-out$sims.list$m[x,y,z] # 1-m
# w*u # N*(1-m)


q <- out$sims.list$m[x,y,z] # 1-m
p <- jags.data.m$matMp[1,1] # probability
m <- log(p/(1-p)) # logit

# Beta distribution----
## OK - after a tonne of work, I realized that the dnorm distribution in JAGS was giving values of maturity <0 and >1 - not possible.  Bernoulli distribution is probability of a single trialbeing 0 or 1.  Binomial is success/trial based on p.  But neither of these give the distribution of a probability.  For that, you need the beta distribution.  And this makes JAGS give sensibile values
plot(density(rbeta(10000, p*100, (1-p)*100)))



# plots of N when multiplited by maturity 
out$sims.list$m[x,y,z]
log(exp(out$sims.list$N[x,y,z])*(1-out$sims.list$m[x,y,z]))
log(exp(out$sims.list$N[x,y,z]))
plot(density(log(exp(out$sims.list$N[,y,z])*(1-out$sims.list$m[,y,z]))))
plot(density(exp(out$sims.list$N[,y,z])*(1-out$sims.list$m[,y,z])))

# doing the same but with the actual data
exp(jags.data.m$matI[y,z])*jags.data.m$matMp[y,z]
log(exp(jags.data.m$matI[y,z])*jags.data.m$matM[y,z])
jags.data.m$matCAA[y,z]
log(jags.data.m$matCAA[y,z])
jags.data.m$matI-log(jags.data.m$matCAA)

# Survival values ---- 
## I don't know if these are "accurate" but at least are within the realm of believability or not completely stupid
plot(density(out$sims.list$s[,,1]))
plot(density(out$sims.list$s[,,2]))

plot(density(out$sims.list$s[,y,z]))
plot(density(out$sims.list$s[,y,2]))

plot(density(out$sims.list$si))  # where the fuck does this come from???????
plot(density(out$sims.list$s[,i,2]))
# plot(density(out$sims.list$sm[,1]))
# plot(density(out$sims.list$si[,2]))
# plot(density(out$sims.list$sm[,2]))

apply(out$sims.list$si, 2, 'median')
apply(out$sims.list$s, 2, 'median')
#apply(out$sims.list$sm, 2, 'median')

# calculations but with assummed mortalities - I think that this was all just to see if I could get values greater than zero
# exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + (exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1]/1000)*.5    
# 
# exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 
# (exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1]/1000)*.5
# log((exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1]/1000)*.5)


# exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + (exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1]/1000)*.5
# 
# exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]*.5
# 
# 
# log(exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + (exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]-jags.data.m$matCAA[,1]/1000)*.5)
# 
# log(exp(jags.data.m$matI[,1])*(1-jags.data.m$matM[,1])*.2 + exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]*.5)
# 
# (exp(jags.data.m$matI[,2])*jags.data.m$matM[,2]-jags.data.m$matCAA[,2]/1000)*.5
# log((exp(jags.data.m$matI[,2])*jags.data.m$matM[,2]-jags.data.m$matCAA[,2]/1000)*.5)


# # look at survival values - all ## Already looked at this above and not using logit_s anymore
# str(out$sims.list)
# #plot(density(head(out$sims.list$s[,,1])))
# plot(density(out$sims.list$s[,,1]))
# head(out$sims.list$s[,,1])
# apply(head(out$sims.list$s[,,1]), 2, 'median')
# apply(out$sims.list$s[,,1], 2, 'median')
# apply(out$sims.list$s[,,1], 2, 'mean')
# plot(density(head(out$sims.list$logit_s[,,1])))
# plot(density(out$sims.list$logit_s[,,1]))
# 
# apply(head(out$sims.list$m[,,1]), 2, 'median')
# apply(head(out$sims.list$m[,,2]), 2, 'median')
# plot(density(out$sims.list$m[,,1]))
# plot(density(out$sims.list$m[,,2]))

#plot(density(head(out$sims.list$s[,,2])))
# plot(density(out$sims.list$s[,,2]))
# apply(head(out$sims.list$s[,,2]), 2, 'median')
# apply(head(out$sims.list$s[,,2]), 2, 'mean')
# plot(density(head(out$sims.list$logit_s[,,2])))
# 
# apply(out$sims.list$alpha, 2, 'median')
# 

# Below is to plot the time series to try and figure out what is going on ito 2010??.
age2 = apply(out$sims.list$N[,,1], 2, 'median')
age3 = apply(out$sims.list$N[,,2], 2, 'median')
age4 = apply(out$sims.list$N[,,3], 2, 'median')

# cohort graphs ----
tmp2 <- as.data.frame(cbind(year = 1985:2024, N2=age2, N3 = lead(age3), N4=lead(age4, 2)))
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

# Mariano exercise - are declines exponential like we expect??
tmp <- as.data.frame(cbind(year = 1985:2024, I2 = jags.data.m$matI[,1], I3 = lead(jags.data.m$matI[,2], 1), I4 = lead(jags.data.m$matI[,3], 2)))

str(tmp)

tmp3_long <- tmp %>% 
   pivot_longer(cols = c("I2", "I3", "I4"))
str(tmp3_long)

tmp3_long$exp <- rep(NA, nrow(tmp3_long))
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

# shows change in cohort
p <- ggplot()
p <- p + geom_bar(data = tmp3_long, aes (x = name, y = value, fill = factor(name), group=1), stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("lightgoldenrod2", "darkgreen", "red"))
#p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
#p <- p + geom_line(aes(x = name, y = exp)) 
p <- p + facet_wrap(~ year)
p

# trying for the above but with immature separated from mature
## not sure if this was every successful
tmp_mat <- as.data.frame(cbind(
   year = 1985:2024, 
   I2 = log(exp(jags.data.m$matI[,1])*jags.data.m$matM[,1]), 
   I3 = lead(log(exp(jags.data.m$matI[,2]),1)*lead(jags.data.m$matM[,1]), 1),    I4 = lead(log(exp(jags.data.m$matI[,3]),1)*lead(jags.data.m$matM[,1], 2))
))

tmp_mat_long <- tmp_mat %>% # tmp_mat was tmp_imm....not sure why bc no sign of tmp_imm so this may have been a typo or I overwrote somehting
   pivot_longer(cols = c("I2", "I3", "I4"))
str(tmp_mat_long)

tmp3_long$mat <- "all"
tmp_mat_long$mat <- "mat"
tmp_mat_long$exp <- NA
tmp_test <- rbind(tmp3_long, tmp_mat_long)
tmp_test$mat <- as.factor(tmp_test$mat)
str(tmp_test)

## facet over all years
p <- ggplot(data = tmp_test, aes (x = name, y = value, fill=mat, colour = mat, alpha=mat))
p <- p + geom_col(position=position_identity())
p <- p + scale_colour_manual(values=c("lightblue", "pink"))
p <- p + scale_fill_manual(values=c("lightblue","pink"))
p <- p + theme_bw()
p <- p + geom_line(data = tmp_test[!is.na(tmp_test$exp),], aes(x = name, y = exp)) 
p <- p + facet_wrap(~ year)
p

## by year
p <- ggplot(data = tmp_test %>%
               filter(year == 1985), aes (x = name, y = value, fill=mat, colour = mat, alpha=mat))
p <- p + geom_col(position=position_identity())
p <- p + scale_colour_manual(values=c("lightblue", "pink"))
p <- p + scale_fill_manual(values=c("lightblue","pink"))
p <- p + theme_bw()
p <- p + geom_line(data = tmp_test[!is.na(tmp_test$exp),], aes(x = name, y = exp)) 
#p <- p + facet_wrap(~ year)
p

#p <- ggplot(data = tmp_test[!is.na(tmp_test$exp),], aes(x = name, y = exp))
#p <- p + geom_line()
#p


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

tmp2 <- as.data.frame(cbind(year = 1985:2024, N2=age2, N3 = lead(age3), N4=lead(age4, 2)))
tmp2$immN <- log(exp(tmp2$N2)*(1-jags.data.m$matM[,1]))
tmp2$matN <- log(exp(tmp2$N2)*jags.data.m$matM[,1])

names_matI <- c()
tmp4 <- as.data.frame(jags.data.m$matI)
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

# I and N for age-2 and -3 with immatures
p <- ggplot(data = tmp6_long, aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "red", "darkgreen", "gray"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p

# I and N for age-2 and -3 with immatures
## pre-collapse
p <- ggplot(data = tmp6_long[1:36,], aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "red", "darkgreen", "gray"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p

# I and N for age-2 and -3 with immatures
## post-collapse
p <- ggplot(data = tmp6_long[37:234,], aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "red", "darkgreen", "gray"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p

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
tmp <- as.data.frame(cbind(year = 1985:2024, N2 = calc$N2, N3 = calc$N3, N4 = calc$N4))
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
tmp <- as.data.frame(cbind(year = 1985:2024, jags.data.m$matI, jags.data.m$matMp))
#tmp <- rename(tmp, I2 = V2, I3 = V3, I4 = V4, m = V5)
tmp$s2 <- lead(exp(tmp$I3), 1)/(exp(tmp$I2) - exp(tmp$I2)*tmp$mat2)
str(tmp)

plot(tmp$year, tmp$s2)
abline(h=1)

# S ----
calc$s
