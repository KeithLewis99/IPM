## Start----
#libraries


rm(list=ls())
# load data----

# Source files
source("IPM_dat.R")
source("IPM_fun.R")

##disaggregated data----

# age breakdown by strata
p01 <- ggplot(data = df_dis_strata, aes(x = year, y = abundance, colour = stratum)) + geom_point() + facet_grid(age ~ .)
p01

# strata breakdown by age
df_dis_strata %>% filter(stratum == "A"|stratum == "B"|stratum == "C") %>%
ggplot(aes(x = year, y = abundance, colour = age)) + geom_point() + facet_grid(stratum ~ .)

df_dis_strata %>% filter(stratum == "D"|stratum == "E"|stratum == "F") %>%
     ggplot(aes(x = year, y = abundance, colour = age)) + geom_point() + facet_grid(stratum ~ .)

df_dis_strata %>% filter(stratum == "H"|stratum == "I"|stratum == "J") %>%
     ggplot(aes(x = year, y = abundance, colour = age)) + geom_point() + facet_grid(stratum ~ .)

df_dis_strata %>% filter(stratum == "K"|stratum == "L"|stratum == "M") %>%
     ggplot(aes(x = year, y = abundance, colour = age)) + geom_point() + facet_grid(stratum ~ .)

# Check amounts and percents of Unknown biomass
# TAKE HOME MESSAGE: Unknown biomass is usually <1% of biomass but in 1999 it was ~ 3%
# known biomass
check1 <- df_dis %>%
     group_by(year) %>%
     filter(age != "Unknown") %>%
     summarize(biomass = sum(biomass))
check1

#Unknown biomass
check2 <- df_dis %>%
     group_by(year) %>%
     filter(age == "Unknown") %>%
     summarize(biomass = sum(biomass))
check2 # 5.7kt of Unknown in 1999

# percent Unknown biomass
check3 <- left_join(check1, check2, by = "year") %>%  
     mutate((biomass.y/biomass.x)*100)
check3 

# these are equivalent
cbind(jd$LD, lag(jd$LD,2), jd$I2)
cbind(jd$LD, jd$LD, lead(jd$I2,2))

# this just shows that there's no real relation between any of the age classes and therefore, the Schaekell approach won't work very well.
p1 <- ggplot(data = df_dis_tab, aes(x = I2, y = lead(I3))) + geom_point()
p1
p2 <- ggplot(data = df_dis_tab, aes(x = I2, y = lead(I4))) + geom_point()
p2
p3 <- ggplot(data = df_dis_tab, aes(x = I2, y = lead(I4, 2))) + geom_point()
p3
p4 <- ggplot(data = df_dis_tab, aes(x = I2, y = I)) + geom_point()
p4
p5 <- ggplot(data = df_dis_tab, aes(x = I2, y = lead(I))) + geom_point()
p5
p8 <- ggplot(data = df_dis_tab, aes(x = I2, y = lead(I, 2))) + geom_point()
p8
p6 <- ggplot(data = df_dis_tab, aes(x = I3, y = lead(I4))) + geom_point()
p6
p7 <- ggplot(data = df_dis_tab, aes(x = I3, y = lead(I))) + geom_point()
p7

# I2 v I3
#p <- ggplot(jd, aes(x = log10(I2) , y = lead(log10(I3), 1), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
#ggplotly(p)

#library(ggcorrplot)
#corr <- round(cor(df_dis_tab[, 2:6], use = "complete.obs", method = c("pearson")), 2)

#ggcorrplot(corr, hc.order = T, type = "lower", lab = T)



# Just want to see the relative contribution of each age class
df_dis_summF <- df_dis_summ %>%
     filter(age !=1)

p21 <- ggplot(df_dis_summF, aes(x = year, y = abun, fill = age))
p21 <- p21 + geom_bar(stat = "identity")
p21

p22 <- ggplot(df_dis_summF, aes(x = year, y = biomass, fill = age))
p22 <- p22 + geom_bar(stat = "identity")
p22

## maturity ----
# plot maturity data
plot(df_mat$year, df_mat$mat)

p60 <- ggplot(jdy, aes(x = year, y = m, text = paste("I2: ", I2, "\n", sep = ""))) + geom_point() + theme_bw()

p61 <- ggplot(jdy, aes(x = I2, y = m, text = paste("year: ", year, "\n", sep = ""))) + geom_point() + theme_bw()
summary(lm(m ~ I2, data = jd))

p62 <- ggplot(jdy, aes(x = I3, y = m, text = paste("year: ", year, "\n", sep = ""))) + geom_point() + theme_bw()
p63 <- ggplot(jdy, aes(x = log(exp(I2) + exp(I3)), y = m, text = paste("year: ", year, "\n", sep = ""))) + geom_point() + theme_bw()

p64 <- ggplot(jdy, aes(x = m, y = lead(I3), text = paste("year: ", year, "\n", sep = ""))) + geom_point() + theme_bw()

summary(lm(lead(I3,1) ~ m, data = jdy))

## larval density ----
# graph larval data time series
ggplot(df_ld, aes(x = year, y = larvae)) + geom_point() + theme_bw()

# calculation of priors for LD
range(df_ld$larvae, na.rm = T)

# 2018 doesn't match but condition the best ever that year, 2021 doesn't match but survey incomplete, 2010 doesn't match but survey known to be low.keith
jd$year <- 1999:2023

p10 <- ggplot(jd, aes(x = lag(LD,2) , y = I2)) + geom_point()

# LD v I2
p11 <- ggplot(jdy, aes(x = lag(log10(LD),2) , y = log10(I2), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p11)

p11a <- ggplot(jd_raw, aes(x = lag(log10(LD),2) , y = log10(I2), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p11a)



# LD v I3
p12 <- ggplot(jdy, aes(x = log10(LD) , y = lead(log10(I3), 1), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p12)



# explore prior dist----

# this is the acoustic survey divided by (1,000) and then transformed to the natural log scale 
range(jdy$I2, na.rm = T)
range(jdy$I3, na.rm = T)

mean(jdy$I2, na.rm = T)
mean(jdy$I3, na.rm = T)

sd(jdy$I2, na.rm = T)
sd(jdy$I3, na.rm = T)

# informative priors
x_I2 <- rnorm(1000, 9.27, 1.03)
plot(density(x_I2)) # density approaces zero at 6 & 12


x_I3 <- rnorm(1000, 8.23, 0.95)
plot(density(x_I3)) # density approaces zero at 5 & 11.5

# uniformative priors - based on the above work, I played around with the sd until I got wider priors - using these for N2 and N3
x_I2u <- rnorm(1000, 9.27, 3.03)
plot(density(x_I2u)) # density approaces zero at 0 & 16
#So, if sd^2 = 9, then precision is 1/9

x_I3u <- rnorm(1000, 8.23, 3)
plot(density(x_I3u)) # density approaces zero at 0 & 16
#So, if sd^2 = 9, then precision is 1/9



