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
jd$year <- 1985:2023

p10 <- ggplot(jd, aes(x = lag(LD,2) , y = I2)) + geom_point()

# LD v I2
p11 <- ggplot(jdy, aes(x = lag(log10(LD),2) , y = log10(I2), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p11)

#p11a <- ggplot(jd_raw, aes(x = lag(log10(LD),2) , y = log10(I2), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
#ggplotly(p11a)



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


# abundance-at-age ----
source("C:/Users/lewiske/Documents/capelin_LRP/analyses/capelinLRP/simpleLRP_FUN.R")
df_lag <- read.csv("C:/Users/lewiske/Documents/capelin_LRP/analyses/capelinLRP/data/lag.csv")
df_lag <- df_lag[, c(2, 6:8, 10)]

# df_aaa <- jags.data.m %>%
#    #select(as.data.frame (matI)) %>%
#    bind_rows() %>%
#    select(matI, LD, TI, CO) %>%
#    mutate_if(is.numeric, round, 2)
# #   rename(a1 = matI[,1])
# df_aaa$year <- as.numeric(1985:2023)

df_aaa <- left_join(df_lag, round(df_dis_tabLog[, 1:5], 2), by = "year")

str(df_aaa)

# Scatter1(df = df_aaa, xaxis = avg_densityt_2, yaxis = I2, colour = year,
#          c1 = "Year: ", c2 = "Larval density: ", c3 = "Abundance: ",
#          xlab = "Year", ylab = "Capelin abundance (billions)",
#          filename = "figs/2-cond-rank-year.pdf", save = "no")

# Cohort-maturity----
# create data frame to show realationships between mature/imm age 2, age 3 and age 4
tmp1 <- as.data.frame(cbind(year = 1985:2021, age2 = jags.data.m$matI[1:37,1], perMat = jags.data.m$matM[1:37,1]))
tmp1$imm <- log(exp(tmp1$age2)*(1-tmp1$perMat))
tmp1$mat <- log(exp(tmp1$age2)*tmp1$perMat)
tmp <- as.data.frame(cbind(tmp1, age3 = lead(jags.data.m$matI[1:37,2]), age4 = lead(jags.data.m$matI[1:37,3], 2)))

tmp_long <- pivot_longer(tmp, cols = c("mat", "imm", "age3", "age4"))

level_order <- c("mat", "imm", "age3", "age4")
tmp_long$name <- factor(tmp_long$name, levels=level_order)

# whole time series
p <- ggplot(data = tmp_long, aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", width = 1, position=position_dodge())
p <- p + scale_fill_manual(values = c("black", "lightgoldenrod2", "darkgreen", "red"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p
#ggsave("figs/recruitment1986_2020.png", width = 7, height = 5)

tmp$check <- tmp$age2-tmp$age3


# post collapse
p <- ggplot(data = tmp_long[25:148,], aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("black", "lightgoldenrod2", "darkgreen", "red"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p
#ggsave("figs/recruitment1990_2020.png", width = 7, height = 5)


# post collapse with Recovery for maturity and q
## even with these adjustments, there are years where age 2 < age 3
tmp1$mat <- log(exp(tmp1$mat)*0.26)
tmp1$imm <- log(exp(tmp1$imm)/0.5)
tmp <- as.data.frame(cbind(tmp1, age3 = lead(jags.data.m$matI[1:37,2]), age4 = lead(jags.data.m$matI[1:37,3], 2)))
tmp_long <- pivot_longer(tmp, cols = c("mat", "imm", "age3", "age4"))
level_order <- c("mat", "imm", "age3", "age4")
tmp_long$name <- factor(tmp_long$name, levels=level_order)

p <- ggplot(data = tmp_long[25:148,], aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("black", "lightgoldenrod2", "darkgreen", "red"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p
#ggsave("figs/recruitment1990_2020_modified.png", width = 7, height = 5)


# redo but with the age2 in as well.
tmp_long <- pivot_longer(tmp, cols = c("mat", "imm", "age2", "age3", "age4"))
str(tmp_long)
tmp_long$name <- as.factor(tmp_long$name)
tmp_long$name1 <- as.factor(tmp_long$name)

level_order <- c("age2", "mat", "imm", "age3", "age4")
tmp_long$name <- factor(tmp_long$name, levels=level_order)

p <- ggplot(data = tmp_long, aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "darkgreen", "red"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p


p <- ggplot(data = tmp_long[31:185,], aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("orange", "black", "lightgoldenrod2", "darkgreen", "red"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw()
p


# Z & M Barents Sea (BS) style----
## See Notation.Rmd in C:\Users\lewiske\Documents\capelin_LRP\analyses\capelinLRP
### Have asked Ale and Hannah about getting the proper equestions here.  Both have responded but do not have the actual equations.

df_dis_tab$Z2 <- -log(lead(df_dis_tab$I3,1)/
                         (df_dis_tab$I2*(1-df_mat1_per$mat2*0.01)))
mean(df_dis_tab$Z2[df_dis_tab$Z2>0], na.rm = T)

write.csv(df_dis_tab$Z2, "M.csv")

df_dis_tab$Z3 <- -log(lead(df_dis_tab$I4,1)/
                         (df_dis_tab$I3*(1-df_mat1_per$mat3*0.01)))
mean(df_dis_tab$Z3[df_dis_tab$Z3>0], na.rm = T)

#Hilborn and Walters 11.2.2
median(-(log(lead(df_dis_tab$I3,1)) - log(df_dis_tab$I2)), na.rm = T)
sd(-(log(lead(df_dis_tab$I3,1)) - log(df_dis_tab$I2)), na.rm = T)
sd(-(log(lead(df_dis_tab$I3[-7],1)) - log(df_dis_tab$I2[-6])), na.rm = T)
median(-(log(lead(df_dis_tab$I4,1)) - log(df_dis_tab$I3)), na.rm = T)
median(-(log(lead(df_dis_tab$I5,1)) - log(df_dis_tab$I4)), na.rm = T)
median(-(log(lead(df_dis_tab$I6,1)) - log(df_dis_tab$I5)), na.rm = T)

# df_dis_tab$M <- -log((lead(df_dis_tab$I3,1) + lead(matCAA[1:37,2],1) + matCAA[1:37,1])/(df_dis_tab$I2*(1-df_mat1_per$mat2*0.01)))                    

# something is wrong here.  The matCAA is supposed to be immatures but too many negatives.  The problem is that there are basically no imature age 3 so I've modified the equation
df_dis_tab$M3 <- -log((lead(df_dis_tab$I4,1) + lead(matCAA[1:37,3],1) + matCAA[1:37,2])/(df_dis_tab$I3*(1-df_mat1_per$mat3*0.01)))  

-log((lead(df_dis_tab$I4,1) + lead(matCAA[1:37,3],1) + matCAA[1:37,2])/(df_dis_tab$I3))  

df_dis_tab$Mi <- -log((lead(df_dis_tab$I3,1) + lead(matCAA[1:37,2],1) + matCAA[1:37,1])/(df_dis_tab$I2))                    

year <- c(1985:2021)

mortTab <- df_dis_tab[, c(1, 9:12)]

# remove outliers as per BS



# mature biomass----
# biomass-at-age-1985-2012 - from FRan -> put these values in ssb_calculator and then 
df_ssb_FM <- read_csv("data/ssb.csv")
plot(df_ssb_FM$ssb, lag(df_ssb_FM$abundance, 2))


# pivot the data - longer to wider with the disaggregated abundance as columns
# get abundance/biomass by age
# these values are close to those in df_cap but not exact.  
# 1999-present
df_ssb_tab <- df_baa_filter[, c(1:2, 6)] %>%
   #  filter(age != 1 & age != 5) %>%
   filter(age != 1) %>%
   #pivot_wider(names_from = age, values_from = biomass) %>%
   pivot_wider(names_from = age, values_from = matbio) %>%
   #rename(a2 = '2', a3 = '3', a4 = '4')
   rename(ssb2 = '2', ssb3 = '3', ssb4 = '4', ssb5 = '5') %>%
   mutate(ssb = sum(c_across(starts_with("s")), na.rm = T)) %>%
   mutate(var = var(c_across(starts_with("s")), na.rm = T))  %>%
   mutate(sd = sd(c_across(starts_with("s")), na.rm = T))
df_ssb_tab

#incredibly, I can't figure out how to get pivot_wider to fill in the missing years with NA!!!  So using this crude but proven method
df_tmp <- df_ssb_tab[1:3,]
df_tmp[, 1:8] <- NA
df_tmp$year[1:3] <- c(2006, 2016, 2020)
df_tmp

# bind the blank years (NAs) with the data
df_ssb_tab <- bind_rows(df_tmp, df_ssb_tab) %>% 
   arrange(year)
df_ssb_tab$ssb6 <- NA
df_ssb_tab <- df_ssb_tab[,c(1:5, 9, 6)]
str(df_ssb_tab, give.attr = F)

srr_1998 <- left_join(df_ssb_tab, df_dis_tab, by = "year") %>%
   rename(abundance = I2) %>%
   mutate_at(vars(abundance), ~ ./ 1000)


# combine the data sets as needed.
if(disaggregated == "1985-present") {
   df_ssb <- rbind(df_ssb_FM[1:14,], srr_1998[,c(1,9,7)])
} else {
   df_ssb
} 


df_ssb$abundance_tp2 <- lead(df_ssb$abundance,2)

write.csv(df_ssb, "data/ssb_all.csv")
