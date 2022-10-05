# Data file for IPM for capelin


# NOTE THAT I AM IN A LOT OF DOUBT ABOUT THE NUMBERS, ESPECIALLY FOR PROPORTIONS WHICH DON'T MATCH THE 2020 SAR.  THE AGE DISAGGREGATED ONLY GOES BACK TO 1999!!!!!!  The age disaggregated data that I downloaded from teh database also does not correspond to the master list in caplein2021.xlsx.....so all of this needs to be vetted by Fran and Aaron before we do ANYTHING WITH IT!!!!!

# Set up a project - see the below link fir directions.
#https://happygitwithr.com/rstudio-git-github.html

# But basically:
# 1.	Set up a Git repo on GitHub.
# 2.	Create the project in R - New Project - VErsion Control - Git
# 3. type "git add -A" in the terminal
# 4.	Create a bunch of directories automatically (see below)
# 5. Copy git -ignore file

#Create a "name_dat.R" file
#put this file in the folder with the project and create the following subfolders
if(!dir.exists("archive"))dir.create("archive")
#if(!dir.exists("data"))dir.create("data") # best to keep the data centralized
if(!dir.exists("figs"))dir.create("figs") #for publication quality only
if(!dir.exists("output"))dir.create("output") # for tables and figures
if(!dir.exists("ms"))dir.create("ms") # manuscript
if(!dir.exists("report"))dir.create("report") #for rmd report
if(!dir.exists("refs"))dir.create("refs") #for rmd report


## Start----
# shouldn't need the above after the first day
#libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(purrr)


rm(list=ls())
options(dplyr.print_max = 1e9)


# Source files
source("IPM_fun.R")

# variables
#save <- "no"
disaggregated <- "1985-present" # "1999-present"
# disaggregated <- "1999-present"

# load data----

# aggregated data----
df_cap <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/capelin-2021.csv")

str(df_cap)
head(df_cap)
# add extra years to end the time series
df_tmp <- df_cap[1:4,]
df_tmp[, 1:8] <- NA
df_tmp$year[1:4] <- c(2020:2023)
df_tmp
df_cap <- rbind(df_cap, df_tmp)

df_cap$var_abun <- (df_cap$abundance_med*(log(df_cap$ab_lci)-log(df_cap$abundance_med))/1.96)^2

# priors for full data set - get mean of the values from 1985-1990 on ln scale
log(mean(c(456,239,149,409,366,464)))
log(347)
log(347000)


##disaggregated data----
# units in millions
df_dis <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/age-disaggregated-2022.csv")
str(df_dis)

# bring in the historical data - ideally, this should all be in one step
df_dis_all <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/capelin_age_disaggregate_abundance.csv")
str(df_dis_all)

# Manipulate 1999-2021 data first, then 1985-2021

# this is just to get a breakdown by age and stratum
df_dis_strata <- df_dis %>%
     group_by(year, age) %>%
     filter(age != "Unknown" & age != 1 & age !=5)

# summarize maturity and abudnance
# Note that this variance is based on variance among strata.
df_dis_summ <- df_dis %>%
     group_by(year, age) %>%
     filter(age != "Unknown") %>%
     #summarise(mat = mean(prop_mat), biomass = sum(biomass))
     summarise(mat = mean(prop_mat, na.rm = T), abun = sum(abundance), biomass = sum(biomass), varA = var(abundance, na.rm = T), varB = var(biomass, na.rm = T))

df_dis_summ  
str(df_dis_summ)

#add missing years
df_tmp <- df_dis_summ[1:3,]
df_tmp[, 1:4] <- NA
df_tmp$year[1:3] <- c(2006, 2016, 2020)
df_tmp

# bind summarized data with missing data
df_dis_summ <- bind_rows(df_tmp, df_dis_summ) %>% 
     arrange(year)


# pivot the data - longer to wider with the disaggregated abundance as columns
# get abundance/biomass by age
# these values are close to those in df_cap but not exact.  
df_dis_tab <- df_dis_summ[, c(1:2,4)] %>%
     #  filter(age != 1 & age != 5) %>%
     filter(age != 1) %>%
     #pivot_wider(names_from = age, values_from = biomass) %>%
     pivot_wider(names_from = age, values_from = abun) %>%
     #rename(a2 = '2', a3 = '3', a4 = '4')
     rename(I2 = '2', I3 = '3', I4 = '4', I5 = '5') %>%
     mutate(I = sum(c_across(starts_with("I")), na.rm = T)) %>%
     mutate(var = var(c_across(starts_with("I")), na.rm = T))  %>%
     mutate(sd = sd(c_across(starts_with("I")), na.rm = T))
df_dis_tab

#incredibly, I can't figure out how to get pivot wider to fill in the missing years with NA!!!  So using this crude but proven method
df_tmp <- df_dis_tab[1:3,]
df_tmp[, 1:4] <- NA
df_tmp$year[1:3] <- c(2006, 2016, 2020)
df_tmp

# bind the blank years (NAs) with the data
df_dis_tab <- bind_rows(df_tmp, df_dis_tab) %>% 
     arrange(year)

df_dis_1998 <- df_dis_all[1:14, c(1, 3:6, 8:9)] #remove age 1 and age 6

# rename the columns and multiply by 1000 to match the scale of the earlier data set
df_dis_1998 <- df_dis_1998 %>%
        rename("mature" = "age2PerMat", "I2" = "age2", "I3" = "age3", "I4" = "age4", "I5" = "age5") %>%
        mutate(I2 = I2*1000) %>%
        mutate(I3 = I3*1000) %>%
        mutate(I4 = I4*1000) %>%
        mutate(I5 = I5*1000)

# create a total I column and blanks for variance and SD        
df_dis_1998$I <- rowSums(df_dis_1998[, 2:5], na.rm = T)
df_dis_1998$var <- NA
df_dis_1998$sd <- NA

# combine the data sets as needed.
if(disaggregated == "1985-present") {
        df_dis_tab <- rbind(df_dis_1998[, c(1:5, 8:10)], df_dis_tab)
} else {
        df_dis_tab
} 

# pivot the data - longer to wider with the disaggregated abundance as columns
# abundance value in natural logarithms
df_dis_tabLog <- df_dis_tab %>%
     #  mutate(loga2 = log(a2), loga3 = log(a3), loga4 = log(a4)) %>%
     # select(year, loga2, loga3, loga4)
     mutate(I2 = log(I2), I3 = log(I3), I4 = log(I4), I = log(I)) %>%
     mutate(var = log(var), na.rm = T) %>%
     mutate(sd = log(sd), na.rm = T) %>%
     select(year, I2, I3, I4, I, var, sd)

df_dis_tabLog
range(df_dis_tabLog$var, na.rm = T)
range(df_dis_tabLog$sd, na.rm = T)


# imputation ----
## the below are two possible approaches to resolving the 2010 missing fish problem
# First approach: calculate the difference between A2-A3 and A3-A4 fish make a dataframe
df_23 <- df_dis_tabLog$I2[15:37] - lead(df_dis_tabLog$I3[15:37])
df_34 <- df_dis_tabLog$I3[15:37] - lead(df_dis_tabLog$I4[15:37])

df_tmp <- as.data.frame(cbind(year = 1999:2021, I2_3 = df_23, I3_4 = df_34))

# What is the mean difference in the above table
df_mean <- df_tmp %>%
        filter(year != 2010) %>%
        summarise(I23 = mean(I2_3, na.rm = T), I34 = mean(I3_4, na.rm = T))
df_mean

# extract the value for the age at year from teh above dataframe

if(disaggregated == "1985-present") {
        I2_2009 <- df_dis_tabLog$I2[26] 
        I3_2011 <- df_dis_tabLog$I3[28]
        I3_2009 <- df_dis_tabLog$I3[26]
        
} else {
        I2_2009 <- df_dis_tabLog$I2[11] 
        I3_2011 <- df_dis_tabLog$I3[13]
        I3_2009 <- df_dis_tabLog$I3[11]
}

# add or subtract the mean values from the above age at year to get estimate of missing fish
I2_2010 <- as.numeric(I3_2011 + df_mean[1])
I3_2010 <- as.numeric(I2_2009 - df_mean[1])
I4_2010 <- as.numeric(I3_2009 - df_mean[2])
       
# total estimate for 2010
I2010 <- log(exp(I3_2010) + exp(I2_2010) + exp(I4_2010))
# create oject to compare to min values below.
cbind(I2 = I2_2010, I3 = I3_2010, I4 = I4_2010, I = I2010)

# Second approach: What is the max or minimum difference in the above table.  The goal here is to mimic a really poor year bc 2010 had terrible ice and terrible condition
df_max <- df_tmp %>%
        filter(year != 2010 & year != 2012) %>%
        summarise(I23 = max(I2_3, na.rm = T), I34 = max(I3_4, na.rm = T), I23min = min(I2_3, na.rm = T))

# add or subtract the max/min values from the above age at year to get estimate of missing fish
I2_2010min <- as.numeric(I3_2011 + df_max[3])  # this is really a min to make for smallest increase possible
I3_2010min <- as.numeric(I2_2009 - df_max[1])
I4_2010min <- as.numeric(I3_2009 - df_max[2])

# total estimate for 2010
I2010min <- log(exp(I3_2010min) + exp(I2_2010min) + exp(I4_2010min))
# create oject to compare to min values below.
cbind(I2m = I2_2010min, I3m = I3_2010min, I4m = I4_2010min, Im = I2010min)


# if(disaggregated == "1985-present") {
#         df_dis_tabLog[26,2]  <- I2_2010 
#         df_dis_tabLog[26,3]  <- I3_2010 
#         df_dis_tabLog[26,4]  <- I4_2010 
# } else {
#         df_dis_tabLog[11,2]  <- I2_2010min 
#         df_dis_tabLog[11,3]  <- I3_2010min 
#         df_dis_tabLog[11,4]  <- I4_2010min 
#         
# }


# USSR data 1981-1992----

## maturity ----
df_mat <- df_dis_summ %>%
     filter(age == 2 | is.na(age))
#df_mat$mat[19] <- 0.3 # this is just a place holder until we figure out what is going on.
str(df_mat)

#impute data - place holder
imp <- mean(df_mat$mat, na.rm = T) 
#df_mat$mat[7] <- imp
df_mat$mat[8] <- imp
df_mat$mat[18] <- imp
df_mat$mat[22] <- imp

# this yields the same as imp but for some reason, you need the group_by() and you don't need it above.........this is the inconsistency with tidyverse that is frustrating

tmp <- df_mat %>%
     group_by() %>%
     summarise(meanMat = mean(mat, na.rm=T))
tmp 

# create an identical data frame to df_mat for years 1985:1998 and append them to more recent data.  


if(disaggregated == "1985-present") {
        df_tmp <- df_mat[1:14,]
        df_tmp[, c(1, 3:7)] <- NA
        df_tmp$year <- c(1985:1998)
        df_tmp$mat <- df_dis_1998$mature/100
        df_mat <- rbind(df_tmp, df_mat)
        # get a mean maturity form 1991:1999
        imp90 <- mean(df_mat$mat[7:15], na.rm = T) 
        #df_mat$mat[7] <- imp
        df_mat$mat[c(9:11, 13:14)] <- imp90
        
} else {
        df_mat
} 

# check the relationships
plot(df_mat$year, df_mat$mat)
plot(df_dis_1998$year, df_dis_1998$perAge2)


df_matM  <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/analyses/capelinLRP/data/springAcoustics-percentMature.csv")
str(df_matM)

# check data
df_matM$age2[is.na(df_matM$age2)]
is.na(df_matM$age2)
df_matM$age3[is.na(df_matM$age3)]
df_matM[,1:3]


# crude imputation for age 2 maturity
for (i in seq_along(df_matM$age2)){
     if(is.na(df_matM$age2[i])){
          df_matM$age2[i] <- mean(df_matM$age2[12:35], na.rm = T)
     }
}
# confirm above works
df_matM[,1:3]

# could also do this with density dependent approach

## larval density ----
df_ld  <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/larvae2001_2022.csv")
str(df_ld)

# add extra years to start the time series
if(disaggregated == "1985-present") {
        df_tmp <- as.data.frame(matrix(NA, 16, 3))
        df_tmp[, 1] <- c(1985:2000)
        names(df_tmp) <- names(df_ld)
        df_ld <- rbind(df_tmp, df_ld)
} else {
        df_tmp <- df_ld[1:2,]
        df_tmp[, 1:3] <- NA
        df_tmp$SurveyYear[1:2] <- c(1999,2000)
        df_tmp
        df_ld <- rbind(df_tmp, df_ld)
} 

# change column names
# df_ld <- df_ld %>% rename(year = SurveyYear,
#                           larvae = `Bellevue_larvae_m-3`,
#                           log_larvae = `log_Bellevue_larvae_m-3`)
# df_ld$lnlarvae <- log(df_ld$larvae)

df_ld <- df_ld %>% rename(year = `Year`,
                          larvae = `Larval densities_ind_m-3`,
                          se_auc = `SE_AUC`) 
df_ld$lnlarvae <- log(df_ld$larvae)
str(df_ld)


## ice ----
#Note that I added in dummy data for 2021
df_ice  <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/ice-m1-2020.csv"
)
str(df_ice)

if(disaggregated == "1985-present") {
        df_ice <- df_ice %>% # get rid of years 1969-1984
            slice(17:54)
} else {
        df_ice <- df_ice %>% # get rid of years 1969-1998
        slice(31:54)
} 



## condition ----
#Note that I added in dummy data for 2021
df_con <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/analyses/capelinLRP/data/condition_ag1_2_MF_out.csv"
)
str(df_con)

df_tmp <- df_con[1:3,]
df_tmp[, 1:2] <- NA
df_tmp$year[1:3] <- c(2019,2020, 2021)
df_tmp
imp <- mean(df_con$meanCond, na.rm = T) 
df_tmp$meanCond[1:3] <- imp
df_con <- rbind(df_con, df_tmp)

if(disaggregated == "1985-present") {
        df_tmp <- as.data.frame(matrix(NA, 10, 2))
        df_tmp[, 1] <- c(1985:1994)
        names(df_tmp) <- names(df_con)
        df_con <- rbind(df_tmp, df_con)
} else {
        df_con <- df_con %>%
                slice(5:27)
        } 




# Bundle data----
num_forecasts = 2 # 2 extra years
jags.data <- ls_jag("yes", "yes", "no")
str(jags.data)
jd <- as.data.frame(jags.data)

# get lengths of jags.data
leng_jd <- rep(NA, 8)
for (i in 1:length(jags.data)){
        len_jd <- length(jags.data[[i]])
        leng_jd[i] <- len_jd
}
leng_jd

# add years for convenience of graphing later
if(disaggregated == "1985-present") {
        yearF <- 1985:2023
        year <- 1985:2021
} else {
        yearF <- 1999:2023
        year <- 1999:2021
}

jdy <- cbind(year = yearF, jd)
str(jdy)

# jd_raw <- ls_jag("no", "no")
# jd_raw <- as.data.frame(jd_raw)
# jd_raw <- cbind(year = year, jd_raw)

#source("IPM_fun.R")
jags.data.m <- ls_jag("yes", "yes", "yes")
str(jags.data.m)




# figure ----
# create data frame to show realationships between mature/imm age 2, age 3 and age 4
tmp1 <- as.data.frame(cbind(year = 1985:2021, age2 = jags.data.m$matI[1:37,1], perMat = jags.data.m$m[1:37]))
tmp1$imm <- log(exp(tmp1$age2)*(1-tmp1$perMat))
tmp1$mat <- log(exp(tmp1$age2)*tmp1$perMat)
tmp <- as.data.frame(cbind(tmp1, age3 = lead(jags.data.m$matI[1:37,2]), age4 = lead(jags.data.m$matI[1:37,3], 2)))

tmp_long <- pivot_longer(tmp, cols = c("mat", "imm", "age3", "age4"))

level_order <- c("mat", "imm", "age3", "age4")
tmp_long$name <- factor(tmp_long$name, levels=level_order)

# whole time series
p <- ggplot(data = tmp_long, aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
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

# trying to zoom in without the 2014 outlier

#install.packages("ggforce")                   # Install & load ggforce package
library("ggforce")

p <- ggplot(data = tmp_long[25:148,], aes (x = year, y = exp(value), fill = factor(name, level_order)))
p <- p + geom_bar(stat="identity", position=position_dodge())
p <- p + scale_fill_manual(values = c("black", "lightgoldenrod2", "darkgreen", "red"))
p <- p + guides(fill=guide_legend(title="Age/Mat"))
p <- p + theme_bw() 
p <- p + facet_zoom(ylim = c(0, 20000))
p
#ggsave("figs/recruitment1986_2020_facet.png", width = 7, height = 5)

