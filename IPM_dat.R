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
#save <- "no"

# load data----

# aggregated data----
df_cap <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/capelin-2021.csv")

str(df_cap)
head(df_cap)
# add extra years to end the time series
df_tmp <- df_cap[1:2,]
df_tmp[, 1:8] <- NA
df_tmp$year[1:2] <- c(2020,2021)
df_tmp
df_cap <- rbind(df_cap, df_tmp)


##disaggregated data----
# units in millions
df_dis <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/age-disaggregated-2022.csv")
str(df_dis)

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
     summarise(mat = mean(prop_mat), abun = sum(abundance), biomass = sum(biomass), varA = var(abundance), varB = var(biomass))
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
df_dis_summ
df_dis_summ$year



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

#incredibly, I can't figure out how to get pivot wider to fill in the missing years with NA!!!  
df_tmp <- df_dis_tab[1:3,]
df_tmp[, 1:4] <- NA
df_tmp$year[1:3] <- c(2006, 2016, 2020)
df_tmp

# bind the blank years (NAs) with the data
df_dis_tab <- bind_rows(df_tmp, df_dis_tab) %>% 
     arrange(year)
df_dis_tab
df_dis_tab$year


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
range(df_dis_tabLog$var)
range(df_dis_tabLog$sd)


## maturity ----
df_mat <- df_dis_summ %>%
     filter(age == 2 | is.na(age))
df_mat$mat[19] <- 0.3 # this is just a place holder until we figure out what is going on.
str(df_mat)

#impute data - place holder
imp <- mean(df_mat$mat, na.rm = T) 
df_mat$mat[7] <- imp
df_mat$mat[8] <- imp
df_mat$mat[18] <- imp
df_mat$mat[22] <- imp

# this yields the same as imp but for some reason, you need the group_by() and you don't need it above.........this is the inconsistency with tidyverse that is frustrating

tmp <- df_mat %>%
     group_by() %>%
     summarise(meanMat = mean(mat, na.rm=T))
tmp 


## larval density ----
df_ld  <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/larvae2001_2021.csv"
)
str(df_ld)

# add extra years to start the time series
df_tmp <- df_ld[1:2,]
df_tmp[, 1:3] <- NA
df_tmp$SurveyYear[1:2] <- c(1999,2000)
df_tmp
df_ld <- rbind(df_tmp, df_ld)

# change column names
df_ld <- df_ld %>% rename(year = SurveyYear,
                          larvae = `Bellevue_larvae_m-3`,
                          log_larvae = `log_Bellevue_larvae_m-3`) 
df_ld$lnlarvae <- log(df_ld$larvae)
#df_ld$lnlarvae[1] <- 6.6
#df_ld$lnlarvae[2] <- 6.6


## ice ----
#Note that I added in dummy data for 2021
df_ice  <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/ice-m1-2020.csv"
)
str(df_ice)

# get rid of years 1969-1998
df_ice <- df_ice %>%
     slice(31:54)


## condition ----
#Note that I added in dummy data for 2021
df_con <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/analyses/capelinLRP/data/condition_ag1_2_MF_out.csv"
)
str(df_con)
df_con <- df_con %>%
        slice(5:24)

df_tmp <- df_con[1:3,]
df_tmp[, 1:2] <- NA
df_tmp$year[1:3] <- c(2019,2020, 2021)
df_tmp
imp <- mean(df_con$meanCond, na.rm = T) 
df_tmp$meanCond[1:3] <- imp
df_con <- rbind(df_con, df_tmp)

# Bundle data----
num_forecasts = 2 # 2 extra years
jags.data <- ls_jag("yes", "yes")
jd <- as.data.frame(jags.data)

