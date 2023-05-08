# Data file for IPM for capelin
## This includes the biomass correcations of 2022 - see A. Adamack for details or 2023 Res Doc.

## Some of below is confusing because there are two data sets.  1985-1998 and 1999-present.  None of them seem to be in the same form.  So a lot of the below is bringing in two data sets and combining them.  Ideally, this should all be in one step but I leave that to Pelagics.


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
if(!dir.exists("data"))dir.create("data") # best to keep the data centralized
if(!dir.exists("data/derived"))dir.create("data/derived") # best to keep the data centralized
if(!dir.exists("figs"))dir.create("figs") #for publication quality only
if(!dir.exists("output"))dir.create("output") # for tables and figures
if(!dir.exists("ms"))dir.create("ms") # manuscript
if(!dir.exists("report"))dir.create("report") #for rmd report
if(!dir.exists("refs"))dir.create("refs") #for rmd report

# Start ----
# shouldn't need the above after the first day
#libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(purrr)

# Source files
source("IPM_fun.R")
options(dplyr.print_max = 1e9)

# variables
#save <- "no"
disaggregated <- "1985-present" # "1999-present"
# disaggregated <- "1999-present"

# load data----

## disaggregated abund/biomass age data----
## 1999-2022 age-disaggregated abundance and biomass but alsom with mature abundance and biomass.  All have lower and upper CIs.  

### Units millions
df_dis <- read_csv("data/abundance and biomass by age and year2.csv")
str(df_dis)

# bring in the historical data - 
## 1985 2017 - # in billions
df_dis_1985 <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/capelin_age_disaggregate_abundance1.csv")
str(df_dis_1985)

# Manipulate 1999-2021 data first, then 1985-2021
## 1999-2021:this is just to get a breakdown by age and stratum and to rename variables
### do not filter out Unknowns at this point
df_dis_summ <- df_dis %>%
   group_by(year, age) %>%
   select(year, age,
          abun = med.abund.age.fran, 
          abun.lci = low.abund.age.fran, abun.uci = up.abund.age.fran, 
          matabun = med.mat.abund.age.fran, 
          matabun.lci = low.mat.abund.age.fran, matabun.uci = up.mat.abund.age.fran
   ) %>% # select and rename 
   mutate_at(vars(abun:matabun.uci), ~ ./ 1000) %>% # convert to billions
   mutate_at(vars(abun:matabun.uci), round, 2) # round

df_dis_summ  
str(df_dis_summ, give.attr = FALSE)

# pivot the data - longer to wider with the disaggregated abundance as columns
# get abundance/biomass by age
# these values are close to those in df_cap but not exact.  
df_dis_tab <- df_dis_summ[, c(1:3)] %>%
     pivot_wider(names_from = age, values_from = abun) %>%
     rename(I1 = '1', I2 = '2', I3 = '3', I4 = '4', I5 = '5') %>%
   mutate(I = rowSums(across(I1:Unknown), na.rm = T)) %>%
   mutate(var = var(c_across(I1:Unknown), na.rm = T))  %>%
   mutate(sd = sd(c_across(I1:Unknown), na.rm = T))
df_dis_tab

# create a df with missing years
df_tmp <- df_dis_tab[1:3,]
df_tmp[, 1:length(df_tmp)] <- NA
df_tmp$year[1:3] <- c(2006, 2016, 2020)
df_tmp

# bind the blank years (NAs) with the data
df_dis_tab <- bind_rows(df_tmp, df_dis_tab) %>% 
     arrange(year)

# rename the columns to match the scale of the earlier data set
## coded this out because we want these to be in billions!!!!
df_dis_1998 <- df_dis_1985 %>%
        rename("mature" = "age2PerMat", I1 = "age1", "I2" = "age2", "I3" = "age3", "I4" = "age4", "I5" = "age5", "I6" = "age6") # %>%

# create a total I column and blanks for variance and SD        
df_dis_1998$I <- rowSums(df_dis_1998[, 2:7], na.rm = T)

for(i in 1:length(df_dis_1998$I)){
   if(df_dis_1998$I[i] == 0)
      df_dis_1998$I[i] <- NA
}

# create blank columns
df_dis_1998$var <- NA
df_dis_1998$sd <- NA

# combine the data sets as needed.
if(disaggregated == "1985-present") {
        df_dis_tab <- rbind(df_dis_1998[, c(1:7, 11:13)], df_dis_tab[, c(1:6, 11, 8:10)])
} else {
        df_dis_tab
} 

write.csv(df_dis_tab, "data/derived/capelin_abundance_1985-2022.csv")

# convert abundance value to natural logarithms
## NOTE THAT this is only I2-I4 because that is what the state space model deals with
df_dis_tabLog <- df_dis_tab %>%
     mutate(I2 = log(I2), I3 = log(I3), I4 = log(I4), I = log(I)) %>%
     mutate(var = log(var), na.rm = T) %>%
     mutate(sd = log(sd), na.rm = T) %>%
     select(year, I2, I3, I4, I, var, sd)

df_dis_tabLog
range(df_dis_tabLog$var, na.rm = T)
range(df_dis_tabLog$sd, na.rm = T)

#geometric mean for low productivity leading to high one.
## https://en.wikipedia.org/wiki/Geometric_mean see this for rationale...even though log numbers, you are taking the product of all numbers exponentiated to 1/n.  Or, you an do sum of ln (n_i)*1/n and exponentiate
### not sure why I did this.
exp(mean(log(df_dis_tabLog[c(15:25, 27:28), ]$I), na.rm=T))

### abun aggreg ----
## aggregate age-disagregated data >= 1999 abundance and CIs
### - note that these are 5th and 95th percentiles which is a bootstrap CI!!!!!
#### note that this step is a bit redundant, i.e., Fran/Aaron have probably already done this but I want my own so that I don't have to continutally ask for data, i.e, I just want to be able to get age-disaggregated data with maturities and run the analysis.
df_ag_1999 <- df_dis_summ[, c(1:5)] %>%
   group_by(year) %>%
   summarize(abundance_med = sum(abun), ab_lci = sum(abun.lci), ab_uci = sum(abun.uci))
 
## aggregate age-disagregated data < 1999
### - note that that I don't have the 5th and 95th percentiles for these and they only go back to 1988 (or is it 1998?)

##### BUT FRAN SEEMS TO HAVE AGGREGATED CIS TO 1988 - SO NEED TO DIG INTO THIS FURTHER***************************************************

df_ag_1985 <- df_dis_1998[, c(1, 11)] 
colnames(df_ag_1985)[2] <- "abundance_med"
df_ag_1985$ab_lci <- NA
df_ag_1985$ab_uci <- NA

# combine the 1985-1998 and 1999-present.
df_agg <- rbind(df_ag_1985, df_ag_1999)

# add years with no data
df_tmp <- df_agg[1:3,]
df_tmp[, 1:length(df_tmp)] <- NA
df_tmp$year[1:3] <- c(2006, 2016, 2020)
df_tmp

# combine wtih age-aggregagated data 1985-present
df_agg <- bind_rows(df_tmp, df_agg) %>% 
   arrange(year)

write_csv(df_agg, "data/derived/capelin_aggregated_abundance_1985-2022.csv")



## biomass disagregated ----
# biomass-at-age-1985-2012 - from FRan
## units (kt)
df_baa_FM <- read_csv("data/baa-1985-2012.csv")
df_baa_FM[,2:7] <- round(df_baa_FM[, 2:7], 1)
str(df_baa_FM, give.attr = F)


## biomass for 1982-1999

vec_bio <- c(446, NA, NA, 3426, 3697, 2576, 4285, 3712, 5783, 138, 138, NA, NA, NA, 47, NA, NA
)
length(vec_bio)
df_bio_agg_1982 <- data.frame(cbind(seq(1982, 1998), vec_bio)) %>%
   rename(year = V1, biomass = vec_bio)


# biomass-at-age 1999-present - from Aaron and teh Shiny App 
## Units - tonnes
df_baa_AA <- read_csv("data/abundance and biomass by age and year2.csv")
str(df_baa_AA, give.attr = F)

# remove Unknowns and Age-1
df_baa_filter <- df_baa_AA %>%
   group_by(year, age) %>%
   #filter(age != "Unknown") %>% 
   select(year, age,
          bio = med.bm.age.fran, 
          bio.lci = low.bm.age.fran, bio.uci = up.bm.age.fran, 
          matbio = med.mat.bm.age.fran, 
          matbio.lci = low.mat.bm.age.fran, matbio.uci = up.mat.bm.age.fran
) %>%
   mutate_at(vars(bio:matbio.uci), ~ ./ 1000) %>% # convert units to kt
   mutate_at(vars(bio:matbio.uci), funs(round(., 2)))



# pivot the data - longer to wider with the disaggregated abundance as columns
# get abundance/biomass by age
# these values are close to those in df_cap but not exact.  
df_baa_tab <- df_baa_filter[, c(1:3)] %>%
   #  filter(age != 1 & age != 5) %>%
   #filter(age != 1) %>%
   #pivot_wider(names_from = age, values_from = biomass) %>%
   pivot_wider(names_from = age, values_from = bio) %>%
   rename(bio1 = '1', bio2 = '2', bio3 = '3', bio4 = '4', bio5 = '5') %>%
   # mutate(biomass = sum(c_across(starts_with("b")), na.rm = T)) %>%
   # mutate(var = var(c_across(starts_with("b")), na.rm = T))  %>%
   # mutate(sd = sd(c_across(starts_with("b")), na.rm = T))
   mutate(biomass = rowSums(across(bio1:Unknown), na.rm = T)) %>%
   mutate(var = var(c_across(bio1:Unknown), na.rm = T))  %>%
   mutate(sd = sd(c_across(bio1:Unknown), na.rm = T))
df_baa_tab

#incredibly, I can't figure out how to get pivot wider to fill in the missing years with NA!!!  So using this crude but proven method
df_tmp <- df_baa_tab[1:3,]
df_tmp[, 1:length(df_tmp)] <- NA
df_tmp$year[1:3] <- c(2006, 2016, 2020)
df_tmp

# bind the blank years (NAs) with the data
df_baa_tab <- bind_rows(df_tmp, df_baa_tab[, c(1:8, 9:10)]) %>% 
   arrange(year)
df_baa_tab$bio6 <- NA
df_baa_tab <- df_baa_tab[,c(1:6, 11, 7:8)]
str(df_baa_tab, give.attr = F)

# create a total biomass column, get relevant rows, and fill in NAs
df_baa_FM$biomass <- NA
df_baa_FM$Unknown <- NA
df_baa_1985 <- df_baa_FM[1:14, c(1:7, 9, 8)]
df_baa_1985 <- left_join(df_baa_1985[, 1:8], df_bio_agg_1982, "year")
#df_baa_1985[c(9:11, 13:14), 8] <- NA


# combine the data sets as needed.
if(disaggregated == "1985-present") {
   df_baa_tab <- rbind(df_baa_1985, df_baa_tab)
} else {
   df_baa_tab
} 

write.csv(df_baa_tab[-9], "data/capelin_biomass_1985-2022.csv", row.names = F)

# pivot the data - longer to wider with the disaggregated abundance as columns
# abundance value in natural logarithms
df_baa_tabLog <- df_baa_tab %>%
   #  mutate(loga2 = log(a2), loga3 = log(a3), loga4 = log(a4)) %>%
   # select(year, loga2, loga3, loga4)
   mutate(B2 = log(bio2), B3 = log(bio3), B4 = log(bio4), B = log(biomass)) %>%
   #mutate(var = log(var), na.rm = T) %>%
   #mutate(sd = log(sd), na.rm = T) %>%
   select(year, B2, B3, B4, B)

df_baa_tabLog
#range(df_baa_tabLog$var, na.rm = T)
#range(df_baa_tabLog$sd, na.rm = T)

#geometric mean for low productivity leading to high one.
## https://en.wikipedia.org/wiki/Geometric_mean see this for rationale...even though log numbers, you are taking the product of all numbers exponentiated to 1/n.  Or, you an do sum of ln (n_i)*1/n and exponentiate
exp(mean(log(df_baa_tabLog[c(15:25, 27:28), ]$B), na.rm=T))

### biomass agg ----

## aggregate age-disagregated data >= 1999
### - note that these are 5th and 95th percentiles and not CIs!!!!!
df_ag_bio_1999 <- df_baa_filter[, c(1:5)] %>%
   group_by(year) %>%
   summarize(biomass_med = sum(bio), bm_lci = sum(bio.lci), bm_uci = sum(bio.uci))

## aggregate age-disagregated data < 1999
### - note that that I don't have the 5th and 95th percentiles for these and they only go back to 1988

# df_ag_bio_1985 <- df_baa_1985[, c(1, 8)] 
# colnames(df_ag_bio_1985)[2] <- "biomass_med"
# df_ag_bio_1985$bm_lci <- NA
# df_ag_bio_1985$bm_uci <- NA

colnames(df_bio_agg_1982)[2] <- "biomass_med"
df_bio_agg_1982$bm_lci <- NA
df_bio_agg_1982$bm_uci <- NA

# combine the data sets as needed.
df_agg_bio <- rbind(df_bio_agg_1982, df_ag_bio_1999)
df_tmp <- df_agg_bio[1:3,]
df_tmp[, 1:length(df_tmp)] <- NA
df_tmp$year[1:3] <- c(2006, 2016, 2020)
df_tmp

df_agg_bio <- bind_rows(df_tmp, df_agg_bio) %>% 
   arrange(year)

#write.csv(df_dis_tab, "capelin_abundance_1985-2021.csv")
write_csv(df_agg_bio, "data/capelin_aggregated_biomass_1985-2022.csv")


## USSR data 1981-1992----

## Trinity Bay ----
### 1999-2019 (update when needed)
df_tb <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/fromAaron/TB_abun_atAge.csv")
str(df_tb, give.attr = F)
head(df_tb)
df_tb$age <- as.factor(df_tb$age)
df_tb$stratum <- as.factor("TB")

df_tb_NAA <- df_tb %>%
     filter(age != "Unknown" & age != "1" & age != "5") %>%
     filter(abundance >0) %>%
     pivot_wider(id_cols = year, names_from = age, values_from = abundance, names_sort = T) %>%
     rename(I2 = '2', I3 = '3', I4 = '4') %>%
     mutate(I2 = log(I2), I3 = log(I3), I4 = log(I4))

str(df_tb_NAA)
df_tb_NAA

# add missing years
df_tmp <- df_tb_NAA[1:8,] 
df_tmp[, 1:4] <- NA
df_tmp$year[1:8] <- c(2006, 2014:2016, 2020:2023)
df_tmp

# bind summarized data with missing data
df_tb_NAA <- bind_rows(df_tmp, df_tb_NAA) %>% 
     arrange(year)
str(df_tb_NAA)

if(disaggregated == "1985-present") {
     df_tmp <- as.data.frame(matrix(NA, 14, 4))
     df_tmp[, 1] <- c(1985:1998)
     names(df_tmp) <- names(df_tb_NAA)
     df_tb_NAA <- rbind(df_tmp, df_tb_NAA)
} else {
     df_con <- df_con %>%
          slice(5:27)
} 

# impute
imp <- colMeans(df_tb_NAA[,2:4], na.rm = T)
df_tb_NAA[c(22,30:32, 36:39), 2:4] <- imp


matITB <- as.matrix(df_tb_NAA[, 2:4])

## maturity ----
## note that variable 'mat' is as a proportion, not a percent - so no need to divide by 100
## this is from the spring acoustic and Shiny app 1999-2022
df_mat <- df_dis_summ %>%
     filter(age == 2 | is.na(age))
#df_mat$mat[19] <- 0.3 # this is just a place holder until we figure out what is going on.
str(df_mat, give.attr=F)

#impute data - place holder
imp <- mean(df_mat$matabun, na.rm = T) 
#df_mat$mat[7] <- imp
# df_mat$matabun[8] <- imp
# df_mat$matabun[18] <- imp
# df_mat$matabun[22] <- imp

# this yields the same as imp but for some reason, you need the group_by() and you don't need it above.........this is the inconsistency with tidyverse that is frustrating

tmp <- df_mat %>%
     group_by() %>%
     summarise(meanMat = mean(matabun, na.rm=T))
tmp 

# create an identical data frame to df_mat for years 1985:1998 and append them to more recent data.  
if(disaggregated == "1985-present") {
        df_tmp <- df_mat[1:14,]
        df_tmp[, c(1, 3:7)] <- NA
        df_tmp$year <- c(1985:1998)
        df_tmp$mat <- df_dis_1998$mature/100 # maturity is a percentage here so divide by 100
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

## mat - diaggregated ----
### USE THIS ONE, NOT CODE UNDER maturity or mat-matrix
## 1999-present - abundance mature based on import from AA file
df_mat_tab <- df_dis_summ[, c(1:2, 6)] %>%
   #  filter(age != 1 & age != 5) %>%
   filter(age != 1) %>%
   #pivot_wider(names_from = age, values_from = biomass) %>%
   pivot_wider(names_from = age, values_from = matabun) %>%
   #rename(a2 = '2', a3 = '3', a4 = '4')
   rename(mat2 = '2', mat3 = '3', mat4 = '4', mat5 = '5') %>%
   mutate(matureAbun = sum(c_across(starts_with("m")), na.rm = T)) %>%
   mutate(var = var(c_across(starts_with("m")), na.rm = T))  %>%
   mutate(sd = sd(c_across(starts_with("m")), na.rm = T))
df_mat_tab
str(df_mat_tab, give.attr=F)

# create empty rows
df_tmp <- df_mat_tab[1:3,]
df_tmp[, 1:length(df_tmp)] <- NA
df_tmp$year[1:3] <- c(2006, 2016, 2020)
df_tmp

# bind the blank years (NAs) with the data
df_mat_tab <- bind_rows(df_tmp, df_mat_tab) %>% 
   arrange(year)

## 1985-2012 - abundance mature
### based on file from the biochar file (BIOCHAR FROM ACOUSTICS_revised to use Monte Carlo abundnace for 1988-1996 in annual page (003).xls; C:\Users\lewiske\Documents\capelin_LRP\IPM\data)  
df_mat_1985 <- read_csv("data/matAbun.csv")
str(df_mat_1985, give.attr = F)



# this is the mature abundance for age-2 to -5 and total (includes unknowns and 6s)
df_mat1 <- rbind(df_mat_1985[1:14, c(-2, -7)], df_mat_tab[,c(1:5, 7)])
str(df_mat1, give.attr = F)


df_mat1_per <- cbind(df_mat1[1], df_mat1[-1]/df_dis_tab[c(3:7)]*100)
str(df_mat1_per, give.attr = F)
write.csv(df_mat1_per, "data/capelin_perMat_1985-2022.csv", row.names = F)

## mat- matrix----
#### 1985-2019
#### These are percentages; this file is used in LRP_dat but only for age2
df_matM  <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/analyses/capelinLRP/data/springAcoustics-percentMature.csv")
str(df_matM, give.attr = F)
head(df_matM) # these are percentages
tail(df_matM)

# check data - just look for NA in one age class
df_matM$age2[is.na(df_matM$age2)]
is.na(df_matM$age2)
df_matM$age3[is.na(df_matM$age3)]
df_matM[,1:3]

df_tmp <- df_matM[1:4,]
df_tmp[, 1:7] <- NA
df_tmp$year[1:4] <- c(2020:2023)

df_matM <- rbind(df_matM, df_tmp)

# just so that its not a zero leading to no age4 fish
df_matM[19,4] <- 99

# crude imputation for age 2 maturity - this is weakest for age-2, a bit better for age 3, and probably a good approximation for age-4.  Note, no NA in pre-collapse period
for (i in seq_along(df_matM$age2)){
     if(is.na(df_matM$age2[i])){
          df_matM$age2[i] <- mean(df_matM$age2[12:35], na.rm = T)
     }
}

for (i in seq_along(df_matM$age3)){
     if(is.na(df_matM$age3[i])){
          df_matM$age3[i] <- mean(df_matM$age3[12:35], na.rm = T)
     }
}

for (i in seq_along(df_matM$age4)){
     if(is.na(df_matM$age4[i])){
          df_matM$age4[i] <- mean(df_matM$age4[12:35], na.rm = T)
     }
}


# confirm above works
df_matM[,1:3]
df_matM[,1:5]

m_matM <- as.matrix(cbind(df_matM$age2/100, df_matM$age3/100, df_matM$age4/100))


# could also do this with density dependent approach

## TB maturity ----
df_tb_matAA <- df_tb %>%
     filter(age != "Unknown" & age != "1" & age != "5") %>%
     pivot_wider(id_cols = year, names_from = age, values_from = prop_mat, names_sort = T) %>%
     rename(m2 = '2', m3 = '3', m4 = '4')

# add missing years
df_tmp <- df_tb_matAA[1:8,] 
df_tmp[, 1:4] <- NA
df_tmp$year[1:8] <- c(2006, 2014:2016, 2020:2023)
df_tmp


# bind summarized data with missing data
df_tb_matAA <- bind_rows(df_tmp, df_tb_matAA) %>% 
     arrange(year)
str(df_tb_matAA)

if(disaggregated == "1985-present") {
     df_tmp <- as.data.frame(matrix(NA, 14, 4))
     df_tmp[, 1] <- c(1985:1998)
     names(df_tmp) <- names(df_tb_matAA)
     df_tb_matAA <- rbind(df_tmp, df_tb_matAA)
} else {
     df_con <- df_con %>%
          slice(5:27)
} 


# impute
imp <- colMeans(df_tb_matAA[,2:4], na.rm = T)
df_tb_matAA[c(22,30:32, 36:39), 2:4] <- imp

# change '1' to high percentage
df_tb_matAA[,3][df_tb_matAA[,3] == 1] <- 0.995
     
m_maaTB <- as.matrix(df_tb_matAA[,2:4])


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
#df_con <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/analyses/capelinLRP/data/condition_ag1_2_MF_out.csv")
df_con <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/fromAaron/condition_2JK_ag1_2_MF_2022.csv")
str(df_con)

# df_tmp <- df_con[1:2,]
# df_tmp[, 1:2] <- NA
# df_tmp$year[1:2] <- c(2022:2023)
# df_tmp
# imp <- mean(df_con$meanCond, na.rm = T) 
# df_tmp$meanCond[1:2] <- imp
# df_con <- rbind(df_con, df_tmp)

if(disaggregated == "1985-present") {
        df_tmp <- as.data.frame(matrix(NA, 10, 2))
        df_tmp[, 1] <- c(1985:1994)
        names(df_tmp) <- names(df_con)
        df_con <- rbind(df_tmp, df_con)
} else {
        df_con <- df_con %>%
                slice(5:27)
        } 
str(df_con)

## catch-at-age----
# CAA 1998-2021
df_caa <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/fromAaron/catchAtAge1998_2021.csv")
str(df_caa, give.attr = FALSE)
df_caa

tmp <- df_caa %>%
     group_by(year, age) %>%
     filter(age != "Unknown") %>%
     summarise(abundance_sum = sum(N_millions), prop_mat_mean = mean(prop_mat))
tmp$age <- as.numeric(tmp$age)
str(tmp, give.attr = FALSE)

df_caa_tab_abun <- tmp[c("year", "age", "abundance_sum")] %>%
     pivot_wider(id_cols = year, names_from = age, values_from = abundance_sum, names_sort = T) %>%
      rename(c1 = '1', c2 = '2', c3 = '3', c4 = '4', c5 = '5', c6 = '6') %>%
     select(year, c2, c3, c4)

df_caa_tab_abun

# insert 2022 and 2023 as average of last 10 years
df_tmp <- as.data.frame(matrix(NA, 2, 4))
df_tmp[, 1] <- c(2022, 2023)
names(df_tmp) <- names(df_caa_tab_abun)
df_caa_tab_abun <- as.data.frame(rbind(df_caa_tab_abun, df_tmp))

imp_post <- colMeans(df_caa_tab_abun[14:24, 2:4], na.rm = T)
df_caa_tab_abun[25, 2:4] <- imp_post
df_caa_tab_abun[26, 2:4] <- imp_post

# tmp <- df_caa %>% 
#    group_by(year, age) %>%
#    filter(age != "Unknown") %>%
#    summarise(biomass_sum = sum(weight_tonnes), prop_mat_mean = mean(prop_mat))
# tmp$age <- as.numeric(tmp$age)
# str(tmp, give.attr = FALSE)
# 
# df_caa_tab_bio <- tmp[c("year", "age", "biomass_sum")] %>%
#    pivot_wider(id_cols = year, names_from = age, values_from = biomass_sum, names_sort = T) %>%
#    rename(c1 = '1', c2 = '2', c3 = '3', c4 = '4', c5 = '5', c6 = '6') %>%
#    select(year, c2, c3, c4)
# 
# df_caa_tab_bio

# insert 2022 and 2023 as average of last 10 years
# df_tmp <- as.data.frame(matrix(NA, 2, 4))
# df_tmp[, 1] <- c(2022, 2023)
# names(df_tmp) <- names(df_caa_tab_bio)
# df_caa_tab_bio <- as.data.frame(rbind(df_caa_tab_bio, df_tmp))
# 
# imp_post <- colMeans(df_caa_tab_bio[14:24, 2:4], na.rm = T)
# df_caa_tab_bio[25, 2:4] <- imp_post
# df_caa_tab_bio[26, 2:4] <- imp_post

# CAA 1982-1997 in landings (tonnes) by age
# df_caa1982_1997 <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/Fran/catchAtAge1982_1997bio.csv")
# str(df_caa1982_1997, give.attr = FALSE)


# summ and clean up
# originally did this for boimass which doesn't match with what I did for the post-collapse nor does it match with the SSM. I think that biomass makes more sense but since I don't have that yet, sticking to abundance. 

# tmp_caa <- df_caa1982_1997 %>%
#    group_by(year) %>%
#    summarize(across(age1:total, sum)) %>%
#    rename(c1 = 'age1', c2 = 'age2', c3 = 'age3', c4 = 'age4', c5 = 'age5', c6 = 'age6') %>%
#    select(year, c2, c3, c4)
# 
# # Note that we don't need 1982 or 1983, 1984 and 1989, 1991, 1992 are missing, and 1995 had no fishery.  Add 2022 and 2023
# 
# df_tmp <- as.data.frame(matrix(NA, 3, 4))
# df_tmp[, 1] <- c(1989, 1991, 1992)
# names(df_tmp) <- names(tmp_caa)
# df_caa_tab_1985_1997 <- (rbind(df_tmp, tmp_caa))
# df_caa_tab_1985_1997 <- df_caa_tab_1985_1997[order(df_caa_tab_1985_1997$year),]
# 
# # average of precollapse years
# imp_pre <- colMeans(df_caa_tab_1985_1997[1:7, 2:4], na.rm = T)
# df_caa_tab_1985_1997[7, 2:4] <- imp_pre
# df_caa_tab_1985_1997[9:10, 2:4] <- 0
# 
# df_caa_all <- rbind(df_caa_tab_1985_1997[3:15,], df_caa_tab_bio)
# 
# matCAA <- as.matrix(df_caa_all[, 2:4])
# str(matCAA)

# note that this actually goes to 1998
df_caa1982_1997 <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/Fran/catchAtAge1982_1997abun.csv")
str(df_caa1982_1997, give.attr = FALSE)


tmp_caa <- df_caa1982_1997 %>%
   group_by(year) %>%
   summarize(across(age1:age6, sum)) %>%
   mutate_at(vars(age1:age6), ~ ./ 1000) %>%
   rename(c1 = 'age1', c2 = 'age2', c3 = 'age3', c4 = 'age4', c5 = 'age5', c6 = 'age6') %>%
   select(year, c2, c3, c4)

# Note that we don't need 1982 or 1983, 1984 and 1989, 1991, 1992 are missing, and 1995 had no fishery.  Add 2022 and 2023

df_tmp <- as.data.frame(matrix(NA, 4, 4))
df_tmp[, 1] <- c(1984, 1989, 1991, 1992)
names(df_tmp) <- names(tmp_caa)
df_caa_tab_1985_1997 <- (rbind(df_tmp, tmp_caa))
df_caa_tab_1985_1997 <- df_caa_tab_1985_1997[order(df_caa_tab_1985_1997$year),]

# average of precollapse years
imp_pre <- colMeans(df_caa_tab_1985_1997[1:7, 2:4], na.rm = T)
df_caa_tab_1985_1997[3, 2:4] <- imp_pre
df_caa_tab_1985_1997[8, 2:4] <- imp_pre
df_caa_tab_1985_1997[c(10:11, 14), 2:4] <- 0


df_caa_all <- rbind(df_caa_tab_1985_1997[3:15,], df_caa_tab_abun)

write.csv(df_caa_all, "data/CAA.csv")

matCAA <- as.matrix(df_caa_all[, 2:4])
str(matCAA)

# df_caa_tab_mat <- tmp[c("year", "age", "prop_mat_mean")] %>%
#      pivot_wider(id_cols = year, names_from = age, values_from = prop_mat_mean, names_sort = T) %>%
#      rename(cm1 = '1', cm2 = '2', cm3 = '3', cm4 = '4', cm5 = '5', cm6 = '6') %>%
#      select(year, cm2, cm3, cm4)
# 
# df_caa_tab_mat
# 
# df_tmp <- as.data.frame(matrix(NA, 15, 4))
# df_tmp[, 1] <- c(1985:1997, 2022:2023)
# names(df_tmp) <- names(df_caa_tab_mat)
# df_caa_tab_mat <- rbind(df_tmp, df_caa_tab_mat)
# df_caa_tab_mat <- df_caa_tab_mat %>%
#      arrange(year) 
# 
# df_caa_tab_mat[, 2:4][df_caa_tab_mat[, 2:4] == 1] <- 0.995
# 
# imp <- colMeans(df_caa_tab_mat[,2:4], na.rm = T)
# df_caa_tab_mat[c(1:13, 38:39), 2:4] <- imp
# 
# 
# matCAA_m <- as.matrix(df_caa_tab_mat[, 2:4])

p <- ggplot(data = df_caa_all, aes(x = year))
p <- p + geom_line(aes(y = log10(c2)), colour = "red")
p <- p + geom_line(aes(y = log10(c3)), colour = "green")
p <- p + geom_line(aes(y = log10(c4)), colour = "blue")
p <- p + geom_line(aes(y = log10(c2+c3+c4)))
p


# Bundle data----
num_forecasts = 2 # 2 extra years
jags.data <- ls_jag("yes", "yes", "no")
str(jags.data)
jd <- as.data.frame(jags.data[2:8])

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
