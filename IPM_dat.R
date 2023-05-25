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
## 1999-2022 abundance- and biomass-at-age but also with mature abundance- and biomass-at-age.  All have lower and upper CIs.  
###from Aaron and the Shiny App 

### Units millions -> convert to billions below
### Units in tonnes -> convert to kilotonnes below
df_dis <- read_csv("data/abundance and biomass by age and year2.csv")
str(df_dis)

# bring in the historical data - 
## 1985 2017 - abundance: # in billions
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
   mutate_at(vars(abun:matabun.uci), round, 2) # round - this may be causing a slight discrepancy between my values and Aarons.

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
df_dis_tab <- as.data.frame(apply(df_dis_tab, 2, function(x) replace(x, is.na(x), 0)))

# create a df with missing years
df_tmp <- df_dis_tab[1:4,]
df_tmp[, 1:length(df_tmp)] <- NA
df_tmp$year[1:nrow(df_tmp[])] <- c(2006, 2016, 2020, 2021)
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
df_dis_tab$I6 <- NA

# combine the data sets 
## filter out Unknowns, mature, perAge2, var, and sd 
if(disaggregated == "1985-present") {
        df_dis_tab <- rbind(df_dis_1998[, c(1:7, 11:13)], df_dis_tab[, c(1:6, 11, 8:10)])
} else {
        df_dis_tab
} 

# abundance-at-age 1985-2022
write.csv(df_dis_tab, "data/derived/capelin_abundance_1985-2022.csv")

# # abundance-at-age 1985-2022 in natural logarithms
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
### - note that that I don't have the 5th and 95th percentiles for df_dis_1998 so I have to import a new data set

# df_ag_1985 <- df_dis_1998[, c(1, 11)] 
# colnames(df_ag_1985)[2] <- "abundance_med"
# df_ag_1985$ab_lci <- NA
# df_ag_1985$ab_uci <- NA
# units in billions and kilotonnes - note that I don't have the 1982 abundance - i'll add the biomass below
df_ag_1985 <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/capelin-2021.csv")
str(df_ag_1985)

# combine the 1985-1998 and 1999-present.
df_agg <- rbind(df_ag_1985[1:14, 1:4], df_ag_1999)

# add years with no data
df_tmp <- df_agg[1:4,]
df_tmp[, 1:length(df_tmp)] <- NA
df_tmp$year[1:nrow(df_tmp)] <- c(2006, 2016, 2020, 2021)
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


# biomass-at-age 1999-present - from Aaron and the Shiny App 
## from the abundance import above df_dis <- read_csv("data/abundance and biomass by age and year2.csv")
## Units - tonnes
###All have lower and upper CIs.  

# Select columns, change names, and convert units to kt and round
## note that rounding may cause some discrepancies between these values and Adamack when they are aggregated
df_baa_filter <- df_dis %>%
   group_by(year, age) %>%
   select(year, age,
          bio = med.bm.age.fran, 
          bio.lci = low.bm.age.fran, bio.uci = up.bm.age.fran, 
          matbio = med.mat.bm.age.fran, 
          matbio.lci = low.mat.bm.age.fran, matbio.uci = up.mat.bm.age.fran
) %>%  # select and rename
   mutate_at(vars(bio:matbio.uci), ~ ./ 1000) %>% # convert units to kt
   mutate_at(vars(bio:matbio.uci), round, 2) # round - this may be causing a slight discrepancy between my values and Aarons.


# pivot the data - longer to wider with the disaggregated abundance as columns
# get abundance/biomass by age
# these values are close to those in df_cap but not exact - perhaps due to rounding?  
df_baa_tab <- df_baa_filter[, c(1:3)] %>%
   pivot_wider(names_from = age, values_from = bio) %>%
   rename(bio1 = '1', bio2 = '2', bio3 = '3', bio4 = '4', bio5 = '5') %>%
   mutate(biomass = rowSums(across(bio1:Unknown), na.rm = T)) %>%
   mutate(var = var(c_across(bio1:Unknown), na.rm = T))  %>%
   mutate(sd = sd(c_across(bio1:Unknown), na.rm = T))
df_baa_tab

# create a df with missing years
df_tmp <- df_baa_tab[1:4,]
df_tmp[, 1:length(df_tmp)] <- NA
df_tmp$year[1:nrow(df_tmp)] <- c(2006, 2016, 2020, 2021)
df_tmp

# bind the blank years (NAs) with the data
df_baa_tab <- bind_rows(df_tmp, df_baa_tab[, c(1:8, 9:10)]) %>% 
   arrange(year)
df_baa_tab$bio6 <- NA
df_baa_tab <- df_baa_tab[,c(1:6, 11, 7:8)]
str(df_baa_tab, give.attr = F)

# create a total biomass column for the 1985-2012 data, get relevant rows, and fill in NAs
df_baa_FM$biomass <- NA
df_baa_FM$Unknown <- NA
df_baa_1985 <- df_baa_FM[1:14, c(1:7, 9, 8)] # just 1985-1998 and relevant columns

## age-aggregated biomass for 1982-1999 - from Mariano I think
## just adding this so that the age-disaggregated biomass has a biomass column and this may be easier than summing the biomass and then expanding the df for the 1982 values.
### see notes below - this may not be a great long term way to deal with this.
vec_bio <- c(446, NA, NA, 3426, 3697, 2576, 4285, 3712, 5783, 138, 138, NA, NA, NA, 47, NA, NA
)
length(vec_bio)

df_bio_agg_1982 <- data.frame(cbind(seq(1982, 1998), vec_bio)) %>%
   rename(year = V1, biomass = vec_bio)

# join the age disaagregated biomass to the aggregated biomass.
df_baa_1985 <- left_join(df_bio_agg_1982, df_baa_1985[, 1:8], "year")


# combine the data sets as needed.
if(disaggregated == "1985-present") {
   df_baa_tab <- rbind(df_baa_1985[,c(1, 3:9, 2)], df_baa_tab)
} else {
   df_baa_tab
} 

# save just the baa for 1985-present because no age disaggregated for 1982-1984.
write.csv(df_baa_tab[4:40, -9], "data/capelin_biomass_1985-2022.csv", row.names = F)

# pivot the data - longer to wider with the disaggregated abundance as columns
# abundance value in natural logarithms
df_baa_tabLog <- df_baa_tab %>%
   mutate(B2 = log(bio2), B3 = log(bio3), B4 = log(bio4), B = log(biomass)) %>%
   select(year, B2, B3, B4, B)

df_baa_tabLog


### biomass agg ----
#### note that I brought in the aggregated data from 1982-present above in a vector - this was a quick and dirty solution for some of the work I did on LRPs right before the RAP - may want to fix.

## aggregate age-disagregated data >= 1999
### - note that these are 5th and 95th percentiles which is a bootstrap CI!!!!!
#### note that this step is a bit redundant, i.e., Fran/Aaron have probably already done this but I want my own so that I don't have to continutally ask for data, i.e, I just want to be able to get age-disaggregated data with maturities and run the analysis.
df_ag_bio_1999 <- df_baa_filter[, c(1:5)] %>%
   group_by(year) %>%
   summarize(biomass_med = sum(bio), bm_lci = sum(bio.lci), bm_uci = sum(bio.uci))

## aggregate age-disagregated data < 1999
### - note that that I don't have the 5th and 95th percentiles for these and they only go back to 1988

# commented out lines are for if I use the 1982 onwards, biomass only
# colnames(df_bio_agg_1982)[2] <- "biomass_med"
# df_bio_agg_1982$bm_lci <- NA
# df_bio_agg_1982$bm_uci <- NA

df_ag_1985[, c(1, 5:8)]

# combine the data sets as needed.
df_agg_bio <- rbind(df_ag_1985[1:14, c(1, 5:7)], df_ag_bio_1999)

# create a df with missing years
df_tmp <- df_agg_bio[1:6,]
df_tmp[1:length(df_tmp$year), ] <- NA
df_tmp$year <- c(1982:1984, 2006, 2016, 2020)
df_tmp[1:3, 2] <- df_bio_agg_1982[1:3, 2]
df_tmp

# bind the blank years (NAs) with the data
df_agg_bio <- bind_rows(df_tmp, df_agg_bio) %>% 
   arrange(year)

write_csv(df_agg_bio, "data/capelin_aggregated_biomass_1985-2022.csv")


## maturity ----
### USE THIS ONE, NOT CODE UNDER maturity or mat-matrix
## 1999-present - abundance mature based on import from AA file
df_mat_tab <- df_dis_summ[, c(1:2, 6)] %>%
   filter(age != 1) %>% # because they aren't mature
   pivot_wider(names_from = age, values_from = matabun) %>%
   rename(mat2 = '2', mat3 = '3', mat4 = '4', mat5 = '5') %>%
   mutate(matureAbun = sum(c_across(starts_with("m")), na.rm = T)) %>%
   mutate(var = var(c_across(starts_with("m")), na.rm = T))  %>%
   mutate(sd = sd(c_across(starts_with("m")), na.rm = T))
df_mat_tab
str(df_mat_tab, give.attr=F)

# # create empty rows
df_tmp <- df_mat_tab[1:4,]
df_tmp[, 1:length(df_tmp)] <- NA
df_tmp$year[1:nrow(df_tmp)] <- c(2006, 2016, 2020, 2021)
df_tmp

# bind the blank years (NAs) with the data
## this is mature-abundance-at-age
df_mat_tab <- bind_rows(df_tmp, df_mat_tab) %>% 
   arrange(year)

## 1985-2012 - abundance mature
### based on file from the biochar file (BIOCHAR FROM ACOUSTICS_revised to use Monte Carlo abundnace for 1988-1996 in annual page (003).xls; C:\Users\lewiske\Documents\capelin_LRP\IPM\data)  
df_mat_1985 <- read_csv("data/matAbun.csv")
str(df_mat_1985, give.attr = F)


# this is the mature abundance for age-2 to -5 and total (includes unknowns and 6s)
df_mat <- rbind(df_mat_1985[1:14, c(-2, -7)], df_mat_tab[,c(1:5, 7)])
str(df_mat, give.attr = F)
matM <- as.matrix(df_mat[, 2:4])


# impute
## Note that this works but do we really want to impute???  No I think.
df_mat[7:38, ] <- lapply(df_mat[7:38, ], function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
matM <- as.matrix(df_mat[, 2:4])

# natural log of maturity
## this may not be needed.
df_mat_tabLog <- df_mat %>%
   mutate(M2 = log(mat2), M3 = log(mat3), M4 = log(mat4), M = log(matureAbun)) %>%
   select(year, M2, M3, M4, M)

df_mat_tabLog

### proprotion ----
# this is just taking the year column separately, then, dividing the mature abundance by the total mature and multiplying by 100 to get a percentage.
df_mat_prop <- cbind(df_mat[1], df_mat[-1]/df_dis_tab[c(3:7)])
str(df_mat_prop, give.attr = F)
write.csv(df_mat_prop, "data/derived/capelin_perMat_1985-2022.csv", row.names = F)

#impute
df_mat_prop[7:38, ] <- lapply(df_mat_prop[7:38, ], function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))


# confirms the means of the years
mean_tmp <- apply(df_mat_prop[7:38, 2:5], 2, function(x) mean(x))

# add NAs for forecasts
matMp <- as.data.frame(apply(df_mat_prop, 2, function(x) c(x, rep(NA, num_forecasts))))
# fill years
matMp[39:40, 1] <- c(2023, 2024)

# subset - remember that this is a loop so the the mean(x) is a vector.  If you do x[7:38,] you get dimension errors
# the c(7:8, 12, 15:21, 23:31, 33:35, 38) is so that we are not calculation averages with imputed values
matMp <- apply(matMp, 2, function(x) replace(x, is.na(x), mean(x[c(7:8, 12, 15:21, 23:31, 33:35, 38)], na.rm = T)))

# note that this could be done more cleanly by simply adding the years 2023 and 2024 to df_mat_prop and then the following 
#df_mat_prop[7:40, ] <- lapply(df_mat_prop[7:40, ], function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
#Keeping the above because it took so long to figure out. 


## Trinity Bay ----
### 1999-2019 (update when needed)
### Note that thie below is for abundance only.  I have not done similar work for biomass altough it would be the same as above.
### Units millions -> convert to billions below
### Units in tonnes -> convert to kilotonnes below
df_tb <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/fromAaron/TB_abun_atAge.csv")
str(df_tb, give.attr = F)
head(df_tb)
df_tb$age <- as.factor(df_tb$age)
df_tb$stratum <- as.factor("TB")

df_tb_NAA <- df_tb %>%
   # filter(age != "Unknown" & age != "1" & age != "5") %>%
   filter(abundance >0) %>%
   pivot_wider(id_cols = year, names_from = age, values_from = abundance, names_sort = T) %>%
   rename(I1 = '1', I2 = '2', I3 = '3', I4 = '4', I5 = '5') %>%
   mutate_at(vars(I1:Unknown), ~ ./ 1000)  
str(df_tb_NAA)
df_tb_NAA

# add missing years
df_tmp <- df_tb_NAA[1:7,] 
df_tmp[, 1:length(df_tb_NAA)] <- NA
df_tmp$year[1:nrow(df_tmp)] <- c(2006, 2014:2016, 2020:2022)
df_tmp

# bind summarized data with missing data
df_tb_NAA <- bind_rows(df_tmp, df_tb_NAA) %>% 
   arrange(year)
str(df_tb_NAA)

if(disaggregated == "1985-present") {
   df_tmp <- as.data.frame(matrix(NA, 14, 7))
   df_tmp[, 1] <- c(1985:1998)
   names(df_tmp) <- names(df_tb_NAA)
   df_tb_NAA <- rbind(df_tmp, df_tb_NAA)
} else {
   df_con <- df_con %>%
      slice(5:27)
} 

# impute
## 2006, 2014:2016, 2020:2022
## Note that this works but do we really want to impute???  No I think.
df_tb_NAA[15:38, 3:5] <- lapply(df_tb_NAA[15:38, 3:5], function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))


# convert to a Matrix
matITB <- as.matrix(df_tb_NAA[, 3:5])


df_tb_NAALog <- df_tb_NAA %>%
   mutate(I1 = log(I1), I2 = log(I2), I3 = log(I3), I4 = log(I4), I5 = log(I5))



## TB maturity ----
df_tb_matAA <- df_tb %>%
     filter(age != "Unknown" & age != "1" & age != "5") %>%
     pivot_wider(id_cols = year, names_from = age, values_from = prop_mat, names_sort = T) %>%
     rename(m2 = '2', m3 = '3', m4 = '4')

# add missing years
df_tmp <- df_tb_matAA[1:7,] 
df_tmp[1:length(df_tb_matAA), ] <- NA
df_tmp$year[1:nrow(df_tmp)] <- c(2006, 2014:2016, 2020:2022)
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
## Note that this works but do we really want to impute???  No I think.
df_tb_matAA[15:38, 2:4] <- lapply(df_tb_matAA[15:38, 2:4], function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))

# change '1' to high percentage
df_tb_matAA[,3][df_tb_matAA[,3] == 1] <- 0.995
     
m_maaTB <- as.matrix(df_tb_matAA[,2:4])

## USSR data 1981-1992----
### This will take a lot of attention from Divya or maybe a M.Sc student.  It is beyond me to try and determine how to compare the spatial and temproal coverages of the various Soviet surveys and how they calculate the biomass and abundance given these discrepancies.


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


df_ld <- df_ld %>% rename(year = `Year`,
                          larvae = `Larval densities_ind_m-3`,
                          se_auc = `SE_AUC`) 
df_ld$lnlarvae <- log(df_ld$larvae)
str(df_ld)


## ice ----
#Note that I added in dummy data for 2021
df_ice  <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/ice-m1-2021.csv"
)
str(df_ice)

if(disaggregated == "1985-present") {
        df_ice <- df_ice %>% # get rid of years 1969-1984
            slice(17:54)
} else {
        df_ice <- df_ice %>% # get rid of years 1969-1998
        slice(31:54)
} 

# add missing years
## Note just doing this because I don't have the 2022 data
df_tmp <- df_ice[1,] 
df_tmp[, 1:length(df_ice)] <- NA
df_tmp$year[1:nrow(df_tmp)] <- c(2022)
df_tmp
df_ice <- rbind(df_tmp, df_ice) %>%
   arrange(year) 

# impute
## Note just doing this because I don't have the 2022 data
df_ice[, 2:4] <- lapply(df_ice[, 2:4], function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))

## condition ----
#Note that I added in dummy data for 2021
#df_con <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/analyses/capelinLRP/data/condition_ag1_2_MF_out.csv")
df_con <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/fromAaron/condition_2JK_ag1_2_MF_2022.csv")
str(df_con)

# add missing year.  Note that 2002 is there just because I don't have the 2023 data
if(disaggregated == "1985-present") {
        df_tmp <- as.data.frame(matrix(NA, 11, 2))
        df_tmp[, 1] <- c(1985:1994, 2022)
        names(df_tmp) <- names(df_con)
        df_con <- rbind(df_tmp, df_con) %>%
                 arrange(year) 
} else {
        df_con <- df_con %>%
                slice(5:27)
        } 
str(df_con)


# impute
## Note just doing this because I don't have the 2022 data
df_con <- as.data.frame(apply(df_con, 2, function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))))


## catch-at-age----
# CAA 1998-2021
df_caa <- read_csv("C:/Users/lewiske/Documents/capelin_LRP/data/fromAaron/catchAtAge1998_2021.csv")
## units millions and tonnes
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
df_tmp <- as.data.frame(matrix(NA, 1, 4))
df_tmp[, 1] <- c(2022)
names(df_tmp) <- names(df_caa_tab_abun)
df_caa_tab_abun <- as.data.frame(rbind(df_caa_tab_abun, df_tmp))

#impute
df_caa_tab_abun[, 2:4] <- lapply(df_caa_tab_abun[, 2:4], function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))


## this is the same as the above but for biomass
### i'm going to leave the commented out code for CAA because I may use this at some point
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

# CAA 1982-1997 biomass
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

# CAA 1982-1997 abundance
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
df_caa_tab_1985_1997 <- rbind(df_tmp, tmp_caa)
df_caa_tab_1985_1997 <- df_caa_tab_1985_1997[order(df_caa_tab_1985_1997$year),]

# average of precollapse years
#impute
df_caa_tab_abun[1:7, 2:4] <- lapply(df_caa_tab_abun[1:7, 2:4], function(x) replace(x, is.na(x), mean(x, na.rm = TRUE)))
# no catches
df_caa_tab_1985_1997[c(10:11, 14), 2:4] <- 0

df_caa_all <- rbind(df_caa_tab_1985_1997[3:15,], df_caa_tab_abun)

write.csv(df_caa_all, "data/CAA.csv")

matCAA <- as.matrix(df_caa_all[, 2:4])
str(matCAA)

## this is the mature CAA - probably not relevant since most fish caught are mature, at least in more recent times.  Perhaps this was not true pre-2000 when lots of discarding.
### i'm going to leave the commented out code for CAA because I may use this at some point

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


#source("IPM_fun.R")
# Bundle data----
num_forecasts = 2 # 2 extra years
jags.data.m <- ls_jag("yes", "yes", "yes")
str(jags.data.m)

# check that the lengths of the lists all match
leng_jd <- rep(NA, 8)
for (i in 1:length(jags.data.m)){
   len_jd <- length(jags.data.m[[i]])
   leng_jd[i] <- len_jd
}
leng_jd

# Note that I haven't done the imputation for df_mat_prop
jd <- as.data.frame(jags.data.m[c(2:5, 7:11)])
