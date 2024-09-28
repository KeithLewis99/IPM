
## TODO
## - expand to 5+
## - get mean weights for age 1
## - get units and scale of capelin abundance at age right
## - get most recent landings estimates
## - get more complete CAA
## - get more complete mat at age for ages 1 to 6
## - why is the sum of index at age * weight at age roughly double the biomass index?


### Capelin ###
#-------------#
library(wham)

library(ggplot2)
library(dplyr)
library(tidyr)
Cap_data <- readRDS("capelin/CapData.rds")

source("capelin/approx_catch_paa_from_figure.r")

years <- as.integer(1985:2022)
ages <- as.integer(1:6)

landings <- read.csv("capelin/capelin_landings_kt.csv") |> 
    subset(Year %in% years)

iaa <- read.csv("capelin/capelin_abundance_1985-2022.csv") |> 
    subset(year %in% years) |> 
    select(year, I1, I2, I3, I4, I5, I6)

# Cap_data$landings <- Cap_data$landings[Cap_data$landings$year < 2022,]

catch <- as.matrix(landings[, 2]) * 1000
catch_cv <- matrix(c(rep(0.5, 7), rep(0.2, 8), rep(.1, 23)), 
                   ncol = 1, nrow = length(years))
catch_Neff <- matrix(50, ncol = 1, nrow = length(years))
catch_paa <- array(approx_paa_mat, dim = c(1, length(years), length(ages)))
use_catch_paa <- matrix(1, nrow = length(years), ncol = 1)
selblock_pointer_fleets <- matrix(1, nrow = length(years))

index_vals <- iaa |> select(I1:I6) |> as.matrix()
index_vals[is.na(index_vals)] <- 0
index <- (rowSums(index_vals, na.rm = TRUE) * 1000000) %>% 
    matrix(nrow = length(years), ncol = 2, byrow = F)
index_cv <- matrix(0.1, ncol = 2, nrow = length(years))
index_Neff <- matrix(50, ncol = 2, nrow = length(years))
index_fracyr <- matrix(5/12, ncol = 2, nrow = length(years))
units_indices <- 2 # 2 is abundance, 1 is biomass indices??
units_index_paa <- 2
use_indices <- matrix(index[, 1] > 0, ncol = 2, nrow = length(years))
use_index_paa <- matrix(index[, 1] > 0, ncol = 2, nrow = length(years))


## juveniles ## BUMP FORWARD ONE YEAR AS THESE LARVAE WILL BE AGE 1 IN JAN
index[,2] <- Hmisc::Lag(Cap_data$LD[seq_along(years)], 1) %>% replace_na(0)
index_cv[,2] <- 0.3
# index_Neff[,2] <- 50
# index_fracyr[,2] <- 5/12
units_indices <- c(2, 2)
units_index_paa <- c(2, 2)
use_indices[,2] <- index[,2] > 0
use_index_paa[,2] <- 0 # PAA NOT NEEDED SINCE ONLY ONE AGE REPRESENTED

ipaa <- prop.table(index_vals, 1)
ipaa[is.nan(ipaa)] <- 0
index_paa <- array(0, dim = c(2,length(years),length(ages)))
for(y in 1:length(years)){
    for(a in 1:length(ages)){
        index_paa[1,y,a] <- ifelse(is.na(ipaa[y,a]), 0, ipaa[y,a])
        index_paa[2,y,a] <- ifelse(a == 1, 1, 0)
    }
}
selblock_pointer_indices <- matrix(rep(c(2, 3), each = length(years)), 
                                   ncol = 2, nrow = length(years))


waa_vals <- xtabs(w ~ year+age, data = Cap_data$WaA, 
                  subset = age %in% ages)
ind <- waa_vals[, "1"] > 0
waa_vals[ , "1"] <- mean(waa_vals[ind, "1"])
# ind <- waa_vals[, "6"] > 0
# waa_vals[ , "6"] <- mean(waa_vals[ind, "6"])
waa_vals <- unname(waa_vals)

waa_ex <- apply(waa_vals, 2, \(x)mean(x, na.rm=T))
waa_full <- rbind(matrix(waa_ex, nrow = 13, ncol = length(ages), byrow = T), 
                  waa_vals,
                  matrix(waa_ex, nrow = 1, ncol = length(ages), byrow = T))
waa <- array(waa_full / 1000, dim = c(1, length(years), length(ages)))

mat <- array(cbind(0, Cap_data$matMp[seq_along(years), ], 1, 1), dim = c(1, length(years), length(ages)))
waa_pointer_indices <- c(1, 1)
waa_pointer_fleets <- 1
waa_pointer_ssb <- 1
fracyr_ssb <-  matrix(0, ncol = 1, nrow = length(years))




## Does the scale make sense? (convert index numbers at age from billions to thousands, and multiply by weight at age in kg)
index_tonnes <- rowSums(index_vals[, c(2:4)] * 1000000 / 2 * waa[1,, c(2:4)])
plotly::plot_ly(y = index_tonnes, x = years) |> plotly::add_lines()
## sum of index at age multiplied by weight at age is double the biomass index. Don't know why.

index[, 1] <- index[, 1] / 2 # TEMP FIX: HALF THE CURRENT INDEX



selectivity <- list(model = c("age-specific", "age-specific", "age-specific"), 
                    n_selblocks = 3,
                    re = c("2dar1", "none", "none"),
                    initial_pars = list(c(0, rep(0.1, length(ages) - 1)), #Catches
                                        c(0.5, rep(1, length(ages) - 1)), #, #Acoustics
                                        c(1, rep(0, length(ages) - 1))), #Juveniles
                    map_pars = list(c(NA, rep(1, length(ages) - 1)), 
                                    c(2, rep(NA, length(ages) - 1)), 
                                    c(rep(NA, length(ages))))
                    )

max_age <- 6
4.899 * max_age ^ (-0.916) # approx = 1; Then et al. (2015)
M_in <- list(initial_MAA = array(1, dim = c(1,1,length(years),length(ages)) ))

q_in <- list(q_upper = c(1, 1)) # survey q should not exceed 1

F_in <- list(
    F = cbind(rep(5, length(years))),
    map_F = cbind(rep(NA, length(years))))

basic_info <- list(
  n_stocks = 1L,
  ages = ages,
  n_seasons = 1L,
  n_fleets = 1L,
  fracyr_SSB = fracyr_ssb,
  maturity = mat,
  years = years,
  waa = waa,
  waa_pointer_ssb = waa_pointer_ssb
)

catch_info <- list(
  n_fleets = NCOL(catch),
  agg_catch = catch,
  agg_catch_cv = catch_cv,
  catch_paa = catch_paa,
  use_catch_paa = use_catch_paa,
  catch_Neff = catch_Neff,
  selblock_pointer_fleets = selblock_pointer_fleets,
  waa_pointer_fleets = waa_pointer_fleets
)

index_info <- list(
  n_indices = NCOL(index),
  agg_indices = index,
  units_indices = units_indices,
  units_index_paa = units_index_paa,
  agg_index_cv = index_cv,
  fracyr_indices = index_fracyr,
  use_indices = use_indices,
  use_index_paa = use_index_paa,
  index_paa = index_paa,
  index_Neff = index_Neff,
  selblock_pointer_indices = selblock_pointer_indices,
  waa_pointer_indices = waa_pointer_indices
)

NAA_in <- list(N1_model = "equilibrium")

# ecov_info <- set_ecov()

input_all <- prepare_wham_input(basic_info = basic_info, 
                                selectivity = selectivity, 
                                catch_info = catch_info, 
                                index_info = index_info, 
                                M = M_in, F = F_in,
                                catchability = q_in, 
                                NAA_re = NAA_in,
                                age_comp = "logistic-normal-miss0") 

fit <- fit_wham(input_all, do.fit = F, do.retro = F, do.brps = F, do.osa = F)
fit$fn()
fit$rep$NAA[1,1,,]

fit <- fit_wham(input_all, do.fit = T, do.retro = F, do.brps = F, do.osa = F, do.sdrep = T)
fit$fn()
fit$opt
fit$sdrep
round(fit$rep$NAA[1,1,,], 2)
matplot(fit$rep$NAA[1,1,,], type = "l")

# plot_wham_output(fit, res = 600, dir.main = file.path(getwd(), "capelin"))

input1 <- input_all
fit1 <- fit

## Recruitment and cohort deviations -------------------------------------------

input2 <- set_NAA(input1, 
                  list(recruitment_model = 2,
                       sigma = "rec+1",
                       cor = "2dar1"))

fit2 <- fit_wham(input2, do.fit = T, do.retro = F, do.brps = F, do.osa = F, do.sdrep = T)
fit2$fn()
fit2$opt
fit2$sdrep
round(fit2$rep$NAA[1,1,,], 2)
matplot(fit2$rep$NAA[1,1,,], type = "l")

# plot_wham_output(fit2, res = 600, dir.main = file.path(getwd(), "capelin"))


## M deviations ----------------------------------------------------------------

input3 <- set_NAA(input1, 
                  list(recruitment_model = 2,
                       sigma = "rec+1",
                       cor = "2dar1")) |> 
    set_M(list(mean_model = "estimate-M",
               re_model = matrix("none")))

fit3 <- fit_wham(input3, do.fit = T, do.retro = F, do.brps = F, do.osa = F, do.sdrep = T)
fit3$opt
fit3$sdrep

# plot_wham_output(fit3, res = 600, dir.main = file.path(getwd(), "capelin"))





