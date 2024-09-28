
## TODO
## - expand to 5+
## - get mean weights for age 1
## - get units and scale of capelin abundance at age right
## - get most recent landings estimates
## - get more complete CAA
## - get more complete mat at age for ages 1 to 6
## - why is the sum of index at age * weight at age roughly double the biomass index?
## - are catch proportions at age based on numbers or biomass?
## - estimate observation error for surveys lacking a CV (all except spring survey)

## - consider replacing zeros in the paa with a small value
## - try to use direchlet-multinomial and supply a large Neff

## - apply q prior to all acoustic surveys
## - try estimating age-specific + ar1_y selectivity for the fishery
## - try to impose F blocks and/or selectivity blocks for pre-collapse, moratorium, then post-collapse
##   to capture major fishery changes


### Capelin ###
#-------------#
library(wham) # using dev branch

library(ggplot2)
library(dplyr)
library(tidyr)
Cap_data <- readRDS("wham/data/CapData.rds")

source("wham/data/approx_catch_paa_from_figure.r")

years <- as.integer(1972:2022)
ages <- as.integer(1:6)

approx_paa_mat <- unname(approx_paa_mat[as.character(years), as.character(ages)])

landings <- read.csv("wham/data/capelin_landings_kt.csv") |> 
    subset(Year %in% years)

iaa <- read.csv("wham/data/capelin_abundance_1985-2022.csv") |> 
    right_join(data.frame(year = years), by = "year") |> 
    select(year, I1, I2, I3, I4, I5, I6) |> 
    arrange(year)

baa <- read.csv("wham/data/capelin_biomass_1985-2022.csv") |> 
    right_join(data.frame(year = years), by = "year") |> 
    select(year, bio1, bio2, bio3, bio4, bio5, bio6) |> 
    arrange(year)

can_spring <- read.csv("wham/data/can_spring_acoustic_biomass.csv", col.names = c("year", "can_spring", "lwr", "upr"))
larval_den <- data.frame(year = 1985:2024, larval_den = Cap_data$LD)
can_fall <- read.csv("wham/data/can_fall_acoustic_biomass.csv", col.names = c("year", "can_fall"))
ussr_fall <- read.csv("wham/data/ussr_fall_acoustic_biomass.csv", col.names = c("year", "ussr_fall"))
ussr_spring <- read.csv("wham/data/ussr_spring_acoustic_biomass.csv", col.names = c("year", "ussr_spring"))

## TODO: drop approximation of standard errors if pelagics are able to provide
##       se values rather than just the 95% confidence intervals
approx_se_from_ci <- function(lower, upper) (upper - lower) / (2 * 1.96)
can_spring$approx_se <- approx_se_from_ci(can_spring$lwr, can_spring$upr)
can_spring$approx_cv <- can_spring$approx_se / can_spring$can_spring

# Cap_data$landings <- Cap_data$landings[Cap_data$landings$year < 2022,]

catch <- as.matrix(landings[, 2]) * 1000
catch_cv <- matrix(ifelse(years < 1991, 0.5, ifelse(years < 2000, 0.2, 0.1)), 
                   ncol = 1, nrow = length(years))
catch_Neff <- matrix(50, ncol = 1, nrow = length(years))
catch_paa <- array(approx_paa_mat, dim = c(1, length(years), length(ages)))
use_catch_paa <-as.matrix(rowSums(approx_paa_mat) > 0) |> unname() 
selblock_pointer_fleets <- matrix(1, nrow = length(years))

index <- can_spring |> 
    select(year, can_spring) |> 
    right_join(data.frame(year = years), by = "year") |> 
    arrange(year) |> 
    left_join(can_fall, by = "year") |> 
    left_join(ussr_fall, by = "year") |> 
    left_join(ussr_spring, by = "year") |> 
    left_join(larval_den, by = "year") |> 
    select(-year) |> 
    as.matrix() |> 
    unname()
index[, 1:4] <- index[, 1:4] * 1000 # convert all biomass indices to tonnes
index[is.na(index)] <- 0

index_cv <- can_spring |> 
    select(year, approx_cv) |> 
    right_join(data.frame(year = years), by = "year") |> 
    arrange(year) |> 
    mutate(can_fall_cv = mean(approx_cv, na.rm = TRUE) * 2,
           ussr_fall_cv = mean(approx_cv, na.rm = TRUE),
           ussr_spring_cv = mean(approx_cv, na.rm = TRUE),
           larval_den_cv = 0.3,
           approx_cv = replace_na(approx_cv, mean(approx_cv, na.rm = TRUE))) |> 
    select(-year) |> 
    as.matrix() |> 
    unname()

initial_index_sd_scale <- rep(1, ncol(index))
map_index_sd_scale <- seq.int(ncol(index)) # use to try and estimate observation error
map_index_sd_scale[] <- NA 

index_Neff <- t(replicate(length(years), rep(50, ncol(index))))
index_fracyr <- t(replicate(length(years), c(5 / 12, 9 / 12, 9 / 12, 5 / 12, 0)))
units_indices <- c(1, 1, 1, 1, 2)
units_index_paa <- rep(2, ncol(index))
use_indices <- index > 0

index_vals <- iaa |> select(I1:I6) |> as.matrix()
index_vals[is.na(index_vals)] <- 0
ipaa <- prop.table(index_vals, 1)
ipaa[is.nan(ipaa)] <- 0
index_paa <- array(0, dim = c(ncol(index), length(years), length(ages)))
for(y in 1:length(years)){
    for(a in 1:length(ages)){
        index_paa[1,y,a] <- ifelse(is.na(ipaa[y,a]), 0, ipaa[y,a])
        index_paa[2,y,a] <- ifelse(a == 1, 1, 0)
        index_paa[3:5,y,a] <- ifelse(a == 1, 0, 1)
    }
}
use_index_paa <- apply(index_paa, 1, rowSums) > 0
use_index_paa[, 2:5] <- FALSE # larval survey represents one age; others lack age data

selblock_pointer_indices <- t(replicate(length(years), c(2, 3, 4, 5, 6)))

## Approximate mean weight at age (kg) via biomass at age (tonnes) / numbers at age (thousands)
waa_vals <- (baa[, -1] * 1000) / (iaa[, -1] * 1000000)
waa_vals <- as.matrix(waa_vals)
waa_vals[waa_vals == 0] <- NA
waa_vals[is.nan(waa_vals)] <- NA
for (a in colnames(waa_vals)) {
    ind <- is.na(waa_vals[, a])
    waa_vals[ind, a] <- mean(waa_vals[, a], na.rm = TRUE) 
}
waa <- array(waa_vals, dim = c(1, length(years), length(ages)))

mat_vals <- data.frame(year = 1985:2024, Cap_data$matMp) |> 
    right_join(data.frame(year = years), by = "year") |> 
    arrange(year) |> 
    select(-year) |> 
    as.matrix() |> 
    unname()
for (i in seq(ncol(mat_vals))) {
    ind <- is.na(mat_vals[, i])
    mat_vals[ind, i] <- mean(mat_vals[, i], na.rm = TRUE)
}

mat <- array(cbind(0, mat_vals, 1, 1), dim = c(1, length(years), length(ages)))
waa_pointer_indices <- rep(1, ncol(index))
waa_pointer_fleets <- 1
waa_pointer_ssb <- 1
fracyr_ssb <-  matrix(0, ncol = 1, nrow = length(years))



## Relates to numbers at age. Now using biomass estimates. Assuming that paa based on numbers is representative.
# ## Does the scale make sense? (convert index numbers at age from billions to thousands, and multiply by weight at age in kg)
# index_tonnes <- rowSums(index_vals * 1000000 / 2 * waa[1,,])
# plotly::plot_ly(y = index_tonnes, x = years) |> plotly::add_lines()
# ## sum of index at age multiplied by weight at age is double the biomass index. Don't know why.
# 
# index[, 1] <- index[, 1] / 2 # TEMP FIX: HALF THE CURRENT INDEX



selectivity <- list(model = rep("age-specific", ncol(index) + 1), 
                    n_selblocks = ncol(index) + 1,
                    re = c("2dar1", rep("none", ncol(index))),
                    initial_pars = list(c(0, rep(0.1, length(ages) - 1)), # Catches
                                        c(0.5, rep(1, length(ages) - 1)), # Acoustics
                                        c(0, rep(1, length(ages) - 1)),
                                        c(0, rep(1, length(ages) - 1)),
                                        c(0, rep(1, length(ages) - 1)),
                                        c(1, rep(0, length(ages) - 1))), # Larval
                    map_pars = list(c(NA, rep(1, length(ages) - 1)), 
                                    c(2, rep(NA, length(ages) - 1)), 
                                    c(rep(NA, length(ages))),
                                    c(rep(NA, length(ages))),
                                    c(rep(NA, length(ages))),
                                    c(rep(NA, length(ages))))
                    )

max_age <- 6
4.899 * max_age ^ (-0.916) # approx = 1; Then et al. (2015)
M_in <- list(initial_MAA = array(1, dim = c(1, 1, length(years), length(ages))))

x <- rnorm(100000, sd = 1.5)
hist(x, breaks = 200, col = "grey", border = "grey")
hist(plogis(x), breaks = 200, col = "grey", border = "grey")
q_in <- list(q_upper = rep(1, ncol(index)),
             initial_q = rep(0.5, ncol(index)),
             prior_sd = c(2, rep(NA, ncol(index) - 1))) # survey q should not exceed 1

F_in <- list(
    F = cbind(rep(2, length(years))),
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
  waa_pointer_indices = waa_pointer_indices,
  initial_index_sd_scale = initial_index_sd_scale,
  map_index_sd_scale = map_index_sd_scale
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
                  list(N1_model = "equilibrium",
                       recruitment_model = 2,
                       sigma = "rec+1",
                       cor = "2dar1"))

fit2 <- fit_wham(input2, do.fit = T, do.retro = F, do.brps = F, do.osa = F, do.sdrep = T)
fit2$fn()
fit2$opt
fit2$sdrep
round(fit2$rep$NAA[1,1,,], 2)
matplot(fit2$rep$NAA[1,1,,], type = "l")

# fit2 <- make_osa_residuals(fit2)
# fit2$peels <- retro(fit2)

# plot_wham_output(fit2, res = 600, dir.main = file.path(getwd(), "capelin", "fit2"))


## Estimate M ----------------------------------------------------------------

input3 <- set_NAA(input1, 
                  list(N1_model = "equilibrium",
                       recruitment_model = 2,
                       sigma = "rec+1",
                       cor = "2dar1")) |> 
    set_M(list(mean_model = "estimate-M",
               re_model = matrix("none")))

fit3 <- fit_wham(input3, do.fit = T, do.retro = F, do.brps = F, do.osa = F, do.sdrep = T)
fit3$opt
fit3$sdrep

# plot_wham_output(fit3, res = 600, dir.main = file.path(getwd(), "capelin", "fit3"))


## S-R ----------------------------------------------------------------

input4 <- set_NAA(input1, 
                  list(N1_model = "equilibrium",
                       recruit_model = 3,
                       sigma = "rec+1",
                       cor = "ar1_y"))

fit4 <- fit_wham(input4, do.fit = T, do.retro = F, do.brps = F, do.osa = F, do.sdrep = T)
fit4$opt
fit4$sdrep

# fit4 <- make_osa_residuals(fit4)
# fit4$peels <- retro(fit4)

# plot_wham_output(fit4, res = 600, dir.main = file.path(getwd(), "capelin", "fit4"))


## Maturity effect -------------------------------------------------------------

input5 <- set_NAA(input1, 
                  list(N1_model = "equilibrium",
                       recruitment_model = 2,
                       sigma = "rec+1",
                       cor = "2dar1")) |> 
    set_ecov(list(
        label = "Maturity effect",
        mean = mat[1,,],
        logsigma = matrix(0.01, nrow = length(years), ncol = length(ages)),
        year = years,
        use_obs = mat[1,,] >= 0,
        process_model = "rw",
        M_how = array("lag-0-linear", dim = c(length(ages), 1, length(ages), 1))
    )) |> 
    set_M(list(
        mean_model = "estimate-M",
        initial_MAA = array(1, dim = c(1, 1, length(years), length(ages)))
    ))
input5$map$Ecov_beta_M <- factor(rep(1, length(input5$map$Ecov_beta_M)))

fit5 <- fit_wham(input5, do.fit = T, do.retro = F, do.brps = F, do.osa = F, do.sdrep = T)
fit5$opt
fit5$sdrep

# plot_wham_output(fit5, res = 600, dir.main = file.path(getwd(), "capelin", "fit5"))

