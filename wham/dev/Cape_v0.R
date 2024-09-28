
## TODO
## - expand to 5+


### Capelin ###
#-------------#
library(wham)

library(ggplot2)
library(dplyr)
library(tidyr)
Cap_data <- readRDS("capelin/CapData.rds")

years <- as.integer(1985:2022)
ages <- as.integer(2:4)

landings <- read.csv("capelin/capelin_landings_kt.csv") |> 
    subset(Year %in% years)

# Cap_data$landings <- Cap_data$landings[Cap_data$landings$year < 2022,]

catch_vals <- Cap_data$matCAA[seq_along(years), ] %>% replace_na(0)
catch <- matrix((catch_vals/1000) %>% apply(1, sum), 
                ncol = 1, nrow = length(years))
catch_cv <- matrix(c(rep(1, 7), rep(0.5, 8), rep(.1, 23)), 
                   ncol = 1, nrow = length(years))
catch_Neff <- matrix(50, ncol = 1, nrow = length(years))
cpaa <- catch_vals %>% 
  as.data.frame() %>% 
  mutate(total = c2+c3+c4,
         c2 = c2/total,
         c3 = c3/total,
         c4 = c4/total) %>% 
  select(-total)
catch_paa <- array(0, dim = c(1, length(years), length(ages)))
use_catch_paa <- matrix(1, nrow = length(years), ncol = 1)
use_catch_paa[c(5),] <- 0
selblock_pointer_fleets <- matrix(1, nrow = length(years))


index_vals <- Cap_data$matI[seq_along(years), ]
index <- (rowSums(exp(index_vals))) |> replace_na(0) %>% 
    matrix(nrow = length(years), ncol = 2, byrow = F)
index_cv <- matrix(.1, ncol = 2, nrow = length(years))
index_Neff <- matrix(50, ncol = 2, nrow = length(years))
index_fracyr <- matrix(5/12, ncol = 2, nrow = length(years))
units_indices <- 2 # 2 is abundance, 1 is biomass indices??
units_index_paa <- 2
use_indices <- matrix(index > 0, ncol = 2, nrow = length(years))
use_index_paa <- matrix(index > 0, ncol = 2, nrow = length(years))

## juveniles ##
index[,2] <- Cap_data$LD[seq_along(years)] %>% replace_na(0)
index_cv[,2] <- 1
# index_Neff[,2] <- 50
# index_fracyr[,2] <- 5/12
units_indices <- c(2, 2)
units_index_paa <- c(2, 2)
use_indices[,2] <- index[,2] > 0
use_index_paa[,2] <- index[,2]  > 0

ipaa <- index_vals %>% 
  exp() %>% 
  as.data.frame() %>% 
  mutate(total = I2+I3+I4,
         I2 = I2/total,
         I3 = I3/total,
         I4 = I4/total) %>% 
  select(-total)
index_paa <- array(0, dim = c(2,length(years),length(ages)))
selblock_pointer_indices <- matrix(rep(c(2,3), each = length(years)), 
                                   ncol = 2, nrow = length(years))

for(y in 1:length(years)){
  for(a in 1:length(ages)){
    catch_paa[1,y,a] <- ifelse(is.na(cpaa[y,a]), 0, cpaa[y,a])
    index_paa[1,y,a] <- ifelse(is.na(ipaa[y,a]), 0, ipaa[y,a])
    index_paa[2,y,a] <- ifelse(a == 1, 1, 0)
  }
}


waa_vals <- xtabs(w ~ year+age, data = Cap_data$WaA, 
              subset = age %in% ages) %>% matrix(ncol = 3)
waa_ex <- apply(waa_vals, 2, \(x)mean(x, na.rm=T))
waa_full <- rbind(matrix(waa_ex, nrow = 13, ncol = 3, byrow = T), 
                  waa_vals,
                  matrix(waa_ex, nrow = 1, ncol = 3, byrow = T))
waa <- array(0, dim = c(3, length(years), length(ages)))
waa[1,,] <- matrix(waa_full / 1000, nrow = length(years), ncol = length(ages))
waa[2,,] <- 1
waa[3,,] <- 1

mat <- array(Cap_data$matMp[seq_along(years), ], dim = c(1, length(years), length(ages)))
waa_pointer_indices <- c(1,3)
waa_pointer_fleets <- 1 # all weights set to 1 b/c we have catch abundances
waa_pointer_ssb <- 1
fracyr_ssb <-  matrix(0, ncol = 1, nrow = length(years))


selectivity <- list(model = c("age-specific", "age-specific"), #, "age-specific"), 
                    n_selblocks = 2, # 3,
                    re = c("2dar1", "none"),
                    initial_pars = list(c(.1, .1, .1), #Catches
                                        c(1, 1, 1)), #, #Acoustics
                                        # c(.1, 0, 0)), #Juveniles
                    map_pars = list(c(1,1,1), c(NA,NA,NA)) #, c(6, NA, NA))
                    )

NAA_in <- list(N1_model = "age-specific",
               N1_pars = c(4.9, 3.2, 1.9),
               recruit_model = 2,
               sigma = "rec")

F_in <- list(
  F = cbind(rep(5, length(years))),
  map_F = cbind(rep(NA, length(years))))


max_age <- 6
4.899 * max_age ^ (-0.916) # Then et al. (2015)
M_in <- list(initial_MAA = array(1, dim = c(1,1,length(years),length(ages)) ))

q_in <- list(q_upper = 1)

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
  agg_catch = as.matrix(landings[, 2]) * 1000,
  agg_catch_cv = catch_cv,
  catch_paa = catch_paa,
  use_catch_paa = use_catch_paa,
  catch_Neff = catch_Neff,
  selblock_pointer_fleets = selblock_pointer_fleets,
  waa_pointer_fleets = waa_pointer_fleets
)

index_info <- list(
  n_indices = 1L, # NCOL(index),
  agg_indices = index[, 1, drop = FALSE],
  units_indices = units_indices[1],
  units_index_paa = units_index_paa[1],
  agg_index_cv = index_cv[, 1, drop = FALSE],
  fracyr_indices = index_fracyr[, 1, drop = FALSE],
  use_indices = use_indices[, 1, drop = FALSE],
  use_index_paa = use_index_paa[, 1, drop = FALSE],
  index_paa = array(index_paa[1,,], dim = c(1, length(years), length(ages))),
  index_Neff = index_Neff[, 1, drop = FALSE],
  selblock_pointer_indices = selblock_pointer_indices[, 1, drop = FALSE],
  waa_pointer_indices = waa_pointer_indices[1]
)


# ecov_info <- set_ecov()

input_all <- prepare_wham_input(basic_info = basic_info, 
                                selectivity = selectivity, 
                                catch_info = catch_info, 
                                index_info = index_info, 
                                M = M_in, F = F_in,
                                catchability = q_in) 

fit <- fit_wham(input_all, do.fit = F, do.retro = F, do.brps = F, do.osa = F)
fit$fn()
fit$rep$NAA[1,1,,]

fit <- fit_wham(input_all, do.fit = T, do.retro = F, do.brps = F, do.osa = F, do.sdrep = T)
fit$fn()
fit$opt
fit$sdrep
round(fit$rep$NAA[1,1,,], 2)
matplot(fit$rep$NAA[1,1,,], type = "l")

plot_wham_output(fit, res = 600, dir.main = file.path(getwd(), "capelin"))

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

plot_wham_output(fit2, res = 600, dir.main = file.path(getwd(), "capelin"))

## M deviations ----------------------------------------------------------------

input3 <- set_NAA(input1, 
                  list(N1_model = "equilibrium",
                       recruitment_model = 2,
                       sigma = "rec+1",
                       cor = "2dar1")) |> 
    set_M(list(mean_model = "estimate-M",
               re_model = matrix("none"))) |> 
    set_q(list(q_upper = 20))

fit3 <- fit_wham(input3, do.fit = T, do.retro = F, do.brps = F, do.osa = F, do.sdrep = T)
fit3$opt
fit3$sdrep

fit3 <- make_osa_residuals(fit3)
fit3$peels <- retro(fit3)

plot_wham_output(fit3, res = 600, dir.main = file.path(getwd(), "capelin"))


