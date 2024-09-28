
library(dplyr)

## Cumulative PAA
approx_paa <- readxl::read_xlsx("wham/data/approx_catch_paa_from_figure.xlsx", sheet = "approx_data")

approx_paa$year <- round(approx_paa$year)
approx_paa <- approx_paa |> 
    group_by(year, age) |> 
    summarise(prop = mean(prop)) |> 
    filter(year <= 2018)

approx_paa_mat <- xtabs(prop ~ year + age, data = approx_paa)
ind5 <- approx_paa_mat[, "5"] == 0
approx_paa_mat[ind5, "5"] <- 1
approx_paa_mat[, "6"] <- 1

# pre_collapse_means <- colMeans(approx_paa_mat[as.character(1980:1990), ])
# post_collapse_means <- colMeans(approx_paa_mat[as.character(2010:2018), ])
pre_collapse_means <- post_collapse_means <- rep(0, ncol(approx_paa_mat))

approx_paa_mat <- rbind(t(replicate(8, pre_collapse_means)),
                        approx_paa_mat,
                        t(replicate(4, post_collapse_means)))
rownames(approx_paa_mat) <- 1972:2022

approx_paa_mat <- cbind(0, t(apply(approx_paa_mat, 1, diff)))
colnames(approx_paa_mat) <- 1:6



