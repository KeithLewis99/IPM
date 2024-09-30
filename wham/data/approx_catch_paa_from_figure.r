
library(dplyr)

## Cumulative PAA extracted from Figure 6 in Research Document 2023/076 
## using WebPlotDigitizer (https://automeris.io/wpd/)
approx_paa <- read.csv("wham/data/approx_catch_paa_from_figure.csv")
approx_paa$prop[approx_paa$prop > 0.997] <- 1

approx_paa$year <- round(approx_paa$year)
approx_paa <- approx_paa |> 
    group_by(year, age) |> 
    summarise(prop = mean(prop)) |> 
    filter(year <= 2019)

approx_paa <- rbind(data.frame(year = sort(unique(approx_paa$year)), age = 1, prop = 0),
                    approx_paa)

approx_paa_mat <- xtabs(prop ~ year + age, data = approx_paa)
ind5 <- approx_paa_mat[, "5"] == 0
approx_paa_mat[ind5, "5"] <- 1
approx_paa_mat[, "6"] <- 1
# approx_paa_mat <- round(approx_paa_mat, 4)

# pre_collapse_means <- colMeans(approx_paa_mat[as.character(1980:1990), ])
# post_collapse_means <- colMeans(approx_paa_mat[as.character(2010:2018), ])
pre_collapse_means <- post_collapse_means <- rep(0, ncol(approx_paa_mat))

approx_paa_mat <- rbind(t(replicate(8, pre_collapse_means)),
                        approx_paa_mat,
                        t(replicate(3, post_collapse_means)))
rownames(approx_paa_mat) <- 1972:2022

approx_paa_mat <- cbind(0, t(apply(approx_paa_mat, 1, diff)))
colnames(approx_paa_mat) <- 1:6



