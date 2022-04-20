## Start----

rm(list=ls())
# load data----

# Source files
source("IPM_dat.R")
source("IPM_fun.R")

##disaggregated data----

# age breakdown by strata
ggplot(data = df_dis_strata, aes(x = year, y = abundance, colour = stratum)) + geom_point() + facet_grid(age ~ .)

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


# this just shows that there's no real relation between any of the age classes and therefore, the Schaekell approach won't work very well.
p <- ggplot(data = df_dis_tab, aes(x = I2, y = lead(I3))) + geom_point()
p
p <- ggplot(data = df_dis_tab, aes(x = I2, y = lead(I4))) + geom_point()
p
p <- ggplot(data = df_dis_tab, aes(x = I2, y = lead(I4, 2))) + geom_point()
p
p <- ggplot(data = df_dis_tab, aes(x = I2, y = I)) + geom_point()
p
p <- ggplot(data = df_dis_tab, aes(x = I2, y = lead(I))) + geom_point()
p
p <- ggplot(data = df_dis_tab, aes(x = I3, y = lead(I4))) + geom_point()
p
p <- ggplot(data = df_dis_tab, aes(x = I3, y = lead(I))) + geom_point()
p

# I2 v I3
p <- ggplot(jd, aes(x = log10(I2) , y = lead(log10(I3), 1), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p)

#library(ggcorrplot)
#corr <- round(cor(df_dis_tab[, 2:6], use = "complete.obs", method = c("pearson")), 2)

#ggcorrplot(corr, hc.order = T, type = "lower", lab = T)



# Just want to see the relative contribution of each age class
df_dis_summF <- df_dis_summ %>%
     filter(age !=1)

p <- ggplot(df_dis_summF, aes(x = year, y = abun, fill = age))
p <- p + geom_bar(stat = "identity")
p

p <- ggplot(df_dis_summF, aes(x = year, y = biomass, fill = age))
p <- p + geom_bar(stat = "identity")
p

## maturity ----
# plot maturity data
plot(df_mat$year, df_mat$mat)
ggplot(df_mat, aes(x = year, y = mat)) + geom_point() + theme_bw()


## larval density ----
# graph larval data time series
ggplot(df_ld, aes(x = year, y = larvae)) + geom_point() + theme_bw()

# calculation of priors for LD
range(df_ld$larvae, na.rm = T)

# 2018 doesn't match but condition the best ever that year, 2021 doesn't match but survey incomplete, 2010 doesn't match but survey known to be low.keith
jd$year <- 1999:2023

ggplot(jd, aes(x = lag(LD,2) , y = I2)) + geom_point()

# LD v I2
p <- ggplot(jd, aes(x = lag(log10(LD),2) , y = log10(I2), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p)

# LD v I3
p <- ggplot(jd, aes(x = log10(LD) , y = lead(log10(I3), 1), text = paste("year: ", year, "\n", sep = ""))) + geom_point()
ggplotly(p)



