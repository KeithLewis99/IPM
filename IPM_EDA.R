## Start----

rm(list=ls())
# load data----

# Source files
#source("IPM_dat.R")

##disaggregated data----


ggplot(data = df_dis_strata, aes(x = year, y = abundance, colour = stratum)) + geom_point() + facet_grid(age ~ .)

# Check amounts and percents of Unknown biomass
check1 <- df_dis %>%
     group_by(year) %>%
     filter(age != "Unknown") %>%
     summarize(biomass = sum(biomass))
check1

check2 <- df_dis %>%
     group_by(year) %>%
     filter(age == "Unknown") %>%
     summarize(biomass = sum(biomass))
check2 # 5.7kt of Unknown in 1999

check3 <- left_join(check1, check2, by = "year") %>%  
     mutate((biomass.y/biomass.x)*100)
check3 # so Unknown biomass is usually <1% of biomass but in 1999 it was ~ 3%


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



#library(ggcorrplot)
#corr <- round(cor(df_dis_tab[, 2:6], use = "complete.obs", method = c("pearson")), 2)
#ggcorrplot(corr, hc.order = T, type = "lower", lab = T)



# Just want to see the relative contribution of each age class
df_dis_summF <- df_dis_summ %>%
     filter(age !=1)

p <- ggplot(df_dis_summF, aes(x = year, y = abun, fill = age))
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


