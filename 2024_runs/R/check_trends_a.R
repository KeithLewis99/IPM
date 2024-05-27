years <- 1985:2024

## Main Trend ##
# png(here('2024_runs/output/Index_trend.png'), width = 1000, height = 800, res = 100)
IPMplot(raw, jags.data.m$matI)
# dev.off()

par(mfrow = c(1,1), mar = c(5, 5, 1, 1))
#All three together
# png(here('2024_runs/output/compare_trends.png'), width = 1000, height = 800, res = 100)
{
matplot(matrix(years, nrow = 40, ncol = 3), 
        ls_all$ls_med$Np, type = 'l',
        xlab = 'Year', ylab = 'Index')
matlines(matrix(years, nrow = 40, ncol = 3),
         ls_all$ls_med$N, type = 'b', pch = 3, lty = 2)
matpoints(matrix(years, nrow = 40, ncol = 3),
          jags.data.m$matI, pch = 24)
legend('topright', 
       legend = c('Processed', 'Flat', 'Observed'), 
       pch = c(NA, 3, 24),
       lty = c(1, 2, NA))
}
# dev.off()


#Flat with Processed
matplot(matrix(years, nrow = 40, ncol = 3), xlab = 'Year',
        ls_all$ls_med$Np, type = 'l', lty = 1)
matlines(matrix(years, nrow = 40, ncol = 3),
         ls_all$ls_med$N, type = 'b', pch = 3, lty = 2)
legend('topright', 
       legend = c('Processed', 'Flat'), 
       pch = c(NA, 3),
       lty = c(1, 2))

#Processed with Observed
matplot(matrix(years, nrow = 40, ncol = 3), xlab = 'Year',
        ls_all$ls_med$Np, type = 'l', lty = 1)
matpoints(matrix(years, nrow = 40, ncol = 3),
          jags.data.m$matI, pch = 24)
legend('topright', 
       legend = c('Processed', 'Observed'), 
       pch = c(NA, 24),
       lty = c(1, NA))

#Flat with Observed
matplot(matrix(years, nrow = 40, ncol = 3), xlab = 'Year',
        ls_all$ls_med$N, type = 'b', pch = 3, lty = 1)
matpoints(matrix(years, nrow = 40, ncol = 3),
          jags.data.m$matI, pch = 24)
legend('topright', 
       legend = c('Flat', 'Observed'), 
       pch = c(3, 24),
       lty = c(2, NA))



#Flat Abundances
# ls_all$ls_med[c('N2', 'N3', 'N4')] |>
#    as.data.frame() |>
#    mutate(Year = years) |>
#    pivot_longer(cols = N2:N4, names_to = 'Age', values_to = 'N') |>
#    ggplot(aes(x = Year)) +
#    geom_line(aes(y = N, col = Age))
# 
# #Process-ed Abundances
# ls_all$ls_med$Np |>
#    as.data.frame() |>
#    mutate(Year = years) |>
#    rename('Age 1' = V1, 'Age 2' = V2, 'Age 3' = V3) |>
#    pivot_longer(cols = 'Age 1':'Age 3', names_to = 'Age', values_to = 'N') |>
#    ggplot(aes(x = Year)) +
#    geom_line(aes(y = N, col = Age))
# 
# #Projected Indices (almost identical to N)
# ls_all$ls_med[c('I2.rep', 'I3.rep', 'I4.rep')] |>
#    as.data.frame() |>
#    mutate(Year = years) |>
#    pivot_longer(cols = I2.rep:I4.rep, names_to = 'Age', values_to = 'N') |>
#    ggplot(aes(x = Year)) +
#    geom_line(aes(y = N, col = Age))


#Comparison of N to I -- strange offset
# ls_all$ls_med$I.rep |>
#    as.data.frame() |>
#    rename('Value' = `ls_all$ls_med$I.rep`) |>
#    mutate(Parm = 'I') |>
#    rbind(
#       ls_all$ls_med$N |>
#       as.data.frame() |>
#       mutate(Parm = 'N',
#              Value = log(exp(V1)+exp(V2)+exp(V3))
#              )|>
#          select(Parm, Value)
#    ) |>
#    mutate(Year = rep(years, 2)) |>
#    mutate(lty = factor(ifelse(Parm == 'N', NA, 3)),
#           pch = factor(ifelse(Parm == 'N', 20, NA))) |>
#    ggplot(aes(x = Year, y = Value, col = Parm)) +
#    geom_point(aes(pch = pch))+
#    geom_line(aes(lty = lty))+
#    guides(lty = 'none', pch = 'none')+
#    theme_bw()


#SPAY
# png(here('2024_runs/output/spay.png'), width = 1000, height = 800, res = 100)
spay_plot(data)
# dev.off()

#Recruitment
# png(here('2024_runs/output/rec.png'), width = 1000, height = 800, res = 100)
recruitment_plot(data)
# dev.off()

#Maturities
# png(here('2024_runs/output/mat.png'), width = 1000, height = 800, res = 100)
maturity_plot(data)
# dev.off()

