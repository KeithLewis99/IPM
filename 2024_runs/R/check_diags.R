data <- all.list

#pAbundance-Index residuals
# png(here('2024_runs/output/eps_residuals.png'), width = 1000, height = 800, res = 100)
ls_all$ls_med$eps |>
   as.data.frame() |>
   rename('Age 1' = V1, 'Age 2' = V2, 'Age 3' = V3) |>
   mutate(Year = 1985:2024) |>
   pivot_longer(cols = 'Age 1':'Age 3', names_to = 'Age', values_to = 'Eps') |>
   ggplot(aes(x = Year, y = Eps)) +
   geom_point() +
   # geom_text(aes(x = 2024, y = 2, label = round(sum(Eps), 2))) +
   geom_hline(yintercept = 0, lty = 2) +
   facet_wrap(~Age, ncol = 1)
# dev.off()

## Parameter Convergence ##
# png(here('2024_runs/output/alpha_beta.png'), width = 1000, height = 800, res = 100)
par(mfrow = c(2,1), mar = c(2, 4, 1, .4))
raw$alpha |> plot(type = 'l', ylab = 'alpha')
raw$beta |> plot(type = 'l', ylab = 'beta')
# dev.off()

# png(here('2024_runs/output/delta.png'), width = 1000, height = 800, res = 100)
par(mfrow = c(3,1), mar = c(2, 4, 1, .4))
raw$delta[,1] |> plot(type = 'l', ylab = '')
raw$delta[,2] |> plot(type = 'l', ylab = 'delta')
raw$delta[,3] |> plot(type = 'l', ylab = '')
# dev.off()

# png(here('2024_runs/output/epsilon.png'), width = 1000, height = 800, res = 100)
raw$epsilon[,1] |> plot(type = 'l', ylab = '')
raw$epsilon[,2] |> plot(type = 'l', ylab = 'epsilon')
raw$epsilon[,3] |> plot(type = 'l', ylab = '')
# dev.off()

# png(here('2024_runs/output/gamma.png'), width = 1000, height = 800, res = 100)
raw$gamma[,1] |> plot(type = 'l', ylab = '')
raw$gamma[,2] |> plot(type = 'l', ylab = 'gamma') # this has a trend
raw$gamma[,3] |> plot(type = 'l', ylab = '')
# dev.off()

# png(here('2024_runs/output/resid_acf.png'), width = 1000, height = 800, res = 100)
acf(ls_all$ls_med$eps)
# dev.off()



## OSA Function Was HERE ##
# osaResidPlot(raw$posa, 'Cohort', byage = F)

png(here('2024_runs/output/posa_resid_cohort.png'), width = 1000, height = 800, res = 100)
osaResidPlot(raw$posa, 'Cohort', byage = T)
dev.off()

png(here('2024_runs/output/osa_resid_cohort.png'), width = 1000, height = 800, res = 100)
osaResidPlot(raw$osa, 'Cohort', byage = T)
dev.off()

png(here('2024_runs/output/posa_resid_year.png'), width = 1000, height = 800, res = 100)
osaResidPlot(raw$posa, 'Year')
dev.off()

png(here('2024_runs/output/osa_resid_year.png'), width = 1000, height = 800, res = 100)
osaResidPlot(raw$osa, 'Year')
dev.off()

   

# png(here('2024_runs/output/Z0.png'), width = 1000, height = 800, res = 100)
par(mfrow = c(1,1), mar = c(6,3,2,2))
c('Z0' = ls_all$ls_med$Z0) |>
   barplot(las = 2)
abline(h = 0, lty = 1)
# dev.off()




# c(Mu = raw$tau.proc |> mean(), 
#   Med = raw$tau.proc |> median(), 
#   Upper = raw$tau.proc |> quantile(.95),
#   Lower = raw$tau.proc |> quantile(.05)
# )
# c(Mu = raw$tau.obs |> mean(), 
#   Med = raw$tau.obs |> median(), 
#   Upper = raw$tau.obs |> quantile(.95),
#   Lower = raw$tau.obs |> quantile(.05)
# )

# png(here('2024_runs/output/compare_taus.png'), width = 1000, height = 800, res = 100)
data.frame('Tau proc' = raw$tau.proc[raw$tau.proc < quantile(raw$tau.proc, .90) & raw$tau.proc > quantile(raw$tau.proc, .1)], 
           'Tau obs' = raw$tau.obs[raw$tau.obs < quantile(raw$tau.obs, .90) & raw$tau.obs > quantile(raw$tau.obs, .1)]) |>
boxplot()
# dev.off()

raw$tau.obs|> boxplot()
raw$tau.proc|> boxplot()

#Abundance correlations
c(ls_all$ls_med$rec, NA) |>
   cbind(ls_all$ls_med$N) |>
   as.data.frame() |>
   rename('Age 0' = V1, 'Age 1' = V2, 'Age 2' = V3, 'Age 3' = V4) |>
   as.data.frame() |>
   GGally::ggpairs()

#OSA correlations
# png(here('2024_runs/output/osa_corrs.png'), width = 1000, height = 800, res = 100)
ls_all$ls_med$osa |>
   as.data.frame() |>
   GGally::ggpairs()
# dev.off()

# png(here('2024_runs/output/osa_qq.png'), width = 1000, height = 800, res = 100)
ls_all$ls_med$osa |>
   as.data.frame() |>
   pivot_longer(everything(), names_to = 'N', values_to = 'osa') |>
   mutate(N = gsub('V', '', N), N = factor(as.numeric(N)+1)) |>
   ggplot(aes(sample = osa, col = N)) +
   stat_qq() + geom_qq_line()
# dev.off()

#POSA residuals
# png(here('2024_runs/output/osa_corrs.png'), width = 1000, height = 800, res = 100)
ls_all$ls_med$posa |>
   as.data.frame() |>
   GGally::ggpairs()
# dev.off()

# png(here('2024_runs/output/osa_qq.png'), width = 1000, height = 800, res = 100)
ls_all$ls_med$posa |>
   as.data.frame() |>
   pivot_longer(everything(), names_to = 'N', values_to = 'osa') |>
   mutate(N = gsub('V', '', N), N = factor(as.numeric(N)+1)) |>
   ggplot(aes(sample = osa, col = N)) +
   stat_qq() + geom_qq_line()
# dev.off()

#Some density plots
old.par <- par()

# png(here('2024_runs/output/density_Z0.png'), width = 1000, height = 800, res = 100)
par(old.par)
plot(density(raw$Z0[,1]), col = 'darkolivegreen3', ylim = c(0, 1.5), main = expression(Z[0]))
lines(density(raw$Z0[,2]), col = 'coral')
lines(density(raw$Z0[,3]), col = 'steelblue')
legend('topright', c('Age 1', 'Age 2', 'Age 3'),
       col = c('darkolivegreen3', 'coral', 'steelblue'),
       lty = c(2,2,2))
# dev.off()

# png(here('2024_runs/output/density_alpha.png'),  width = 1000, height = 800, res = 100)
plot(density(data$raw$alpha), main = expression(alpha))
abline(v = 0, lty = 2)
# dev.off()
# png(here('2024_runs/output/density_beta.png'), width = 1000, height = 800, res = 100)
plot(density(data$raw$beta), main = expression(beta))
abline(v = 0, lty = 2)
# dev.off()
# png(here('2024_runs/output/density_gamma.png'), width = 1000, height = 800, res = 100)
plot(density(data$raw$gamma[,1]), col = 'darkolivegreen3', ylim = c(0, 1.5), 
     main = expression(gamma))
lines(density(data$raw$gamma[,2]), col = 'coral')
lines(density(data$raw$gamma[,3]), col = 'steelblue')
abline(v = 0, lty = 2)
legend('topright', c('Age 2', 'Age 3', 'Age 4'),
       col = c('darkolivegreen3', 'coral', 'steelblue'),
       lty = c(2,2,2), cex = .8)
# dev.off()
# png(here('2024_runs/output/density_delta.png'), width = 1000, height = 800, res = 100)
plot(density(data$raw$delta[,1]), col = 'darkolivegreen3', ylim = c(0, 1.0), 
     main = expression(delta))
lines(density(data$raw$delta[,2]), col = 'coral')
lines(density(data$raw$delta[,3]), col = 'steelblue')
abline(v = 0, lty = 2)
legend('topright', c('Age 2', 'Age 3', 'Age 4'),
       col = c('darkolivegreen3', 'coral', 'steelblue'),
       lty = c(2,2,2))
# dev.off()
# png(here('2024_runs/output/density_epsilon.png'), width = 1000, height = 800, res = 100)
plot(density(data$raw$epsilon[,1]), col = 'darkolivegreen3', 
     ylim = c(0, 2.5), xlim = c(-1, 1), 
     main = expression(epsilon))
lines(density(data$raw$epsilon[,2]), col = 'coral')
lines(density(data$raw$epsilon[,3]), col = 'steelblue')
abline(v = 0, lty = 2)
legend('topright', c('Age 2', 'Age 3', 'Age 4'),
       col = c('darkolivegreen3', 'coral', 'steelblue'),
       lty = c(2,2,2), cex = .8)
# dev.off()

# png(here('2024_runs/output/density_esp8504.png'), width = 600, height = 1400, res = 100)
# par(mfrow = c(7, 3), mar = c(3 ,2 ,3,1))
# for(t in 1:20){ #i.e. residuals
#    plot(density(raw$eps[,t,1]), col = 'darkolivegreen3', 
#         ylab = '', xlab = '', main = 1984+t)
#    lines(density(raw$eps[,t,2]), col = 'coral')
#    lines(density(raw$eps[,t,3]), col = 'steelblue')
# }
# plot(NULL)
# legend('center', c('Age 1', 'Age 2', 'Age 3'),
#        col = c('darkolivegreen3', 'coral', 'steelblue'),
#        lty = c(2,2,2), cex = 2)
# dev.off()
# 
# png(here('2024_runs/output/density_eps0524.png'), width = 600, height = 1400, res = 100)
# par(mfrow = c(7, 3), mar = c(3 ,2 ,3,1))
# for(t in 21:40){ #i.e. residuals
#    plot(density(raw$eps[,t,1]), col = 'darkolivegreen3', 
#         ylab = '', xlab = '', main = 1984+t)
#    lines(density(raw$eps[,t,2]), col = 'coral')
#    lines(density(raw$eps[,t,3]), col = 'steelblue')
# }
# plot(NULL)
# legend('center', c('Age 1', 'Age 2', 'Age 3'),
#        col = c('darkolivegreen3', 'coral', 'steelblue'),
#        lty = c(2,2,2), cex = 2)
# dev.off()

# png(here('2024_runs/output/density_sigmaZ.png'), width = 1000, height = 800, res = 100)
par(old.par)
plot(density(raw$sigma.Z))
# dev.off()
# png(here('2024_runs/output/density_sigmaRec.png'), width = 1000, height = 800, res = 100)
plot(density(raw$sigma.rec))
# dev.off()
# png(here('2024_runs/output/density_tauObs.png'), width = 1000, height = 800, res = 100)
plot(density(raw$tau.obs))
plot(density(1/sqrt(raw$tau.obs)))
# dev.off()
# png(here('2024_runs/output/density_tauProc.png'), width = 1000, height = 800, res = 100)
plot(density(all.list$raw$sigma.procR))
plot(density(1/sqrt(all.list$raw$tau.proc)))
# dev.off()

