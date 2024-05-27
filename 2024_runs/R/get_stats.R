
get_stats <- function(jagsout){
   
   samples <- jags.samples(jagsout$model,
                            c('WAIC', 'deviance', 'nll'),
                            type = 'mean',
                            n.iter = 10000,
                            n.burnin = 1000,
                            nthin = 10)

   pWAIC <- samples$WAIC
   WAIC <- samples$deviance + pWAIC
   
   AICs <- c(deviance = sum(samples$deviance),
             pD = jagsout$BUGSoutput$pD,
             DIC = jagsout$BUGSoutput$DIC,
             pWAIC = sum(pWAIC),
             WAIC = sum(WAIC),
             nll = sum(jagsout$BUGSoutput$median$nll)) |>
      round(1)
   
   AICs

}

# stats <- get_stats(jagsout)

# jags.samples(jagsout$model,
#              c('WAIC', 'deviance', 'nll'),
#              type = 'mean',
#              n.iter = 1000,
#              n.burnin = 100,
#              nthin = 1)
