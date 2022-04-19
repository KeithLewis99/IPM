
#' JAGS data format
#'
#' @param log 
#' @param forecast 
#'
#' @return
#' @export
#'
#' @examples
ls_jag <- function(log, forecast){
        #browser()
        if(log == "yes" & forecast == "no"){
                jags.data <- list(#year = df_dis_tab$year,
                        #n.occasions = length(df_dis_tab$year) + num_forecasts,
                        n.occasions = length(df_dis_tab$year),
                        I2 =  df_dis_tabLog$I2,
                        I3 =  df_dis_tabLog$I3,
                        m = df_mat$mat,
                        LD = df_ld$lnlarvae,
                        TI = df_ice$tice,
                        CO = df_ice$tice
                )
        } else if (log == "no" & forecast == "no"){
                jags.data <- list(#year = df_dis_tab$year,
                        n.occasions = length(df_dis_tab$year),
                        I2 =  df_dis_tab$I2,
                        I3 =  df_dis_tab$I3,
                        m = df_mat$mat,
                        LD = df_ld$larvae,
                        TI = df_ice$tice,
                        CO = df_ice$tice
                )
        } else if (log == "no" & forecast == "yes"){
                jags.data <- list(#year = df_dis_tab$year,
                        n.occasions = length(df_dis_tab$year) + num_forecasts,
                        I2 =  c(df_dis_tab$I2, rep(NA, num_forecasts)),
                        I3 =  c(df_dis_tab$I3, rep(NA, num_forecasts)),
                        m = c(df_mat$mat, rep(mean(df_mat$mat), num_forecasts)),
                        LD = c(df_ld$larvae, rep(mean(df_ld$larvae), num_forecasts)),
                        TI = c(df_ice$tice, rep(mean(df_ice$tice), num_forecasts)),
                        CO = c(df_ice$tice, rep(mean(df_ice$tice), num_forecasts))
                )
        } else if(log == "yes" & forecast == "yes"){
                jags.data <- list(#year = df_dis_tab$year,
                        n.occasions = length(df_dis_tab$year) + num_forecasts,
                        I2 =  c(df_dis_tabLog$I2, rep(NA, num_forecasts)),
                        I3 =  c(df_dis_tabLog$I3, rep(NA, num_forecasts)),
                        m = c(df_mat$mat, rep(mean(df_mat$mat), num_forecasts)),
                        LD = c(df_ld$larvae, rep(mean(df_ld$larvae), num_forecasts)),
                        TI = c(df_ice$tice, rep(mean(df_ice$tice), num_forecasts)),
                        CO = c(df_ice$tice, rep(mean(df_ice$tice), num_forecasts))
                )
        }
        return(jags.data)
}


#' ls_out
#' The purpose here is to extract the output from the jags.data
#' @param out 
#'
#' @return
#' @export
#'
#' @examples
ls_out <- function(out){
     
     I2 = out$sims.list$I2
     I3 = out$sims.list$I3
     N3 = out$sims.list$N3
     N2 = out$sims.list$N2
     ls_raw <- list(I2 = I2, I3=I3, N3=N3, N2=N2)
     return(ls_raw)
}

#' ls_med
#' The purpose here is to extract medians, credible intervals, and prediction intervals
#' @return
#' @export
#'
#' @examples
ls_med <- function(ls){
     # convert to simplify code
     I2 <- ls$I2
     I3 <- ls$I3
     N2 <- ls$N2
     N3 <- ls$N3
     
     #summarize output as medians
     I2_med = apply(I2,2,'median') # median values of y_pred
     I3_med = apply(I3,2,'median') # median values of y_pred
     N3_med = apply(N3,2,'median')
     N2_med = apply(N2,2,'median')
     I_med <- log(apply(exp(I3)+exp(I2), 2, 'median'))
     N_med <- log(apply(exp(N3)+exp(N2), 2, 'median'))
     
     # calculate credible and prediction intervals
     I_ci <- log(apply(exp(I3)+exp(I2), 2, 'quantile', c(0.025, 0.975)))
     N_ci <- log(apply(exp(N3)+exp(N2), 2, 'quantile', c(0.025, 0.975)))
     Pr_ci <- log(apply(exp(N3)+exp(N2), 2, 'quantile', c(0.1, 0.9)))
     I = log(exp(jd$I2) + exp(jd$I3))
     ls_med <- list(I2_med = I2_med, I3_med = I3_med, N2_med = N2_med, N3_med = N3_med, I_med = I_med, N_med = N_med, I_ci = I_ci, N_ci = N_ci, Pr_ci = Pr_ci, I = I)
     return(ls_med)
}


#cbind(N2_med, N3_med, N2_med+N3_med)





#' IPM plot
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
ipm_plot <- function(x){
     p <- ggplot()  
     #plot credible interval and median
     p <- p + geom_ribbon(aes(x = year, 
                              ymax = x$N_ci[2, 1:c(ly)], 
                              ymin = x$N_ci[1, 1:c(ly)]),
                          alpha = 0.5, fill = "grey60")
     p <- p + geom_line(aes(x = year, y = x$N_med[1:c(ly)]))
     
     #plot prediction interval and median
     p <- p + geom_ribbon(aes(x = forecast, 
                              ymax = tail(x$Pr_ci[2,], lf), 
                              ymin = tail(x$Pr_ci[1,], lf)),
                          alpha = 0.5, fill = "orange")
     p <- p + geom_line(aes(x = forecast, y = tail(x$N_med, lf)))
     
     # plot points - observation abundance for I2 and I3
     p <- p + geom_point(aes(y = x$I_med, x = c(year, forecast)),
                         shape = 16, 
                         size = 1.5,
                         colour = "red")
     
     # plot error bars from aggregated survey
     p <- p + geom_errorbar(data = df_cap[15:37,], aes(x = year, min=log(ab_lci*1000), ymax=log(ab_uci*1000)), width = 0.3, colour = "black")
     
     # plot points from disaggregated survey
     p <- p + geom_point(aes(y = x$I, x = c(year, forecast)),
                         shape = 16, 
                         size = 1.5)
     p <- p + theme_bw() + ylab("ln(Capelin abundance x 1,000)")
     return(p)
}
