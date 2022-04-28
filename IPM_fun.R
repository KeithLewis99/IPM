
#' JAGS data format
#' 
#' @param log - put the capelin abundance on natural log (yes) or real scale (no)
#' @param forecast - add 2 forecast years - yes or no
#' 
#' @return a list with the the various input data sets; 
#' @export
#'
#' @examples see IPM_fun.R
#' num_forecasts = 2 
#' # 2 extra years jags.data <- ls_jag("yes", "yes")
#' jd <- as.data.frame(jags.data)
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
                        CO = df_con$meanCond
                )
        } else if (log == "no" & forecast == "no"){
                jags.data <- list(#year = df_dis_tab$year,
                        n.occasions = length(df_dis_tab$year),
                        I2 =  df_dis_tab$I2,
                        I3 =  df_dis_tab$I3,
                        m = df_mat$mat,
                        LD = df_ld$larvae,
                        TI = df_ice$tice,
                        CO = df_con$meanCond
                )
        } else if (log == "no" & forecast == "yes"){
                jags.data <- list(#year = df_dis_tab$year,
                        n.occasions = length(df_dis_tab$year) + num_forecasts,
                        I2 =  c(df_dis_tab$I2, rep(NA, num_forecasts)),
                        I3 =  c(df_dis_tab$I3, rep(NA, num_forecasts)),
                        m = c(df_mat$mat, rep(mean(df_mat$mat), num_forecasts)),
                        LD = c(df_ld$larvae, rep(mean(df_ld$larvae), num_forecasts)),
                        TI = c(df_ice$tice, rep(mean(df_ice$tice), num_forecasts)),
                        CO = c(df_con$meanCond, rep(mean(df_ice$tice), num_forecasts))
                )
        } else if(log == "yes" & forecast == "yes"){
                jags.data <- list(#year = df_dis_tab$year,
                        n.occasions = length(df_dis_tab$year) + num_forecasts,
                        I2 =  c(df_dis_tabLog$I2, rep(NA, num_forecasts)),
                        I3 =  c(df_dis_tabLog$I3, rep(NA, num_forecasts)),
                        m = c(df_mat$mat, rep(mean(df_mat$mat), num_forecasts)),
                        LD = as.vector(scale(c(df_ld$larvae, rep(NA, num_forecasts)),10)),
                        TI = as.vector(scale(c(df_ice$tice, rep(mean(df_ice$tice), num_forecasts)),10)),
                        CO = as.vector(scale(c(df_con$meanCond, rep(mean(df_con$meanCond), num_forecasts))))
                )
        }
        return(jags.data)
}


#' ls_out
#' Extract a list from the jags.data list with important variables - feeds function "ls_med" - really just simplifies by renaming
#' @param x - the list of jags output from function ls_jag
#'
#' @return a list
#' @export
#'
#' @examples raw <- ls_out(out) - see IPM_out.R
ls_out <- function(x){
     I2 = x$sims.list$I2
     I3 = x$sims.list$I3
     N3 = x$sims.list$N3
     N2 = x$sims.list$N2
     mu = x$sims.list$mu
     sigmaJ = x$sims.list$sigma
     ls_raw <- list(I2 = I2, I3=I3, N3=N3, N2=N2, mu=mu, sigmaJ=sigmaJ)
     return(ls_raw)
}


#' ls_med
#' The purpose here is to extract medians, credible intervals, and prediction intervals - feeds function "ipm_plot"
#' @param ls - the list of jags output from function ls_out
#'
#' @return A list 
#' @export
#'
#' @examples calc <- ls_med(raw) - see IPM_out.R
ls_med <- function(ls){
     # convert to simplify code
     I2 <- ls$I2
     I3 <- ls$I3
     N2 <- ls$N2
     N3 <- ls$N3
     mu <- ls$mu
     
     #summarize output as medians
     # out$median$N2 produces the same code as below
     I2_med = apply(I2,2,'median') # median values of y_pred
     I3_med = apply(I3,2,'median') # median values of y_pred
     N3_med = apply(N3,2,'median')
     N2_med = apply(N2,2,'median')
     mu_med = apply(mu,2,'median')
     #log(exp(calc$I3_med)+exp(calc$I2_med)) same code except for NA values
     I_med <- log(apply(exp(I3)+exp(I2), 2, 'median'))
     N_med <- log(apply(exp(N3)+exp(N2), 2, 'median'))
     
     # calculate credible and prediction intervals
     I_ci <- log(apply(exp(I3)+exp(I2), 2, 'quantile', c(0.025, 0.975)))
     I2_ci <- apply(I2, 2, 'quantile', c(0.025, 0.975))
     N_ci <- log(apply(exp(N3)+exp(N2), 2, 'quantile', c(0.025, 0.975)))
     Pr_ci <- log(apply(exp(N3)+exp(N2), 2, 'quantile', c(0.1, 0.9)))
     I = log(exp(jd$I2) + exp(jd$I3))
     ls_med <- list(I2_med = I2_med, I3_med = I3_med, N2_med = N2_med, N3_med = N3_med, I_med = I_med, N_med = N_med, mu_med = mu_med, I_ci = I_ci, N_ci = N_ci, Pr_ci = Pr_ci, I = I, I2_ci = I2_ci)
     return(ls_med)
}


#cbind(N2_med, N3_med, N2_med+N3_med)


#' IPM plot
#' Produce a plot of the jags output inlucding median values, credible/prediction intervals and the raw data
#' @param x the output from function ls_med
#'
#' @return
#' @export
#'
#' @examples tmp_plot <- ipm_plot(calc) - see IPM_out.R
#' requires specification of values for ly and lf 
#' year <- 1999:2021
#' ly <- length(year)
#' forecast <- 2022:2023
#' lf <- length(forecast)

ipm_plot <- function(df1 = x, df2 = y){
     p <- ggplot()  
     #plot credible interval and median for process
     p <- p + geom_ribbon(aes(x = year, 
                              ymax = df1$N_ci[2, 1:c(ly)], 
                              ymin = df1$N_ci[1, 1:c(ly)]),
                          alpha = 0.5, fill = "grey60")
     p <- p + geom_line(aes(x = year, y = df1$N_med[1:c(ly)]))
     
     #plot prediction interval and median for process
     p <- p + geom_ribbon(aes(x = forecast, 
                              ymax = tail(df1$Pr_ci[2,], lf), 
                              ymin = tail(df1$Pr_ci[1,], lf)),
                          alpha = 0.5, fill = "orange")
     p <- p + geom_line(aes(x = forecast, y = tail(df1$N_med, lf)))
     #browser()
     # plot points - observation median for I2 and I3
     p <- p + geom_point(aes(y = df1$I_med, x = c(year, forecast)),
                         shape = 16, 
                         size = 1.5,
                         colour = "red")
     
     # plot error bars from aggregated survey
     p <- p + geom_errorbar(data = df_cap[15:37,], aes(x = year, min=log(ab_lci*1000), ymax=log(ab_uci*1000)), width = 0.3, colour = "black")
     
     # plot points from disaggregated survey
     p <- p + geom_point(aes(y = df1$I, x = c(year, forecast)),
                         shape = 16, 
                         size = 1.5)
     # points from teh aggregated survey
     p <- p + geom_point(aes(y = log(df2$abundance_med*1000), x = c(year, forecast)),
                         shape = 18, 
                         size = 1.5,
                         colour = "pink")
     p <- p + theme_bw() + ylab("ln(Capelin abundance x 1,000)")
     return(p)
}



#' Title
#'
#' @param x 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
keep_rows <- function(x,y){
  #browser()
  tmp_vec <- rep(NA, y)
  for(i in 1:y){
    row <- paste0(x,"[",i, "]")
    tmp_vec[i] <- row
  }
  return(tmp_vec)
}




#' Title - DEPRECATED - WENT ABOUT THIS IN THE TOTALLY WRONG WAY.  MUCH EASIER WAY BELOW.  But, I learned some valuable lessons on this including how to work with [[]] and text.
#'
#' #' @param df1 
#' #' @param df2 
#' #' @param df3 
#' #' @param tabName 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' tabParm <- function(df1 = x, df2 = y, df3 = z, tabName){
#'   browser()
#'   ls_tab <- rep(list(list()), length(tabName))
#'   #if()
#'   for(i in 1:length(tabName)){
#'     med = df1[[tabName[i]]]
#'     LL = paste0(df2[[tabName[i]]][1,])
#'     UL = paste0(df2[[tabName[i]]][2,])
#'     Rhat = df3[1:25, "Rhat"]
#'     Neff = df3[1:25, "n.eff"]
#'     tab <- data.frame(cbind(year = 1999:2023, med = med, LL = LL, UL = UL, Rhat = Rhat, Neff = Neff))
#'     ls_tab[[i]] <- tab
#'   }
#'   return(ls_tab)
#' }



#' Title
#' Complete the skeleton.  Perhaps loop through tabNmes would help simplify code.  
#' @param df1 
#' @param tabName 
#'
#' @return
#' @export
#'
#' @examples
tabParm <- function(df1 = x, rows){
  #browser()
  tab <- data.frame(matrix(NA, nrow = 25, ncol = 5))
  #if()
  #for(i in 1:length(rows)){
    tab <- df1[rownames(df1) %in% df_keep_rows, c("mean", "2.5%", "97.5%", "Rhat", "n.eff")]
   # }
    tab1 <- data.frame(cbind(year = 1999:2023, tab))
  return(tab1)
}


# END ----
  
  
  