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
  ls_jag <- function(log, forecast, matrix = NULL){
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
        } else if(log == "yes" & forecast == "yes" & matrix == "no"){
                jags.data <- list(#year = df_dis_tab$year,
            n.occasions = length(df_dis_tab$year) + num_forecasts,
            I2 =  c(df_dis_tabLog$I2, rep(NA, num_forecasts)),
            I3 =  c(df_dis_tabLog$I3, rep(NA, num_forecasts)),
            I4 =  c(df_dis_tabLog$I4, rep(NA, num_forecasts)),
            m = c(df_mat$mat, rep(mean(df_mat$mat), num_forecasts)),
            LD = as.vector(scale(c(df_ld$larvae, rep(NA, num_forecasts-1)),10)),
            TI = as.vector(scale(c(df_ice$tice, rep(mean(df_ice$tice), num_forecasts)),10)),
            CO = as.vector(scale(c(df_con$meanCond, rep(mean(df_con$meanCond, na.rm = T), num_forecasts))))
                )
        } else if (log == "yes" & forecast == "yes" & matrix == "yes"){
          #browser()
          matI <- matrix(NA, nrow=39, ncol = 3)
          matI[,1] <- c(df_dis_tabLog$I2, rep(NA, num_forecasts))
          matI[,2] <- c(df_dis_tabLog$I3, rep(NA, num_forecasts))
          matI[,3] <- c(df_dis_tabLog$I4, rep(NA, num_forecasts))
          
          jags.data <- list(#year = df_dis_tab$year,
            n.occasions = length(df_dis_tab$year) + num_forecasts,
            matI = matI,
            m = c(df_mat$mat, rep(mean(df_mat$mat), num_forecasts)),
            LD = as.vector(scale(c(df_ld$larvae, rep(NA, num_forecasts-1)),10)), # the minus 1 is just a patch for now until I update the other data sets.
            TI = as.vector(scale(c(df_ice$tice, rep(mean(df_ice$tice), num_forecasts)),10)),
            CO = as.vector(scale(c(df_con$meanCond, rep(mean(df_con$meanCond, na.rm = T), num_forecasts))))
            )
            }
        return(jags.data)
}


#' ls_out
#' Extract a list from the jags.data list with important variables - feeds function "ls_med" - really just simplifies by renaming
#' confirmed that this renames arrays as well
#' @param x - the list of jags output from function ls_jag
#'
#' @return a list
#' @export
#'
#' @examples raw <- ls_out(out) - see IPM_out.R
ls_out <- function(x){
  #browser()
  #ls_names <- names(x$sims.list)
  ls_raw <- rep(list(list()), length(x$sims.list)) # make a list of appropriate length
    for(i in 1:length(x$sims.list)){
      df = x$sims.list[names(x$sims.list[i])] # creates a list called df
      ls_raw[i] <- df # put df in the proper slot of list ls_raw
      names(ls_raw)[i] <- names(x$sims.list[i]) # names the slot in the list
    }
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
  #browser()
  ls_med <- rep(list(list()), length(ls))
  
  # for(i in 1:length(ls)){
  #   med = apply(ls[[i]],2,'median')
  #   ls_med[[i]] <- med
  #   names(ls_med)[i] <- names(ls)[i]
  # }

  for(i in 1:length(ls)){
       if(is.matrix(ls[[i]])){ # is it a matrix?
            
            med = apply(ls[[i]],2,'median') # if yes, get the median
            ls_med[[i]] <- med # put the median values in the list slot
            names(ls_med)[i] <- names(ls)[i] # name the list slot
       } else {
            med = apply(ls[[i]],c(2,3),'median') # if its an array, get teh median of the column and the third dimension
            #med_array <- rep(list(list()), ncol(med))
            #for(a in ncol(med)){
            ls_med[[i]] <- med
            names(ls_med)[i] <- names(ls)[i]
       }
  }  
  
  #browser()
  ls_cri <- rep(list(list()), length(ls))

  for(i in 1:length(ls)){
       if(is.matrix(ls[[i]])){ # for matrices
            cri = apply(ls[[i]],2,'quantile', c(0.025, 0.975))
            ls_cri[[i]] <- cri
            names(ls_cri)[i] <- paste0(names(ls)[i], "_cri")
            
       } else { # for arrays
            cri = apply(ls[[i]], c(2,3),'quantile', c(0.025, 0.975))
            ls_cri[[i]] <- cri
            names(ls_cri)[i] <- paste0(names(ls)[i], "_cri")
       }
  }

  
  ls_pri <- rep(list(list()), length(ls))

  for(i in 1:length(ls)){
       if(is.matrix(ls[[i]])){
            pri = apply(ls[[i]],2,'quantile', c(0.1, 0.9))
            ls_pri[[i]] <- pri
            names(ls_pri)[i] <- paste0(names(ls)[i], "_pri")
       } else {
            pri = apply(ls[[i]], c(2,3),'quantile', c(0.1, 0.9))
            ls_pri[[i]] <- pri
            names(ls_pri)[i] <- paste0(names(ls)[i], "_pri")
       }
  }
  #browser()
  
  # calculate the median value for all ages by year
  N_med <- log(apply(exp(ls$N3)+exp(ls$N2) + exp(ls$N4), 2, 'median'))
  # calculate the credible interval for all ages by year
  N_ci <- apply(log(exp(ls$N3)+exp(ls$N2) + exp(ls$N4)), 2, 'quantile', c(0.025, 0.975))
  # calculate the prediction interval based on simulated data using the approach in Schuab and Kerry
  Pr_ci <- log(apply(exp(ls$I3.rep)+exp(ls$I2.rep) + exp(ls$I4.rep), 2, 'quantile', c(0.1, 0.9)))
  
  ls_calc <- list(ls_med = ls_med, ls_cri = ls_cri, ls_pri = ls_pri, N_med = N_med, N_ci = N_ci, Pr_ci = Pr_ci)
     return(ls_calc)
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

ipm_plot <- function(df_med = x, df_cri = y, df_pri = z, df_dat = w){
     p <- ggplot()  
     #plot credible interval and median for process
     #browser()
     p <- p + geom_ribbon(aes(x = year, 
                              ymax = df_cri[2, 1:c(ly)], 
                              ymin = df_cri[1, 1:c(ly)]),
                          alpha = 0.5, fill = "grey60")
     p <- p + geom_line(aes(x = year, y = df_med[1:c(ly)]))
     
     #plot prediction interval and median for process
     p <- p + geom_ribbon(aes(x = forecast, 
                              ymax = tail(df_pri[2,], lf), 
                              ymin = tail(df_pri[1,], lf)),
                          alpha = 0.5, fill = "orange")
     p <- p + geom_line(aes(x = forecast, y = tail(df_med, lf)))
     #browser()
     # plot points - observation median for I2 and I3
     # p <- p + geom_point(aes(y = df1$I_med, x = c(year, forecast)),
     #                     shape = 16, 
     #                     size = 1.5,
     #                     colour = "red")
     # 
     # plot error bars from aggregated survey
     p <- p + geom_errorbar(data = df_dat, aes(x = year, min=log(ab_lci*1000), ymax=log(ab_uci*1000)), width = 0.3, colour = "black")
     
     # plot points from disaggregated survey
     # p <- p + geom_point(aes(y = df1$I, x = c(year, forecast)),
     #                     shape = 16, 
     #                     size = 1.5)
     # # points from teh aggregated survey
     p <- p + geom_point(aes(y = log(df_dat$abundance_med*1000), x = c(year, forecast)),
                         shape = 18, 
                         size = 2,
                         colour = "red")
     p <- p + ylab("ln(Capelin abundance x 1,000)") + xlab("Year")
     p <- p + theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25)) + theme_bw() 
     return(p)
}

pe_plot <- function(df_med = x, df_cri = y, df_pri = z){
     #' PE plot
     #' Produce a plot of the process error output inlucding median values, credible/prediction intervals and the raw data
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
     
          p <- ggplot()  
          #plot credible interval and median for process
          #browser()
          p <- p + geom_ribbon(aes(x = year, 
                                   ymax = df_cri[2, 1:c(ly)], 
                                   ymin = df_cri[1, 1:c(ly)]),
                               alpha = 0.5, fill = "grey60")
          p <- p + geom_line(aes(x = year, y = df_med[1:c(ly)]))
          
          #plot prediction interval and median for process
          p <- p + geom_ribbon(aes(x = forecast, 
                                   ymax = tail(df_pri[2,], lf), 
                                   ymin = tail(df_pri[1,], lf)),
                               alpha = 0.5, fill = "orange")
          p <- p + geom_line(aes(x = forecast, y = tail(df_med, lf)))
          #browser()
          # plot points - observation median for I2 and I3
          p <- p + ylab("Process error") + xlab("Year")
          p <- p + theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25)) + theme_bw() 
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
keep_rows <- function(x,y, yr = NULL){
  #browser()
  tmp_vec <- rep(NA, y)
if(matrix == "yes"){
     for(i in 1:y){
          row <- paste0(x, "[", i, ",", yr, "]")
          tmp_vec[i] <- row
     }
     } else {
     for(i in 1:y){
          row <- paste0(x,"[",i, "]")
          tmp_vec[i] <- row
          }
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
     if(disaggregated == "1985-present"){
          nrow = length(1985:2023)
     } else {
          nrow = length(1999: 2023)
     }
  tab <- data.frame(matrix(NA, nrow = nrow, ncol = 5))
  #if()
  #for(i in 1:length(rows)){
    tab <- df1[rownames(df1) %in% df_keep_rows, c("mean", "2.5%", "97.5%", "Rhat", "n.eff")]
   # }
    if(disaggregated == "1985-present"){
         tab1 <- data.frame(cbind(year = 1985:2023, tab))
    } else {
         tab1 <- data.frame(cbind(year = 1999:2023, tab))         
    }

  return(tab1)
}


#' postPriors()----
#'
#' @param df - the chain of values of the parameter
#' @param df2 - the distribution of the prior based on randomly generated data
#' @param df3 - the chain values within the credible intervals of the parameter
#' @param limits - an ojbect to set the limits on the grid 
#' @param x_label - label the graph with the name of the parm
#' @param priormean -  the mean value of the prior
#' @param priorsd - the SD of the prior
#' @param by_bin - 
#'
#' @return - figure of the posterior distribution, the prior distribution and the 95% credible interval as a rug.  Abline = 0
#' @export
#'
#' @example p1 <- postPriors(df = alpha$df, df2 = prior, df3 = alpha$df_cred, limits, x_label, priormean, priorsd, by_bin = bin_1)
postPriors <- function(df, df2, df3, limits = limits, x_label = x_label, priormean, priorsd, by_bin = 1, vline="yes"){
  #browser()
  p <- ggplot() + theme_bw(base_size = 20)
  p <- p + coord_cartesian(xlim = c(limits[1], limits[2])) 
  p <- p + geom_histogram(data = as.data.frame(df), 
                          #aes(x = V1, y=..density..),
                          aes(x = df, y=..density..), 
                          col="grey", fill="grey",
                          #alpha=.6,
                          breaks = seq(limits[1], limits[2], by=by_bin)) + xlab(x_label) + ylab("Density")
  p <- p + geom_density(data = as.data.frame(df2), aes(x=df2), size = 1, alpha = 0.3) 
  p <- p + geom_rug(data=data.frame(y=df3), aes(x=y), colour = "red", size=3, alpha=1) 
  if(vline == "yes"){
    p <- p + geom_vline(aes(xintercept = 0), colour = "red", size = 2)
  #  p <- p + annotate(geom = "text", x = 4.6, y = 0.5, label = paste("=", round(pB, 2)), colour = "black")
  }
  return(p)
}


#' post_param
#' The purpose is to bundle up all the variables that go into the postPriors figures.
#' @return
#' @export
#'
#' @examples
post_param <- function(param = w, priormean = NULL, priorsd = NULL, jags = z, x_label = NULL){
  #browser()
if(param == "alpha"|param == "beta"|param=="alpha[1]"|param=="beta[1]"){
  param = param
  priormean <- priormean
  priorsd <- priorsd
  jags <- jags
  prior <- rnorm(n = 10000, mean = priormean, sd = priorsd)
  limits <- c(min(jags)-0.3, max(jags) + 0.3)
  x_label <- x_label
  bin_1 <- abs(mean(jags)/100)
  df_quant <- quantile(jags, c(0.025, 0.975))
  df_cred <- subset(jags, jags > df_quant[1] & jags < df_quant[2])

} else if (param == "gamma"|param == "gamma[1]"){
  #browser()
  param = param
  jags <- jags
  prior <- runif(n = 10000, min = 0, max = 100)
  limits <- c(0, max(jags) + 0.1)
  x_label <- expression(paste(t[italic(ice)], "-MRI-gamma"))
  bin_1 <- mean(jags)/100
  df_quant <- quantile(jags, c(0.025, 0.975))
  df_cred <- subset(jags, jags > df_quant[1] & jags < df_quant[2])

} else if (param == "delta"|param == "delta[1]"){
  #browser()
  param = param
  jags <- jags
  prior <- rgamma(n = 10000, shape = priormean, rate = priorsd)*100
  limits <- c(min(jags)-0.3, max(jags) + 0.3)
  x_label <- expression(paste(t[italic(ice)], "width - delta"))
  bin_1 <- mean(jags)/100
  df_quant <- quantile(jags, c(0.025, 0.975))
  df_cred <- subset(jags, jags > df_quant[1] & jags < df_quant[2])
} 
  list(param=param,priormean = priormean, priorsd=priorsd, jags=jags, prior=prior, limits=limits, x_label=x_label, bin_1=bin_1, df_quant = df_quant, df_cred=df_cred)
}

#' rsq_bayes()----
#'
#' @param ypred - vector values of y from Bayesian chain 
#' @param out - list of output values from Bayesian analysis
#'
#' @return - Bayesian R-squared
#' @export
#'
#' @examples R1_r2 <- rsq_bayes(ypred = y_pred, out = run_recruit)
# see jags_bayes_R2 for proof that this matches the Gelmen code
rsq_bayes <- function(ypred = ypred, out=out){
  #browser()
  # variance of predicted values
  y_var = apply(ypred, 1, 'var')
  
  # variance of residuals
  res_pred = out$BUGSoutput$sims.list$Res
  #res_med = apply(res_pred,2,'median')
  res_var = apply(res_pred,1,'var')
  
  # Distribution of Bayesian R-squared
  bay_rsq <- y_var/(res_var + y_var)
  return(bay_rsq)
}


# END ----
  
  
  