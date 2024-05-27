
IPMplot <- function(raw, inputI){
   
   Ntemp <- raw$Np |>
      apply(c(1,2), sum) |>
      as.data.frame() |>
      pivot_longer(V1:V40, names_to = 'Year', values_to = 'N') |>
      mutate(Year = as.numeric(gsub('^V', '', Year))+1984) |>
      reframe(median = median(N), 
              mean = mean(N),
              p05 = quantile(N, .05),
              p10 = quantile(N, .10),
              p90 = quantile(N, .90),
              p95 = quantile(N, .95),
              .by = 'Year') 
   Itemp <- raw$I3.rep |>
      as.data.frame() |>
      pivot_longer(V1:V40, names_to = 'Year', values_to = 'I') |>
      mutate(Year = as.numeric(gsub('^V', '', Year))+1984) |>
      reframe(median = median(I), 
              mean = mean(I),
              p05 = quantile(I, .05),
              p10 = quantile(I, .10),
              p90 = quantile(I, .90),
              p95 = quantile(I, .95),
              .by = 'Year')
   
   ggplot()+
      geom_line(data = Ntemp, aes(x = Year, y = median), col = 'steelblue')+
      geom_ribbon(data = Ntemp, aes(x = Year, ymin = p05, ymax = p95), 
                  fill = 'steelblue', alpha = .4)+
      geom_line(data = Itemp, aes(x = Year, y = median), col = 'coral')+
      geom_ribbon(data = Itemp, aes(x = Year, ymin = p10, ymax = p90),
                  fill = 'coral', alpha = .2)+
      geom_point(data = inputI |>
                    as.data.frame() |>
                    mutate(Year = 1985:2024,
                           I = log(exp(I2)+exp(I3)+exp(I4))),
                 aes(x = Year, y = I))+
      geom_line(data = data.frame(x=2000,y=1, col = c('coral', 'steelblue', 'black')), 
                aes(x = x, y = y, col = col),lty = NA)+ #This is a dummy df used for legend
      geom_point(data = data.frame(x=2000,y=1, col = c('coral', 'steelblue', 'black')), 
                 aes(x = x, y = y, col = col), pch = NA)+ #This is a dummy df used for legend
      theme_minimal()+
      ylab('Abundance (Index)')+
      scale_color_manual(values = c('steelblue', 'coral','black'),
                         labels = c('N', 'Irep', 'Iobs'),
                         name = "")+
      guides(color = guide_legend(override.aes = list(pch = c(NA, NA, 20), 
                                                      size = c(NA, NA, 3),
                                                      lty = c(1, 1, NA)))
      )+
      theme(legend.position = 'right')
   
}



osaResidPlot <- function(data, .by = 'Year', byage = T, ylim = c(-10, 10)){
   
   df1 <- data[,,1] |> 
      as.data.frame() |>
      rename_with(\(x)gsub('V', '', x), everything()) |>
      mutate(Age = '1') |>
      pivot_longer('1':'39', names_to = 'Year', values_to = 'value')
   df2 <- data[,,1] |>
      as.data.frame() |>
      rename_with(\(x)gsub('V', '', x), everything()) |>
      mutate(Age = '2') |>
      pivot_longer('1':'39', names_to = 'Year', values_to = 'value')
   df3 <- data[,,1] |> 
      as.data.frame() |>
      rename_with(\(x)gsub('V', '', x), everything()) |>
      mutate(Age = '3') |>
      pivot_longer('1':'39', names_to = 'Year', values_to = 'value')
   data <- rbind(df1, df2, df3) |>
      mutate(Year = as.numeric(Year)+1984,
             Cohort = Year-as.numeric(Age),
             Age = gsub('^', 'Age ', Age))
   # return(data)
   
   fig <- data |>
      group_by(all_of(.by), Age) |>
      mutate(med = median(value)) |>
      ungroup() |>
      ggplot(aes(x = get(.by), y = value)) +
      geom_point() + 
      geom_hline(yintercept = 0, linetype = 2) +
      geom_line(aes(y = med), col = 'coral') + 
      ylim(ylim[1], ylim[2])+
      xlab(.by)+ylab("Res")
   if(byage){
      fig <- fig +
         facet_wrap(~Age, nrow = 3)
   }
   
   return(fig)
   
}


# index_plot <- function(data, .byage = T){
#    
#    data <- data$med[c('N2', 'N3', 'N4')] |>
#       as.data.frame() |>
#       mutate(Year = years) |>
#       rename('Age 2' = N2, 'Age 3' = N3, 'Age 4' = N4) |>
#       pivot_longer(cols = 'Age 2':'Age 4', names_to = 'Age', values_to = 'Value') |>
#       mutate(Parm = 'N') |>
#       rbind(
#          data$med[c('I2.rep', 'I3.rep', 'I4.rep')] |>
#             as.data.frame() |>
#             mutate(Year = years) |>
#             rename('Age 2' = I2.rep, 'Age 3' = I3.rep, 'Age 4' = I4.rep) |>
#             pivot_longer(cols = 'Age 2':'Age 4', names_to = 'Age', values_to = 'Value') |>
#             mutate(Parm = 'I')
#       )
#       
#    if(!.byage){
#       data <- data |>
#          reframe(Value = sum(Value),
#                  .by = c("Year", 'Parm')) |>
#          pivot_wider(id_cols = c('Year'), 
#                      values_from = Value, names_from = Parm) 
#    }else{
#       data <- data |>
#          pivot_wider(id_cols = c('Age', 'Year'), 
#                      values_from = Value, names_from = Parm) 
#    }
#    
#    fig <- data |>
#       ggplot(aes(x = Year)) +
#       geom_line(aes(y = N))+
#       geom_point(aes(y = I))+
#       theme_bw()
#    
#    if(.byage){
#       return(fig + facet_wrap(~Age, nrow = 3))
#    }
#    fig
#    
# }
# index_plot(data, T)

maturity_plot <- function(data){
   
   data.frame(M1 = data$med$m[,1],
              M1.lo = data$cri$m_cri[1,,1],
              M1.hi = data$cri$m_cri[2,,1],
              M2 = data$med$m[,2],
              M2.lo = data$cri$m_cri[1,,2],
              M2.hi = data$cri$m_cri[2,,2],
              M3 = data$med$m[,3],
              M3.lo = data$cri$m_cri[1,,3],
              M3.hi = data$cri$m_cri[2,,3],
              Year = 1985:2024) |>
      pivot_longer(M1:M3.hi, values_to = 'Val', names_to = 'Parm') |>
      mutate(Age = as.factor(substr(Parm, 2, 2)),
             Stat = sapply(Parm, \(x){
                y <- substr(x, 4, 5)
                y <- ifelse(y == '', 'med', y)
                y
             })
      ) |>
      select(-Parm) |>
      pivot_wider(id_cols = c(Year, Age), values_from = Val, names_from = Stat) |>
      ggplot(aes(x = Year)) +
      geom_line(aes(y = med, col = Age))+
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = Age), alpha = .2)+
      theme_bw()+
      ylab("Proportion Mature")
   
}

recruitment_plot <- function(data){
   
   data.frame(Rec = data$med$rec,
              Rec.lo = data$cri$rec_cri[1,],
              Rec.hi = data$cri$rec_cri[2,],
              Recp = data$med$recp,
              Recp.lo = data$cri$recp_cri[1,],
              Recp.hi = data$cri$recp_cri[2,],
              Year = 1985:2023) |>
      ggplot(aes(x = Year)) + 
      geom_line(aes(y = Rec, col = 'black'))+
      geom_ribbon(aes(ymin = Rec.lo, ymax = Rec.hi), fill = 'grey', alpha = .2)+
      geom_line(aes(y = Recp, col = 'red'))+
      geom_ribbon(aes(ymin = Recp.lo, ymax = Recp.hi), fill = 'pink', alpha = .2)+
      theme_bw()+
      scale_color_manual(name = 'With process error', labels = c('Yes', 'No'), 
                         values = c('red', 'black'))+
      theme(legend.position = 'top', legend.direction = 'horizontal')
   
}

spay_plot <- function(data){
   
   c(data$med$rec, NA) |>
      cbind(data$med$N) |>
      as.data.frame() |>
      mutate(Year = 1985:2024) |>
      rename('Age 1' = V1, 'Age 2' = V2, 'Age 3' = V3, 'Age 4' = V4) |>
      pivot_longer(cols = 'Age 1':'Age 4', names_to = 'Age', values_to = 'N') |>
      group_by(Age) |>
      mutate(N = (N-mean(N, na.rm = T))/sd(N, na.rm = T)) |>
      ggplot(aes(x = Year, y = Age, size = N, color = ifelse(N > 0, 'steelblue', 'coral')),
             alpha = .4) +
      geom_abline(slope = 1, intercept = -(1982:2025) , linetype = 2, col = 'gray') +
      geom_point() +
      guides(color = 'none', size = 'none')+
      theme_bw()
   
}
# spay_plot(data)



plot_condition <- function(data, input){
   
   input$CO <- input$CO[11:37]
   plot(1995:2021, input$CO, ylab = "Condition", xlab = "Year", 
        type = "b", xlim = c(1985, 2024),
        ylim = c(-8, 8))
   points(c(1985:1995, 2022:2024), data$med$eCO, col = "navy")
   arrows(x0 = 1:11+1984, x1 = 1:11+1984, 
          y0 = data$pri$eCO_pri[1,1:11], 
          y1 = data$pri$eCO_pri[2,1:11],
          length=0.05, angle=90, code=3, lwd = 2, col = "blue")
   arrows(x0 = 2022:2024, x1 = 2022:2024, 
          y0 = data$pri$eCO_pri[1,12:14], 
          y1 = data$pri$eCO_pri[2,12:14],
          length=0.05, angle=90, code=3, lwd = 2, col = "blue")
   
}


plotLaA <- function(data, input){
   LaA <- function(a){
      Linf <- data$med$Linf
      K <- data$med$K
      return(Linf * (1-exp(-K * a)))
   }
   A <- input$lAmax
   
   plot(input$lengths ~ input$lages,
        xlab = "Age", ylab = "Length", 
        xlim = c(0, A), ylim = c(0, 200))
   curve(LaA, 0, A, add = T)
   abline(h = data$med$Linf, lty = 2)
   arrows(x0 = 1:A, x1 = 1:A, 
          y0 = data$med$LaA+1.96*data$med$sigma.laa,# y0 = data$pri$LaA_pri[1,], 
          y1 = data$med$LaA-1.96*data$med$sigma.laa, # y1 = data$pri$LaA_pri[2,],
          length=0.1, angle=90, code=3, lwd = 2, col = "blue")
   text(3, 10, paste0("RMSE = ", data$med$RMSE %>% round(1)))

}

# plot(jags.data.m$lengths ~ jags.data.m$lages,
#      xlab = "Age", ylab = "Length", 
#      xlim = c(0, 7), ylim = c(0, 200))
# curve(ls_all$ls_med$Linf * (1-exp(-ls_all$ls_med$K * x)), 0, 6, add = T)
# abline(h = ls_all$ls_med$Linf, lty = 2)
# arrows(x0 = 1:6, x1 = 1:6, 
#        y0 = ls_all$ls_med$LaA+1.96*ls_all$ls_med$sigma.laa,# y0 = ls_all$ls_pri$LaA_pri[1,], 
#        y1 = ls_all$ls_med$LaA-1.96*ls_all$ls_med$sigma.laa, # y1 = ls_all$ls_pri$LaA_pri[2,],
#        length=0.1, angle=90, code=3, lwd = 2, col = "blue")
# text(x=6, y=0, labels = paste0("RMSE = ", ls_all$ls_med$RMSE %>% round(1)))
