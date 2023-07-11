## Post processing functions for frequentist simulations

## Set up color pallette using color blind friendly colors ##
## Resource: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with grey:
cbPalette <- c(  "#0072B2","#E69F00", "#D55E00", "#CC79A7","#999999", "#E69F00","#F0E442", "#56B4E9", "#009E73",  "#CC79A7")

# To use for fills, add scale_fill_manual(values=cbPalette)

# To use for line and point colors, add scale_colour_manual(values=cbPalette)


## Get percent bias of counterfactuals
# Uses all time points
get_perc_bias_Y0_gsc <- function(pop_trt,fit_gsc,Y0_matrix,ind){
  raw_bias_list <- list()
  perc_bias_list <- list()
  avg_perc_bias <- list()
  mean_perc_bias <- list()
  for(i in 1:length(fit_gsc)) {
    true_Y0_rate <- 100000 * (Y0_matrix[,i]/pop_trt)
    raw_bias_list[[i]] <- (as.vector(fit_gsc[[i]]$Y.ct) - (100000 * (Y0_matrix[,i])/pop_trt)) #fit_gsc returns rate
    perc_bias_list[[i]] <- 100*(raw_bias_list[[i]])/true_Y0_rate #(median trt effect - true trt effect)/true trt
    avg_perc_bias[i] <- median(abs(perc_bias_list[[i]][ind]))
    mean_perc_bias[i] <- mean(abs(perc_bias_list[[i]][ind]))
  }
  
  return(list(avg_perc_bias=avg_perc_bias,mean_perc_bias = mean_perc_bias))
}
#---------------------------------------------------------------------------------
## Abs Percent Bias Y0 Function GSC
medianAbsPercBiasY0_year_gsc <- function(pop_trt,fit_gsc,Y0_matrix,ind){
  medianAbsPercBias <- c()
  for(i in 1:length(fit_gsc)){
    # Find true treatment effect rate
    true_Y0_rate <- 100000 * (Y0_matrix[ind,i]/pop_trt[ind])# Find Rate
    
    # Take average of counties for each time point
    true_Y0_mean_rate <- colMeans(matrix(true_Y0_rate , ncol = 7, byrow = T)) # take median across counties
    
    # extract the treated indices
    # Convert to rate
    mean_vec_rate <- as.vector(fit_gsc[[i]]$Y.ct)[ind]
    
    # For each simulation, take average counterfactual for each year over all counties 
    # should be a vector of length 7
    est_Y0 <- colMeans(matrix(mean_vec_rate, ncol = 7, byrow = T))
    
    # find absolute percent bias Y0 
    absPercBias <- abs(((est_Y0  - true_Y0_mean_rate)/true_Y0_mean_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  
  # median of the simulations
  median <- median(medianAbsPercBias)
  LB <- quantile(medianAbsPercBias, probs = 0.05, na.rm = T)
  UB <- quantile(medianAbsPercBias, probs = 0.95, na.rm = T)
  
  return(list(median = median, LB = LB, UB = UB))
}

## Abs Percent Bias Y0 ASC
medianAbsPercBiasY0_year_asc <- function(pop_trt,fit_asc,Y0_matrix,ind){
  medianAbsPercBias <- c()
  for(i in 1:length(fit_asc)){
    # Find true treatment effect rate
    true_Y0_rate <- 100000 * (Y0_matrix[ind,i]/pop_trt[ind])# Find Rate
    
    # Take average of counties for each time point
    true_Y0_mean_rate <- colMeans(matrix(true_Y0_rate , ncol = 7, byrow = T)) # take median across counties
    
    # extract the treated indices
    # Convert to rate
    mean_vec_rate <- predict(fit_asc[[i]], att = F)[9:15]
    
    # For each simulation, take average counterfactual for each year over all counties 
    # should be a vector of length 7
    est_Y0 <- colMeans(matrix(mean_vec_rate, ncol = 7, byrow = T))
    
    # find absolute percent bias Y0 
    absPercBias <- abs(((est_Y0  - true_Y0_mean_rate)/true_Y0_mean_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  
  # median of the simulations
  median <- median(medianAbsPercBias)
  LB <- quantile(medianAbsPercBias, probs = 0.05, na.rm = T)
  UB <- quantile(medianAbsPercBias, probs = 0.95, na.rm = T)
  
  return(list(median = median, LB = LB, UB = UB))
}

#------------------------------------------------------------------------------------

## Function for GSC
medianAbsPercBiasATT_year_gsc <- function(fit_gsc,Y0_matrix,Y1_matrix,true_trt_effect){
  medianAbsPercBias <- c()
  for(i in 1:length(fit_gsc)){
    # Find true treatment effect rate
    true_trt_effect_rate <- (true_trt_effect/pop_trt[ind] * 100000)[,i]# Find Rate
    # Take average of counties for each time point
    true_trt_effect_median_rate <- colMeans(matrix(true_trt_effect_rate , ncol = 7, byrow = T)) # take median across counties
    
    # For each simulation, get median of 1000 MCMC samples. Then only extract the treated indices.
    # Convert to rate
    Y0_vec_rate <- as.vector(fit_gsc[[i]]$Y.ct)[ind]
    
    # For each simulation, convert observed Y1 to rate
    # should be a vector of length 6
    Y1_matrix_rate <- (Y1_matrix[,i]/pop_trt[ind]) * 100000
    
    # Find ATT for each time point
    TE <- Y1_matrix_rate - Y0_vec_rate
    ATT_vec <- colMeans(matrix(TE, ncol = 7, byrow = T)) # one ATT for each time
    
    # find absolute percent bias ATT 
    absPercBias <- abs(((ATT_vec  - true_trt_effect_median_rate)/true_trt_effect_median_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  # median of the simulations
  median <- median(medianAbsPercBias)
  LB <- quantile(medianAbsPercBias, probs = 0.05, na.rm = T)
  UB <- quantile(medianAbsPercBias, probs = 0.95, na.rm = T)
  
  return(list(median = median, LB = LB, UB = UB))
}

## Function for ASC
medianAbsPercBiasATT_year_asc <- function(fit_asc,Y0_matrix,Y1_matrix,true_trt_effect){
  medianAbsPercBias <- c()
  for(i in 1:length(fit_gsc)){
    # Find true treatment effect rate
    true_trt_effect_rate <- (true_trt_effect/pop_trt[ind] * 100000)[,i]# Find Rate
    # Take average of counties for each time point
    true_trt_effect_median_rate <- colMeans(matrix(true_trt_effect_rate , ncol = 7, byrow = T)) # take median across counties
    
    # Get median estimated treatment effect rate
    est_trt_effect <- predict(fit_asc[[i]], att = T)[9:15]
    
    # find absolute percent bias ATT 
    absPercBias <- abs(((est_trt_effect  - true_trt_effect_median_rate)/true_trt_effect_median_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }

  # median of the simulations
  median <- median(medianAbsPercBias)
  LB <- quantile(medianAbsPercBias, probs = 0.05, na.rm = T)
  UB <- quantile(medianAbsPercBias, probs = 0.95, na.rm = T)
  
  return(list(median = median, LB = LB, UB = UB))
}

#----------------------------------------------------------------------------------------------------
medianAbsPercBias_time_plot <- function(Y0_matrix_1, fit_gsc_mc_1,
                                        Y0_matrix_3, fit_gsc_mc_3,
                                        Y0_matrix_7, fit_gsc_mc_7,
                                        fit_gsc, fit_asc,
                                        pop_trt, ylim, ind){
  # Matrix completion in gsynth
  gsc_mc_1 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_1,Y0_matrix_1,ind)
  gsc_mc_3 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_3,Y0_matrix_3,ind)
  gsc_mc_7 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_7,Y0_matrix_7,ind)
  
  # GSC in gsynth
  gsc <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc,Y0_matrix_1,ind)
  
  # ASC in augsynth
  asc <- medianAbsPercBiasY0_year_asc(pop_trt,fit_asc,Y0_matrix_1,ind)
  
  k <- c("k=1","k=3","k=7", "NA", "NA")
  Model <- c(rep("MC in gsynth", 3), "GSC", "ASC")
  Median <- c(gsc_mc_1, gsc_mc_3, gsc_mc_7, gsc, asc)
  df <- data.frame(Model,Median,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1", "NA"),
           Model = fct_relevel(Model, "MC in gsynth", "GSC", "ASC"))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Median, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7', 'NA')) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(Median,2)),position = position_dodge(1), hjust = -.2, size=2.8)+
    #ggtitle("Absolute Percent Bias of the Counterfactuals (Non-Rare)") +
    coord_flip(ylim = ylim ) + 
    scale_x_discrete(limits = rev) +
    #scale_y_continuous(name="Median", limits=c(0, 90000)) +
    theme_bw() + #adds plot outline vs. theme_minimal
    theme(axis.text=element_text(color="black"),
          legend.title = element_blank()
    ) 
}

#----------------------------------------------------------------------------------------------------
medianAbsPercBiasATT_time_plot <- function(Y0_matrix_1, Y1_matrix_1, true_trt_effect_1, fit_gsc_mc_1,
                                        Y0_matrix_3, Y1_matrix_3, true_trt_effect_3, fit_gsc_mc_3,
                                        Y0_matrix_7, Y1_matrix_7, true_trt_effect_7, fit_gsc_mc_7,
                                        fit_gsc, fit_asc,
                                        pop_trt, ylim, ind){
  # Matrix completion in gsynth
  gsc_mc_1 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  gsc_mc_3 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  gsc_mc_7 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  
  # GSC in gsynth
  gsc <-medianAbsPercBiasATT_year_gsc(fit_gsc,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  
  # ASC in augsynth
  asc <-medianAbsPercBiasATT_year_asc(fit_asc,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  
  k <- c("k=1","k=3","k=7", "NA", "NA")
  Model <- c(rep("MC in gsynth", 3), "GSC", "ASC")
  Median <- c(gsc_mc_1, gsc_mc_3, gsc_mc_7, gsc, asc)
  df <- data.frame(Model,Median,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1", "NA"),
           Model = fct_relevel(Model, "MC in gsynth", "GSC", "ASC"))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Median, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7', 'NA')) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(Median,2)),position = position_dodge(1), hjust = -.2, size=2.8)+
    #ggtitle("Absolute Percent Bias of the Counterfactuals (Non-Rare)") +
    coord_flip(ylim = ylim ) + 
    scale_x_discrete(limits = rev) +
    #scale_y_continuous(name="Median", limits=c(0, 90000)) +
    theme_bw() + #adds plot outline vs. theme_minimal
    theme(axis.text=element_text(color="black"),
          legend.title = element_blank()
    ) 
}

freq_table <- function(Y0_matrix_1, Y1_matrix_1, true_trt_effect_1, fit_gsc_mc_1,
                       Y0_matrix_3, Y1_matrix_3, true_trt_effect_3, fit_gsc_mc_3,
                       Y0_matrix_7, Y1_matrix_7, true_trt_effect_7, fit_gsc_mc_7,
                       fit_gsc, fit_asc, pop_trt, ind){
  # Matrix completion in gsynth
  Y0_gsc_mc_1 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_1,Y0_matrix_1,ind)$median
  Y0_gsc_mc_3 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_3,Y0_matrix_3,ind)$median
  Y0_gsc_mc_7 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_7,Y0_matrix_7,ind)$median
  
  # GSC in gsynth
  Y0_gsc <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc,Y0_matrix_1,ind)$median
  
  # ASC in augsynth
  Y0_asc <- medianAbsPercBiasY0_year_asc(pop_trt,fit_asc,Y0_matrix_1,ind)$median
  
  # Matrix completion in gsynth
  ATT_gsc_mc_1 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)$median
  ATT_gsc_mc_3 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)$median
  ATT_gsc_mc_7 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)$median
  
  # GSC in gsynth
  ATT_gsc <-medianAbsPercBiasATT_year_gsc(fit_gsc,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)$median
  
  # ASC in augsynth
  ATT_asc <-medianAbsPercBiasATT_year_asc(fit_asc,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)$median
  
  # Matrix completion in gsynth
  UB_Y0_gsc_mc_1 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_1,Y0_matrix_1,ind)$UB
  UB_Y0_gsc_mc_3 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_3,Y0_matrix_3,ind)$UB
  UB_Y0_gsc_mc_7 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_7,Y0_matrix_7,ind)$UB
  
  # GSC in gsynth
  UB_Y0_gsc <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc,Y0_matrix_1,ind)$UB
  
  # ASC in augsynth
  UB_Y0_asc <- medianAbsPercBiasY0_year_asc(pop_trt,fit_asc,Y0_matrix_1,ind)$UB
  
  # Matrix completion in gsynth
  UB_ATT_gsc_mc_1 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)$UB
  UB_ATT_gsc_mc_3 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)$UB
  UB_ATT_gsc_mc_7 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)$UB
  
  # GSC in gsynth
  UB_ATT_gsc <-medianAbsPercBiasATT_year_gsc(fit_gsc,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)$UB
  
  # ASC in augsynth
  UB_ATT_asc <-medianAbsPercBiasATT_year_asc(fit_asc,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)$UB
  
  # Matrix completion in gsynth
  LB_Y0_gsc_mc_1 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_1,Y0_matrix_1,ind)$LB
  LB_Y0_gsc_mc_3 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_3,Y0_matrix_3,ind)$LB
  LB_Y0_gsc_mc_7 <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc_mc_7,Y0_matrix_7,ind)$LB
  
  # GSC in gsynth
  LB_Y0_gsc <- medianAbsPercBiasY0_year_gsc(pop_trt,fit_gsc,Y0_matrix_1,ind)$LB
  
  # ASC in augsynth
  LB_Y0_asc <- medianAbsPercBiasY0_year_asc(pop_trt,fit_asc,Y0_matrix_1,ind)$LB
  
  # Matrix completion in gsynth
  LB_ATT_gsc_mc_1 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)$LB
  LB_ATT_gsc_mc_3 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)$LB
  LB_ATT_gsc_mc_7 <-medianAbsPercBiasATT_year_gsc(fit_gsc_mc_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)$LB
  
  # GSC in gsynth
  LB_ATT_gsc <-medianAbsPercBiasATT_year_gsc(fit_gsc,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)$LB
  
  # ASC in augsynth
  LB_ATT_asc <-medianAbsPercBiasATT_year_asc(fit_asc,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)$LB
  
  
  # Set up table
  Model <- c(rep("MC in gsynth", 3), "GSC", "ASC")
  k <- c("1","3","7", "NA", "NA")
  absPercBias_Y0 <- round(c(Y0_gsc_mc_1, Y0_gsc_mc_3, Y0_gsc_mc_7, Y0_gsc, Y0_asc),2)
  LB_absPercBias_Y0 <- round(c(LB_Y0_gsc_mc_1, LB_Y0_gsc_mc_3, LB_Y0_gsc_mc_7, LB_Y0_gsc, LB_Y0_asc),2)
  UB_absPercBias_Y0 <- round(c(UB_Y0_gsc_mc_1, UB_Y0_gsc_mc_3, UB_Y0_gsc_mc_7, UB_Y0_gsc, UB_Y0_asc),2)
  
  absPercBias_ATT <- round(c(ATT_gsc_mc_1, ATT_gsc_mc_3, ATT_gsc_mc_7, ATT_gsc, ATT_asc),2)
  LB_absPercBias_ATT <- round(c(LB_ATT_gsc_mc_1, LB_ATT_gsc_mc_3, LB_ATT_gsc_mc_7,LB_ATT_gsc, LB_ATT_asc),2)
  UB_absPercBias_ATT <- round(c(UB_ATT_gsc_mc_1, UB_ATT_gsc_mc_3, UB_ATT_gsc_mc_7,UB_ATT_gsc, UB_ATT_asc),2)
  
  df <- cbind(Model,k,absPercBias_Y0,LB_absPercBias_Y0,UB_absPercBias_Y0,absPercBias_ATT,LB_absPercBias_ATT,UB_absPercBias_ATT)
  
  return(df)
}

