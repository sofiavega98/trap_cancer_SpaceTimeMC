## Set up color pallette using color blind friendly colors ##
## Resource: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with grey:
cbPalette <- c(  "#0072B2","#E69F00", "#D55E00", "#CC79A7","#999999", "#E69F00","#F0E442", "#56B4E9", "#009E73",  "#CC79A7")

# To use for fills, add scale_fill_manual(values=cbPalette)

# To use for line and point colors, add scale_colour_manual(values=cbPalette)


# This is a function to set up the data after simulation results
## data_big is the combined simulatin results list
setup_sim_data <- function(data_big){
  num_sim <- length(data_big) 
  ## 1 - Y0_matrix
  Y0_matrix_datalist = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y0_matrix"]
    Y0_matrix_datalist[[i]] <- dat
  }
  
  Y0_matrix <- dplyr::bind_cols(Y0_matrix_datalist) %>% as.matrix() %>% unname()
  
  ## 2 - Y1_matrix
  Y1_matrix_datalist = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y1_matrix"]
    Y1_matrix_datalist[[i]] <- dat
  }
  
  Y1_matrix <- dplyr::bind_cols(Y1_matrix_datalist) %>% as.matrix() %>% unname()
  
  ## 3 - Mu_trt_ori
  Mu_trt_ori = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Mu_trt_ori"]
    Mu_trt_ori = c(Mu_trt_ori, dat)
  }
  
  ## 4 - Mu_trt_space
  Mu_trt_space = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Mu_trt_space"]
    Mu_trt_space = c(Mu_trt_space, dat)
  }
  
  
  ## 5 - Mu_trt_spacetime ICAR(1)
  Mu_trt_spacetime_ICAR = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Mu_trt_spacetime_ICAR"]
    Mu_trt_spacetime_ICAR = c(Mu_trt_spacetime_ICAR, dat)
  }
  
  
  ## 6 - Mu_trt_spacetime AR(1)
  Mu_trt_spacetime_AR = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Mu_trt_spacetime_AR"]
    Mu_trt_spacetime_AR = c(Mu_trt_spacetime_AR, dat)
  }
  
  ## 7 - Mu_trt_space_comb_time_shrink
  
  Mu_trt_space_comb_time_shrink = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Mu_trt_space_comb_time_shrink"]
    Mu_trt_space_comb_time_shrink = c(Mu_trt_space_comb_time_shrink, dat)
  }
  
  ## 8 - Mu_trt_lasso
  
  Mu_trt_lasso = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Mu_trt_lasso"]
    Mu_trt_lasso = c(Mu_trt_lasso, dat)
  }
  
  ## 2 - Mu_trt_ori
  Y_pred_ori = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y_pred_ori"]
    Y_pred_ori = c(Y_pred_ori, dat)
  }
  
  ## 4 - Mu_trt_space
  Y_pred_space = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y_pred_space"]
    Y_pred_space = c(Y_pred_space, dat)
  }
  
  
  ## 5 - Mu_trt_spacetime ICAR(1)
  Y_pred_spacetime_ICAR = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y_pred_spacetime_ICAR"]
    Y_pred_spacetime_ICAR = c(Y_pred_spacetime_ICAR, dat)
  }
  
  
  ## 6 - Mu_trt_spacetime AR(1)
  Y_pred_spacetime_AR = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y_pred_spacetime_AR"]
    Y_pred_spacetime_AR = c(Y_pred_spacetime_AR, dat)
  }
  
  ## 7 - Mu_trt_space_comb_time_shrink
  
  Y_pred_space_comb_time_shrink = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y_pred_space_comb_time_shrink"]
    Y_pred_space_comb_time_shrink = c(Y_pred_space_comb_time_shrink, dat)
  }
  
  ## 8 - Mu_trt_lasso
  
  Y_pred_lasso = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y_pred_lasso"]
    Y_pred_lasso = c(Y_pred_lasso, dat)
  }
  
  # 7 - Mu_trt_space_comb_time_shrink
  
  Y1_obs_rate = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y1_obs_rate"]
    Y1_obs_rate = c(Y1_obs_rate, dat)
  }
  
  ## 8 - Mu_trt_lasso
  
  Y0_full_rate = list()
  
  for (i in 1:num_sim){
    dat <- data_big[[i]]["Y0_full_rate"]
    Y0_full_rate = c(Y0_full_rate, dat)
  }
 
  
  
  ### quick clean
  
  rm(i, dat, Y0_matrix_datalist, Y1_matrix_datalist)
  
  ### return lists
  return(list(Y0_matrix,Y1_matrix,Mu_trt_ori,Mu_trt_space,Mu_trt_spacetime_ICAR,Mu_trt_spacetime_AR,
           Mu_trt_space_comb_time_shrink,Mu_trt_lasso,
           Y_pred_ori,Y_pred_space,Y_pred_spacetime_ICAR,Y_pred_spacetime_AR,
           Y_pred_space_comb_time_shrink,Y_pred_lasso,
           Y0_full_rate,Y1_obs_rate))
  
}
#------------------------------------------------------------------------------------------------------
# take difference in  estimated Y(0) and true Y(0) and convert to rate
# Note: this works under the assumption that for the simulations we tested two possible 
# values of k 

diff_post_obs_rate <- function(Y0_matrix,Mu_trt,pop_trt){
  ## k = 3
  full_diff_list_1_3 <- list() #initialize list of differences for 45 samples
  full_rate_list_1_3 <- list() #initialize list of rates of differences for 45 samples
  for(i in 1:ncol(Y0_matrix)) {
    full_diff_list_1_3[[i]] <- sweep(Mu_trt[[i]],2,Y0_matrix[,i]) #posterior samples - observed Y0 values
    full_rate_list_1_3[[i]] <- 100000 * sweep(full_diff_list_1_3[[i]],2,pop_trt,"/")  #difference/population each is 1000x90
  }
  
  ## k = 7
  #full_diff_list_1_7 <- list() #initialize list of differences for 45 samples
  #full_rate_list_1_7 <- list() #initialize list of rates of differences for 45 samples
  #for(i in 1:ncol(Y0_matrix)) {
  #  full_diff_list_1_7[[i]] <- sweep(Mu_trt[[i]][[2]],2,Y0_matrix[,i]) #posterior samples - observed Y0 values
  # full_rate_list_1_7[[i]] <- 100000 * sweep(full_diff_list_1_7[[i]],2,pop_trt,"/")  #difference/population each is 1000x90
  #}
  
#return(list(full_rate_list_1_3,full_rate_list_1_7))
  return(full_rate_list_1_3)

}

diff_post_obs <- function(Y0_matrix,Mu_trt,pop_trt){
  full_diff_list_1_3 <- list() #initialize list of differences for 45 samples
  full_rate_list_1_3 <- list() #initialize list of rates of differences for 45 samples
  for(i in 1:ncol(Y0_matrix)) {
    full_diff_list_1_3[[i]] <- sweep(Mu_trt[[i]],2,Y0_matrix[,i]) #posterior samples - observed Y0 values
    full_rate_list_1_3[[i]] <- 100000 * sweep(full_diff_list_1_3[[i]],2,pop_trt,"/")  #difference/population each is 1000x90
  }
  
  return(full_diff_list_1_3)
  
}

#-----------------------------------------------------------------------------------------------------
# Get median of 1000 MCMC posterior samples for each Y(0) difference (length 90) and 
# Take Median of these medians to get one median est Y(0) - true Y(0) for each simulation



get_med_raw_bias <- function(full_rate_list){
  median_rate_list <- list()
  for(i in 1:length(full_rate_list)) {
    
    median_rate_list[[i]] <- colMedians(full_rate_list[[i]]) #gets median of 1000 post samples
    
  }
  
  med <- c()
  for(i in 1:length(median_rate_list)){
    med[i] <- median(median_rate_list[[i]])
  }
  
  return(med)
}

get_med_raw_bias_gsc <- function(pop_trt,fit_gsc,Y0_matrix){
  raw_bias_list <- list()
  for(i in 1:length(fit_gsc)) {
    raw_bias_list[[i]] <- 100000 * (as.vector(fit_gsc[[i]]$fit_gsc$Y.ct) - Y0_matrix[,i])/pop_trt
  }
  
  med <- c()
  for(i in 1:length(raw_bias_list)){
    med[i] <- median(raw_bias_list[[i]])
  }
  
  return(med)
}



get_mean_raw_bias <- function(full_rate_list){
  mean_rate_list <- list()
  for(i in 1:length(full_rate_list)) {
    
    mean_rate_list[[i]] <- colMeans(full_rate_list[[i]]) #gets mean of 1000 post samples
    
  }
  
  mean <- c()
  for(i in 1:length(mean_rate_list)){
    mean[i] <- mean(mean_rate_list[[i]])
  }
  
  return(mean)
}

get_mean_raw_bias_gsc <- function(pop_trt,fit_gsc,Y0_matrix){
  raw_bias_list <- list()
  for(i in 1:length(fit_gsc)) {
    raw_bias_list[[i]] <- 100000 * (as.vector(fit_gsc[[i]]$fit_gsc$Y.ct) - Y0_matrix[,i])/pop_trt
  }
  
  mean <- c()
  for(i in 1:length(raw_bias_list)){
    mean[i] <- mean(raw_bias_list[[i]])
  }
  
  return(mean)
}
#---------------------------------------------------------------------------------------
# Percent Bias of the counterfactuals

get_perc_bias_Y0 <- function(full_rate_list,Y0_matrix,pop_trt,ind){
  # full rate list is the raw bias in rates for each of the posterior samples
  median_rate_list <- list()
  perc_bias_list <- list()
  avg_perc_bias <- list()
  mean_perc_bias <- list()
  for(i in 1:length(full_rate_list)) {
    true_Y0_rate <- 100000 * (Y0_matrix[,i]/pop_trt)
    median_rate_list[[i]] <- colMedians(full_rate_list[[i]]) #gets median of 1000 post samples
    perc_bias_list[[i]] <- 100*(median_rate_list[[i]])/true_Y0_rate #(median trt effect - true trt effect)/true trt
    avg_perc_bias[i] <- median(abs(perc_bias_list[[i]][ind]))
    mean_perc_bias[i] <- mean(abs(perc_bias_list[[i]][ind]))
  }
  
  
  return(list(avg_perc_bias=avg_perc_bias, mean_perc_bias = mean_perc_bias))
}



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

plot_avg_perc_bias_Y0_gsc <- function(ori_3,space_3,spacetime_ICAR_3,spacetime_AR_3,spacetime_shrink,lasso, gsc){
  num_sim <- length(ori_3[[1]])
  #avg percent bias
  percent_bias_df <- abs(data.frame(ori_3 = unlist(ori_3$avg_perc_bias),space_3 = unlist(space_3$avg_perc_bias),
                                    spacetime_ICAR_3 = unlist(spacetime_ICAR_3$avg_perc_bias),
                                    spacetime_AR_3 = unlist(spacetime_AR_3$avg_perc_bias),
                                    spacetime_shrink = unlist(spacetime_shrink$avg_perc_bias),
                                    lasso = unlist(lasso$avg_perc_bias),
                                    gsc = unlist(gsc$avg_perc_bias)))
  percent_bias_df.m <- melt(percent_bias_df)
  k <- c(rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),
         rep("k_max=7",num_sim),rep("k=3",num_sim),rep("k=3",num_sim))
  model <- c(rep("Original",num_sim),rep("Space",num_sim),rep("Space-Time \n ICAR",num_sim),
             rep("Space-Time \n AR(1)",num_sim),rep("Space-Time \n Shrinkage",num_sim),
             rep("Shrinkage \n Lasso",num_sim),
             rep("Gsynth", num_sim))
  
  
  # Plot
  df_bias <- cbind(percent_bias_df.m ,k,model)
  p <- df_bias  %>%
    # Add a column called 'type': do we want to highlight the group or not?
    #mutate( model=ifelse(variable=="spacetime_shrink","Shrinkage Spatio-Temporal","Spatio-Temporal")) %>%
    ggplot( aes(x=model, y=value, color = k)) + 
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = quantile(df_bias$value, c(0.1, .999), na.rm = T)) + 
    xlab("Model") + ylab("Absolute Percent Bias") +
    scale_color_manual(values=c( "#00AFBB","#E7B800","#FC4E07","#E69F00", "#56B4E9")) + 
    ggtitle("Comparing Absolute Percent Biases of Y(0)")
  
  # Table
  
  temp <- as.data.frame(do.call(cbind, lapply(percent_bias_df, function(x) round(summary(x), digits = 3))))
  
  df <- cbind(temp[,1:7])
  df[7,1] <- 0
  df[7,2] <- 0
  df[7,3] <- 0
  df[7,4] <- 0
  df[7,5] <- 0
  df[7,7] <- 0
  
  colnames(df) <- c("k=3","k=3","k=3","k=3","k_max=7","k=3","k=3")
  t <- df %>%
    kbl(caption = "Summary Statistics for Absolute Percent Bias") %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c(" " = 1, "Original Model" = 1, "Space Model" = 1, "Space-Time ICAR Model" = 1, "Space-Time AR(1) Model" = 1,"Spatio-Temporal \n Shrinkage Model" = 1, "Shrinkage \n Lasso" = 1, "Gsynth Model" = 1))
  
  return(list(p,t))
}

plot_avg_perc_bias_Y0_gsc_log <- function(ori_3,space_3,spacetime_ICAR_3,spacetime_AR_3,spacetime_shrink, gsc){
  num_sim <- length(ori_3[[1]])
  #avg percent bias
  percent_bias_df <- abs(data.frame(ori_3 = unlist(ori_3$avg_perc_bias),space_3 = unlist(space_3$avg_perc_bias),
                                    spacetime_ICAR_3 = unlist(spacetime_ICAR_3$avg_perc_bias),
                                    spacetime_AR_3 = unlist(spacetime_AR_3$avg_perc_bias),
                                    spacetime_shrink = unlist(spacetime_shrink$avg_perc_bias),
                                    lasso = unlist(lasso$avg_perc_bias),
                                    gsc = unlist(gsc$avg_perc_bias)))
  percent_bias_df.m <- melt(percent_bias_df)
  k <- c(rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k_max=7",num_sim),rep("k=3",num_sim),rep("k=3",num_sim))
  model <- c(rep("Original",num_sim),rep("Space",num_sim),rep("Space-Time \n ICAR",num_sim),
             rep("Space-Time \n AR(1)",num_sim),rep("Space-Time \n Shrinkage",num_sim),
             rep("Shrinkage \n Lasso",num_sim),
             rep("Gsynth", num_sim))
  
  
  # Plot
  df_bias <- cbind(percent_bias_df.m ,k,model)
  p <- df_bias  %>%
    # Add a column called 'type': do we want to highlight the group or not?
    #mutate( model=ifelse(variable=="spacetime_shrink","Shrinkage Spatio-Temporal","Spatio-Temporal")) %>%
    ggplot( aes(x=model, y=log10(value), color = k)) + 
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = quantile(df_bias$value, c(0.1, .97), na.rm = T)) + 
    xlab("Model") + ylab("Absolute Percent Bias") +
    scale_color_manual(values=c( "#00AFBB","#E7B800","#FC4E07","#E69F00", "#56B4E9")) + 
    ggtitle("Comparing Absolute Percent Biases of Y(0)")
  
  # Table
  
  temp <- as.data.frame(do.call(cbind, lapply(percent_bias_df, function(x) round(summary(x), digits = 3))))
  
  df <- cbind(temp[,1:7])
  #df[7,6] <- 0
  
  colnames(df) <- c("k=3","k=3","k=3","k=3","k_max=7","k=3","k=3")
  t <- df %>%
    kbl(caption = "Summary Statistics for Absolute Percent Bias") %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c(" " = 1, "Original Model" = 1, "Space Model" = 1, "Space-Time ICAR Model" = 1, "Space-Time AR(1) Model" = 1,"Spatio-Temporal \n Shrinkage Model" = 1, "Shrinkage \n Lasso" = 1, "Gsynth Model" = 1))
  
  return(list(p,t))
}

#----------------------------------------------------------------------------------------------------

# Plot raw biases box plot and table
raw_bias_box_2k <- function(med_1_3,med_1_7,med_2_3,med_2_7,med_3_3,med_3_7,med_4_3,med_4_7,med_5){
  num_sim <- length(med_1_3)
  med_df <- abs(data.frame(original_3 = med_1_3,original_7 = med_1_7, space_3 = med_2_3, space_7 = med_2_7, spacetime_ICAR_3 = med_3_3,spacetime_ICAR_7 = med_3_7,spacetime_AR_3 = med_4_3,spacetime_AR_7 = med_4_7,spacetime_shrink=med_5))
  med.m <- melt(med_df)
  k <- c(rep("k=3",num_sim),rep("k=7",num_sim),rep("k=3",num_sim),rep("k=7",num_sim),rep("k=3",num_sim),rep("k=7",num_sim),rep("k=3",num_sim),rep("k=7",num_sim),rep("k_max=7",num_sim))
  model <- c(rep("Original",num_sim*2),rep("Space",num_sim*2),rep("Space-Time \n ICAR",num_sim*2),rep("Space-Time \n AR(1)",num_sim*2),rep("Space-Time \n Shrinkage",num_sim))
  # Plot
  df_med <- cbind(med.m,k,model)
  p <- df_med  %>%
    # Add a column called 'type': do we want to highlight the group or not?
    #mutate( model=ifelse(variable=="spacetime_shrink","Space Time Shrinkage","Spatio-Temporal")) %>%
    ggplot( aes(x=model, y=value, color = k)) + geom_boxplot() + xlab("Model") + ylab("Median Est. Y(0) - True Y(0)") +scale_color_manual(values=c( "#00AFBB","#E7B800","#FC4E07")) + ggtitle("Comparing Raw Biases of Y(0)")
  
  temp <- as.data.frame(do.call(cbind, lapply(med_df, function(x) round(summary(x), digits = 3))))
  
  df <- cbind(temp[,1:9])
  #df[7,1] <- 0
  #df[7,6] <- 0
  #df[7,7] <- 0
  #df[7,8] <- 0
  #df[7,9] <- 0
  #df[4,3] <- NA
  #df[4,5] <- NA
  colnames(df) <- c("k=3","k=7","k=3","k=7","k=3","k=7","k=3","k=7","k_max=7")
  t <- df %>%
    kbl(caption = "Summary Statistics for Estimated Y(0) - True Y(0)") %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c(" " = 1, "Original Model" = 2, "Space Model" = 2, "Space-Time ICAR Model" = 2, "Space-Time AR(1) Model" = 2,"Space-Time Shrinkage Model" = 1))
  
  return(list(p,t))
}

# Plot raw biases box plot and table
raw_bias_box <- function(med_1_3,med_2_3,med_3_3,med_4_3,med_5){
  num_sim <- length(med_1_3)
  med_df <- abs(data.frame(original_3 = med_1_3, space_3 = med_2_3,  spacetime_ICAR_3 = med_3_3,spacetime_AR_3 = med_4_3,spacetime_shrink=med_5))
  med.m <- melt(med_df)
  k <- c(rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k_max=7",num_sim))
  model <- c(rep("Original",num_sim),rep("Space",num_sim),rep("Space-Time \n ICAR",num_sim),rep("Space-Time \n AR(1)",num_sim),rep("Space-Time \n Shrinkage",num_sim))
  # Plot
  df_med <- cbind(med.m,k,model)
  p <- df_med  %>%
    # Add a column called 'type': do we want to highlight the group or not?
    #mutate( model=ifelse(variable=="spacetime_shrink","Space Time Shrinkage","Spatio-Temporal")) %>%
    ggplot( aes(x=model, y=value, color = k)) + geom_boxplot() + xlab("Model") + ylab("Median Est. Y(0) - True Y(0)") +scale_color_manual(values=c( "#00AFBB","#E7B800","#FC4E07")) + ggtitle("Comparing Raw Biases of Y(0)")
  
  temp <- as.data.frame(do.call(cbind, lapply(med_df, function(x) round(summary(x), digits = 3))))
  
  df <- cbind(temp[,1:5])
  #df[7,1] <- 0
  #df[7,6] <- 0
  #df[7,7] <- 0
  #df[7,8] <- 0
  #df[7,9] <- 0
  #df[4,3] <- NA
  #df[4,5] <- NA
  colnames(df) <- c("k=3","k=3","k=3","k=3","k_max=7")
  t <- df %>%
    kbl(caption = "Summary Statistics for Estimated Y(0) - True Y(0)") %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c(" " = 1, "Original Model" = 1, "Space Model" = 1, "Space-Time ICAR Model" = 1, "Space-Time AR(1) Model" = 1,"Space-Time Shrinkage Model" = 1))
  
  return(list(p,t))
}

# Plot raw biases box plot and table
raw_bias_box_gsc <- function(med_1_3,med_2_3,med_3_3,med_4_3,med_5,med_6,med_gsc){
  num_sim <- length(med_1_3)
  med_df <- abs(data.frame(original_3 = med_1_3, space_3 = med_2_3,  spacetime_ICAR_3 = med_3_3,spacetime_AR_3 = med_4_3,spacetime_shrink=med_5, med_lasso = med_6, gsc = med_gsc))
  med.m <- melt(med_df)
  k <- c(rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k_max=7",num_sim),rep("k=3",num_sim),rep("k=3",num_sim))
  model <- c(rep("Original",num_sim),rep("Space",num_sim),rep("Space-Time \n ICAR",num_sim),rep("Space-Time \n AR(1)",num_sim),rep("Space-Time \n Shrinkage",num_sim),rep("Shrinkage \n Lasso", num_sim),rep("Gsynth",num_sim))
  # Plot
  df_med <- cbind(med.m,k,model)
  p <- df_med  %>%
    # Add a column called 'type': do we want to highlight the group or not?
    #mutate( model=ifelse(variable=="spacetime_shrink","Space Time Shrinkage","Spatio-Temporal")) %>%
    ggplot( aes(x=model, y=value, color = k)) + geom_boxplot() + xlab("Model") + ylab("Median Est. Y(0) - True Y(0)") +scale_color_manual(values=c( "#00AFBB","#E7B800","#FC4E07")) + ggtitle("Comparing Raw Biases of Y(0)")
  
  temp <- as.data.frame(do.call(cbind, lapply(med_df, function(x) round(summary(x), digits = 3))))
  
  df <- cbind(temp[,1:7])
  df[7,1] <- 0
  df[7,2] <- 0
  df[7,3] <- 0
  df[7,4] <- 0
  df[7,5] <- 0
  df[7,7] <- 0
  #df[7,6] <- 0
  #df[7,7] <- 0
  #df[7,8] <- 0
  #df[7,9] <- 0
  #df[4,3] <- NA
  #df[4,5] <- NA
  #df <- cbind(temp[,6],temp[,1:5])
  #df[7,7] <- 0
  colnames(df) <- c("k=3","k=3","k=3","k=3","k_max=7", "k=3", "k=3")
  t <- df %>%
    kbl(caption = "Summary Statistics for Estimated Y(0) - True Y(0)") %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c(" " = 1, "Original Model" = 1, "Space Model" = 1, "Space-Time ICAR Model" = 1, "Space-Time AR(1) Model" = 1,"Space-Time Shrinkage Model" = 1, "Shrinkage Lasso" = 1, "Gsynth Model" = 1))
  
  return(list(p,t))
}
#---------------------------------------------------------------------------------------------------
# percent bias function
# average across timepoints
perc_bias_function_2 <- function(Mu_trt,Y1_matrix,Y0_matrix,pop_trt,true_trt_effect,ind) {
  med_list <- list()
  perc_bias_list <- list()
  avg_perc_bias <- list()
  true_trt_rate_list <- list()
  mean_perc_bias <- list()
  for(i in 1:ncol(Y1_matrix)) {
    full_diff <- -1*sweep(Mu_trt[[i]][,ind],2,Y1_matrix[,i], FUN="-") #-1(estimated Y(0) - observed Y(1)) (42 x 1000)
    full_rate <- 100000 * sweep(full_diff,2,pop_trt[ind],"/") #difference/population*100000
    med_list[[i]] <- colMedians(full_rate) #get avg of all MCMC in each time point can I use median here? (42 x 1)
    true_trt_rate_list[[i]] <- 100000*(true_trt_effect[,i]/pop_trt[ind]) #(42 x 1)
    perc_bias_list[[i]] <- 100*(med_list[[i]]-true_trt_rate_list[[i]])/true_trt_rate_list[[i]] #(median trt effect - true trt effect)/true trt
    avg_perc_bias[i] <- median(abs(perc_bias_list[[i]]))
    mean_perc_bias[i] <- mean(abs(perc_bias_list[[i]]))
  }
  
  return(list(med_list=med_list, true_trt_rate_list=true_trt_rate_list, perc_bias_list=perc_bias_list,avg_perc_bias=avg_perc_bias,mean_perc_bias = mean_perc_bias))
  # med list returns estimate of the treatment effect
  
}


perc_bias_function_gsc <- function(fit_gsc,Y1_matrix,Y0_matrix,pop_trt,true_trt_effect,ind) {
  med_list <- list()
  perc_bias_list <- list()
  avg_perc_bias <- list()
  true_trt_rate_list <- list()
  mean_perc_bias <- list()
  for(i in 1:ncol(Y1_matrix)) {
    med_list[[i]] <- (100000 * Y1_matrix[,i]/pop_trt[ind]) - as.vector(fit_gsc[[i]]$fit_gsc$Y.ct)[ind]
    true_trt_rate_list[[i]] <- 100000*(true_trt_effect[,i]/pop_trt[ind]) #(42 x 1)
    perc_bias_list[[i]] <- 100*(med_list[[i]]-true_trt_rate_list[[i]])/true_trt_rate_list[[i]] #(median trt effect - true trt effect)/true trt
    avg_perc_bias[i] <- median(abs(perc_bias_list[[i]]))
    mean_perc_bias[i] <- mean(abs(perc_bias_list[[i]]))
  }
  
  return(list(med_list=med_list, true_trt_rate_list=true_trt_rate_list, perc_bias_list=perc_bias_list,avg_perc_bias=avg_perc_bias, mean_perc_bias = mean_perc_bias))
  # med list returns estimate of the treatment effect
  
}

#---------------------------------------------------------------------------------------------------
# percent bias function
# average across counties
#this doesnt do anything yet
perc_bias_function_counties <- function(Mu_trt,Y1_matrix,Y0_matrix,pop_trt,true_trt_effect,ind) {
  med_list <- list()
  perc_bias_list <- list()
  avg_perc_bias <- list()
  for(i in 1:ncol(Y1_matrix)) {
    full_diff <- -1*sweep(Mu_trt[[i]][,49:90],2,Y1_matrix[,i]) #-1(estimated Y(0) - observed Y(1))
    full_rate <- 100000 * sweep(full_diff,2,pop_trt[49:90],"/") #difference/population*100000
    med_list[[i]] <- colMedians(full_rate) #get avg of all MCMC in each time point can I use median here? (42 x 1)
    true_trt_rate <- 100000*(true_trt_effect[,i]/pop_trt[49:90]) #(42 x 1)
    perc_bias_list[[i]] <- 100*(med_list[[i]]-true_trt_rate)/true_trt_rate #(median trt effect - true trt effect)/true trt
    avg_perc_bias[i] <- median(perc_bias_list[[i]])
  }
  
  return(list(med_list=med_list,perc_bias_list=perc_bias_list,avg_perc_bias=avg_perc_bias))
}

#----------------------------------------------------------------------------------------------------
# plot average percent bias
plot_avg_perc_bias_2k <- function(ori_3,ori_7,space_3,space_7,spacetime_ICAR_3,spacetime_ICAR_7,spacetime_AR_3,spacetime_AR_7,spacetime_shrink){
  num_sim <- length(ori_3[[1]])
  #avg percent bias
  percent_bias_df <- abs(data.frame(ori_3 = unlist(ori_3$avg_perc_bias),ori_7 = unlist(ori_7$avg_perc_bias),space_3 = unlist(space_3$avg_perc_bias),space_7 = unlist(space_7$avg_perc_bias),spacetime_ICAR_3 = unlist(spacetime_ICAR_3$avg_perc_bias),spacetime_ICAR_7 = unlist(spacetime_ICAR_7$avg_perc_bias),spacetime_AR_3 = unlist(spacetime_AR_3$avg_perc_bias),spacetime_AR_7 = unlist(spacetime_AR_7$avg_perc_bias),spacetime_shrink = unlist(spacetime_shrink$avg_perc_bias)))
  percent_bias_df.m <- melt(percent_bias_df)
  k <- c(rep("k=3",num_sim),rep("k=7",num_sim),rep("k=3",num_sim),rep("k=7",num_sim),rep("k=3",num_sim),rep("k=7",num_sim),rep("k=3",num_sim),rep("k=7",num_sim),rep("k_max=7",num_sim))
  model <- c(rep("Original",num_sim*2),rep("Space",num_sim*2),rep("Space-Time \n ICAR",num_sim*2),rep("Space-Time \n AR(1)",num_sim*2),rep("Space-Time \n Shrinkage",num_sim))
  
  
  # Plot
  df_bias <- cbind(percent_bias_df.m ,k,model)
  p <- df_bias  %>%
    # Add a column called 'type': do we want to highlight the group or not?
    #mutate( model=ifelse(variable=="spacetime_shrink","Shrinkage Spatio-Temporal","Spatio-Temporal")) %>%
    ggplot( aes(x=model, y=value, color = k)) + geom_boxplot() + xlab("Model") + ylab("Absolute Percent Bias") +scale_color_manual(values=c( "#00AFBB","#E7B800","#FC4E07","#E69F00", "#56B4E9")) + ggtitle("Comparing Absolute Percent Biases")
  
  # Table
  
  temp <- as.data.frame(do.call(cbind, lapply(percent_bias_df, function(x) round(summary(x), digits = 3))))
  
  df <- cbind(temp[,1:9])
  
  colnames(df) <- c("k=3","k=7","k=3","k=7","k=3","k=7","k=3","k=7","k_max=7")
  t <- df %>%
    kbl(caption = "Summary Statistics for Absolute Percent Bias") %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c(" " = 1, "Original Model" = 2, "Space Model" = 2, "Space-Time ICAR Model" = 2, "Space-Time AR(1) Model" = 2,"Spatio-Temporal \n Shrinkage Model" = 1))

  return(list(p,t))
}

plot_avg_perc_bias <- function(ori_3,space_3,spacetime_ICAR_3,spacetime_AR_3,spacetime_shrink){
  num_sim <- length(ori_3[[1]])
  #avg percent bias
  percent_bias_df <- abs(data.frame(ori_3 = unlist(ori_3$avg_perc_bias),space_3 = unlist(space_3$avg_perc_bias),
                                    spacetime_ICAR_3 = unlist(spacetime_ICAR_3$avg_perc_bias),
                                    spacetime_AR_3 = unlist(spacetime_AR_3$avg_perc_bias),
                                   spacetime_shrink = unlist(spacetime_shrink$avg_perc_bias)))
  percent_bias_df.m <- melt(percent_bias_df)
  k <- c(rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k_max=7",num_sim))
  model <- c(rep("Original",num_sim),rep("Space",num_sim),rep("Space-Time \n ICAR",num_sim),rep("Space-Time \n AR(1)",num_sim),rep("Space-Time \n Shrinkage",num_sim))
  
  
  # Plot
  df_bias <- cbind(percent_bias_df.m ,k,model)
  p <- df_bias  %>%
    # Add a column called 'type': do we want to highlight the group or not?
    #mutate( model=ifelse(variable=="spacetime_shrink","Shrinkage Spatio-Temporal","Spatio-Temporal")) %>%
    ggplot( aes(x=model, y=value, color = k)) + geom_boxplot(outlier.shape = NA) + xlab("Model") + ylab("Absolute Percent Bias") +scale_color_manual(values=c( "#00AFBB","#E7B800","#FC4E07","#E69F00", "#56B4E9")) + ggtitle("Comparing Absolute Percent Biases")
  
  # Table
  
  temp <- as.data.frame(do.call(cbind, lapply(percent_bias_df, function(x) round(summary(x), digits = 3))))
  
  df <- cbind(temp[,1:5])
  
  colnames(df) <- c("k=3","k=3","k=3","k=3","k_max=7")
  t <- df %>%
    kbl(caption = "Summary Statistics for Absolute Percent Bias") %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c(" " = 1, "Original Model" = 1, "Space Model" = 1, "Space-Time ICAR Model" = 1, "Space-Time AR(1) Model" = 1,"Spatio-Temporal \n Shrinkage Model" = 1))
  
  return(list(p,t))
}

plot_avg_perc_bias_gsc <- function(ori_3,space_3,spacetime_ICAR_3,spacetime_AR_3,spacetime_shrink,lasso_3, gsc){
  num_sim <- length(ori_3[[1]])
  #avg percent bias
  percent_bias_df <- abs(data.frame(ori_3 = unlist(ori_3$avg_perc_bias),space_3 = unlist(space_3$avg_perc_bias),
                                    spacetime_ICAR_3 = unlist(spacetime_ICAR_3$avg_perc_bias),
                                    spacetime_AR_3 = unlist(spacetime_AR_3$avg_perc_bias),
                                    spacetime_shrink = unlist(spacetime_shrink$avg_perc_bias),
                                    lasso_3 = unlist(lasso_3$avg_perc_bias),
                                    gsc = unlist(gsc$avg_perc_bias)))
  percent_bias_df.m <- melt(percent_bias_df)
  k <- c(rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k_max=7",num_sim),rep("k=3",num_sim),rep("k=3",num_sim))
  model <- c(rep("Original",num_sim),rep("Space",num_sim),rep("Space-Time \n ICAR",num_sim),
             rep("Space-Time \n AR(1)",num_sim),rep("Space-Time \n Shrinkage",num_sim),
             rep("Shrinkage \n Lasso", num_sim),
             rep("Gsynth", num_sim))
  
  
  # Plot
  df_bias <- cbind(percent_bias_df.m ,k,model)
  p <- df_bias  %>%
    # Add a column called 'type': do we want to highlight the group or not?
    #mutate( model=ifelse(variable=="spacetime_shrink","Shrinkage Spatio-Temporal","Spatio-Temporal")) %>%
    ggplot( aes(x=model, y=value, color = k)) + 
    geom_boxplot(outlier.shape = NA) + 
    coord_cartesian(ylim = quantile(df_bias$value, c(0.1, .9999), na.rm = T)) + 
    xlab("Model") + ylab("Absolute Percent Bias") +
    scale_color_manual(values=c( "#00AFBB","#E7B800","#FC4E07","#E69F00", "#56B4E9")) #+ 
    #ggtitle("Comparing Absolute Percent Biases of the ATT")
  
  # Table
  
  temp <- as.data.frame(do.call(cbind, lapply(percent_bias_df, function(x) round(summary(x), digits = 3))))
  
  df <- cbind(temp[,1:7])
  #df[7,6] <- 0
  df[7,1] <- 0
  df[7,2] <- 0
  df[7,3] <- 0
  df[7,4] <- 0
  df[7,5] <- 0
  df[7,7] <- 0
  
  colnames(df) <- c("k=3","k=3","k=3","k=3","k_max=7","k=3","k=3")
  t <- df %>%
    kbl(caption = "Summary Statistics for Absolute Percent Bias") %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c(" " = 1, "Original Model" = 1, "Space Model" = 1, "Space-Time ICAR Model" = 1, "Space-Time AR(1) Model" = 1,"Spatio-Temporal \n Shrinkage Model" = 1, "Shrinkage Lasso" = 1,"Gsynth Model"=1))
  
  return(list(p,t))
}


plot_avg_perc_bias_gsc_log <- function(ori_3,space_3,spacetime_ICAR_3,spacetime_AR_3,spacetime_shrink, gsc){
  num_sim <- length(ori_3[[1]])
  #avg percent bias
  percent_bias_df <- abs(data.frame(ori_3 = unlist(ori_3$avg_perc_bias),space_3 = unlist(space_3$avg_perc_bias),
                                    spacetime_ICAR_3 = unlist(spacetime_ICAR_3$avg_perc_bias),
                                    spacetime_AR_3 = unlist(spacetime_AR_3$avg_perc_bias),
                                    spacetime_shrink = unlist(spacetime_shrink$avg_perc_bias),
                                    gsc = unlist(gsc$avg_perc_bias)))
  percent_bias_df.m <- melt(percent_bias_df)
  k <- c(rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k=3",num_sim),rep("k_max=7",num_sim),rep("k=3",num_sim))
  model <- c(rep("Original",num_sim),rep("Space",num_sim),rep("Space-Time \n ICAR",num_sim),
             rep("Space-Time \n AR(1)",num_sim),rep("Space-Time \n Shrinkage",num_sim),
             rep("Gsynth", num_sim))
  
  
  # Plot
  df_bias <- cbind(percent_bias_df.m ,k,model)
  p <- df_bias  %>%
    # Add a column called 'type': do we want to highlight the group or not?
    #mutate( model=ifelse(variable=="spacetime_shrink","Shrinkage Spatio-Temporal","Spatio-Temporal")) %>%
    ggplot( aes(x=model, y=log10(value), color = k)) + 
    geom_boxplot() + 
    xlab("Model") + ylab(bquote(~log[10]~"(Absolute Percent Bias)")) +
    scale_color_manual(values=c( "#00AFBB","#E7B800","#FC4E07","#E69F00", "#56B4E9")) #+ 
    #ggtitle("Comparing Absolute Percent Biases of the ATT")
  
  # Table
  
  temp <- as.data.frame(do.call(cbind, lapply(percent_bias_df, function(x) round(summary(x), digits = 3))))
  
  df <- cbind(temp[,1:6])
  df[7,6] <- 0
  
  colnames(df) <- c("k=3","k=3","k=3","k=3","k_max=7","k=3")
  t <- df %>%
    kbl(caption = "Summary Statistics for Absolute Percent Bias") %>%
    kable_classic_2(full_width = F) %>%
    add_header_above(c(" " = 1, "Original Model" = 1, "Space Model" = 1, "Space-Time ICAR Model" = 1, "Space-Time AR(1) Model" = 1,"Spatio-Temporal \n Shrinkage Model" = 1, "Gsynth Model"))
  
  return(list(p,t))
}
#------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
raw_v_est_plot_sim <- function(data=NM_new,pop_trt,Mu_trt,Y0_matrix,n_exp) {
  full_time <- unique(data$YEAR_DX)
  
  ### Posterior
  post_Mu_trt_rate_f <- 100000 * sweep(Mu_trt,2,pop_trt,"/")
  post_Mu_trt_mat <- matrix(colMeans(post_Mu_trt_rate_f),nrow=6,byrow=T)
  post_Mu_trt_avg_rate_f <- colMeans(post_Mu_trt_mat)
  post_Mu_trt_avg_lower_f <- apply(post_Mu_trt_mat, 2, 
                                                   function(x) {quantile(x, probs = 0.05)})
  post_Mu_trt_avg_upper_f <- apply(post_Mu_trt_mat, 2, 
                                                   function(x) {quantile(x, probs = 0.95)})   
  
  #posterior sample df
  post_trt_f <- data.frame(time = full_time, mean = post_Mu_trt_avg_rate_f, 
                           lower = post_Mu_trt_avg_lower_f, upper = post_Mu_trt_avg_upper_f, 
                           type = factor("post"))
  ### Raw
  raw_rate <- 100000 * (Y0_matrix/pop_trt)
  raw_trt_avg_rate_f <- colMeans(matrix(raw_rate,nrow=6,byrow = T))
  raw_trt_avg_lower_f <- apply(matrix(raw_rate,nrow=6,byrow = T), 2, function(x) {quantile(x, probs = 0.05)})
  raw_trt_avg_upper_f <- apply(matrix(raw_rate,nrow=6,byrow = T),2,  function(x) {quantile(x, probs = 0.95)})  
  
  #observed sample df
  raw_trt_f <- data.frame(time = full_time, mean = raw_trt_avg_rate_f , 
                          lower =raw_trt_avg_lower_f , upper = raw_trt_avg_upper_f  , 
                          type = factor("raw"))
  
  
  to_plot <- rbind(post_trt_f, raw_trt_f)
  
  p <- ggplot(data = to_plot, aes(x = time, y = mean, group = type, color = type)) + geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.1) + 
    ylab( 'Number of Cases per 100000') + geom_vline(xintercept = 1996, color = 'black') +
    ggtitle("Comparing Raw Treated Average Rate with \n Estimated Posterior Average Rate for All Counties")  
  
  return(p)
  
}

raw_v_est_plot_sim_gsc <- function(data=NM_new,pop_trt,gsc_Mu,Y0_matrix,n_exp) {
  full_time <- unique(data$YEAR_DX)
  
  ### Posterior
  post_Mu_trt_rate <- 100000 * as.vector(gsc_Mu)/pop_trt
  post_Mu_trt_mat <- matrix(post_Mu_trt_rate,nrow=6,byrow=T)
  post_Mu_trt_avg_rate_f <- colMeans(post_Mu_trt_mat)
  post_Mu_trt_avg_lower_f <- apply(post_Mu_trt_mat, 2, 
                                                   function(x) {quantile(x, probs = 0.05)})
                                             
  post_Mu_trt_avg_upper_f <- apply(post_Mu_trt_mat, 2, 
                                   function(x) {quantile(x, probs = 0.95)})
  
  #posterior sample df
  post_trt_f <- data.frame(time = full_time, mean = post_Mu_trt_avg_rate_f, 
                           lower = post_Mu_trt_avg_lower_f, upper = post_Mu_trt_avg_upper_f, 
                           type = factor("post"))
  ### Raw
  raw_rate <- 100000 * (Y0_matrix/pop_trt)
  raw_trt_avg_rate_f <- colMeans(matrix(raw_rate,nrow=6,byrow = T))
  #raw_trt_avg_lower_f <- lapply(raw_trt_avg_rate_f,  function(x) {quantile(x, probs = 0.05)})
  #raw_trt_avg_upper_f <- lapply(raw_trt_avg_rate_f,  function(x) {quantile(x, probs = 0.95)})  
  
  #observed sample df
  raw_trt_f <- data.frame(time = full_time, mean = raw_trt_avg_rate_f , 
                          lower =raw_trt_avg_rate_f , upper = raw_trt_avg_rate_f , 
                          type = factor("raw"))
  
  
  to_plot <- rbind(post_trt_f, raw_trt_f)
  
  p <- ggplot(data = to_plot, aes(x = time, y = mean, group = type, color = type)) + geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.1) + 
    ylab( 'Number of Cases per 100000') + geom_vline(xintercept = 1996, color = 'black') +
    ggtitle("Comparing Raw Treated Average Rate with \n Estimated Posterior Average Rate for All Counties")  
  
  return(p)
  
}

##---------------------------------------------------------------------------------------------------
## Check Convergence of Models
# k in 1,3,7
# our six models
conv_vec <- function(Mu_trt_ori_1,Mu_trt_space_1,Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1,Mu_trt_lasso_1,
                     Mu_trt_ori_3,Mu_trt_space_3,Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3,Mu_trt_lasso_3,
                     Mu_trt_ori_7,Mu_trt_space_7,Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7,Mu_trt_lasso_7,
                     Mu_trt_space_comb_time_shrink_7){
  
  ## k = 1
  ori_1 <- c()
  for(i in 1:length(Mu_trt_ori_1)){
    if(length(unique(which(is.na(Mu_trt_ori_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      ori_1[i] <- 1
    }
  }
  ori_conv_1 <- sum(ori_1,na.rm = T)
  
  space_1 <- c()
  for(i in 1:length(Mu_trt_space_1)){
    if(length(unique(which(is.na(Mu_trt_space_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      space_1[i] <- 1
    }
  }
  space_conv_1 <- sum(space_1,na.rm = T)
  
  spacetime_ICAR_1 <- c()
  for(i in 1:length(Mu_trt_spacetime_ICAR_1)){
    if(length(unique(which(is.na(Mu_trt_spacetime_ICAR_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_ICAR_1[i] <- 1
    }
  }
  spacetime_ICAR_conv_1 <- sum(spacetime_ICAR_1,na.rm = T)
  
  spacetime_AR_1 <- c()
  for(i in 1:length(Mu_trt_spacetime_AR_1)){
    if(length(unique(which(is.na(Mu_trt_spacetime_AR_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_AR_1[i] <- 1
    }
  }
  spacetime_AR_conv_1 <- sum(spacetime_AR_1,na.rm = T)
  
  lasso_1 <- c()
  for(i in 1:length(Mu_trt_lasso_1)){
    if(length(unique(which(is.na(Mu_trt_lasso_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      lasso_1[i] <- 1
    }
  }
  lasso_conv_1 <- sum(lasso_1,na.rm = T)
  
  ## k = 3
  ori_3 <- c()
  for(i in 1:length(Mu_trt_ori_3)){
    if(length(unique(which(is.na(Mu_trt_ori_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      ori_3[i] <- 1
    }
  }
  ori_conv_3 <- sum(ori_3,na.rm = T)
  
  space_3 <- c()
  for(i in 1:length(Mu_trt_space_3)){
    if(length(unique(which(is.na(Mu_trt_space_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      space_3[i] <- 1
    }
  }
  space_conv_3 <- sum(space_3,na.rm = T)
  
  spacetime_ICAR_3 <- c()
  for(i in 1:length(Mu_trt_spacetime_ICAR_3)){
    if(length(unique(which(is.na(Mu_trt_spacetime_ICAR_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_ICAR_3[i] <- 1
    }
  }
  spacetime_ICAR_conv_3 <- sum(spacetime_ICAR_3,na.rm = T)
  
  spacetime_AR_3 <- c()
  for(i in 1:length(Mu_trt_spacetime_AR_3)){
    if(length(unique(which(is.na(Mu_trt_spacetime_AR_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_AR_3[i] <- 1
    }
  }
  spacetime_AR_conv_3 <- sum(spacetime_AR_3,na.rm = T)
  
  lasso_3 <- c()
  for(i in 1:length(Mu_trt_lasso_3)){
    if(length(unique(which(is.na(Mu_trt_lasso_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      lasso_3[i] <- 1
    }
  }
  lasso_conv_3 <- sum(lasso_3,na.rm = T)
  
  ## k = 7
  ori_7 <- c()
  for(i in 1:length(Mu_trt_ori_7)){
    if(length(unique(which(is.na(Mu_trt_ori_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      ori_7[i] <- 1
    }
  }
  ori_conv_7 <- sum(ori_7,na.rm = T)
  
  space_7 <- c()
  for(i in 1:length(Mu_trt_space_7)){
    if(length(unique(which(is.na(Mu_trt_space_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      space_7[i] <- 1
    }
  }
  space_conv_7 <- sum(space_7,na.rm = T)
  
  spacetime_ICAR_7 <- c()
  for(i in 1:length(Mu_trt_spacetime_ICAR_7)){
    if(length(unique(which(is.na(Mu_trt_spacetime_ICAR_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_ICAR_7[i] <- 1
    }
  }
  spacetime_ICAR_conv_7 <- sum(spacetime_ICAR_7,na.rm = T)
  
  spacetime_AR_7 <- c()
  for(i in 1:length(Mu_trt_spacetime_AR_7)){
    if(length(unique(which(is.na(Mu_trt_spacetime_AR_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_AR_7[i] <- 1
    }
  }
  spacetime_AR_conv_7 <- sum(spacetime_AR_7,na.rm = T)
  
  lasso_7 <- c()
  for(i in 1:length(Mu_trt_lasso_7)){
    if(length(unique(which(is.na(Mu_trt_lasso_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      lasso_7[i] <- 1
    }
  }
  lasso_conv_7 <- sum(lasso_7,na.rm = T)
  
  spacetime_shrink <- c()
  for(i in 1:length(Mu_trt_space_comb_time_shrink_7)){
    if(length(unique(which(is.na(Mu_trt_space_comb_time_shrink_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_shrink[i] <- 1
    }
  }
  spacetime_shrink_conv <- sum(spacetime_shrink,na.rm = T)
  
  ## Make Table
  k <- c(rep("k=1",5),rep("k=3",5),rep("k=7",6))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Space-Time \n Shrinkage")
  conv <- c(ori_conv_1,space_conv_1,spacetime_ICAR_conv_1,spacetime_AR_conv_1,lasso_conv_1,
            ori_conv_3,space_conv_3,spacetime_ICAR_conv_3,spacetime_AR_conv_3,lasso_conv_3,
            ori_conv_7,space_conv_7,spacetime_ICAR_conv_7,spacetime_AR_conv_7,lasso_conv_7,spacetime_shrink_conv)/100*100
  
  conv_vec <- cbind(k,Model,conv)
  
  return(conv_vec)
  
}

## Check Convergence of Models
# k in 1,3,7
# our six models
conv_vec <- function(Mu_trt_ori_1,Mu_trt_space_1,Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1,Mu_trt_lasso_1,
                     Mu_trt_ori_3,Mu_trt_space_3,Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3,Mu_trt_lasso_3,
                     Mu_trt_ori_7,Mu_trt_space_7,Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7,Mu_trt_lasso_7,
                     Mu_trt_space_comb_time_shrink_7){
  
  ## k = 1
  ori_1 <- c()
  for(i in 1:length(Mu_trt_ori_1)){
    if(length(unique(which(is.na(Mu_trt_ori_1[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      ori_1[i] <- 1
    }
  }
  ori_conv_1 <- sum(ori_1,na.rm = T)
  
  space_1 <- c()
  for(i in 1:length(Mu_trt_space_1)){
    if(length(unique(which(is.na(Mu_trt_space_1[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      space_1[i] <- 1
    }
  }
  space_conv_1 <- sum(space_1,na.rm = T)
  
  spacetime_ICAR_1 <- c()
  for(i in 1:length(Mu_trt_spacetime_ICAR_1)){
    if(length(unique(which(is.na(Mu_trt_spacetime_ICAR_1[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      spacetime_ICAR_1[i] <- 1
    }
  }
  spacetime_ICAR_conv_1 <- sum(spacetime_ICAR_1,na.rm = T)
  
  spacetime_AR_1 <- c()
  for(i in 1:length(Mu_trt_spacetime_AR_1)){
    if(length(unique(which(is.na(Mu_trt_spacetime_AR_1[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      spacetime_AR_1[i] <- 1
    }
  }
  spacetime_AR_conv_1 <- sum(spacetime_AR_1,na.rm = T)
  
  lasso_1 <- c()
  for(i in 1:length(Mu_trt_lasso_1)){
    if(length(unique(which(is.na(Mu_trt_lasso_1[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      lasso_1[i] <- 1
    }
  }
  lasso_conv_1 <- sum(lasso_1,na.rm = T)
  
  ## k = 3
  ori_3 <- c()
  for(i in 1:length(Mu_trt_ori_3)){
    if(length(unique(which(is.na(Mu_trt_ori_3[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      ori_3[i] <- 1
    }
  }
  ori_conv_3 <- sum(ori_3,na.rm = T)
  
  space_3 <- c()
  for(i in 1:length(Mu_trt_space_3)){
    if(length(unique(which(is.na(Mu_trt_space_3[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      space_3[i] <- 1
    }
  }
  space_conv_3 <- sum(space_3,na.rm = T)
  
  spacetime_ICAR_3 <- c()
  for(i in 1:length(Mu_trt_spacetime_ICAR_3)){
    if(length(unique(which(is.na(Mu_trt_spacetime_ICAR_3[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      spacetime_ICAR_3[i] <- 1
    }
  }
  spacetime_ICAR_conv_3 <- sum(spacetime_ICAR_3,na.rm = T)
  
  spacetime_AR_3 <- c()
  for(i in 1:length(Mu_trt_spacetime_AR_3)){
    if(length(unique(which(is.na(Mu_trt_spacetime_AR_3[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      spacetime_AR_3[i] <- 1
    }
  }
  spacetime_AR_conv_3 <- sum(spacetime_AR_3,na.rm = T)
  
  lasso_3 <- c()
  for(i in 1:length(Mu_trt_lasso_3)){
    if(length(unique(which(is.na(Mu_trt_lasso_3[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      lasso_3[i] <- 1
    }
  }
  lasso_conv_3 <- sum(lasso_3,na.rm = T)
  
  ## k = 7
  ori_7 <- c()
  for(i in 1:length(Mu_trt_ori_7)){
    if(length(unique(which(is.na(Mu_trt_ori_7[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      ori_7[i] <- 1
    }
  }
  ori_conv_7 <- sum(ori_7,na.rm = T)
  
  space_7 <- c()
  for(i in 1:length(Mu_trt_space_7)){
    if(length(unique(which(is.na(Mu_trt_space_7[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      space_7[i] <- 1
    }
  }
  space_conv_7 <- sum(space_7,na.rm = T)
  
  spacetime_ICAR_7 <- c()
  for(i in 1:length(Mu_trt_spacetime_ICAR_7)){
    if(length(unique(which(is.na(Mu_trt_spacetime_ICAR_7[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      spacetime_ICAR_7[i] <- 1
    }
  }
  spacetime_ICAR_conv_7 <- sum(spacetime_ICAR_7,na.rm = T)
  
  spacetime_AR_7 <- c()
  for(i in 1:length(Mu_trt_spacetime_AR_7)){
    if(length(unique(which(is.na(Mu_trt_spacetime_AR_7[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      spacetime_AR_7[i] <- 1
    }
  }
  spacetime_AR_conv_7 <- sum(spacetime_AR_7,na.rm = T)
  
  lasso_7 <- c()
  for(i in 1:length(Mu_trt_lasso_7)){
    if(length(unique(which(is.na(Mu_trt_lasso_7[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      lasso_7[i] <- 1
    }
  }
  lasso_conv_7 <- sum(lasso_7,na.rm = T)
  
  spacetime_shrink <- c()
  for(i in 1:length(Mu_trt_space_comb_time_shrink_7)){
    if(length(unique(which(is.na(Mu_trt_space_comb_time_shrink_7[[i]])*1 == 1, arr.ind = TRUE))) > 0){
      spacetime_shrink[i] <- 1
    }
  }
  spacetime_shrink_conv <- sum(spacetime_shrink,na.rm = T)
  
  ## Make Table
  k <- c(rep("k=1",5),rep("k=3",5),rep("k=7",6))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Space-Time \n Shrinkage")
  conv <- c(ori_conv_1,space_conv_1,spacetime_ICAR_conv_1,spacetime_AR_conv_1,lasso_conv_1,
            ori_conv_3,space_conv_3,spacetime_ICAR_conv_3,spacetime_AR_conv_3,lasso_conv_3,
            ori_conv_7,space_conv_7,spacetime_ICAR_conv_7,spacetime_AR_conv_7,lasso_conv_7,spacetime_shrink_conv)/100*100
  
  conv_vec <- cbind(k,Model,conv)
  
  return(conv_vec)
  
}


conv_table <- function(Mu_trt_ori_1,Mu_trt_space_1,Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1,Mu_trt_lasso_1,
                       Mu_trt_ori_3,Mu_trt_space_3,Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3,Mu_trt_lasso_3,
                       Mu_trt_ori_7,Mu_trt_space_7,Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7,Mu_trt_lasso_7,
                       Mu_trt_space_comb_time_shrink_7){
  ## k = 1
  ori_1 <- c()
  for(i in 1:length(Mu_trt_ori_1)){
    if(length(unique(which(is.na(Mu_trt_ori_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      ori_1[i] <- 1
    }
  }
  ori_conv_1 <- sum(ori_1,na.rm = T)
  
  space_1 <- c()
  for(i in 1:length(Mu_trt_space_1)){
    if(length(unique(which(is.na(Mu_trt_space_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      space_1[i] <- 1
    }
  }
  space_conv_1 <- sum(space_1,na.rm = T)
  
  spacetime_ICAR_1 <- c()
  for(i in 1:length(Mu_trt_spacetime_ICAR_1)){
    if(length(unique(which(is.na(Mu_trt_spacetime_ICAR_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_ICAR_1[i] <- 1
    }
  }
  spacetime_ICAR_conv_1 <- sum(spacetime_ICAR_1,na.rm = T)
  
  spacetime_AR_1 <- c()
  for(i in 1:length(Mu_trt_spacetime_AR_1)){
    if(length(unique(which(is.na(Mu_trt_spacetime_AR_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_AR_1[i] <- 1
    }
  }
  spacetime_AR_conv_1 <- sum(spacetime_AR_1,na.rm = T)
  
  lasso_1 <- c()
  for(i in 1:length(Mu_trt_lasso_1)){
    if(length(unique(which(is.na(Mu_trt_lasso_1[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      lasso_1[i] <- 1
    }
  }
  lasso_conv_1 <- sum(lasso_1,na.rm = T)
  
  ## k = 3
  ori_3 <- c()
  for(i in 1:length(Mu_trt_ori_3)){
    if(length(unique(which(is.na(Mu_trt_ori_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      ori_3[i] <- 1
    }
  }
  ori_conv_3 <- sum(ori_3,na.rm = T)
  
  space_3 <- c()
  for(i in 1:length(Mu_trt_space_3)){
    if(length(unique(which(is.na(Mu_trt_space_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      space_3[i] <- 1
    }
  }
  space_conv_3 <- sum(space_3,na.rm = T)
  
  spacetime_ICAR_3 <- c()
  for(i in 1:length(Mu_trt_spacetime_ICAR_3)){
    if(length(unique(which(is.na(Mu_trt_spacetime_ICAR_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_ICAR_3[i] <- 1
    }
  }
  spacetime_ICAR_conv_3 <- sum(spacetime_ICAR_3,na.rm = T)
  
  spacetime_AR_3 <- c()
  for(i in 1:length(Mu_trt_spacetime_AR_3)){
    if(length(unique(which(is.na(Mu_trt_spacetime_AR_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_AR_3[i] <- 1
    }
  }
  spacetime_AR_conv_3 <- sum(spacetime_AR_3,na.rm = T)
  
  lasso_3 <- c()
  for(i in 1:length(Mu_trt_lasso_3)){
    if(length(unique(which(is.na(Mu_trt_lasso_3[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      lasso_3[i] <- 1
    }
  }
  lasso_conv_3 <- sum(lasso_3,na.rm = T)
  
  ## k = 7
  ori_7 <- c()
  for(i in 1:length(Mu_trt_ori_7)){
    if(length(unique(which(is.na(Mu_trt_ori_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      ori_7[i] <- 1
    }
  }
  ori_conv_7 <- sum(ori_7,na.rm = T)
  
  space_7 <- c()
  for(i in 1:length(Mu_trt_space_7)){
    if(length(unique(which(is.na(Mu_trt_space_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      space_7[i] <- 1
    }
  }
  space_conv_7 <- sum(space_7,na.rm = T)
  
  spacetime_ICAR_7 <- c()
  for(i in 1:length(Mu_trt_spacetime_ICAR_7)){
    if(length(unique(which(is.na(Mu_trt_spacetime_ICAR_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_ICAR_7[i] <- 1
    }
  }
  spacetime_ICAR_conv_7 <- sum(spacetime_ICAR_7,na.rm = T)
  
  spacetime_AR_7 <- c()
  for(i in 1:length(Mu_trt_spacetime_AR_7)){
    if(length(unique(which(is.na(Mu_trt_spacetime_AR_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_AR_7[i] <- 1
    }
    }
  spacetime_AR_conv_7 <- sum(spacetime_AR_7,na.rm = T)
  
  lasso_7 <- c()
  for(i in 1:length(Mu_trt_lasso_7)){
    if(length(unique(which(is.na(Mu_trt_lasso_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      lasso_7[i] <- 1
    }
  }
  lasso_conv_7 <- sum(lasso_7,na.rm = T)
  
  spacetime_shrink <- c()
  for(i in 1:length(Mu_trt_space_comb_time_shrink_7)){
    if(length(unique(which(is.na(Mu_trt_space_comb_time_shrink_7[[i]])*1 == 1, arr.ind = TRUE)[,1])) > 0){
      spacetime_shrink[i] <- 1
    }
  }
  spacetime_shrink_conv <- sum(spacetime_shrink,na.rm = T)
  
  ## Make Table
  k <- c(rep("k=1",6),rep("k=3",6),rep("k=7",7))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Space-Time \n Shrinkage","Shrinkage \n Lasso","Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth")
  df <- cbind(c(ori_conv_1,space_conv_1,spacetime_ICAR_conv_1,spacetime_AR_conv_1,lasso_conv_1,NA),
              c(ori_conv_3,space_conv_3,spacetime_ICAR_conv_3,spacetime_AR_conv_3,lasso_conv_3,NA),
              c(ori_conv_7,space_conv_7,spacetime_ICAR_conv_7,spacetime_AR_conv_7,lasso_conv_7,spacetime_shrink_conv))/100*100
  rownames(df) <- c("Original","Space","Space-Time ICAR","Space-Time AR","Lasso","Space-Time Shrinkage")
  colnames(df) <- c("k=1","k=3","k=7")
  
  # Make Table
  df %>%
    kbl(caption = "Percent NAs for 100 Simulations") %>%
    kable_classic_2(full_width = F)
}

#######################################################################################################

## ------------------------- Percent Bias Y0 (Medians) ------------------------------ ##
PercBiasY0_plot_median <- function(Y0_matrix_1,Mu_trt_ori_1,Mu_trt_space_1, Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1, Mu_trt_space_comb_time_shrink_1, Mu_trt_lasso_1,fit_gsc_1,
                                   Y0_matrix_3,Mu_trt_ori_3,Mu_trt_space_3, Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3, Mu_trt_space_comb_time_shrink_3, Mu_trt_lasso_3,fit_gsc_3,
                                   Y0_matrix_7,Mu_trt_ori_7,Mu_trt_space_7, Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7, Mu_trt_space_comb_time_shrink_7, Mu_trt_lasso_7,fit_gsc_7,
                                   ylim){
  ##ylim is a vector of limits for the axis

  
  # k = 1
  ## original model
  full_rate_list_1_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_ori_1,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  ## space model
  full_rate_list_2_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_space_1,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  
  ## spacetime ICAR model
  full_rate_list_3_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_spacetime_ICAR_1,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_3_3","full_rate_list_3_7"))), envir = parent.frame())
  
  ## spacetime AR model
  full_rate_list_4_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_spacetime_AR_1,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_4_3","full_rate_list_4_7"))), envir = parent.frame())
  
  ## spacetime shrink model
  full_rate_list_5_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_space_comb_time_shrink_1,pop_trt)
  
  ## shrinkage lasso model
  full_rate_list_6_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_lasso_1,pop_trt)
  
  
  # k = 3
  ## original model
  full_rate_list_1_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_ori_3,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  ## space model
  full_rate_list_2_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_space_3,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  
  ## spacetime ICAR model
  full_rate_list_3_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_spacetime_ICAR_3,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_3_3","full_rate_list_3_7"))), envir = parent.frame())
  
  ## spacetime AR model
  full_rate_list_4_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_spacetime_AR_3,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_4_3","full_rate_list_4_7"))), envir = parent.frame())
  
  ## spacetime shrink model
  full_rate_list_5_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_space_comb_time_shrink_3,pop_trt)
  
  ## shrinkage lasso model
  full_rate_list_6_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_lasso_3,pop_trt)
  
  # k = 7
  ## original model
  full_rate_list_1_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_ori_7,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  ## space model
  full_rate_list_2_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_space_7,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  
  ## spacetime ICAR model
  full_rate_list_3_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_spacetime_ICAR_7,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_3_3","full_rate_list_3_7"))), envir = parent.frame())
  
  ## spacetime AR model
  full_rate_list_4_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_spacetime_AR_7,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_4_3","full_rate_list_4_7"))), envir = parent.frame())
  
  ## spacetime shrink model
  full_rate_list_5_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_space_comb_time_shrink_7,pop_trt)
  
  ## shrinkage lasso model
  full_rate_list_6_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_lasso_7,pop_trt)
  
  # k = 1
  ori_1 <- get_perc_bias_Y0(full_rate_list_1_1,Y0_matrix_1,pop_trt,ind)
  
  space_1 <- get_perc_bias_Y0(full_rate_list_2_1,Y0_matrix_1,pop_trt,ind)
  
  spacetime_ICAR_1 <- get_perc_bias_Y0(full_rate_list_3_1,Y0_matrix_1,pop_trt,ind)
  
  spacetime_AR_1 <- get_perc_bias_Y0(full_rate_list_4_1,Y0_matrix_1,pop_trt,ind)
  
  spacetime_shrink_1 <- get_perc_bias_Y0(full_rate_list_5_1,Y0_matrix_1,pop_trt,ind)
  
  lasso_1 <- get_perc_bias_Y0(full_rate_list_6_1,Y0_matrix_1,pop_trt,ind)
  
  gsc_1 <- get_perc_bias_Y0_gsc(pop_trt,fit_gsc_1,Y0_matrix_1,ind)
  
  
  # k = 3
  ori_3 <- get_perc_bias_Y0(full_rate_list_1_3,Y0_matrix_3,pop_trt,ind)
  
  space_3 <- get_perc_bias_Y0(full_rate_list_2_3,Y0_matrix_3,pop_trt,ind)
  
  spacetime_ICAR_3 <- get_perc_bias_Y0(full_rate_list_3_3,Y0_matrix_3,pop_trt,ind)
  
  spacetime_AR_3 <- get_perc_bias_Y0(full_rate_list_4_3,Y0_matrix_3,pop_trt,ind)
  
  spacetime_shrink_3 <- get_perc_bias_Y0(full_rate_list_5_3,Y0_matrix_3,pop_trt,ind)
  
  lasso_3 <- get_perc_bias_Y0(full_rate_list_6_3,Y0_matrix_3,pop_trt,ind)
  
  gsc_3 <- get_perc_bias_Y0_gsc(pop_trt,fit_gsc_3,Y0_matrix_3,ind)
  
  # k = 7
  ori_7 <- get_perc_bias_Y0(full_rate_list_1_7,Y0_matrix_7,pop_trt,ind)
  
  space_7 <- get_perc_bias_Y0(full_rate_list_2_7,Y0_matrix_7,pop_trt,ind)
  
  spacetime_ICAR_7 <- get_perc_bias_Y0(full_rate_list_3_7,Y0_matrix_7,pop_trt,ind)
  
  spacetime_AR_7 <- get_perc_bias_Y0(full_rate_list_4_7,Y0_matrix_7,pop_trt,ind)
  
  spacetime_shrink_7 <- get_perc_bias_Y0(full_rate_list_5_7,Y0_matrix_7,pop_trt,ind)
  
  lasso_7 <- get_perc_bias_Y0(full_rate_list_6_7,Y0_matrix_7,pop_trt,ind)
  
  gsc_7 <- get_perc_bias_Y0_gsc(pop_trt,fit_gsc_7,Y0_matrix_7,ind)
  
  # Barchart of Mean Percent Bias of Y0
  # Get means
  ori_mean_1 <- median(abs(unlist(ori_1$avg_perc_bias)),na.rm = T)
  space_mean_1 <-  median(abs(unlist(space_1$avg_perc_bias)),na.rm = T)
  spacetime_ICAR_mean_1 <-  median(abs(unlist(spacetime_ICAR_1$avg_perc_bias)),na.rm = T)
  spacetime_AR_mean_1 <-  median(abs(unlist(spacetime_AR_1$avg_perc_bias)),na.rm = T)
  spacetime_shrink_mean_1 <-  median(abs(unlist(spacetime_shrink_1$avg_perc_bias)),na.rm = T)
  lasso_mean_1 <-  median(abs(unlist(lasso_1$avg_perc_bias)),na.rm = T)
  gsc_mean_1 <-  median(abs(unlist(gsc_1$avg_perc_bias)),na.rm = T)
  
  ori_mean_3 <- median(abs(unlist(ori_3$avg_perc_bias)),na.rm = T)
  space_mean_3 <-  median(abs(unlist(space_3$avg_perc_bias)),na.rm = T)
  spacetime_ICAR_mean_3 <-  median(abs(unlist(spacetime_ICAR_3$avg_perc_bias)),na.rm = T)
  spacetime_AR_mean_3 <-  median(abs(unlist(spacetime_AR_3$avg_perc_bias)),na.rm = T)
  spacetime_shrink_mean_3 <-  median(abs(unlist(spacetime_shrink_3$avg_perc_bias)),na.rm = T)
  lasso_mean_3 <-  median(abs(unlist(lasso_3$avg_perc_bias)),na.rm = T)
  gsc_mean_3 <-  median(abs(unlist(gsc_3$avg_perc_bias)),na.rm = T)
  
  ori_mean_7 <- median(abs(unlist(ori_7$avg_perc_bias)),na.rm = T)
  space_mean_7 <-  median(abs(unlist(space_7$avg_perc_bias)),na.rm = T)
  spacetime_ICAR_mean_7 <-  median(abs(unlist(spacetime_ICAR_7$avg_perc_bias)),na.rm = T)
  spacetime_AR_mean_7 <-  median(abs(unlist(spacetime_AR_7$avg_perc_bias)),na.rm = T)
  spacetime_shrink_mean_7 <-  median(abs(unlist(spacetime_shrink_7$avg_perc_bias)),na.rm = T)
  lasso_mean_7 <-  median(abs(unlist(lasso_7$avg_perc_bias)),na.rm = T)
  gsc_mean_7 <-  median(abs(unlist(gsc_7$avg_perc_bias)),na.rm = T)
  
  # make df
  #k <- c(rep("k=1",7),rep("k=3",7),rep("k=7",7))
  #Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth", "Space-Time \n Shrinkage",
  #           "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth","Space-Time \n Shrinkage",
  #           "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
  #           "Space-Time \n Shrinkage")
  #Mean <- c(ori_mean_1,space_mean_1,spacetime_ICAR_mean_1,spacetime_AR_mean_1,
  #          lasso_mean_1,gsc_mean_1,spacetime_shrink_mean_1, ori_mean_3,space_mean_3,spacetime_ICAR_mean_3,spacetime_AR_mean_3,
  #          lasso_mean_3,gsc_mean_3,spacetime_shrink_mean_3,ori_mean_7,space_mean_7,spacetime_ICAR_mean_7,spacetime_AR_mean_7,
  #          lasso_mean_7,gsc_mean_7,spacetime_shrink_mean_7  )
  #df <- data.frame(Model,Mean,k) %>%
  #  mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
  #         Model = fct_relevel(Model,
  #                            "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  
  # Below is the way I want to construct the df
  # Right now, I'm using above because I am testing the space-time shrinkage with different k
  k <- c(rep("k=1",6),rep("k=3",6),rep("k=7",7))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Space-Time \n Shrinkage")
  Median <- c(ori_mean_1,space_mean_1,spacetime_ICAR_mean_1,spacetime_AR_mean_1,
            lasso_mean_1,gsc_mean_1, ori_mean_3,space_mean_3,spacetime_ICAR_mean_3,spacetime_AR_mean_3,
            lasso_mean_3,gsc_mean_3,ori_mean_7,space_mean_7,spacetime_ICAR_mean_7,spacetime_AR_mean_7,
            lasso_mean_7,gsc_mean_7,spacetime_shrink_mean_7  )
  df <- data.frame(Model,Median,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
           Model = fct_relevel(Model,
                               "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  
  # Thin ggplot environment
  rm(list=setdiff(ls(), c("df", "ylim")))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Median, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7')) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(Median,2)),position = position_dodge(1), hjust = -.2, size=2.8)+
    #ggtitle("Absolute Percent Bias of the Counterfactuals (Non-Rare)") +
    coord_flip(ylim = ylim ) + 
    scale_x_discrete(limits = rev) +
    #scale_y_continuous(name="Mean", limits=c(0, 90)) +
    theme_bw() + #adds plot outline vs. theme_minimal
    #theme(legend.position="none") #removes legend
    theme(axis.text=element_text(color="black"),
          legend.title = element_blank()) 
  
  
}

#################################################################################################

PercBiasATT_plot_median <- function(Y0_matrix_1,Mu_trt_ori_1,Mu_trt_space_1, Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1, Mu_trt_space_comb_time_shrink_1, Mu_trt_lasso_1,fit_gsc_1,true_trt_effect_1,Y1_matrix_1,
                                    Y0_matrix_3,Mu_trt_ori_3,Mu_trt_space_3, Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3, Mu_trt_space_comb_time_shrink_3, Mu_trt_lasso_3,fit_gsc_3,true_trt_effect_3,Y1_matrix_3,
                                    Y0_matrix_7,Mu_trt_ori_7,Mu_trt_space_7, Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7, Mu_trt_space_comb_time_shrink_7, Mu_trt_lasso_7,fit_gsc_7,true_trt_effect_7,Y1_matrix_7,
                                    ylim, pop_trt){
  Mu_trt_ori_1_ <- list()
  for(i in 1:length(Mu_trt_ori_1)) {
    Mu_trt_ori_1_[[i]] <- Mu_trt_ori_1[[i]]
  }
  
  
  Mu_trt_space_1_ <- list()
  for(i in 1:length(Mu_trt_space_1)) {
    Mu_trt_space_1_[[i]] <- Mu_trt_space_1[[i]]
  }
  
  
  
  Mu_trt_spacetime_ICAR_1_ <- list()
  for(i in 1:length(Mu_trt_spacetime_ICAR_1)) {
    Mu_trt_spacetime_ICAR_1_[[i]] <- Mu_trt_spacetime_ICAR_1[[i]]
  }
  
  
  
  Mu_trt_spacetime_AR_1_ <- list()
  for(i in 1:length(Mu_trt_spacetime_AR_1)) {
    Mu_trt_spacetime_AR_1_[[i]] <- Mu_trt_spacetime_AR_1[[i]]
  }
  
  
  Mu_trt_lasso_1_ <- list()
  for(i in 1:length(Mu_trt_lasso_1)) {
    Mu_trt_lasso_1_[[i]] <- Mu_trt_lasso_1[[i]]
  }
  
  Mu_trt_spacetime_shrink_1_ <- list()
  for(i in 1:length(Mu_trt_space_comb_time_shrink_1)) {
    Mu_trt_spacetime_shrink_1_[[i]] <- Mu_trt_space_comb_time_shrink_1[[i]]
  }
  
  Mu_trt_ori_3_ <- list()
  for(i in 1:length(Mu_trt_ori_3)) {
    Mu_trt_ori_3_[[i]] <- Mu_trt_ori_3[[i]]
  }
  
  
  Mu_trt_space_3_ <- list()
  for(i in 1:length(Mu_trt_space_3)) {
    Mu_trt_space_3_[[i]] <- Mu_trt_space_3[[i]]
  }
  
  
  
  Mu_trt_spacetime_ICAR_3_ <- list()
  for(i in 1:length(Mu_trt_spacetime_ICAR_3)) {
    Mu_trt_spacetime_ICAR_3_[[i]] <- Mu_trt_spacetime_ICAR_3[[i]]
  }
  
  
  
  Mu_trt_spacetime_AR_3_ <- list()
  for(i in 1:length(Mu_trt_spacetime_AR_3)) {
    Mu_trt_spacetime_AR_3_[[i]] <- Mu_trt_spacetime_AR_3[[i]]
  }
  
  
  Mu_trt_spacetime_shrink_3_ <- list()
  for(i in 1:length(Mu_trt_space_comb_time_shrink_3)) {
    Mu_trt_spacetime_shrink_3_[[i]] <- Mu_trt_space_comb_time_shrink_3[[i]]
  }
  
  Mu_trt_lasso_3_ <- list()
  for(i in 1:length(Mu_trt_lasso_3)) {
    Mu_trt_lasso_3_[[i]] <- Mu_trt_lasso_3[[i]]
  }
  
  Mu_trt_ori_7_ <- list()
  for(i in 1:length(Mu_trt_ori_7)) {
    Mu_trt_ori_7_[[i]] <- Mu_trt_ori_7[[i]]
  }
  
  
  Mu_trt_space_7_ <- list()
  for(i in 1:length(Mu_trt_space_7)) {
    Mu_trt_space_7_[[i]] <- Mu_trt_space_7[[i]]
  }
  
  
  
  Mu_trt_spacetime_ICAR_7_ <- list()
  for(i in 1:length(Mu_trt_spacetime_ICAR_7)) {
    Mu_trt_spacetime_ICAR_7_[[i]] <- Mu_trt_spacetime_ICAR_7[[i]]
  }
  
  
  
  Mu_trt_spacetime_AR_7_ <- list()
  for(i in 1:length(Mu_trt_spacetime_AR_7)) {
    Mu_trt_spacetime_AR_7_[[i]] <- Mu_trt_spacetime_AR_7[[i]]
  }
  
  
  Mu_trt_lasso_7_ <- list()
  for(i in 1:length(Mu_trt_lasso_7)) {
    Mu_trt_lasso_7_[[i]] <- Mu_trt_lasso_7[[i]]
  }
  
  Mu_trt_spacetime_shrink_7_ <- list()
  for(i in 1:length(Mu_trt_space_comb_time_shrink_7)) {
    Mu_trt_spacetime_shrink_7_[[i]] <- Mu_trt_space_comb_time_shrink_7[[i]]
  }
  
  ## k = 1
  
  ori_1 <- perc_bias_function_2(Mu_trt_ori_1_,Y1_matrix_1,Y0_matrix_1,pop_trt,true_trt_effect_1, ind)
  
  space_1 <- perc_bias_function_2(Mu_trt_space_1_,Y1_matrix_1,Y0_matrix,pop_trt,true_trt_effect_1, ind)
  
  spacetime_ICAR_1 <- perc_bias_function_2(Mu_trt_spacetime_ICAR_1_,Y1_matrix_1,Y0_matrix,pop_trt,true_trt_effect_1, ind)
  
  spacetime_AR_1 <- perc_bias_function_2(Mu_trt_spacetime_AR_1_,Y1_matrix_1,Y0_matrix,pop_trt,true_trt_effect_1, ind)
  
  lasso_1 <- perc_bias_function_2(Mu_trt_lasso_1_,Y1_matrix_1,Y0_matrix_1,pop_trt,true_trt_effect_1, ind)
  
  gsc_1 <- perc_bias_function_gsc(fit_gsc_1,Y1_matrix_1,Y0_matrix_1,pop_trt,true_trt_effect_1, ind)
  
  spacetime_shrink_1 <- perc_bias_function_2(Mu_trt_spacetime_shrink_1_,Y1_matrix_1,Y0_matrix_1,pop_trt,true_trt_effect_1, ind)
  
  
  ## k = 3
  ori_3 <- perc_bias_function_2(Mu_trt_ori_3_,Y1_matrix_3,Y0_matrix_3,pop_trt,true_trt_effect_3, ind)
  
  space_3 <- perc_bias_function_2(Mu_trt_space_3_,Y1_matrix_3,Y0_matrix,pop_trt,true_trt_effect_3, ind)
  
  spacetime_ICAR_3 <- perc_bias_function_2(Mu_trt_spacetime_ICAR_3_,Y1_matrix_3,Y0_matrix,pop_trt,true_trt_effect_3, ind)
  
  spacetime_AR_3 <- perc_bias_function_2(Mu_trt_spacetime_AR_3_,Y1_matrix_3,Y0_matrix,pop_trt,true_trt_effect_3, ind)
  
  lasso_3 <- perc_bias_function_2(Mu_trt_lasso_3_,Y1_matrix_3,Y0_matrix_3,pop_trt,true_trt_effect_3, ind)
  
  gsc_3 <- perc_bias_function_gsc(fit_gsc_3,Y1_matrix_3,Y0_matrix_3,pop_trt,true_trt_effect_3, ind)
  
  spacetime_shrink_3 <- perc_bias_function_2(Mu_trt_spacetime_shrink_3_,Y1_matrix_3,Y0_matrix_3,pop_trt,true_trt_effect_3, ind)
  
  
  
  ## k = 7
  ori_7 <- perc_bias_function_2(Mu_trt_ori_7_,Y1_matrix_7,Y0_matrix_7,pop_trt,true_trt_effect_7, ind)
  
  space_7 <- perc_bias_function_2(Mu_trt_space_7_,Y1_matrix_7,Y0_matrix,pop_trt,true_trt_effect_7, ind)
  
  spacetime_ICAR_7 <- perc_bias_function_2(Mu_trt_spacetime_ICAR_7_,Y1_matrix_7,Y0_matrix,pop_trt,true_trt_effect_7, ind)
  
  spacetime_AR_7 <- perc_bias_function_2(Mu_trt_spacetime_AR_7_,Y1_matrix_7,Y0_matrix,pop_trt,true_trt_effect_7, ind)
  
  lasso_7 <- perc_bias_function_2(Mu_trt_lasso_7_,Y1_matrix_7,Y0_matrix_7,pop_trt,true_trt_effect_7, ind)
  
  gsc_7 <- perc_bias_function_gsc(fit_gsc_7,Y1_matrix_7,Y0_matrix_7,pop_trt,true_trt_effect_7, ind)
  
  spacetime_shrink_7 <- perc_bias_function_2(Mu_trt_spacetime_shrink_7_,Y1_matrix_7,Y0_matrix_7,pop_trt,true_trt_effect_7, ind)
  
  ## Bar Chart ##
  # Get means
  ori_mean_1 <- median(abs(unlist(ori_1$avg_perc_bias)),na.rm = T)
  space_mean_1 <-  median(abs(unlist(space_1$avg_perc_bias)),na.rm = T)
  spacetime_ICAR_mean_1 <-  median(abs(unlist(spacetime_ICAR_1$avg_perc_bias)),na.rm = T)
  spacetime_AR_mean_1 <-  median(abs(unlist(spacetime_AR_1$avg_perc_bias)),na.rm = T)
  spacetime_shrink_mean_1 <-  median(abs(unlist(spacetime_shrink_1$avg_perc_bias)),na.rm = T)
  lasso_mean_1 <-  median(abs(unlist(lasso_1$avg_perc_bias)),na.rm = T)
  gsc_mean_1 <-  median(abs(unlist(gsc_1$avg_perc_bias)),na.rm = T)
  
  ori_mean_3 <- median(abs(unlist(ori_3$avg_perc_bias)),na.rm = T)
  space_mean_3 <-  median(abs(unlist(space_3$avg_perc_bias)),na.rm = T)
  spacetime_ICAR_mean_3 <-  median(abs(unlist(spacetime_ICAR_3$avg_perc_bias)),na.rm = T)
  spacetime_AR_mean_3 <-  median(abs(unlist(spacetime_AR_3$avg_perc_bias)),na.rm = T)
  spacetime_shrink_mean_3 <-  median(abs(unlist(spacetime_shrink_3$avg_perc_bias)),na.rm = T)
  lasso_mean_3 <-  median(abs(unlist(lasso_3$avg_perc_bias)),na.rm = T)
  gsc_mean_3 <-  median(abs(unlist(gsc_3$avg_perc_bias)),na.rm = T)
  
  ori_mean_7 <- median(abs(unlist(ori_7$avg_perc_bias)),na.rm = T)
  space_mean_7 <-  median(abs(unlist(space_7$avg_perc_bias)),na.rm = T)
  spacetime_ICAR_mean_7 <-  median(abs(unlist(spacetime_ICAR_7$avg_perc_bias)),na.rm = T)
  spacetime_AR_mean_7 <-  median(abs(unlist(spacetime_AR_7$avg_perc_bias)),na.rm = T)
  spacetime_shrink_mean_7 <-  median(abs(unlist(spacetime_shrink_7$avg_perc_bias)),na.rm = T)
  lasso_mean_7 <-  median(abs(unlist(lasso_7$avg_perc_bias)),na.rm = T)
  gsc_mean_7 <-  median(abs(unlist(gsc_7$avg_perc_bias)),na.rm = T)
  
  # make df
  #k <- c(rep("k=1",7),rep("k=3",7),rep("k=7",7))
  #Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth", "Space-Time \n Shrinkage",
  #           "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth","Space-Time \n Shrinkage",
  #           "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
  #           "Space-Time \n Shrinkage")
  #Mean <- c(ori_mean_1,space_mean_1,spacetime_ICAR_mean_1,spacetime_AR_mean_1,
  #          lasso_mean_1,gsc_mean_1,spacetime_shrink_mean_1, ori_mean_3,space_mean_3,spacetime_ICAR_mean_3,spacetime_AR_mean_3,
  #          lasso_mean_3,gsc_mean_3,spacetime_shrink_mean_3,ori_mean_7,space_mean_7,spacetime_ICAR_mean_7,spacetime_AR_mean_7,
  #          lasso_mean_7,gsc_mean_7,spacetime_shrink_mean_7  )
  #df <- data.frame(Model,Mean,k) %>%
  #  mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
  #         Model = fct_relevel(Model,
  #                            "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  
  # Below is the way I want to construct the df
  # Right now, I'm using above because I am testing the space-time shrinkage with different k
  k <- c(rep("k=1",6),rep("k=3",6),rep("k=7",7))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Space-Time \n Shrinkage")
  Median <- c(ori_mean_1,space_mean_1,spacetime_ICAR_mean_1,spacetime_AR_mean_1,
            lasso_mean_1,gsc_mean_1, ori_mean_3,space_mean_3,spacetime_ICAR_mean_3,spacetime_AR_mean_3,
            lasso_mean_3,gsc_mean_3,ori_mean_7,space_mean_7,spacetime_ICAR_mean_7,spacetime_AR_mean_7,
            lasso_mean_7,gsc_mean_7,spacetime_shrink_mean_7  )
  df <- data.frame(Model,Median,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
           Model = fct_relevel(Model,
                               "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  # Thin ggplot environment
  rm(list=setdiff(ls(), c("df", "ylim")))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Median, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7')) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(Median,2)),position = position_dodge(1), hjust = -.2, size=2.8)+
    #ggtitle("Absolute Percent Bias of the Counterfactuals (Non-Rare)") +
    coord_flip(ylim = ylim ) + 
    scale_x_discrete(limits = rev) +
    #scale_y_continuous(name="Mean", limits=c(0, 90000)) +
    theme_bw() + #adds plot outline vs. theme_minimal
    theme(axis.text=element_text(color="black"),
          legend.title = element_blank()
    ) 
}

##############################################################################################
# absolute percent bias ATT mean

PercBiasATT_plot_mean <- function(Y0_matrix_1,Mu_trt_ori_1,Mu_trt_space_1, Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1, Mu_trt_space_comb_time_shrink_1, Mu_trt_lasso_1,fit_gsc_1,true_trt_effect_1,Y1_matrix_1,
                                  Y0_matrix_3,Mu_trt_ori_3,Mu_trt_space_3, Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3, Mu_trt_space_comb_time_shrink_3, Mu_trt_lasso_3,fit_gsc_3,true_trt_effect_3,Y1_matrix_3,
                                  Y0_matrix_7,Mu_trt_ori_7,Mu_trt_space_7, Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7, Mu_trt_space_comb_time_shrink_7, Mu_trt_lasso_7,fit_gsc_7,true_trt_effect_7,Y1_matrix_7,
                                  ylim, pop_trt){
  Mu_trt_ori_1_ <- list()
  for(i in 1:length(Mu_trt_ori_1)) {
    Mu_trt_ori_1_[[i]] <- Mu_trt_ori_1[[i]][[1]]
  }
  
  
  Mu_trt_space_1_ <- list()
  for(i in 1:length(Mu_trt_space_1)) {
    Mu_trt_space_1_[[i]] <- Mu_trt_space_1[[i]][[1]]
  }
  
  
  
  Mu_trt_spacetime_ICAR_1_ <- list()
  for(i in 1:length(Mu_trt_spacetime_ICAR_1)) {
    Mu_trt_spacetime_ICAR_1_[[i]] <- Mu_trt_spacetime_ICAR_1[[i]][[1]]
  }
  
  
  
  Mu_trt_spacetime_AR_1_ <- list()
  for(i in 1:length(Mu_trt_spacetime_AR_1)) {
    Mu_trt_spacetime_AR_1_[[i]] <- Mu_trt_spacetime_AR_1[[i]][[1]]
  }
  
  
  Mu_trt_lasso_1_ <- list()
  for(i in 1:length(Mu_trt_lasso_1)) {
    Mu_trt_lasso_1_[[i]] <- Mu_trt_lasso_1[[i]][[1]]
  }
  
  Mu_trt_spacetime_shrink_1_ <- list()
  for(i in 1:length(Mu_trt_space_comb_time_shrink_1)) {
    Mu_trt_spacetime_shrink_1_[[i]] <- Mu_trt_space_comb_time_shrink_1[[i]][[1]]
  }
  
  Mu_trt_ori_3_ <- list()
  for(i in 1:length(Mu_trt_ori_3)) {
    Mu_trt_ori_3_[[i]] <- Mu_trt_ori_3[[i]][[1]]
  }
  
  
  Mu_trt_space_3_ <- list()
  for(i in 1:length(Mu_trt_space_3)) {
    Mu_trt_space_3_[[i]] <- Mu_trt_space_3[[i]][[1]]
  }
  
  
  
  Mu_trt_spacetime_ICAR_3_ <- list()
  for(i in 1:length(Mu_trt_spacetime_ICAR_3)) {
    Mu_trt_spacetime_ICAR_3_[[i]] <- Mu_trt_spacetime_ICAR_3[[i]][[1]]
  }
  
  
  
  Mu_trt_spacetime_AR_3_ <- list()
  for(i in 1:length(Mu_trt_spacetime_AR_3)) {
    Mu_trt_spacetime_AR_3_[[i]] <- Mu_trt_spacetime_AR_3[[i]][[1]]
  }
  
  
  Mu_trt_spacetime_shrink_3_ <- list()
  for(i in 1:length(Mu_trt_space_comb_time_shrink_3)) {
    Mu_trt_spacetime_shrink_3_[[i]] <- Mu_trt_space_comb_time_shrink_3[[i]][[1]]
  }
  
  Mu_trt_lasso_3_ <- list()
  for(i in 1:length(Mu_trt_lasso_3)) {
    Mu_trt_lasso_3_[[i]] <- Mu_trt_lasso_3[[i]][[1]]
  }
  
  Mu_trt_ori_7_ <- list()
  for(i in 1:length(Mu_trt_ori_7)) {
    Mu_trt_ori_7_[[i]] <- Mu_trt_ori_7[[i]][[1]]
  }
  
  
  Mu_trt_space_7_ <- list()
  for(i in 1:length(Mu_trt_space_7)) {
    Mu_trt_space_7_[[i]] <- Mu_trt_space_7[[i]][[1]]
  }
  
  
  
  Mu_trt_spacetime_ICAR_7_ <- list()
  for(i in 1:length(Mu_trt_spacetime_ICAR_7)) {
    Mu_trt_spacetime_ICAR_7_[[i]] <- Mu_trt_spacetime_ICAR_7[[i]][[1]]
  }
  
  
  
  Mu_trt_spacetime_AR_7_ <- list()
  for(i in 1:length(Mu_trt_spacetime_AR_7)) {
    Mu_trt_spacetime_AR_7_[[i]] <- Mu_trt_spacetime_AR_7[[i]][[1]]
  }
  
  
  Mu_trt_lasso_7_ <- list()
  for(i in 1:length(Mu_trt_lasso_7)) {
    Mu_trt_lasso_7_[[i]] <- Mu_trt_lasso_7[[i]][[1]]
  }
  
  Mu_trt_spacetime_shrink_7_ <- list()
  for(i in 1:length(Mu_trt_space_comb_time_shrink_7)) {
    Mu_trt_spacetime_shrink_7_[[i]] <- Mu_trt_space_comb_time_shrink_7[[i]][[1]]
  }
  
  ## k = 1
  
  ori_1 <- perc_bias_function_2(Mu_trt_ori_1_,Y1_matrix_1,Y0_matrix_1,pop_trt,true_trt_effect_1, ind)
  
  space_1 <- perc_bias_function_2(Mu_trt_space_1_,Y1_matrix_1,Y0_matrix,pop_trt,true_trt_effect_1, ind)
  
  spacetime_ICAR_1 <- perc_bias_function_2(Mu_trt_spacetime_ICAR_1_,Y1_matrix_1,Y0_matrix,pop_trt,true_trt_effect_1, ind)
  
  spacetime_AR_1 <- perc_bias_function_2(Mu_trt_spacetime_AR_1_,Y1_matrix_1,Y0_matrix,pop_trt,true_trt_effect_1, ind)
  
  lasso_1 <- perc_bias_function_2(Mu_trt_lasso_1_,Y1_matrix_1,Y0_matrix_1,pop_trt,true_trt_effect_1, ind)
  
  gsc_1 <- perc_bias_function_gsc(fit_gsc_1,Y1_matrix_1,Y0_matrix_1,pop_trt,true_trt_effect_1, ind)
  
  spacetime_shrink_1 <- perc_bias_function_2(Mu_trt_spacetime_shrink_1_,Y1_matrix_1,Y0_matrix_1,pop_trt,true_trt_effect_1, ind)
  
  
  ## k = 3
  ori_3 <- perc_bias_function_2(Mu_trt_ori_3_,Y1_matrix_3,Y0_matrix_3,pop_trt,true_trt_effect_3, ind)
  
  space_3 <- perc_bias_function_2(Mu_trt_space_3_,Y1_matrix_3,Y0_matrix,pop_trt,true_trt_effect_3, ind)
  
  spacetime_ICAR_3 <- perc_bias_function_2(Mu_trt_spacetime_ICAR_3_,Y1_matrix_3,Y0_matrix,pop_trt,true_trt_effect_3, ind)
  
  spacetime_AR_3 <- perc_bias_function_2(Mu_trt_spacetime_AR_3_,Y1_matrix_3,Y0_matrix,pop_trt,true_trt_effect_3, ind)
  
  lasso_3 <- perc_bias_function_2(Mu_trt_lasso_3_,Y1_matrix_3,Y0_matrix_3,pop_trt,true_trt_effect_3, ind)
  
  gsc_3 <- perc_bias_function_gsc(fit_gsc_3,Y1_matrix_3,Y0_matrix_3,pop_trt,true_trt_effect_3, ind)
  
  spacetime_shrink_3 <- perc_bias_function_2(Mu_trt_spacetime_shrink_3_,Y1_matrix_3,Y0_matrix_3,pop_trt,true_trt_effect_3, ind)
  
  
  
  ## k = 7
  ori_7 <- perc_bias_function_2(Mu_trt_ori_7_,Y1_matrix_7,Y0_matrix_7,pop_trt,true_trt_effect_7, ind)
  
  space_7 <- perc_bias_function_2(Mu_trt_space_7_,Y1_matrix_7,Y0_matrix,pop_trt,true_trt_effect_7, ind)
  
  spacetime_ICAR_7 <- perc_bias_function_2(Mu_trt_spacetime_ICAR_7_,Y1_matrix_7,Y0_matrix,pop_trt,true_trt_effect_7, ind)
  
  spacetime_AR_7 <- perc_bias_function_2(Mu_trt_spacetime_AR_7_,Y1_matrix_7,Y0_matrix,pop_trt,true_trt_effect_7, ind)
  
  lasso_7 <- perc_bias_function_2(Mu_trt_lasso_7_,Y1_matrix_7,Y0_matrix_7,pop_trt,true_trt_effect_7, ind)
  
  gsc_7 <- perc_bias_function_gsc(fit_gsc_7,Y1_matrix_7,Y0_matrix_7,pop_trt,true_trt_effect_7, ind)
  
  spacetime_shrink_7 <- perc_bias_function_2(Mu_trt_spacetime_shrink_7_,Y1_matrix_7,Y0_matrix_7,pop_trt,true_trt_effect_7, ind)
  
  
  
  ## Bar Chart ##
  # Get means
  ori_mean_1 <- mean(abs(unlist(ori_1$mean_perc_bias)),na.rm = T)
  space_mean_1 <-  mean(abs(unlist(space_1$mean_perc_bias)),na.rm = T)
  spacetime_ICAR_mean_1 <-  mean(abs(unlist(spacetime_ICAR_1$mean_perc_bias)),na.rm = T)
  spacetime_AR_mean_1 <-  mean(abs(unlist(spacetime_AR_1$mean_perc_bias)),na.rm = T)
  spacetime_shrink_mean_1 <-  mean(abs(unlist(spacetime_shrink_1$mean_perc_bias)),na.rm = T)
  lasso_mean_1 <-  mean(abs(unlist(lasso_1$mean_perc_bias)),na.rm = T)
  gsc_mean_1 <-  mean(abs(unlist(gsc_1$mean_perc_bias)),na.rm = T)
  
  ori_mean_3 <- mean(abs(unlist(ori_3$mean_perc_bias)),na.rm = T)
  space_mean_3 <-  mean(abs(unlist(space_3$mean_perc_bias)),na.rm = T)
  spacetime_ICAR_mean_3 <-  mean(abs(unlist(spacetime_ICAR_3$mean_perc_bias)),na.rm = T)
  spacetime_AR_mean_3 <-  mean(abs(unlist(spacetime_AR_3$mean_perc_bias)),na.rm = T)
  spacetime_shrink_mean_3 <-  mean(abs(unlist(spacetime_shrink_3$mean_perc_bias)),na.rm = T)
  lasso_mean_3 <-  mean(abs(unlist(lasso_3$mean_perc_bias)),na.rm = T)
  gsc_mean_3 <-  mean(abs(unlist(gsc_3$mean_perc_bias)),na.rm = T)
  
  ori_mean_7 <- mean(abs(unlist(ori_7$mean_perc_bias)),na.rm = T)
  space_mean_7 <-  mean(abs(unlist(space_7$mean_perc_bias)),na.rm = T)
  spacetime_ICAR_mean_7 <-  mean(abs(unlist(spacetime_ICAR_7$mean_perc_bias)),na.rm = T)
  spacetime_AR_mean_7 <-  mean(abs(unlist(spacetime_AR_7$mean_perc_bias)),na.rm = T)
  spacetime_shrink_mean_7 <-  mean(abs(unlist(spacetime_shrink_7$mean_perc_bias)),na.rm = T)
  lasso_mean_7 <-  mean(abs(unlist(lasso_7$mean_perc_bias)),na.rm = T)
  gsc_mean_7 <-  mean(abs(unlist(gsc_7$mean_perc_bias)),na.rm = T)
  
  # make df
  #k <- c(rep("k=1",7),rep("k=3",7),rep("k=7",7))
  #Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth", "Space-Time \n Shrinkage",
  #           "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth","Space-Time \n Shrinkage",
  #           "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
  #           "Space-Time \n Shrinkage")
  #Mean <- c(ori_mean_1,space_mean_1,spacetime_ICAR_mean_1,spacetime_AR_mean_1,
  #          lasso_mean_1,gsc_mean_1,spacetime_shrink_mean_1, ori_mean_3,space_mean_3,spacetime_ICAR_mean_3,spacetime_AR_mean_3,
  #          lasso_mean_3,gsc_mean_3,spacetime_shrink_mean_3,ori_mean_7,space_mean_7,spacetime_ICAR_mean_7,spacetime_AR_mean_7,
  #          lasso_mean_7,gsc_mean_7,spacetime_shrink_mean_7  )
  #df <- data.frame(Model,Mean,k) %>%
  #  mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
  #         Model = fct_relevel(Model,
  #                            "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  
  # Below is the way I want to construct the df
  # Right now, I'm using above because I am testing the space-time shrinkage with different k
  k <- c(rep("k=1",6),rep("k=3",6),rep("k=7",7))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Space-Time \n Shrinkage")
  Mean <- c(ori_mean_1,space_mean_1,spacetime_ICAR_mean_1,spacetime_AR_mean_1,
            lasso_mean_1,gsc_mean_1, ori_mean_3,space_mean_3,spacetime_ICAR_mean_3,spacetime_AR_mean_3,
            lasso_mean_3,gsc_mean_3,ori_mean_7,space_mean_7,spacetime_ICAR_mean_7,spacetime_AR_mean_7,
            lasso_mean_7,gsc_mean_7,spacetime_shrink_mean_7  )
  df <- data.frame(Model,Mean,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
           Model = fct_relevel(Model,
                               "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  # Thin ggplot environment
  rm(list=setdiff(ls(), c("df", "ylim")))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Mean, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7')) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(Mean,2)),position = position_dodge(1), hjust = -.2, size=2.8)+
    #ggtitle("Absolute Percent Bias of the Counterfactuals (Non-Rare)") +
    coord_flip(ylim = ylim ) + 
    scale_x_discrete(limits = rev) +
    #scale_y_continuous(name="Mean", limits=c(0, 90000)) +
    theme_bw() + #adds plot outline vs. theme_minimal
    theme(axis.text=element_text(color="black"),
          legend.title = element_blank()
    ) 
}
##############################################################################################

## Plot mean raw bias per 100,000

RawBiasRate_plot_mean <- function(Y0_matrix_1,Mu_trt_ori_1,Mu_trt_space_1, Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1, Mu_trt_space_comb_time_shrink_1, Mu_trt_lasso_1,fit_gsc_1,
                                  Y0_matrix_3,Mu_trt_ori_3,Mu_trt_space_3, Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3, Mu_trt_space_comb_time_shrink_3, Mu_trt_lasso_3,fit_gsc_3,
                                  Y0_matrix_7,Mu_trt_ori_7,Mu_trt_space_7, Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7, Mu_trt_space_comb_time_shrink_7, Mu_trt_lasso_7,fit_gsc_7,
                                  ylim,ind){
  # k = 1
  ## original model
  full_rate_list_1_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_ori_1,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  ## space model
  full_rate_list_2_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_space_1,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  
  ## spacetime ICAR model
  full_rate_list_3_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_spacetime_ICAR_1,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_3_3","full_rate_list_3_7"))), envir = parent.frame())
  
  ## spacetime AR model
  full_rate_list_4_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_spacetime_AR_1,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_4_3","full_rate_list_4_7"))), envir = parent.frame())
  
  ## spacetime shrink model
  full_rate_list_5_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_space_comb_time_shrink_1,pop_trt)
  
  ## shrinkage lasso model
  full_rate_list_6_1 <- diff_post_obs_rate(Y0_matrix_1,Mu_trt_lasso_1,pop_trt)
  
  
  # k = 3
  ## original model
  full_rate_list_1_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_ori_3,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  ## space model
  full_rate_list_2_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_space_3,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  
  ## spacetime ICAR model
  full_rate_list_3_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_spacetime_ICAR_3,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_3_3","full_rate_list_3_7"))), envir = parent.frame())
  
  ## spacetime AR model
  full_rate_list_4_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_spacetime_AR_3,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_4_3","full_rate_list_4_7"))), envir = parent.frame())
  
  ## spacetime shrink model
  full_rate_list_5_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_space_comb_time_shrink_3,pop_trt)
  
  ## shrinkage lasso model
  full_rate_list_6_3 <- diff_post_obs_rate(Y0_matrix_3,Mu_trt_lasso_3,pop_trt)
  
  # k = 7
  ## original model
  full_rate_list_1_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_ori_7,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  ## space model
  full_rate_list_2_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_space_7,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_1_3","full_rate_list_1_7"))), envir = parent.frame()) 
  
  
  ## spacetime ICAR model
  full_rate_list_3_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_spacetime_ICAR_7,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_3_3","full_rate_list_3_7"))), envir = parent.frame())
  
  ## spacetime AR model
  full_rate_list_4_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_spacetime_AR_7,pop_trt)
  #list2env(setNames(res_list,c(list("full_rate_list_4_3","full_rate_list_4_7"))), envir = parent.frame())
  
  ## spacetime shrink model
  full_rate_list_5_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_space_comb_time_shrink_7,pop_trt)
  
  ## shrinkage lasso model
  full_rate_list_6_7 <- diff_post_obs_rate(Y0_matrix_7,Mu_trt_lasso_7,pop_trt)
  
  
  ## k = 1
  ori_1 <- mean(get_mean_raw_bias(full_rate_list_1_1)[ind],na.rm = T)
  
  space_1 <- mean(get_mean_raw_bias(full_rate_list_2_1)[ind],na.rm = T)
  
  spacetime_ICAR_1 <- mean(get_mean_raw_bias(full_rate_list_3_1)[ind],na.rm = T)
  
  spacetime_AR_1 <- mean(get_mean_raw_bias(full_rate_list_4_1)[ind],na.rm = T)
  
  spacetime_shrink_1 <- mean(get_mean_raw_bias(full_rate_list_5_1)[ind],na.rm = T)
  
  lasso_1 <- mean(get_mean_raw_bias(full_rate_list_6_1)[ind],na.rm = T)
  
  gsc_1 <- mean(get_mean_raw_bias_gsc(pop_trt,fit_gsc_1,Y0_matrix_1)[ind],na.rm = T)
  
  ## k = 3
  ori_3 <- mean(get_mean_raw_bias(full_rate_list_1_3)[ind],na.rm = T)
  
  space_3 <- mean(get_mean_raw_bias(full_rate_list_2_3)[ind],na.rm = T)
  
  spacetime_ICAR_3 <- mean(get_mean_raw_bias(full_rate_list_3_3)[ind],na.rm = T)
  
  spacetime_AR_3 <- mean(get_mean_raw_bias(full_rate_list_4_3)[ind],na.rm = T)
  
  spacetime_shrink_3 <- mean(get_mean_raw_bias(full_rate_list_5_3)[ind],na.rm = T)
  
  lasso_3 <- mean(get_mean_raw_bias(full_rate_list_6_3)[ind],na.rm = T)
  
  gsc_3 <- mean(get_mean_raw_bias_gsc(pop_trt,fit_gsc_3,Y0_matrix_3)[ind],na.rm = T)
  
  ## k = 7
  ori_7 <- mean(get_mean_raw_bias(full_rate_list_1_7)[ind],na.rm = T)
  
  space_7 <- mean(get_mean_raw_bias(full_rate_list_2_7)[ind],na.rm = T)
  
  spacetime_ICAR_7 <- mean(get_mean_raw_bias(full_rate_list_3_7)[ind],na.rm = T)
  
  spacetime_AR_7<- mean(get_mean_raw_bias(full_rate_list_4_7)[ind],na.rm = T)
  
  spacetime_shrink_7 <- mean(get_mean_raw_bias(full_rate_list_5_7)[ind],na.rm = T)
  
  lasso_7 <- mean(get_mean_raw_bias(full_rate_list_6_7)[ind],na.rm = T)
  
  gsc_7 <- mean(get_mean_raw_bias_gsc(pop_trt,fit_gsc_7,Y0_matrix_7)[ind],na.rm = T)
  
  k <- c(rep("k=1",6),rep("k=3",6),rep("k=7",7))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Space-Time \n Shrinkage")
  Mean <- c(ori_1,space_1,spacetime_ICAR_1,spacetime_AR_1,
            lasso_1,gsc_1, ori_3,space_3,spacetime_ICAR_3,spacetime_AR_3,
            lasso_3,gsc_3,ori_7,space_7,spacetime_ICAR_7,spacetime_AR_7,
            lasso_7,gsc_7,spacetime_shrink_7  )
  df <- data.frame(Model,Mean,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
           Model = fct_relevel(Model,
                               "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  # Thin ggplot environment
  rm(list=setdiff(ls(), c("df", "ylim")))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Mean, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7')) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(Mean,2)),position = position_dodge(1), hjust = -.2, size=2.8)+
    #ggtitle("Absolute Percent Bias of the Counterfactuals (Non-Rare)") +
    coord_flip(ylim = ylim ) + 
    scale_x_discrete(limits = rev) +
    #scale_y_continuous(name="Mean", limits=c(0, 90000)) +
    theme_bw() + #adds plot outline vs. theme_minimal
    theme(axis.text=element_text(color="black"),
          legend.title = element_blank()
    ) 
}

#######################################################################################################

## Function for my models
meanAbsPercBias_county <- function(Mu_trt,Y0_matrix){
  meanAbsPercBias <- c()
  for(i in 1:length(Mu_trt)){
    # Find true county mean rate
    Y0_rate <- (Y0_matrix[ind,i]/pop_trt[ind]) * 100000 # Find Rate
    true_Y0_mean_rate <- rowMeans(matrix(Y0_rate , ncol = 7, byrow = T)) # take mean across counties
    
    # For each simulation, get median of 1000 MCMC samples. Then only extract the treated indices.
    med_vec <- colMedians(Mu_trt[[i]][[1]],na.rm = T)[ind]
    
    # Convert to rate
    med_vec_rate <- (med_vec/pop_trt[ind]) * 100000
    
    # For each simulation, take average counterfactual for each county over all time points 
    # should be a vector of length 6
    county_means <- rowMeans(matrix(med_vec_rate, ncol = 7, byrow = T))
    
    absPercBias <- abs(((true_Y0_mean_rate - county_means)/true_Y0_mean_rate) * 100)
    
    meanAbsPercBias[i] <- mean(absPercBias,na.rm = T)
  }
  return(mean(meanAbsPercBias))
}


## Function for GSC
meanAbsPercBias_county_gsc <- function(fit_gsc,Y0_matrix){
  meanAbsPercBias <- c()
  for(i in 1:length(fit_gsc)){
    # Find true county mean rate
    Y0_rate <- (Y0_matrix[ind,i]/pop_trt[ind]) * 100000 # Find Rate
    true_Y0_mean_rate <- rowMeans(matrix(Y0_rate , ncol = 7, byrow = T)) # take mean across counties
    
    # For each simulation, get median of 1000 MCMC samples. Then only extract the treated indices.
    # Convert to rate
    med_vec_rate <- as.vector(fit_gsc[[i]]$fit_gsc$Y.ct)[ind]
    
    # For each simulation, take average counterfactual for each county over all time points 
    # should be a vector of length 6
    county_means <- rowMeans(matrix(med_vec_rate, ncol = 7, byrow = T))
    
    absPercBias <- abs(((true_Y0_mean_rate - county_means)/true_Y0_mean_rate) * 100)
    
    meanAbsPercBias[i] <- mean(absPercBias)
  }
  return(mean(meanAbsPercBias))
}

#######################################################################################################

## Function for my models
medianAbsPercBias_county <- function(Mu_trt,Y0_matrix){
  medianAbsPercBias <- c()
  for(i in 1:length(Mu_trt)){
    # Find true county median rate
    Y0_rate <- (Y0_matrix[ind,i]/pop_trt[ind]) * 100000 # Find Rate
    true_Y0_median_rate <- rowMedians(matrix(Y0_rate , ncol = 7, byrow = T)) # take median across counties
    
    # For each simulation, get median of 1000 MCMC samples. Then only extract the treated indices.
    med_vec <- colMedians(Mu_trt[[i]][[1]],na.rm = T)[ind]
    
    # Convert to rate
    med_vec_rate <- (med_vec/pop_trt[ind]) * 100000
    
    # For each simulation, take average counterfactual for each county over all time points 
    # should be a vector of length 6
    county_medians <- rowMedians(matrix(med_vec_rate, ncol = 7, byrow = T))
    
    absPercBias <- abs(((true_Y0_median_rate - county_medians)/true_Y0_median_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  return(median(medianAbsPercBias))
}


## Function for GSC
medianAbsPercBias_county_gsc <- function(fit_gsc,Y0_matrix){
  medianAbsPercBias <- c()
  for(i in 1:length(fit_gsc)){
    # Find true county median rate
    Y0_rate <- (Y0_matrix[ind,i]/pop_trt[ind]) * 100000 # Find Rate
    true_Y0_median_rate <- rowMedians(matrix(Y0_rate , ncol = 7, byrow = T)) # take median across counties
    
    # For each simulation, get median of 1000 MCMC samples. Then only extract the treated indices.
    # Convert to rate
    med_vec_rate <- as.vector(fit_gsc[[i]]$fit_gsc$Y.ct)[ind]
    
    # For each simulation, take average counterfactual for each county over all time points 
    # should be a vector of length 6
    county_medians <- rowMedians(matrix(med_vec_rate, ncol = 7, byrow = T))
    
    absPercBias <- abs(((true_Y0_median_rate - county_medians)/true_Y0_median_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias)
  }
  return(median(medianAbsPercBias))
}
###############################################################################################
meanAbsPercBias_county_plot <- function(Y0_matrix_1,Mu_trt_ori_1,Mu_trt_space_1, Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1, Mu_trt_space_comb_time_shrink_1, Mu_trt_lasso_1,fit_gsc_1,
                                        Y0_matrix_3,Mu_trt_ori_3,Mu_trt_space_3, Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3, Mu_trt_space_comb_time_shrink_3, Mu_trt_lasso_3,fit_gsc_3,
                                        Y0_matrix_7,Mu_trt_ori_7,Mu_trt_space_7, Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7, Mu_trt_space_comb_time_shrink_7, Mu_trt_lasso_7,fit_gsc_7,
                                        ylim,ind){
  ## k = 1
  ori_mean_1 <- meanAbsPercBias_county(Mu_trt_ori_1,Y0_matrix_1)
  space_mean_1 <- meanAbsPercBias_county(Mu_trt_space_1,Y0_matrix_1)
  spacetime_ICAR_mean_1 <- meanAbsPercBias_county(Mu_trt_spacetime_ICAR_1,Y0_matrix_1)
  spacetime_AR_mean_1 <- meanAbsPercBias_county(Mu_trt_spacetime_AR_1,Y0_matrix_1)
  spacetime_shrink_mean_1 <- meanAbsPercBias_county(Mu_trt_space_comb_time_shrink_1,Y0_matrix_1)
  lasso_mean_1 <- meanAbsPercBias_county(Mu_trt_lasso_1,Y0_matrix_1)
  gsc_mean_1 <- meanAbsPercBias_county_gsc(fit_gsc_1,Y0_matrix_1)
  
  ## k = 3
  ori_mean_3 <- meanAbsPercBias_county(Mu_trt_ori_3,Y0_matrix_3)
  space_mean_3 <- meanAbsPercBias_county(Mu_trt_space_3,Y0_matrix_3)
  spacetime_ICAR_mean_3 <- meanAbsPercBias_county(Mu_trt_spacetime_ICAR_3,Y0_matrix_3)
  spacetime_AR_mean_3 <- meanAbsPercBias_county(Mu_trt_spacetime_AR_3,Y0_matrix_3)
  spacetime_shrink_mean_3 <- meanAbsPercBias_county(Mu_trt_space_comb_time_shrink_3,Y0_matrix_3)
  lasso_mean_3 <- meanAbsPercBias_county(Mu_trt_lasso_3,Y0_matrix_3)
  gsc_mean_3 <- meanAbsPercBias_county_gsc(fit_gsc_3,Y0_matrix_3)
  
  ## k = 7
  ori_mean_7 <- meanAbsPercBias_county(Mu_trt_ori_7,Y0_matrix_7)
  space_mean_7 <- meanAbsPercBias_county(Mu_trt_space_7,Y0_matrix_7)
  spacetime_ICAR_mean_7 <- meanAbsPercBias_county(Mu_trt_spacetime_ICAR_7,Y0_matrix_7)
  spacetime_AR_mean_7 <- meanAbsPercBias_county(Mu_trt_spacetime_AR_7,Y0_matrix_7)
  spacetime_shrink_mean_7 <- meanAbsPercBias_county(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7)
  lasso_mean_7 <- meanAbsPercBias_county(Mu_trt_lasso_7,Y0_matrix_7)
  gsc_mean_7 <- meanAbsPercBias_county_gsc(fit_gsc_7,Y0_matrix_7)
  
  k <- c(rep("k=1",6),rep("k=3",6),rep("k=7",7))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Space-Time \n Shrinkage")
  Mean <- c(ori_mean_1,space_mean_1,spacetime_ICAR_mean_1,spacetime_AR_mean_1,
            lasso_mean_1,gsc_mean_1, ori_mean_3,space_mean_3,spacetime_ICAR_mean_3,spacetime_AR_mean_3,
            lasso_mean_3,gsc_mean_3,ori_mean_7,space_mean_7,spacetime_ICAR_mean_7,spacetime_AR_mean_7,
            lasso_mean_7,gsc_mean_7,spacetime_shrink_mean_7  )
  df <- data.frame(Model,Mean,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
           Model = fct_relevel(Model,
                             "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  # Thin ggplot environment
  rm(list=setdiff(ls(), c("df", "ylim")))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Mean, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7')) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(Mean,2)),position = position_dodge(1), hjust = -.2, size=2.8)+
    #ggtitle("Absolute Percent Bias of the Counterfactuals (Non-Rare)") +
    coord_flip(ylim = ylim ) + 
    scale_x_discrete(limits = rev) +
    #scale_y_continuous(name="Mean", limits=c(0, 90000)) +
    theme_bw() + #adds plot outline vs. theme_minimal
    theme(axis.text=element_text(color="black"),
          legend.title = element_blank()
    ) 
}

###############################################################################################
medianAbsPercBias_county_plot <- function(Y0_matrix_1,Mu_trt_ori_1,Mu_trt_space_1, Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1, Mu_trt_space_comb_time_shrink_1, Mu_trt_lasso_1,fit_gsc_1,
                                        Y0_matrix_3,Mu_trt_ori_3,Mu_trt_space_3, Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3, Mu_trt_space_comb_time_shrink_3, Mu_trt_lasso_3,fit_gsc_3,
                                        Y0_matrix_7,Mu_trt_ori_7,Mu_trt_space_7, Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7, Mu_trt_space_comb_time_shrink_7, Mu_trt_lasso_7,fit_gsc_7,
                                        ylim,ind){
  ## k = 1
  ori_median_1 <- medianAbsPercBias_county(Mu_trt_ori_1,Y0_matrix_1)
  space_median_1 <- medianAbsPercBias_county(Mu_trt_space_1,Y0_matrix_1)
  spacetime_ICAR_median_1 <- medianAbsPercBias_county(Mu_trt_spacetime_ICAR_1,Y0_matrix_1)
  spacetime_AR_median_1 <- medianAbsPercBias_county(Mu_trt_spacetime_AR_1,Y0_matrix_1)
  spacetime_shrink_median_1 <- medianAbsPercBias_county(Mu_trt_space_comb_time_shrink_1,Y0_matrix_1)
  lasso_median_1 <- medianAbsPercBias_county(Mu_trt_lasso_1,Y0_matrix_1)
  gsc_median_1 <- medianAbsPercBias_county_gsc(fit_gsc_1,Y0_matrix_1)
  
  ## k = 3
  ori_median_3 <- medianAbsPercBias_county(Mu_trt_ori_3,Y0_matrix_3)
  space_median_3 <- medianAbsPercBias_county(Mu_trt_space_3,Y0_matrix_3)
  spacetime_ICAR_median_3 <- medianAbsPercBias_county(Mu_trt_spacetime_ICAR_3,Y0_matrix_3)
  spacetime_AR_median_3 <- medianAbsPercBias_county(Mu_trt_spacetime_AR_3,Y0_matrix_3)
  spacetime_shrink_median_3 <- medianAbsPercBias_county(Mu_trt_space_comb_time_shrink_3,Y0_matrix_3)
  lasso_median_3 <- medianAbsPercBias_county(Mu_trt_lasso_3,Y0_matrix_3)
  gsc_median_3 <- medianAbsPercBias_county_gsc(fit_gsc_3,Y0_matrix_3)
  
  ## k = 7
  ori_median_7 <- medianAbsPercBias_county(Mu_trt_ori_7,Y0_matrix_7)
  space_median_7 <- medianAbsPercBias_county(Mu_trt_space_7,Y0_matrix_7)
  spacetime_ICAR_median_7 <- medianAbsPercBias_county(Mu_trt_spacetime_ICAR_7,Y0_matrix_7)
  spacetime_AR_median_7 <- medianAbsPercBias_county(Mu_trt_spacetime_AR_7,Y0_matrix_7)
  spacetime_shrink_median_7 <- medianAbsPercBias_county(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7)
  lasso_median_7 <- medianAbsPercBias_county(Mu_trt_lasso_7,Y0_matrix_7)
  gsc_median_7 <- medianAbsPercBias_county_gsc(fit_gsc_7,Y0_matrix_7)
  
  k <- c(rep("k=1",6),rep("k=3",6),rep("k=7",7))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Space-Time \n Shrinkage")
  Median <- c(ori_median_1,space_median_1,spacetime_ICAR_median_1,spacetime_AR_median_1,
            lasso_median_1,gsc_median_1, ori_median_3,space_median_3,spacetime_ICAR_median_3,spacetime_AR_median_3,
            lasso_median_3,gsc_median_3,ori_median_7,space_median_7,spacetime_ICAR_median_7,spacetime_AR_median_7,
            lasso_median_7,gsc_median_7,spacetime_shrink_median_7  )
  df <- data.frame(Model,Median,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
           Model = fct_relevel(Model,
                               "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  
  # Thin ggplot environment
  rm(list=setdiff(ls(), c("df", "ylim")))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Median, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7')) +
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



#######################################################################################################
#ATT
## Function for my models
medianAbsPercBiasATT_county <- function(Mu_trt,Y0_matrix,Y1_matrix, true_trt_effect){
  medianAbsPercBias <- c()
  for(i in 1:length(Mu_trt)){
    # Find true treatment effect rate
    true_trt_effect_rate <- (true_trt_effect/pop_trt[ind] * 100000)[,i]# Find Rate
    true_trt_effect_median_rate <- rowMedians(matrix(true_trt_effect_rate , ncol = 7, byrow = T)) # take median across counties
    
    # For each simulation, get median of 1000 MCMC samples. Then only extract the treated indices.
    med_vec <- colMedians(Mu_trt[[i]][[1]],na.rm = T)[ind]
    
    # Convert to rate
    med_vec_rate <- (med_vec/pop_trt[ind]) * 100000
    
    # For each simulation, take average counterfactual for each county over all time points 
    # should be a vector of length 6
    county_medians_Y0 <- rowMedians(matrix(med_vec_rate, ncol = 7, byrow = T))
    
    # For each simulation, take average Y1 for each county over all time points 
    # should be a vector of length 6
    Y1_matrix_rate <- (Y1_matrix[,i]/pop_trt[ind]) * 100000
    county_medians_Y1 <- rowMedians(matrix(Y1_matrix_rate, ncol = 7, byrow = T))
    
    # Get median estimated treatment effect rate
    est_trt_effect <- county_medians_Y1 - county_medians_Y0
  
    # find absolute percent bias ATT 
    absPercBias <- abs(((est_trt_effect  - true_trt_effect_median_rate)/true_trt_effect_median_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  return(median(medianAbsPercBias))
}


## Function for GSC
medianAbsPercBiasATT_county_gsc <- function(fit_gsc,Y0_matrix,Y1_matrix,true_trt_effect){
  medianAbsPercBias <- c()
  for(i in 1:length(fit_gsc)){
    # Find true treatment effect rate
    true_trt_effect_rate <- (true_trt_effect/pop_trt[ind] * 100000)[,i]# Find Rate
    true_trt_effect_median_rate <- rowMedians(matrix(true_trt_effect_rate , ncol = 7, byrow = T)) # take median across counties
    
    # For each simulation, get median of 1000 MCMC samples. Then only extract the treated indices.
    # Convert to rate
    med_vec_rate <- as.vector(fit_gsc[[i]]$fit_gsc$Y.ct)[ind]
    
    # For each simulation, take average counterfactual for each county over all time points 
    # should be a vector of length 6
    county_medians_Y0 <- rowMedians(matrix(med_vec_rate, ncol = 7, byrow = T))
    
    # For each simulation, take average Y1 for each county over all time points 
    # should be a vector of length 6
    Y1_matrix_rate <- (Y1_matrix[,i]/pop_trt[ind]) * 100000
    county_medians_Y1 <- rowMedians(matrix(Y1_matrix_rate, ncol = 7, byrow = T))
    
    # Get median estimated treatment effect rate
    est_trt_effect <- county_medians_Y1 - county_medians_Y0
    
    # find absolute percent bias ATT 
    absPercBias <- abs(((est_trt_effect  - true_trt_effect_median_rate)/true_trt_effect_median_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  return(median(medianAbsPercBias))
}


###############################################################################################
medianAbsPercBiasATT_county_plot <- function(Y0_matrix_1,Mu_trt_ori_1,Mu_trt_space_1, Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1, Mu_trt_space_comb_time_shrink_1, Mu_trt_lasso_1,fit_gsc_1,true_trt_effect_1,Y1_matrix_1,
                                             Y0_matrix_3,Mu_trt_ori_3,Mu_trt_space_3, Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3, Mu_trt_space_comb_time_shrink_3, Mu_trt_lasso_3,fit_gsc_3,true_trt_effect_3,Y1_matrix_3,
                                             Y0_matrix_7,Mu_trt_ori_7,Mu_trt_space_7, Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7, Mu_trt_space_comb_time_shrink_7, Mu_trt_lasso_7,fit_gsc_7,true_trt_effect_7,Y1_matrix_7,
                                             ylim = c(0,90000), pop_trt, ind){
  ## k = 1
  ori_median_1 <- medianAbsPercBiasATT_county(Mu_trt_ori_1,Y0_matrix_1, Y1_matrix_1,true_trt_effect_1)
  space_median_1 <- medianAbsPercBiasATT_county(Mu_trt_space_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  spacetime_ICAR_median_1 <- medianAbsPercBiasATT_county(Mu_trt_spacetime_ICAR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  spacetime_AR_median_1 <- medianAbsPercBiasATT_county(Mu_trt_spacetime_AR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  spacetime_shrink_median_1 <- medianAbsPercBiasATT_county(Mu_trt_space_comb_time_shrink_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  lasso_median_1 <- medianAbsPercBiasATT_county(Mu_trt_lasso_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  gsc_median_1 <- medianAbsPercBiasATT_county_gsc(fit_gsc_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  
  ## k = 3
  ori_median_3 <- medianAbsPercBiasATT_county(Mu_trt_ori_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  space_median_3 <- medianAbsPercBiasATT_county(Mu_trt_space_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  spacetime_ICAR_median_3 <- medianAbsPercBiasATT_county(Mu_trt_spacetime_ICAR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  spacetime_AR_median_3 <- medianAbsPercBiasATT_county(Mu_trt_spacetime_AR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  spacetime_shrink_median_3 <- medianAbsPercBiasATT_county(Mu_trt_space_comb_time_shrink_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  lasso_median_3 <- medianAbsPercBiasATT_county(Mu_trt_lasso_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  gsc_median_3 <- medianAbsPercBiasATT_county_gsc(fit_gsc_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  
  ## k = 7
  ori_median_7 <- medianAbsPercBiasATT_county(Mu_trt_ori_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  space_median_7 <- medianAbsPercBiasATT_county(Mu_trt_space_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  spacetime_ICAR_median_7 <- medianAbsPercBiasATT_county(Mu_trt_spacetime_ICAR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  spacetime_AR_median_7 <- medianAbsPercBiasATT_county(Mu_trt_spacetime_AR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  spacetime_shrink_median_7 <- medianAbsPercBiasATT_county(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  lasso_median_7 <- medianAbsPercBiasATT_county(Mu_trt_lasso_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  gsc_median_7 <- medianAbsPercBiasATT_county_gsc(fit_gsc_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  
  
  k <- c(rep("k=1",6),rep("k=3",6),rep("k=7",7))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Space-Time \n Shrinkage")
  Median <- c(ori_median_1,space_median_1,spacetime_ICAR_median_1,spacetime_AR_median_1,
              lasso_median_1,gsc_median_1, ori_median_3,space_median_3,spacetime_ICAR_median_3,spacetime_AR_median_3,
              lasso_median_3,gsc_median_3,ori_median_7,space_median_7,spacetime_ICAR_median_7,spacetime_AR_median_7,
              lasso_median_7,gsc_median_7,spacetime_shrink_median_7  )
  df <- data.frame(Model,Median,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
           Model = fct_relevel(Model,
                               "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  # Thin ggplot environment
  rm(list=setdiff(ls(), c("df", "ylim")))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Median, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7')) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(Median,2)),position = position_dodge(1), hjust = -.2, size=2.8)+
    #ggtitle("Absolute Percent Bias of the ATT (Non-Rare, Rate = .5)") +
    coord_flip(ylim = ylim ) + 
    scale_x_discrete(limits = rev) +
    #scale_y_continuous(name="Median", limits=c(0, 90000)) +
    theme_bw() + #adds plot outline vs. theme_minimal
    theme(axis.text=element_text(color="black"),
          legend.title = element_blank()
    ) 
}

#######################################################################################################
#ATT
## Function for my models
meanAbsPercBiasATT_county <- function(Mu_trt,Y0_matrix,Y1_matrix, true_trt_effect){
  meanAbsPercBias <- c()
  for(i in 1:length(Mu_trt)){
    # Find true treatment effect rate
    true_trt_effect_rate <- (true_trt_effect/pop_trt[ind] * 100000)[,i]# Find Rate
    true_trt_effect_mean_rate <- rowMeans(matrix(true_trt_effect_rate , ncol = 7, byrow = T)) # take mean across counties
    
    # For each simulation, get mean of 1000 MCMC samples. Then only extract the treated indices.
    med_vec <- colMeans(Mu_trt[[i]][[1]],na.rm = T)[ind]
    
    # Convert to rate
    med_vec_rate <- (med_vec/pop_trt[ind]) * 100000
    
    # For each simulation, take average counterfactual for each county over all time points 
    # should be a vector of length 6
    county_means_Y0 <- rowMeans(matrix(med_vec_rate, ncol = 7, byrow = T))
    
    # For each simulation, take average Y1 for each county over all time points 
    # should be a vector of length 6
    Y1_matrix_rate <- (Y1_matrix[,i]/pop_trt[ind]) * 100000
    county_means_Y1 <- rowMeans(matrix(Y1_matrix_rate, ncol = 7, byrow = T))
    
    # Get mean estimated treatment effect rate
    est_trt_effect <- county_means_Y1 - county_means_Y0
    
    # find absolute percent bias ATT 
    absPercBias <- abs(((est_trt_effect  - true_trt_effect_mean_rate)/true_trt_effect_mean_rate) * 100)
    
    meanAbsPercBias[i] <- mean(absPercBias,na.rm = T)
  }
  return(mean(meanAbsPercBias))
}


## Function for GSC
meanAbsPercBiasATT_county_gsc <- function(fit_gsc,Y0_matrix,Y1_matrix,true_trt_effect){
  meanAbsPercBias <- c()
  for(i in 1:length(fit_gsc)){
    # Find true treatment effect rate
    true_trt_effect_rate <- (true_trt_effect/pop_trt[ind] * 100000)[,i]# Find Rate
    true_trt_effect_mean_rate <- rowMeans(matrix(true_trt_effect_rate , ncol = 7, byrow = T)) # take mean across counties
    
    # For each simulation, get mean of 1000 MCMC samples. Then only extract the treated indices.
    # Convert to rate
    med_vec_rate <- as.vector(fit_gsc[[i]]$fit_gsc$Y.ct)[ind]
    
    # For each simulation, take average counterfactual for each county over all time points 
    # should be a vector of length 6
    county_means_Y0 <- rowMeans(matrix(med_vec_rate, ncol = 7, byrow = T))
    
    # For each simulation, take average Y1 for each county over all time points 
    # should be a vector of length 6
    Y1_matrix_rate <- (Y1_matrix[,i]/pop_trt[ind]) * 100000
    county_means_Y1 <- rowMeans(matrix(Y1_matrix_rate, ncol = 7, byrow = T))
    
    # Get mean estimated treatment effect rate
    est_trt_effect <- county_means_Y1 - county_means_Y0
    
    # find absolute percent bias ATT 
    absPercBias <- abs(((est_trt_effect  - true_trt_effect_mean_rate)/true_trt_effect_mean_rate) * 100)
    
    meanAbsPercBias[i] <- mean(absPercBias,na.rm = T)
  }
  return(mean(meanAbsPercBias))
}


###############################################################################################
meanAbsPercBiasATT_county_plot <- function(Y0_matrix_1,Mu_trt_ori_1,Mu_trt_space_1, Mu_trt_spacetime_ICAR_1,Mu_trt_spacetime_AR_1, Mu_trt_space_comb_time_shrink_1, Mu_trt_lasso_1,fit_gsc_1,true_trt_effect_1,Y1_matrix_1,
                                           Y0_matrix_3,Mu_trt_ori_3,Mu_trt_space_3, Mu_trt_spacetime_ICAR_3,Mu_trt_spacetime_AR_3, Mu_trt_space_comb_time_shrink_3, Mu_trt_lasso_3,fit_gsc_3,true_trt_effect_3,Y1_matrix_3,
                                           Y0_matrix_7,Mu_trt_ori_7,Mu_trt_space_7, Mu_trt_spacetime_ICAR_7,Mu_trt_spacetime_AR_7, Mu_trt_space_comb_time_shrink_7, Mu_trt_lasso_7,fit_gsc_7,true_trt_effect_7,Y1_matrix_7,
                                           ylim = c(0,90000), pop_trt, ind){
  ## k = 1
  ori_mean_1 <- meanAbsPercBiasATT_county(Mu_trt_ori_1,Y0_matrix_1, Y1_matrix_1,true_trt_effect_1)
  space_mean_1 <- meanAbsPercBiasATT_county(Mu_trt_space_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  spacetime_ICAR_mean_1 <- meanAbsPercBiasATT_county(Mu_trt_spacetime_ICAR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  spacetime_AR_mean_1 <- meanAbsPercBiasATT_county(Mu_trt_spacetime_AR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  spacetime_shrink_mean_1 <- meanAbsPercBiasATT_county(Mu_trt_space_comb_time_shrink_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  lasso_mean_1 <- meanAbsPercBiasATT_county(Mu_trt_lasso_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  gsc_mean_1 <- meanAbsPercBiasATT_county_gsc(fit_gsc_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1)
  
  ## k = 3
  ori_mean_3 <- meanAbsPercBiasATT_county(Mu_trt_ori_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  space_mean_3 <- meanAbsPercBiasATT_county(Mu_trt_space_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  spacetime_ICAR_mean_3 <- meanAbsPercBiasATT_county(Mu_trt_spacetime_ICAR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  spacetime_AR_mean_3 <- meanAbsPercBiasATT_county(Mu_trt_spacetime_AR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  spacetime_shrink_mean_3 <- meanAbsPercBiasATT_county(Mu_trt_space_comb_time_shrink_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  lasso_mean_3 <- meanAbsPercBiasATT_county(Mu_trt_lasso_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  gsc_mean_3 <- meanAbsPercBiasATT_county_gsc(fit_gsc_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3)
  
  ## k = 7
  ori_mean_7 <- meanAbsPercBiasATT_county(Mu_trt_ori_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  space_mean_7 <- meanAbsPercBiasATT_county(Mu_trt_space_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  spacetime_ICAR_mean_7 <- meanAbsPercBiasATT_county(Mu_trt_spacetime_ICAR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  spacetime_AR_mean_7 <- meanAbsPercBiasATT_county(Mu_trt_spacetime_AR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  spacetime_shrink_mean_7 <- meanAbsPercBiasATT_county(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  lasso_mean_7 <- meanAbsPercBiasATT_county(Mu_trt_lasso_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  gsc_mean_7 <- meanAbsPercBiasATT_county_gsc(fit_gsc_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7)
  
  k <- c(rep("k=1",6),rep("k=3",6),rep("k=7",7))
  Model <- c("Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso","Gsynth",
             "Space-Time \n Shrinkage")
  Mean <- c(ori_mean_1,space_mean_1,spacetime_ICAR_mean_1,spacetime_AR_mean_1,
              lasso_mean_1,gsc_mean_1, ori_mean_3,space_mean_3,spacetime_ICAR_mean_3,spacetime_AR_mean_3,
              lasso_mean_3,gsc_mean_3,ori_mean_7,space_mean_7,spacetime_ICAR_mean_7,spacetime_AR_mean_7,
              lasso_mean_7,gsc_mean_7,spacetime_shrink_mean_7  )
  df <- data.frame(Model,Mean,k) %>%
    mutate(k = fct_relevel(k, "k=7", "k=3", "k=1"),
           Model = fct_relevel(Model,
                               "Gsynth","Original","Space","Space-Time \n ICAR","Space-Time \n AR(1)","Shrinkage \n Lasso"))
  # Thin ggplot environment
  rm(list=setdiff(ls(), c("df", "ylim")))
  
  # plot
  # Minimal theme + blue fill color
  ggplot(data=df, aes(x=Model, y=Mean, fill = k)) +
    scale_fill_manual(values=cbPalette,
                      breaks=c('k=1', 'k=3', 'k=7')) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=round(Mean,2)),position = position_dodge(1), hjust = -.2, size=2.8)+
    #ggtitle("Absolute Percent Bias of the ATT (Non-Rare, Rate = .5)") +
    coord_flip(ylim = ylim ) + 
    scale_x_discrete(limits = rev) +
    #scale_y_continuous(name="Mean", limits=c(0, 90000)) +
    theme_bw() + #adds plot outline vs. theme_minimal
    theme(axis.text=element_text(color="black"),
          legend.title = element_blank()
    ) 
}


####################################################################################################
# Function to get median absolute percent bias of the counterfactuals 
# This function averages all counties for each time point to compare averages
# This is done to compare appropriately to ASC by Ben-Michael et al.
medianAbsPercBiasY0_year <- function(Mu_trt,Y0_matrix,ind){
  medianAbsPercBias <- c()
  for(i in 1:length(Mu_trt)){
  #for(i in c(seq(1:83),seq(85,99))){
    #print(i)
    # Find true year mean rate
    Y0_rate <- (Y0_matrix[ind,i]/pop_trt[ind]) * 100000 # Find Rate
    true_Y0_median_rate <- colMeans(matrix(Y0_rate , ncol = 7, byrow = T)) # take mean across years
    
    # For each simulation, get median of 1000 MCMC samples. Then only extract the treated indices.
    med_vec <- colMedians(Mu_trt[[i]],na.rm = T)[ind]
    
    # Convert to rate
    med_vec_rate <- (med_vec/pop_trt[ind]) * 100000
    
    # For each simulation, take average counterfactual for each time over all counties
    # should be a vector of length 
    county_medians <- colMeans(matrix(med_vec_rate, ncol = 7, byrow = T))
    
    absPercBias <- abs(((true_Y0_median_rate - county_medians)/true_Y0_median_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  
  median <- median(medianAbsPercBias)
  LB <- quantile(medianAbsPercBias, probs = 0.25, na.rm = T)
  UB <- quantile(medianAbsPercBias, probs = 0.75, na.rm = T)
  IQR <- IQR(medianAbsPercBias)
  
  return(list(median = median, LB = LB, UB = UB, IQR = IQR))
}

#--------------------------------------------------------------------------------------------

# Function to get median absolute percent bias of the ATT
# This function averages all counties for each time point to compare averages
# This is done to compare appropriately to ASC by Ben-Michael et al.
medianAbsPercBiasATT_year <- function(Mu_trt,Y0_matrix,Y1_matrix, true_trt_effect, ind){
  medianAbsPercBias <- c()
  for(i in 1:length(Mu_trt)){
  #for(i in c(seq(1:83),seq(85,99))){
    # Find true treatment effect rate
    true_trt_effect_rate <- (true_trt_effect/pop_trt[ind] * 100000)[,i]# Find Rate
    # Take average across counties for each time point
    true_trt_effect_median_rate <- colMeans(matrix(true_trt_effect_rate , ncol = 7, byrow = T)) 
    
    # For each simulation, get median of 1000 MCMC samples. Then only extract the treated indices.
    med_vec <- colMedians(Mu_trt[[i]],na.rm = T)[ind]
    
    # Convert to rate
    Y0_vec_rate <- (med_vec/pop_trt[ind]) * 100000
    Y1_matrix_rate <- (Y1_matrix[,i]/pop_trt[ind]) * 100000
    
    # Calculate ATT for each time point
    TE <- Y1_matrix_rate - Y0_vec_rate
    ATT_vec <- colMeans(matrix(TE, ncol = 7, byrow = T)) # one ATT for each time
    
    # find absolute percent bias ATT 
    absPercBias <- abs(((ATT_vec  - true_trt_effect_median_rate)/true_trt_effect_median_rate) * 100)
    
    medianAbsPercBias[i] <- median(absPercBias,na.rm = T)
  }
  median <- median(medianAbsPercBias)
  LB <- quantile(medianAbsPercBias, probs = 0.25, na.rm = T)
  UB <- quantile(medianAbsPercBias, probs = 0.75, na.rm = T)
  IQR <- IQR(medianAbsPercBias)
  
  return(list(median = median, LB = LB, UB = UB, IQR = IQR))
}

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

## Function to filter out -1s
filtering <- function(Y_pred){
  for(i in 1:length(Y_pred)){
    if(length(which(Y_pred[[i]] == -1)) > 0){
      
      # Identify the rows -1s (identify which MCMC iteration has the overflow problem)
      rows <- which(Y_pred[[i]] == -1, arr.ind = TRUE)[,1] #[,1] since first column of the which statement identify rows
      
      # Filter out
      print("filtering")
      Y_pred[[i]] <- Y_pred[[i]][-rows,]
      
    }
  }
  return(Y_pred)
}


