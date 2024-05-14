# Functions for application section

#---------------------------------------------------------------------------------------------------------

# Run gsynth on applied data 
# Convert cases into rates for gsynth

run_gsynth <- function(data_full_hisp, k = 2){
  ######################################################################################################
  # Input                                                                                              #
  # data_full_hisp: data set containing cases (CL_CASES), indicator for treated (C), population (POP)  #
  # k: number of factors for matrix completion                                                         #
  # Output                                                                                             #
  # fit_gsc: gsynth fit                                                                                #
  ######################################################################################################
  
  ## Convert
  data_full_hisp$rate <- 100000 * (data_full_hisp$CL_CASES/data_full_hisp$POP)
  
  ## Run Gysnth
  fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_hisp, estimator = "mc", 
                    index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                    se = TRUE, k=k, nboots = 1000, parallel = FALSE)
  return(fit_gsc)
}

###########################################################################################################

# Diagnostics functions

# Function to create a convergence table

conv_table <- function(Mu_trt_ori, Mu_trt_space, Mu_trt_spacetime_ICAR, Mu_trt_spacetime_AR, 
                       Mu_trt_space_comb_time_shrink, Mu_trt_lasso){
  ######################################################################################################
  # Input                                                                                              #
  # The STAN results for each of our models                                                            #
  # Output                                                                                             #
  # convergence table                                                                                  #
  ######################################################################################################
  
  # Vanilla Model
  ori_conv <- length(unique(which(is.na(Mu_trt_ori)*1 == 1, arr.ind = TRUE)[,1])) 
  
  # Space Model
  space_conv <- length(unique(which(is.na(Mu_trt_space)*1 == 1, arr.ind = TRUE)[,1])) 
  
  # Space-Time ICAR Model
  spacetime_ICAR_conv <- length(unique(which(is.na(Mu_trt_spacetime_ICAR)*1 == 1, arr.ind = TRUE)[,1])) 

  # Space-Time AR Model
  spacetime_AR_conv <- length(unique(which(is.na(Mu_trt_spacetime_AR)*1 == 1, arr.ind = TRUE)[,1])) 
 
  # Space-Time Shrinkage Model
  space_comb_time_shrink_conv <- length(unique(which(is.na(Mu_trt_space_comb_time_shrink)*1 == 1, arr.ind = TRUE)[,1]))
  
  # Lasso Model
  lasso_conv <- length(unique(which(is.na(Mu_trt_lasso)*1 == 1, arr.ind = TRUE)[,1]))
  
  model <- c(rep("Vanilla",5),rep("Space",5),rep("Space-Time ICAR",5),rep("Space-Time AR",5),rep("Space-Time Shrinkage",5), rep("Lasso",5))
  k <- rep(c("k=2","k=2","k=2","k=2","k=3","k=2"),5)
  df <- rbind(ori_conv,space_conv,spacetime_ICAR_conv,spacetime_AR_conv,space_comb_time_shrink_conv, lasso_conv)/1000*100
  rownames(df) <- c("Vanilla","Space","Space-Time ICAR","Space-Time AR","Space-Time Shrinkage", "Bayesian Lasso")

  # Make Table
  table <- df %>%
    kbl(caption = "Percent NAs for 1000 Iterations") %>%
    kable_classic_2(full_width = F)
  
  return(table)
}
#---------------------------------------------------------------------------------------------------------

# Function to create a convergence table

rhat_table_YPred <- function(Mu_trt_ori, Mu_trt_space, Mu_trt_spacetime_ICAR, Mu_trt_spacetime_AR, 
                       Mu_trt_space_comb_time_shrink, Mu_trt_lasso, ind){
  ######################################################################################################
  # Input                                                                                              #
  # The STAN results for each of our models                                                            #
  # ind: indices for treated values to be included                                                     #                                                  
  # Output                                                                                             #
  # rhat table                                                                                         #
  ######################################################################################################
  library(rstan)
  # Vanilla Model
  ori_summary <- summary(fit_sim_ori, pars = c("Y_pred"))$summary[ind,]
  rhat_ori <- mean(ori_summary[,"Rhat"])
  
  # Space Model
  space_summary <- summary(fit_sim_space, pars = c("Y_pred"))$summary[ind,]
  rhat_space <- mean(space_summary[,"Rhat"])
 
  # ICAR Model
  ICAR_summary <- summary(fit_sim_space_time_ICAR, pars = c("Y_pred"))$summary[ind,]
  rhat_ICAR <- mean(ICAR_summary[,"Rhat"])
  
  # AR Model
  AR_summary <- summary(fit_sim_space_time_AR, pars = c("Y_pred"))$summary[ind,]
  rhat_AR <- mean(AR_summary[,"Rhat"])
  
  # Shrinkage Model
  shrink_summary <- summary(fit_sim_space_comb_time_shrink, pars = c("Y_pred"))$summary[ind,]
  rhat_shrink <- mean(shrink_summary[,"Rhat"])
  #rhat_shrink <- NA
  
  # Lasso Model
  lasso_summary <- summary(fit_sim_lasso, pars = c("Y_pred"))$summary[ind,]
  rhat_lasso <- mean(lasso_summary[,"Rhat"])

  # Table of Results
  df <- rbind(rhat_ori,rhat_space,rhat_ICAR,rhat_AR,rhat_shrink, rhat_lasso)
  rownames(df) <- c("Vanilla","Space","Space-Time ICAR","Space-Time AR","Space-Time Shrinkage", "Bayesian Lasso")
  
  # Make Table
  table <- df %>%
    kbl(caption = "Mean R-Hat Values for Y(0)") %>%
    kable_classic_2(full_width = F)
  
  return(table)
}

rhat_table_YPred_v2 <- function(fit_sim_ori, fit_sim_space, fit_sim_space_time_ICAR, fit_sim_space_time_AR, 
                                fit_sim_space_comb_time_shrink,fit_sim_lasso, ind){
  ######################################################################################################
  # Input                                                                                              #
  # The STAN results for each of our models                                                            #
  # ind: indices for treated values to be included                                                     #                                                  
  # Output                                                                                             #
  # rhat df.                                                                                           #
  ######################################################################################################
  
  library(rstan)
  
  # Vanilla Model
  ori_summary <- summary(fit_sim_ori, pars = c("Y_pred"))$summary[ind,]
  rhat_ori <- mean(ori_summary[,"Rhat"])
  
  # Space Model
  space_summary <- summary(fit_sim_space, pars = c("Y_pred"))$summary[ind,]
  rhat_space <- mean(space_summary[,"Rhat"])
  
  # ICAR Model
  ICAR_summary <- summary(fit_sim_space_time_ICAR, pars = c("Y_pred"))$summary[ind,]
  rhat_ICAR <- mean(ICAR_summary[,"Rhat"])
  
  # AR Model
  AR_summary <- summary(fit_sim_space_time_AR, pars = c("Y_pred"))$summary[ind,]
  rhat_AR <- mean(AR_summary[,"Rhat"])
  #rhat_AR <- NA
  
  # Shrinkage Model
  shrink_summary <- summary(fit_sim_space_comb_time_shrink, pars = c("Y_pred"))$summary[ind,]
  rhat_shrink <- mean(shrink_summary[,"Rhat"])
  #rhat_shrink <- NA
  
  # Lasso Model
  lasso_summary <- summary(fit_sim_lasso, pars = c("Y_pred"))$summary[ind,]
  rhat_lasso <- mean(lasso_summary[,"Rhat"])
  
  # Table of Results
  df <- rbind(rhat_ori,rhat_space,rhat_ICAR,rhat_AR, rhat_lasso, rhat_shrink)
  rownames(df) <- c("Vanilla","Space","Space-Time ICAR","Space-Time AR", "Space-Time Lasso", "Space-Time Shrinkage")
  
  
  return(df)
}


# Function to create a convergence table for smooth and non-smooth data

rhat_table_YPred_wSmooth <- function(fit_sim_ori, fit_sim_space, fit_sim_spacetime_ICAR, fit_sim_spacetime_AR, 
                                     fit_sim_space_comb_time_shrink, fit_sim_lasso, 
                                     fit_sim_ori_smooth, fit_sim_space_smooth, fit_sim_spacetime_ICAR_smooth, fit_sim_spacetime_AR_smooth, 
                                     fit_sim_space_comb_time_shrink_smooth, fit_sim_lasso_smooth, ind){
  ######################################################################################################
  # Input                                                                                              #
  # The STAN results for each of our models                                                            #
  # ind: indices for treated values to be included                                                     #                                                  
  # Output                                                                                             #
  # rhat table                                                                                         #
  ######################################################################################################
  library(rstan)
  
  # Vanilla Model
  ori_summary <- summary(fit_sim_ori, pars = c("Y_pred"))$summary[ind,]
  rhat_ori <- mean(ori_summary[,"Rhat"])
  
  # Space Model
  space_summary <- summary(fit_sim_space, pars = c("Y_pred"))$summary[ind,]
  rhat_space <- mean(space_summary[,"Rhat"])
  
  # ICAR Model
  ICAR_summary <- summary(fit_sim_space_time_ICAR, pars = c("Y_pred"))$summary[ind,]
  rhat_ICAR <- mean(ICAR_summary[,"Rhat"])
  
  # AR Model
  AR_summary <- summary(fit_sim_space_time_AR, pars = c("Y_pred"))$summary[ind,]
  rhat_AR <- mean(AR_summary[,"Rhat"])
  
  # Shrinkage Model
  shrink_summary <- summary(fit_sim_space_comb_time_shrink, pars = c("Y_pred"))$summary[ind,]
  rhat_shrink <- mean(shrink_summary[,"Rhat"])
  #rhat_shrink <- NA
  
  # Lasso Model
  lasso_summary <- summary(fit_sim_lasso, pars = c("Y_pred"))$summary[ind,]
  rhat_lasso <- mean(lasso_summary[,"Rhat"])
  
  ## Smooth
  # Vanilla Model
  ori_summary_smooth <- summary(fit_sim_ori_smooth, pars = c("Y_pred"))$summary[ind,]
  rhat_ori_smooth <- mean(ori_summary_smooth[,"Rhat"])
  
  # Space Model
  space_summary_smooth <- summary(fit_sim_space_smooth, pars = c("Y_pred"))$summary[ind,]
  rhat_space_smooth <- mean(space_summary_smooth[,"Rhat"])
  
  # ICAR Model
  ICAR_summary_smooth <- summary(fit_sim_space_time_ICAR_smooth, pars = c("Y_pred"))$summary[ind,]
  rhat_ICAR_smooth <- mean(ICAR_summary_smooth[,"Rhat"])
  
  # AR Model
  AR_summary_smooth <- summary(fit_sim_space_time_AR_smooth, pars = c("Y_pred"))$summary[ind,]
  rhat_AR_smooth <- mean(AR_summary_smooth[,"Rhat"])
  
  # Shrinkage Model
  shrink_summary_smooth <- summary(fit_sim_space_comb_time_shrink_smooth, pars = c("Y_pred"))$summary[ind,]
  rhat_shrink_smooth <- mean(shrink_summary_smooth[,"Rhat"])
  #rhat_shrink <- NA
  
  # Lasso Model
  lasso_summary_smooth <- summary(fit_sim_lasso_smooth, pars = c("Y_pred"))$summary[ind,]
  rhat_lasso_smooth <- mean(lasso_summary_smooth[,"Rhat"])
  
  # Table of Results
  df <- cbind(c(rhat_ori,rhat_space,rhat_ICAR,rhat_AR, rhat_lasso, rhat_shrink),
              c(rhat_ori_smooth,rhat_space_smooth,rhat_ICAR_smooth,rhat_AR_smooth, rhat_lasso_smooth, rhat_shrink_smooth))
  rownames(df) <- c("Vanilla","Space","Space-Time ICAR","Space-Time AR", "Space-Time Lasso","Space-Time Shrinkage")
  
  # Make Table
  table <- df %>%
    kbl(caption = "Mean R-Hat Values for Y(0)") %>%
    kable_classic_2(full_width = F)
  
  return(table)
}

#---------------------------------------------------------------------------------------------------------
# Outputs raw data vs estimated Y(0) averaging all counties
raw_v_est_plot <- function(data,pop_trt,mu_matrix_l,treated,n_trt,m,Y1) {
  full_Y1_obs1 <- Y1
  
  full_time <- unique(data$YEAR_DX)
  
  
  ### Posterior
  post_Mu_trt_rate_l <- 100000 * sweep(mu_matrix_l,2,pop_trt,"/") #get rate (1000 x 195)
  post_Mu_trt_avg_rate_l <- colMeans(post_Mu_trt_rate_l) #take median of 1000 iterations (length 195)
  post_Mu_trt_matrix <- matrix(post_Mu_trt_avg_rate_l,nrow=n_trt,ncol = m,byrow=T) #(13x15)
  
  
  post_Mu_trt_avg_rate_l <- colMeans(post_Mu_trt_matrix, na.rm = T) #take median of counties in each time point (length 15)
  post_Mu_trt_avg_lower_l <- apply(post_Mu_trt_matrix, 2, function(x) {quantile(x, probs = 0.05, na.rm = T)})
  post_Mu_trt_avg_upper_l <- apply(post_Mu_trt_matrix, 2, function(x) {quantile(x, probs = 0.95, na.rm = T)})    
  
  #posterior sample df
  post_trt_l <- data.frame(time = full_time, mean = post_Mu_trt_avg_rate_l, 
                           lower = post_Mu_trt_avg_lower_l, upper = post_Mu_trt_avg_upper_l, 
                           type = factor("post"))
  ### Raw
  raw_trt_rate_l <- 100000 * full_Y1_obs1/pop_trt #get rate
  full_raw_trt_matrix <- matrix(raw_trt_rate_l,nrow=n_trt,ncol = m,byrow=T)
  
  #raw_trt_avg_rate_l <- colMeans(raw_trt_rate_l)
  raw_trt_avg_rate_l <- colMeans(full_raw_trt_matrix, na.rm = T)
  raw_trt_avg_lower_l <- apply(full_raw_trt_matrix, 2, function(x) {quantile(x, probs = 0.05, na.rm = T)})
  raw_trt_avg_upper_l <- apply(full_raw_trt_matrix, 2, function(x) {quantile(x, probs = 0.95, na.rm = T)})  
  
  #observed sample df
  raw_trt_l <- data.frame(time = full_time, mean = raw_trt_avg_rate_l , 
                          lower = raw_trt_avg_lower_l, upper = raw_trt_avg_upper_l, 
                          type = factor("raw"))
  
  
  to_plot <- rbind(post_trt_l, raw_trt_l)
  
  p1<- ggplot(data = to_plot, aes(x = time, y = mean, group = type, color = type)) + geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.1) + 
    ylab( 'Number of Cases per 100000') + geom_vline(xintercept = 1995, color = 'black') +
    ggtitle("Comparing Raw Treated Average Rate with \n Estimated Posterior Average Rate")  
  
  plot(p1)
  
}
#####################################################################################
# Raw vs. Estimated Plot for SVT
raw_v_est_plot_svt <- function(data,pop_trt,mu_matrix_l,treated,n_trt,m,raw){
  
  full_Y1_obs1 <- raw
  
  full_time <- unique(data$YEAR_DX)
  
  ### Posterior
  post_Mu_trt_matrix <- mu_matrix_l # # counties x # of years
  
  
  post_Mu_trt_avg_rate_l <- colMeans(post_Mu_trt_matrix, na.rm = T) #take median of counties in each time point (length 15)
  post_Mu_trt_avg_lower_l <- apply(post_Mu_trt_matrix, 2, function(x) {quantile(x, probs = 0.05, na.rm = T)}) # get quantiles since averaging over counties
  post_Mu_trt_avg_upper_l <- apply(post_Mu_trt_matrix, 2, function(x) {quantile(x, probs = 0.95, na.rm = T)})    
  
  #posterior sample df
  post_trt_l <- data.frame(time = full_time, mean = post_Mu_trt_avg_rate_l, 
                           lower = post_Mu_trt_avg_lower_l, upper = post_Mu_trt_avg_upper_l, 
                           type = factor("estimated"))
  
  ### Raw
  raw_trt_rate_l <- 100000 * full_Y1_obs1/pop_trt #get rate
  full_raw_trt_matrix <- matrix(raw_trt_rate_l,nrow=n_trt,ncol = m,byrow=T)
  
  #raw_trt_avg_rate_l <- colMeans(raw_trt_rate_l)
  raw_trt_avg_rate_l <- colMeans(full_raw_trt_matrix, na.rm = T)
  raw_trt_avg_lower_l <- apply(full_raw_trt_matrix, 2, function(x) {quantile(x, probs = 0.05, na.rm = T)})
  raw_trt_avg_upper_l <- apply(full_raw_trt_matrix, 2, function(x) {quantile(x, probs = 0.95, na.rm = T)})  
  
  #observed sample df
  raw_trt_l <- data.frame(time = full_time, mean = raw_trt_avg_rate_l , 
                          lower = raw_trt_avg_lower_l, upper = raw_trt_avg_upper_l, 
                          type = factor("raw"))
  
  
  to_plot <- rbind(post_trt_l, raw_trt_l)
  
  p1<- ggplot(data = to_plot, aes(x = time, y = mean, group = type, color = type)) + geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.1) + 
    ylab( 'Number of Cases per 100000') + geom_vline(xintercept = 1995, color = 'black') +
    ggtitle("Comparing Raw Treated Average Rate with \n Estimated Average Rate")  
  
  plot(p1)
  
}

#####################################################################################
# Raw vs. Estimated Plot for GSC
raw_v_est_plot_gsc  <- function(data,pop_trt,mu_matrix_l,treated,n_trt,m,raw){
  
  full_Y1_obs1 <- raw
  
  full_time <- unique(data$YEAR_DX)
  
  ### Posterior
  post_Mu_trt_matrix <- t(mu_matrix_l$Y.ct) # # counties x # of years
  
  
  post_Mu_trt_avg_rate_l <- colMeans(post_Mu_trt_matrix, na.rm = T) #take median of counties in each time point (length 15)
  post_Mu_trt_avg_lower_l <- apply(post_Mu_trt_matrix, 2, function(x) {quantile(x, probs = 0.05, na.rm = T)}) # get quantiles since averaging over counties
  post_Mu_trt_avg_upper_l <- apply(post_Mu_trt_matrix, 2, function(x) {quantile(x, probs = 0.95, na.rm = T)})    
  
  #posterior sample df
  post_trt_l <- data.frame(time = full_time, mean = post_Mu_trt_avg_rate_l, 
                           lower = post_Mu_trt_avg_lower_l, upper = post_Mu_trt_avg_upper_l, 
                           type = factor("estimated"))
  
  ### Raw
  raw_trt_rate_l <- 100000 * full_Y1_obs1/pop_trt #get rate
  full_raw_trt_matrix <- matrix(raw_trt_rate_l,nrow=n_trt,ncol = m,byrow=T)
  
  #raw_trt_avg_rate_l <- colMeans(raw_trt_rate_l)
  raw_trt_avg_rate_l <- colMeans(full_raw_trt_matrix, na.rm = T)
  raw_trt_avg_lower_l <- apply(full_raw_trt_matrix, 2, function(x) {quantile(x, probs = 0.05, na.rm = T)})
  raw_trt_avg_upper_l <- apply(full_raw_trt_matrix, 2, function(x) {quantile(x, probs = 0.95, na.rm = T)})  
  
  #observed sample df
  raw_trt_l <- data.frame(time = full_time, mean = raw_trt_avg_rate_l , 
                          lower = raw_trt_avg_lower_l, upper = raw_trt_avg_upper_l, 
                          type = factor("raw"))
  
  
  to_plot <- rbind(post_trt_l, raw_trt_l)
  
  p1<- ggplot(data = to_plot, aes(x = time, y = mean, group = type, color = type)) + geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.1) + 
    ylab( 'Number of Cases per 100000') + geom_vline(xintercept = 1995, color = 'black') +
    ggtitle("Comparing Raw Treated Average Rate with \n Estimated Average Rate")  
  
  plot(p1)
  
}


###################################################################################################
# get 
raw_v_est_county_plot <- function(data,pop_trt,mu_matrix_l,treated,n_trt,m,Y1,county_id) {
  full_Y1_obs1 <- Y1
  
  full_time <- unique(data$YEAR_DX)
  
  ### Get treated indices for county
  county_ind <- seq((county_id-1)*m + 1, county_id*m)
  
  ### Posterior
  post_Mu_trt_rate_l <- 100000 * sweep(mu_matrix_l,2,pop_trt,"/") #get rate (1000 x 195)
  post_Mu_trt_avg_rate_l <- colMeans(post_Mu_trt_rate_l) #take median of 1000 iterations (length 195)
  post_Mu_trt_matrix <- matrix(post_Mu_trt_avg_rate_l,nrow=n_trt,ncol = m,byrow=T) #(13x15)
  
  
  post_Mu_trt_avg_rate_l <- post_Mu_trt_matrix[county_id,]  #(length 15)
  post_Mu_trt_avg_lower_l <- apply(post_Mu_trt_rate_l[,county_ind], 2, function(x) {quantile(x, probs = 0.05, na.rm = T)})
  post_Mu_trt_avg_upper_l <- apply(post_Mu_trt_rate_l[,county_ind], 2, function(x) {quantile(x, probs = 0.95, na.rm = T)})    
  
  #posterior sample df
  post_trt_l <- data.frame(time = full_time, mean = post_Mu_trt_avg_rate_l, 
                           lower = post_Mu_trt_avg_lower_l, upper = post_Mu_trt_avg_upper_l, 
                           type = factor("post"))
  ### Raw
  raw_trt_rate_l <- 100000 * full_Y1_obs1/pop_trt #get rate
  full_raw_trt_matrix <- matrix(raw_trt_rate_l,nrow=n_trt,ncol = m,byrow=T)
  
  #raw_trt_avg_rate_l <- colMeans(raw_trt_rate_l)
  raw_trt_avg_rate_l <- full_raw_trt_matrix[county_id,]
  raw_trt_avg_lower_l <- NA
  raw_trt_avg_upper_l <- NA
  
  #observed sample df
  raw_trt_l <- data.frame(time = full_time, mean = raw_trt_avg_rate_l , 
                          lower = raw_trt_avg_lower_l, upper = raw_trt_avg_upper_l, 
                          type = factor("raw"))
  
  
  to_plot <- rbind(post_trt_l, raw_trt_l)
  
  p1<- ggplot(data = to_plot, aes(x = time, y = mean, group = type, color = type)) + geom_line() + geom_point() + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.1) + 
    ylab( 'Number of Cases per 100000') + geom_vline(xintercept = 1995, color = 'black') +
    ggtitle("Comparing Raw Treated Average Rate with \n Estimated Posterior Average Rate")  
  
  plot(p1)
  
}
#######################################################################################################

# Results

# ATT function
get_ATT <- function(Mu_trt,Y1,pop_trt,n_trt,m_trt,m) {
  ATT_mat <- matrix(NA,nrow=m_trt, ncol = nrow(Mu_trt))
  Y1_temp <- matrix(Y1, nrow=n_trt, byrow = T) #put y1 in matrix (13x7)
  pop_trt_temp <- matrix(pop_trt, nrow = n_trt, byrow = T)[,(m-m_trt+1):m] #(13 x 7)
  
  for(i in 1:nrow(Mu_trt)){ #for all MCMC samples
    Mu_trt_temp <- Mu_trt[i,] #(1 x 208) 195 = 13 counties x m time points take one iteration
    Mu_trt_temp_mat <- matrix(Mu_trt_temp, nrow = n_trt, byrow = T) #turn into (13 x 15) matrix
    Y0_est <- Mu_trt_temp_mat[,(m-m_trt+1):m] #estimated counterfactuals (13x7) matrix (trt counties x trt time)
    
    diff <- Y1_temp - Y0_est # observed Y(1) - estimated Y(0)
    rate <- 100000 * (diff / pop_trt_temp) # get rate (13 x 7)
    ATT_mat[,i] <- colMeans(rate)  #get avg of all counties in each time point  (1 x 7) and store into column
  }
  
  return(ATT_mat)
}


get_ATT_full <- function(Mu_trt,Y1,pop_trt,n_trt,m_trt,m) {
  ATT_mat <- matrix(NA,nrow=m, ncol = nrow(Mu_trt))
  Y1_temp <- matrix(Y1, nrow=n_trt, byrow = T) #put y1 in matrix (13x7)
  pop_trt_temp <- matrix(pop_trt, nrow = n_trt, byrow = T)[,1:m] #(13 x 7)
  
  for(i in 1:nrow(Mu_trt)){ #for all MCMC samples
    
    Mu_trt_temp <- Mu_trt[i,] #(1 x 208) 195 = 13 counties x m time points take one iteration
    Mu_trt_temp_mat <- matrix(Mu_trt_temp, nrow = n_trt, byrow = T) #turn into (13 x 15) matrix
    Y0_est <- Mu_trt_temp_mat #estimated counterfactuals (13x7) matrix (trt counties x trt time)
    
    diff <- Y1_temp - Y0_est # observed Y(1) - estimated Y(0)
    rate <- 100000 * (diff / pop_trt_temp) # get rate (13 x 7)
    ATT_mat[,i] <- colMeans(rate)  #get avg of all counties in each time point  (1 x 7) and store into column
  }
  
  return(ATT_mat)
}

# ATT function gsc
get_ATT_gsc <- function(Mu_trt_gsc,Y1,pop_trt,ind,n_trt) {
  med_list <- c()
  ATT <- c()
  med_list <- (100000 * (Y1)/pop_trt) - Mu_trt_gsc
  #perc_bias_list[[i]] <- 100*(med_list[[i]]-true_trt_rate_list[[i]])/true_trt_rate_list[[i]] #(median trt effect - true trt effect)/true trt
  ATT <- colMeans(matrix(med_list,nrow = n_trt, byrow=T))
  
  return(ATT)
  # med list returns estimate of the treatment effect
  
}

#-------------------------------------------------------------------------------------------------------



ATT_CACT_plot<- function(data_full_hisp, Mu_trt_ori, Mu_trt_space, Mu_trt_spacetime_ICAR, Mu_trt_spacetime_AR, 
                        Mu_trt_space_comb_time_shrink, Mu_trt_lasso, fit_gsc){
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- get_Y1_obs1(treated1)
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  ind <- c()
  for(i in 0:(n_trt_CA-1)){
    ind <- append(ind,seq((m_trt-1),m) + (i*m))
  }
  
  # Subset STAN output (Cali is first 5 counties)
  #Mu_trt_ori_CA <- Mu_trt_ori[,1:(length(unique(treated1$FIPS))*m_trt)]
  #Mu_trt_space_CA <- Mu_trt_space[,1:(length(unique(treated1$FIPS))*m_trt)]
  #Mu_trt_spacetime_ICAR_CA <- Mu_trt_spacetime_ICAR[,1:(length(unique(treated1$FIPS))*m_trt)]
  #Mu_trt_spacetime_AR_CA <- Mu_trt_spacetime_AR[,1:(length(unique(treated1$FIPS))*m_trt)]
  #Mu_trt_lasso_CA <- Mu_trt_lasso[,1:(length(unique(treated1$FIPS))*m_trt)]
  #Mu_trt_space_comb_time_shrink_CA <- Mu_trt_space_comb_time_shrink[,1:(length(unique(treated1$FIPS))*m_trt)]
  #Mu_trt_gsc <- as.vector(t(fit_gsc$Y.ct)[1:(length(unique(treated1$FIPS))),(m_trt-1):m])
  
  Mu_trt_ori_CA <- Mu_trt_ori[,ind]
  Mu_trt_space_CA <- Mu_trt_space[,ind]
  Mu_trt_spacetime_ICAR_CA <- Mu_trt_spacetime_ICAR[,ind]
  Mu_trt_spacetime_AR_CA <- Mu_trt_spacetime_AR[,ind]
  Mu_trt_lasso_CA <- Mu_trt_lasso[,ind]
  Mu_trt_space_comb_time_shrink_CA <- Mu_trt_space_comb_time_shrink[,ind]
  Mu_trt_gsc <- as.vector(t(fit_gsc$Y.ct)[1:(length(unique(treated1$FIPS))),(m_trt-1):m])
  
  
  # Subset Population
  pop_trt_CA <- treated1 %>% filter(YEAR_DX >= "1995")
  pop_trt_CA <- pop_trt_CA$POP
  
  
  # Calculate ATT
  # Vanilla Model
  ATT_ori_CA <- get_ATT(Mu_trt_ori_CA,Y1_CA,pop_trt_CA,n_trt_CA,m_trt,m_trt)
  
  # Space Model
  ATT_space_CA <- get_ATT(Mu_trt_space_CA,Y1_CA,pop_trt_CA,n_trt_CA,m_trt,m_trt)
  #ATT_space_CA <- NA
  
  # Space-Time Model ICAR
  ATT_spacetime_ICAR_CA <- get_ATT(Mu_trt_spacetime_ICAR_CA,Y1_CA,pop_trt_CA,n_trt_CA,m_trt,m_trt)
  #ATT_spacetime_ICAR_CA <- NA
  
  # Space-Time Model AR
  ATT_spacetime_AR_CA <- get_ATT(Mu_trt_spacetime_AR_CA,Y1_CA,pop_trt_CA,n_trt_CA,m_trt,m_trt)
  
  # Shrinkage Model
  ATT_space_comb_time_shrink_CA <- get_ATT(Mu_trt_space_comb_time_shrink_CA,Y1_CA,pop_trt_CA,n_trt_CA,m_trt,m_trt)
  
  # Lasso Model
  ATT_lasso_CA <- get_ATT(Mu_trt_lasso_CA,Y1_CA,pop_trt_CA,n_trt_CA,m_trt,m_trt)
  
  # Gsynth Model
  ATT_gsc <- get_ATT_gsc(Mu_trt_gsc,Y1_CA,pop_trt_CA,ind,n_trt_CA)
  #ATT_gsc <- fit_gsc$att[(m-m_trt+1):m] # just treated years
  #gsc_bounds <- fit_gsc$est.att[(m-m_trt+1):m,c("CI.lower","CI.upper")] # inference for att
  gsc_bounds <- rbind(rep(0,m_trt),rep(0,m_trt))
  
  # Summarize ATT
  # Vanilla Model
  med_ori <- rowMedians(ATT_ori_CA, na.rm = TRUE)
  bounds_ori <- apply(ATT_ori_CA, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  # Space Model
  med_space <- rowMedians(ATT_space_CA, na.rm = TRUE)
  bounds_space <- apply(ATT_space_CA, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  #med_space <- rep(0,m_trt)
  #bounds_space <- rbind(rep(0,m_trt),rep(0,m_trt))
  
  # Space-Time ICAR Model 
  med_spacetime_ICAR <- rowMedians(ATT_spacetime_ICAR_CA, na.rm = TRUE)
  bounds_spacetime_ICAR <- apply(ATT_spacetime_ICAR_CA, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  #med_spacetime_ICAR <- rep(0,m_trt)
  #bounds_spacetime_ICAR <- rbind(rep(0,m_trt),rep(0,m_trt))
  
  # Space-Time AR Model 
  med_spacetime_AR <- rowMedians(ATT_spacetime_AR_CA, na.rm = TRUE)
  bounds_spacetime_AR <- apply(ATT_spacetime_AR_CA, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  # Space-Time Shrinkage Model 
  med_space_comb_time_shrink <- rowMedians(ATT_space_comb_time_shrink_CA, na.rm = TRUE)
  bounds_space_comb_time_shrink <- apply(ATT_space_comb_time_shrink_CA, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  # for now since this hasnt finished running
  #med_space_comb_time_shrink <- rep(0,m_trt)
  #bounds_space_comb_time_shrink <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  # Lasso Model
  med_lasso <- rowMedians(ATT_lasso_CA, na.rm = TRUE)
  bounds_lasso <- apply(ATT_lasso_CA, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  #med_lasso <- rep(0,m_trt)
  #bounds_lasso <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  ## Plot
  model <- c(rep("Vanilla",m_trt),rep("Space",m_trt),rep("Space-Time \n ICAR",m_trt),rep("Space-Time \n AR",m_trt),rep("Space-Time \n Shrinkage",m_trt), rep("Bayesian \n Lasso", m_trt), rep("Gsynth",m_trt))
  year <- rep(c("1995","1996","1997","1998","1999","2000","2001","2002","2003"),7)
  #UB <- c(med_ori,med_space,med_spacetime_ICAR,med_spacetime_AR,med_space_comb_time_shrink,med_lasso,ATT_gsc)
  LB <-  c(bounds_ori[1,],bounds_space[1,],bounds_spacetime_ICAR[1,],bounds_spacetime_AR[1,],bounds_space_comb_time_shrink[1,], bounds_lasso[1,], gsc_bounds[1,])
  UB <-  c(bounds_ori[2,],bounds_space[2,],bounds_spacetime_ICAR[2,],bounds_spacetime_AR[2,],bounds_space_comb_time_shrink[2,],bounds_lasso[2,],gsc_bounds[2,])
  ATT <- c(med_ori,med_space,med_spacetime_ICAR,med_spacetime_AR,med_space_comb_time_shrink,med_lasso,ATT_gsc)
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  ATT_df = ATT_df %>%
    mutate(Model = fct_relevel(Model, "Gsynth", "Vanilla", "Space", "Space-Time \n ICAR", "Space-Time \n AR", "Space-Time \n Shrinkage", "Bayesian \n Lasso"))
  
  plot_CA <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = Model))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year, fill = Model), alpha = 0.1)+
    geom_hline(yintercept=0)+
    #scale_colour_manual("",values="blue")+
    #scale_fill_manual("",values="grey12")
    #theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+
    facet_wrap(~Model, nrow = 1) + 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Estimated ATT for California")
  # Observed Median Rate ") + geom_vline(xintercept = 1995, color = 'red')
  
  # Get ATT JUST for Connecticut
  Y1_CT <- Y1[(length(Y1_CA)+1):length(Y1)]
  
  n_trt_CT <- length(unique(unique(data_full_hisp %>% filter(str_starts(FIPS,"09")))$FIPS))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq(m_trt-1,m) + (j*m))
  }
  
  # Subset STAN output (Cali is first 5 counties)
  Mu_trt_ori_CT <- Mu_trt_ori[,ind_CT]
  Mu_trt_space_CT <- Mu_trt_space[,ind_CT]
  Mu_trt_spacetime_ICAR_CT <- Mu_trt_spacetime_ICAR[,ind_CT]
  Mu_trt_spacetime_AR_CT <- Mu_trt_spacetime_AR[,ind_CT]
  Mu_trt_space_comb_time_shrink_CT <- Mu_trt_space_comb_time_shrink[,ind_CT]
  Mu_trt_lasso_CT <- Mu_trt_lasso[,ind_CT]
  Mu_trt_gsc <- as.vector(t(fit_gsc$Y.ct)[(n_trt_CA+1):(n_trt_CA + n_trt_CT),(m_trt-1):m])
  
  
  # Subset Population
  treated2 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  pop_trt_CT <- treated2 %>% filter(YEAR_DX >= "1995")
  pop_trt_CT <- pop_trt_CT$POP
  
  # Calculate ATT
  # Vanilla Model
  ATT_ori_CT <- get_ATT(Mu_trt_ori_CT,Y1_CT,pop_trt_CT,n_trt_CT,m_trt,m_trt)
  
  # Space Model
  ATT_space_CT <- get_ATT(Mu_trt_space_CT,Y1_CT,pop_trt_CT,n_trt_CT,m_trt,m_trt)
  #ATT_space_CT <- NA
  
  # Space-Time Model ICAR
  ATT_spacetime_ICAR_CT <- get_ATT(Mu_trt_spacetime_ICAR_CT,Y1_CT,pop_trt_CT,n_trt_CT,m_trt,m_trt)
  
  # Space-Time Model AR
  ATT_spacetime_AR_CT <- get_ATT(Mu_trt_spacetime_AR_CT,Y1_CT,pop_trt_CT,n_trt_CT,m_trt,m_trt)
  
  # Space-Time Model AR
  ATT_space_comb_time_shrink_CT <- get_ATT(Mu_trt_space_comb_time_shrink_CT,Y1_CT,pop_trt_CT,n_trt_CT,m_trt,m_trt)
  #ATT_space_comb_time_shrink <- NA
  
  # Lasso Model
  ATT_lasso_CT <- get_ATT(Mu_trt_lasso_CT,Y1_CT,pop_trt_CT,n_trt_CT,m_trt,m_trt)
  #ATT_lasso_CT <- NA
  
  # Gsynth Model
  ATT_gsc <- get_ATT_gsc(Mu_trt_gsc,Y1_CT,pop_trt_CT,ind,n_trt_CT)
  gsc_bounds <- rbind(rep(0,m_trt),rep(0,m_trt))
  
  # Summarize ATT
  # Vanilla Model
  med_ori <- rowMedians(ATT_ori_CT,na.rm = TRUE)
  bounds_ori <- apply(ATT_ori_CT, 1, function(x) quantile(x, probs = c(.025,.975),na.rm = TRUE))
  
  # Space Model
  med_space <- rowMedians(ATT_space_CT,na.rm = TRUE)
  bounds_space <- apply(ATT_space_CT, 1, function(x) quantile(x, probs = c(.025,.975),na.rm = TRUE))
  #med_space <- rep(0,m_trt)
  #bounds_space <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  # Space-Time ICAR Model 
  med_spacetime_ICAR <- rowMedians(ATT_spacetime_ICAR_CT,na.rm = TRUE)
  bounds_spacetime_ICAR <- apply(ATT_spacetime_ICAR_CT, 1, function(x) quantile(x, probs = c(.025,.975),na.rm = TRUE))
  #med_spacetime_ICAR <- rep(0,m_trt)
  #bounds_spacetime_ICAR <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  # Space-Time AR Model 
  med_spacetime_AR <- rowMedians(ATT_spacetime_AR_CT,na.rm = TRUE)
  bounds_spacetime_AR <- apply(ATT_spacetime_AR_CT, 1, function(x) quantile(x, probs = c(.025,.975),na.rm = TRUE))
  
  # Space-Time Shrinkage Model 
  med_space_comb_time_shrink <- rowMedians(ATT_space_comb_time_shrink_CT)
  bounds_space_comb_time_shrink <- apply(ATT_space_comb_time_shrink_CT, 1, function(x) quantile(x, probs = c(.025,.975)))
  # for now since this hasnt finished running
  #med_space_comb_time_shrink <- rep(0,m_trt)
  #bounds_space_comb_time_shrink <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  # Lasso Model
  med_lasso <- rowMedians(ATT_lasso_CT, na.rm = TRUE)
  bounds_lasso <- apply(ATT_lasso_CT, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  #med_lasso <- rep(0,m_trt)
  #bounds_lasso <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  ## Plot
  model <- c(rep("Vanilla",m_trt),rep("Space",m_trt),rep("Space-Time \n ICAR",m_trt),rep("Space-Time \n AR",m_trt),rep("Space-Time \n Shrinkage",m_trt), rep("Bayesian \n Lasso", m_trt), rep("Gsynth",m_trt))
  year <- rep(c("1995","1996","1997","1998","1999","2000","2001","2002","2003"),7)
  ATT <- c(med_ori,med_space,med_spacetime_ICAR,med_spacetime_AR,med_space_comb_time_shrink,med_lasso,ATT_gsc)
  UB <- c(med_ori,med_space,med_spacetime_ICAR,med_spacetime_AR,med_space_comb_time_shrink,med_lasso,ATT_gsc)
  LB <-  c(bounds_ori[1,],bounds_space[1,],bounds_spacetime_ICAR[1,],bounds_spacetime_AR[1,],bounds_space_comb_time_shrink[1,], bounds_lasso[1,], gsc_bounds[1,])
  UB <-  c(bounds_ori[2,],bounds_space[2,],bounds_spacetime_ICAR[2,],bounds_spacetime_AR[2,],bounds_space_comb_time_shrink[2,],bounds_lasso[2,],gsc_bounds[2,])
  
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  ATT_df = ATT_df %>%
    mutate(Model = fct_relevel(Model, "Gsynth", "Vanilla", "Space", "Space-Time \n ICAR", "Space-Time \n AR", "Space-Time \n Shrinkage", "Bayesian \n Lasso"))
  
  plot_CT <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = Model))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year, fill = Model), alpha = 0.1)+
    geom_hline(yintercept=0)+
    #scale_colour_manual("",values="blue")+
    #scale_fill_manual("",values="grey12")
    #theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+
    facet_wrap(~Model, nrow = 1) + 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Estimated ATT for Connecticut")
  # Observed Median Rate ") + geom_vline(xintercept = 1995, color = 'red')
  
  return(list(plot_CA, plot_CT))
}

#-------------------------------------------------------------------------------------------------------


# Function for aggregate ATT
get_agg_ATT <- function(Mu_trt,Y1,pop_trt) {
  
  # Convert to rates
  Mu_trt_rate <- 100000 * sweep(Mu_trt, 2, pop_trt, "/" )
  Y1_rate <- 100000 * Y1/pop_trt
  
  # Get an ATT for each MCMC iteration
  ATT_vec <- rowMeans(-1*sweep(Mu_trt_rate, 2, Y1_rate, "-")) #ATT= Y1-Y0 = -1*(Y0-Y1)
  
  # Get median of MCMC ATT
  ATT <- median(ATT_vec)
  
  # Get 95% CI
  CI <- quantile(ATT_vec, probs = c(.025,.975),na.rm = TRUE)
  
  return(list(ATT_vec = ATT_vec, ATT = ATT, CI = c(CI)))
}


# Function for aggregate ATT for gsynth
get_agg_ATT_gsc <- function(Y0_est_rate,Y1,pop_trt,n_trt,fit_gsc) {
  
  # Convert to rates
  #Y0_est_rate <- as.vector(t(fit_gsc$Y.ct)[,(m_trt-1):m])
  
  # subset population to be treated counties at treated times
  #pop_trt_temp <- matrix(pop_trt,nrow = n_trt, byrow = T)[,(m-m_trt+1):m] #(13 x 7)
  pop_trt_temp <- matrix(pop_trt,nrow = n_trt, byrow = T)
    
  # convert Y1 into rate
  Y1_rate <- 100000 * Y1/as.vector(t(pop_trt_temp))
  
  
  # Get an ATT 
  ATT <- mean(as.vector(Y1_rate) - Y0_est_rate) # this matches gsynth
  
  
  # Get 95% CI (Note: I can only get CIs for all treated years)
  CI <- fit_gsc$est.avg[3:4]
  
  return(list(ATT = ATT, CI = c(CI)))
}

# Function for aggregate ATT for svt
get_agg_ATT_svt <- function(Y0_est_rate,Y1,pop_trt,n_trt) {
  
  # Convert to rates
  #Y0_est_rate <- as.vector(t(fit_gsc$Y.ct)[,(m_trt-1):m])
  
  # subset population to be treated counties at treated times
  #pop_trt_temp <- matrix(pop_trt,nrow = n_trt, byrow = T)[,(m-m_trt+1):m] #(13 x 7)
  pop_trt_temp <- matrix(pop_trt,nrow = n_trt, byrow = T)
  
  # convert Y1 into rate
  Y1_rate <- 100000 * Y1/as.vector(t(pop_trt_temp))
  
  
  # Get an ATT 
  ATT <- mean(as.vector(Y1_rate) - Y0_est_rate) 
  
  
  
  
  return(list(ATT = ATT))
}


#---------------------------------------------------------------------------------------------
# Function to create a table with the total ATT for CA and CT

ATT_CACT_table <- function(data_full_hisp, Mu_trt_ori, Mu_trt_space, Mu_trt_spacetime_ICAR, Mu_trt_spacetime_AR, 
                           Mu_trt_space_comb_time_shrink, Mu_trt_lasso, fit_gsc,trt_year){
  
  # Split CA and CT
  
  # Get CA
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- get_Y1_obs1(treated1)
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  # Find which indices are CA
  ind <- c()
  for(i in 0:(n_trt_CA-1)){
    ind <- append(ind,seq((m_trt-1),m) + (i*m))
  }
  
  # Subset CA
  Mu_trt_ori_CA <- Mu_trt_ori[,ind]
  Mu_trt_space_CA <- Mu_trt_space[,ind]
  Mu_trt_spacetime_ICAR_CA <- Mu_trt_spacetime_ICAR[,ind]
  Mu_trt_spacetime_AR_CA <- Mu_trt_spacetime_AR[,ind]
  Mu_trt_lasso_CA <- Mu_trt_lasso[,ind]
  #Mu_trt_space_comb_time_shrink_CA <- NA
  Mu_trt_space_comb_time_shrink_CA <- Mu_trt_space_comb_time_shrink[,ind]
  Mu_trt_gsc_CA <- as.vector(t(fit_gsc$Y.ct)[1:(length(unique(treated1$FIPS))),(m_trt-1):m])
  
  # Subset Population for CA
  pop_trt_CA <- treated1 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CA <- pop_trt_CA$POP
  
  # Get CT
  Y1_CT <- Y1[(length(Y1_CA)+1):length(Y1)]
  
  n_trt_CT <- length(unique(unique(data_full_hisp %>% filter(str_starts(FIPS,"09")))$FIPS))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq(m_trt-1,m) + (j*m))
  }
  
  # Subset STAN output (Cali is first 5 counties)
  Mu_trt_ori_CT <- Mu_trt_ori[,ind_CT]
  Mu_trt_space_CT <- Mu_trt_space[,ind_CT]
  Mu_trt_spacetime_ICAR_CT <- Mu_trt_spacetime_ICAR[,ind_CT]
  Mu_trt_spacetime_AR_CT <- Mu_trt_spacetime_AR[,ind_CT]
  #Mu_trt_space_comb_time_shrink_CT <- NA
  Mu_trt_space_comb_time_shrink_CT <- Mu_trt_space_comb_time_shrink[,ind_CT]
  Mu_trt_lasso_CT <- Mu_trt_lasso[,ind_CT]
  Mu_trt_gsc_CT <- as.vector(t(fit_gsc$Y.ct)[(n_trt_CA+1):(n_trt_CA + n_trt_CT),(m_trt-1):m])
  
  # Subset Population
  treated2 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  pop_trt_CT <- treated2 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CT <- pop_trt_CT$POP
  
  # California
  print("ori")
  agg_ATT_ori_CA_sub <- get_agg_ATT(Mu_trt_ori_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_ori_CA_sub <-  get_agg_ATT(Mu_trt_ori_CA,Y1_CA,pop_trt_CA)$CI
  
  print("space")
  agg_ATT_space_CA_sub <- get_agg_ATT(Mu_trt_space_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_space_CA_sub <-  get_agg_ATT(Mu_trt_space_CA,Y1_CA,pop_trt_CA)$CI
  
  print("ICAR")
  agg_ATT_ICAR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_ICAR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CA,Y1_CA,pop_trt_CA)$CI
  
  print("AR")
  agg_ATT_AR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_AR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CA,Y1_CA,pop_trt_CA)$CI
  
  print("shrink")
  #agg_ATT_shrink_CA_sub <- NA
  agg_ATT_shrink_CA_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_shrink_CA_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CA,Y1_CA,pop_trt_CA)$CI
  
  print("lasso")
  agg_ATT_lasso_CA_sub <- get_agg_ATT(Mu_trt_lasso_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_lasso_CA_sub <- get_agg_ATT(Mu_trt_lasso_CA,Y1_CA,pop_trt_CA)$CI
  
  agg_ATT_gsc_CA_sub <- get_agg_ATT_gsc(Mu_trt_gsc_CA,Y1_CA,pop_trt_CA,n_trt_CA,fit_gsc)$ATT
  print("to CT")
  # Connecticut
  agg_ATT_ori_CT_sub <- get_agg_ATT(Mu_trt_ori_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_ori_CT_sub <-  get_agg_ATT(Mu_trt_ori_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_space_CT_sub <- get_agg_ATT(Mu_trt_space_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_space_CT_sub <-  get_agg_ATT(Mu_trt_space_CT,Y1_CT,pop_trt_CT)$CI

  agg_ATT_ICAR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_ICAR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CT,Y1_CT,pop_trt_CT)$CI

  agg_ATT_AR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_AR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CT,Y1_CT,pop_trt_CT)$CI
  
  #agg_ATT_shrink_CT_sub <- NA
  agg_ATT_shrink_CT_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_shrink_CT_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_lasso_CT_sub <- get_agg_ATT(Mu_trt_lasso_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_lasso_CT_sub <- get_agg_ATT(Mu_trt_lasso_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_gsc_CT_sub <- get_agg_ATT_gsc(Mu_trt_gsc_CT,Y1_CT,pop_trt_CT,n_trt_CT,fit_gsc)$ATT
  
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI_ori_CA <- CI_fun(agg_bounds_ori_CA_sub)
  CI_space_CA <- CI_fun(agg_bounds_space_CA_sub)
  CI_ICAR_CA <- CI_fun(agg_bounds_ICAR_CA_sub)
  CI_AR_CA <- CI_fun(agg_bounds_AR_CA_sub)
  CI_shrink_CA <- CI_fun(agg_bounds_shrink_CA_sub)
  #CI_shrink_CA <- NA
  CI_lasso_CA <- CI_fun(agg_bounds_lasso_CA_sub)
  
  CI_ori_CT <- CI_fun(agg_bounds_ori_CT_sub)
  CI_space_CT <- CI_fun(agg_bounds_space_CT_sub)
  CI_ICAR_CT <- CI_fun(agg_bounds_ICAR_CT_sub)
  CI_AR_CT <- CI_fun(agg_bounds_AR_CT_sub)
  CI_shrink_CT <- CI_fun(agg_bounds_shrink_CT_sub)
  #CI_shrink_CT <- NA
  CI_lasso_CT <- CI_fun(agg_bounds_lasso_CT_sub)
  
  # Make kable table
  Model <- c("Gsynth",rep("Vanilla",1),rep("Space",1),rep("Space-Time \n ICAR",1),rep("Space-Time \n AR",1),rep("Space-Time \n Shrinkage",1), rep("Bayesian \n Lasso", 1))
  df <- cbind(Model, 
              round(c(agg_ATT_gsc_CA_sub,agg_ATT_ori_CA_sub,agg_ATT_space_CA_sub,agg_ATT_ICAR_CA_sub,
                      agg_ATT_AR_CA_sub,agg_ATT_shrink_CA_sub,agg_ATT_lasso_CA_sub), 2),
              c("NA",CI_ori_CA,CI_space_CA,CI_ICAR_CA,CI_AR_CA,CI_shrink_CA,CI_lasso_CA),
              round(c(agg_ATT_gsc_CT_sub,agg_ATT_ori_CT_sub,agg_ATT_space_CT_sub,agg_ATT_ICAR_CT_sub,
                      agg_ATT_AR_CT_sub,agg_ATT_shrink_CT_sub,agg_ATT_lasso_CT_sub), 2),
              c("NA",CI_ori_CT,CI_space_CT,CI_ICAR_CT,CI_AR_CT,CI_shrink_CT,CI_lasso_CT))
  
  rownames(df) <- NULL
  colnames(df) <- c("Model", "ATT", "95% CI", "ATT", "95% CI")
  #table <- kbl(df) %>%
  #  kable_classic() %>%
  #  add_header_above(c(" " , "California" = 2, "Connecticut" = 2))
  
  return(df)
}

#---------------------------------------------------------------------------------------------
# Function to create a table with the total ATT for CA and CT and Overall

ATT_CACT_overall_table <- function(data_full_hisp, Mu_trt_ori, Mu_trt_space, Mu_trt_spacetime_ICAR, Mu_trt_spacetime_AR, 
                           Mu_trt_space_comb_time_shrink, Mu_trt_lasso, fit_gsc, fit_svt, trt_year){
  library(kableExtra)
  # Split CA and CT
  
  # Get CA
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- get_Y1_obs1(treated1)
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  # Find which indices are CA
  ind <- c()
  for(i in 0:(n_trt_CA-1)){
    ind <- append(ind,seq((m_trt-1),m) + (i*m))
  }
  
  # Subset CA
  Mu_trt_ori_CA <- Mu_trt_ori[,ind]
  Mu_trt_space_CA <- Mu_trt_space[,ind]
  Mu_trt_spacetime_ICAR_CA <- Mu_trt_spacetime_ICAR[,ind]
  Mu_trt_spacetime_AR_CA <- Mu_trt_spacetime_AR[,ind]
  Mu_trt_lasso_CA <- Mu_trt_lasso[,ind]
  Mu_trt_space_comb_time_shrink_CA <- Mu_trt_space_comb_time_shrink[,ind]
  Mu_trt_gsc_CA <- as.vector(t(fit_gsc$Y.ct)[1:(length(unique(treated1$FIPS))),(m_trt-1):m])
  Mu_trt_svt_CA <-  as.vector(t(fit_svt$X[1:(length(unique(treated1$FIPS))),(m_trt-1):m]))
  
  # Subset Population for CA
  pop_trt_CA <- treated1 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CA <- pop_trt_CA$POP
  
  # Get CT
  Y1_CT <- Y1[(length(Y1_CA)+1):length(Y1)]
  
  n_trt_CT <- length(unique(unique(data_full_hisp %>% filter(str_starts(FIPS,"09")))$FIPS))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq(m_trt-1,m) + (j*m))
  }
  
  # Subset STAN output (Cali is first 5 counties)
  Mu_trt_ori_CT <- Mu_trt_ori[,ind_CT]
  Mu_trt_space_CT <- Mu_trt_space[,ind_CT]
  Mu_trt_spacetime_ICAR_CT <- Mu_trt_spacetime_ICAR[,ind_CT]
  Mu_trt_spacetime_AR_CT <- Mu_trt_spacetime_AR[,ind_CT]
  Mu_trt_space_comb_time_shrink_CT <- Mu_trt_space_comb_time_shrink[,ind_CT]
  Mu_trt_lasso_CT <- Mu_trt_lasso[,ind_CT]
  Mu_trt_gsc_CT <- as.vector(t(fit_gsc$Y.ct)[(n_trt_CA+1):(n_trt_CA + n_trt_CT),(m_trt-1):m])
  Mu_trt_svt_CT <-  as.vector(t(fit_svt$X[(n_trt_CA+1):(n_trt_CA + n_trt_CT),(m_trt-1):m]))
  
  # Subset Population
  treated2 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  pop_trt_CT <- treated2 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CT <- pop_trt_CT$POP
  
  # California
  
  agg_ATT_ori_CA_sub <- get_agg_ATT(Mu_trt_ori_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_ori_CA_sub <-  get_agg_ATT(Mu_trt_ori_CA,Y1_CA,pop_trt_CA)$CI
  
  
  agg_ATT_space_CA_sub <- get_agg_ATT(Mu_trt_space_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_space_CA_sub <-  get_agg_ATT(Mu_trt_space_CA,Y1_CA,pop_trt_CA)$CI
  
  
  agg_ATT_ICAR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_ICAR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CA,Y1_CA,pop_trt_CA)$CI
  
  
  agg_ATT_AR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_AR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CA,Y1_CA,pop_trt_CA)$CI

  agg_ATT_shrink_CA_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_shrink_CA_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CA,Y1_CA,pop_trt_CA)$CI
  
 
  agg_ATT_lasso_CA_sub <- get_agg_ATT(Mu_trt_lasso_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_lasso_CA_sub <- get_agg_ATT(Mu_trt_lasso_CA,Y1_CA,pop_trt_CA)$CI
  
  agg_ATT_gsc_CA_sub <- get_agg_ATT_gsc(Mu_trt_gsc_CA,Y1_CA,pop_trt_CA,n_trt_CA,fit_gsc)$ATT
  
  agg_ATT_svt_CA_sub <- get_agg_ATT_svt(Mu_trt_svt_CA,Y1_CA,pop_trt_CA,n_trt_CA)$ATT
  
  # Connecticut
  agg_ATT_ori_CT_sub <- get_agg_ATT(Mu_trt_ori_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_ori_CT_sub <-  get_agg_ATT(Mu_trt_ori_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_space_CT_sub <- get_agg_ATT(Mu_trt_space_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_space_CT_sub <-  get_agg_ATT(Mu_trt_space_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_ICAR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_ICAR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_AR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_AR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CT,Y1_CT,pop_trt_CT)$CI

  agg_ATT_shrink_CT_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_shrink_CT_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_lasso_CT_sub <- get_agg_ATT(Mu_trt_lasso_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_lasso_CT_sub <- get_agg_ATT(Mu_trt_lasso_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_gsc_CT_sub <- get_agg_ATT_gsc(Mu_trt_gsc_CT,Y1_CT,pop_trt_CT,n_trt_CT,fit_gsc)$ATT
  
  agg_ATT_svt_CT_sub <- get_agg_ATT_svt(Mu_trt_svt_CT,Y1_CT,pop_trt_CT,n_trt_CT)$ATT
  
    
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI_ori_CA <- CI_fun(agg_bounds_ori_CA_sub)
  CI_space_CA <- CI_fun(agg_bounds_space_CA_sub)
  CI_ICAR_CA <- CI_fun(agg_bounds_ICAR_CA_sub)
  CI_AR_CA <- CI_fun(agg_bounds_AR_CA_sub)
  CI_shrink_CA <- CI_fun(agg_bounds_shrink_CA_sub)
  CI_lasso_CA <- CI_fun(agg_bounds_lasso_CA_sub)
  
  CI_ori_CT <- CI_fun(agg_bounds_ori_CT_sub)
  CI_space_CT <- CI_fun(agg_bounds_space_CT_sub)
  CI_ICAR_CT <- CI_fun(agg_bounds_ICAR_CT_sub)
  CI_AR_CT <- CI_fun(agg_bounds_AR_CT_sub)
  CI_shrink_CT <- CI_fun(agg_bounds_shrink_CT_sub)
  CI_lasso_CT <- CI_fun(agg_bounds_lasso_CT_sub)
  
  # Overall ATT
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq((m_trt-1),m) + (i*m))
  }
  
  Mu_trt_gsc <- as.vector(t(fit_gsc$Y.ct)[,(m_trt-1):m])
  
  Mu_trt_svt <- as.vector(t(fit_svt$X[1:(n_trt_CA + n_trt_CT),(m_trt-1):m]))
  
  agg_ATT_ori_sub <- get_agg_ATT(Mu_trt_ori[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_ori_sub <-  get_agg_ATT(Mu_trt_ori[,ind],Y1,pop_trt[ind])$CI
  
  agg_ATT_space_sub <- get_agg_ATT(Mu_trt_space[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_space_sub <-  get_agg_ATT(Mu_trt_space[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_space_sub <- NA
  #agg_bounds_space_sub <- NA
  
  agg_ATT_ICAR_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_ICAR_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_ICAR_sub <- NA
  #agg_bounds_ICAR_sub <- NA
  
  agg_ATT_AR_sub <- get_agg_ATT(Mu_trt_spacetime_AR[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_AR_sub <- get_agg_ATT(Mu_trt_spacetime_AR[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_AR_sub <- NA
  #agg_bounds_AR_sub <- NA
  
  agg_ATT_shrink_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_shrink_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_shrink_sub <- NA
  #agg_bounds_shrink_sub <- NA
  
  agg_ATT_lasso_sub <- get_agg_ATT(Mu_trt_lasso[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_lasso_sub <- get_agg_ATT(Mu_trt_lasso[,ind],Y1,pop_trt[ind])$CI
  
  agg_ATT_gsc_sub <- get_agg_ATT_gsc(Mu_trt_gsc,Y1,pop_trt[ind],n_trt,fit_gsc)$ATT
  agg_bounds_gsc_sub <- get_agg_ATT_gsc(Mu_trt_gsc,Y1,pop_trt[ind],n_trt,fit_gsc)$CI
  
  agg_ATT_svt_sub <- get_agg_ATT_svt(Mu_trt_svt,Y1,pop_trt[ind],n_trt)$ATT
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI_ori <- CI_fun(agg_bounds_ori_sub)
  CI_space <- CI_fun(agg_bounds_space_sub)
  CI_ICAR <- CI_fun(agg_bounds_ICAR_sub)
  CI_AR <- CI_fun(agg_bounds_AR_sub)
  CI_shrink <- CI_fun(agg_bounds_shrink_sub)
  CI_lasso <- CI_fun(agg_bounds_lasso_sub)
  CI_gsc <- CI_fun(agg_bounds_gsc_sub)
  
  # Make kable table
  Model <- c("Gsynth",rep("Vanilla",1),rep("Space",1),rep("Space-Time \n ICAR",1),rep("Space-Time \n AR",1),rep("Space-Time \n Shrinkage",1), rep("Bayesian \n Lasso", 1))
  df <- cbind(Model, 
              round(c(agg_ATT_gsc_sub,agg_ATT_ori_sub,agg_ATT_space_sub,agg_ATT_ICAR_sub,
                      agg_ATT_AR_sub,agg_ATT_shrink_sub,agg_ATT_lasso_sub), 2),
              c(CI_gsc,CI_ori,CI_space,CI_ICAR,CI_AR,CI_shrink,CI_lasso))
  
  rownames(df) <- NULL
  colnames(df) <- c("Model", "ATT", "95% CI")
  
  
  # Make kable table
  Model <- c("Gsynth","SVT MC",rep("Vanilla",1),rep("Space",1),rep("Space-Time \n ICAR",1),rep("Space-Time \n AR",1), rep("Bayesian \n Lasso", 1),rep("Space-Time \n Shrinkage",1))
  df <- cbind(Model, 
              round(c(agg_ATT_gsc_CA_sub,agg_ATT_svt_CA_sub,agg_ATT_ori_CA_sub,agg_ATT_space_CA_sub,agg_ATT_ICAR_CA_sub,
                      agg_ATT_AR_CA_sub,agg_ATT_lasso_CA_sub,agg_ATT_shrink_CA_sub), 2),
              c("NA","NA",CI_ori_CA,CI_space_CA,CI_ICAR_CA,CI_AR_CA,CI_lasso_CA,CI_shrink_CA),
              round(c(agg_ATT_gsc_CT_sub,agg_ATT_svt_CT_sub,agg_ATT_ori_CT_sub,agg_ATT_space_CT_sub,agg_ATT_ICAR_CT_sub,
                      agg_ATT_AR_CT_sub,agg_ATT_lasso_CT_sub,agg_ATT_shrink_CT_sub), 2),
              c("NA","NA",CI_ori_CT,CI_space_CT,CI_ICAR_CT,CI_AR_CT,CI_lasso_CT,CI_shrink_CT), 
              round(c(agg_ATT_gsc_sub,agg_ATT_svt_sub,agg_ATT_ori_sub,agg_ATT_space_sub,agg_ATT_ICAR_sub,
                      agg_ATT_AR_sub,agg_ATT_lasso_sub,agg_ATT_shrink_sub), 2),
              c(CI_gsc,"NA",CI_ori,CI_space,CI_ICAR,CI_AR,CI_lasso,CI_shrink))
  
  rownames(df) <- NULL
  colnames(df) <- c("Model", "ATT", "95% CI", "ATT", "95% CI", "ATT", "95% CI")
  table <- kbl(df, format = "latex") %>%
    kable_classic() %>%
    add_header_above(c(" " , "California" = 2, "Connecticut" = 2, "Overall" = 2))
  
  return(table)
}

ATT_CACT_overall_df <- function(data_full_hisp, Mu_trt_ori, Mu_trt_space, Mu_trt_spacetime_ICAR, Mu_trt_spacetime_AR, 
                                   Mu_trt_space_comb_time_shrink, Mu_trt_lasso, fit_gsc,fit_svt,trt_year){
  library(kableExtra)
  # Split CA and CT
  
  # Get CA
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- get_Y1_obs1(treated1)
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  # Find which indices are CA
  ind <- c()
  for(i in 0:(n_trt_CA-1)){
    ind <- append(ind,seq((m_trt-1),m) + (i*m))
  }
  
  # Subset CA
  Mu_trt_ori_CA <- Mu_trt_ori[,ind]
  Mu_trt_space_CA <- Mu_trt_space[,ind]
  Mu_trt_spacetime_ICAR_CA <- Mu_trt_spacetime_ICAR[,ind]
  Mu_trt_spacetime_AR_CA <- Mu_trt_spacetime_AR[,ind]
  Mu_trt_lasso_CA <- Mu_trt_lasso[,ind]
  Mu_trt_space_comb_time_shrink_CA <- Mu_trt_space_comb_time_shrink[,ind]
  Mu_trt_gsc_CA <- as.vector(t(fit_gsc$Y.ct)[1:(length(unique(treated1$FIPS))),(m_trt-1):m])
  Mu_trt_svt_CA <-  as.vector(t(fit_svt$X[1:(length(unique(treated1$FIPS))),(m_trt-1):m]))
  
  # Subset Population for CA
  pop_trt_CA <- treated1 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CA <- pop_trt_CA$POP
  
  # Get CT
  Y1_CT <- Y1[(length(Y1_CA)+1):length(Y1)]
  
  n_trt_CT <- length(unique(unique(data_full_hisp %>% filter(str_starts(FIPS,"09")))$FIPS))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq(m_trt-1,m) + (j*m))
  }
  
  # Subset STAN output (Cali is first 5 counties)
  Mu_trt_ori_CT <- Mu_trt_ori[,ind_CT]
  Mu_trt_space_CT <- Mu_trt_space[,ind_CT]
  Mu_trt_spacetime_ICAR_CT <- Mu_trt_spacetime_ICAR[,ind_CT]
  Mu_trt_spacetime_AR_CT <- Mu_trt_spacetime_AR[,ind_CT]
  Mu_trt_space_comb_time_shrink_CT <- Mu_trt_space_comb_time_shrink[,ind_CT]
  Mu_trt_lasso_CT <- Mu_trt_lasso[,ind_CT]
  Mu_trt_gsc_CT <- as.vector(t(fit_gsc$Y.ct)[(n_trt_CA+1):(n_trt_CA + n_trt_CT),(m_trt-1):m])
  Mu_trt_svt_CT <-  as.vector(t(fit_svt$X[(n_trt_CA+1):(n_trt_CA + n_trt_CT),(m_trt-1):m]))
  
  # Subset Population
  treated2 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  pop_trt_CT <- treated2 %>% filter(YEAR_DX >= trt_year)
  pop_trt_CT <- pop_trt_CT$POP
  
  # California
  
  agg_ATT_ori_CA_sub <- get_agg_ATT(Mu_trt_ori_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_ori_CA_sub <-  get_agg_ATT(Mu_trt_ori_CA,Y1_CA,pop_trt_CA)$CI
  
  
  agg_ATT_space_CA_sub <- get_agg_ATT(Mu_trt_space_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_space_CA_sub <-  get_agg_ATT(Mu_trt_space_CA,Y1_CA,pop_trt_CA)$CI
  
  
  agg_ATT_ICAR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_ICAR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CA,Y1_CA,pop_trt_CA)$CI
  
  
  agg_ATT_AR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_AR_CA_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CA,Y1_CA,pop_trt_CA)$CI
  
  agg_ATT_shrink_CA_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_shrink_CA_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CA,Y1_CA,pop_trt_CA)$CI
  
  
  agg_ATT_lasso_CA_sub <- get_agg_ATT(Mu_trt_lasso_CA,Y1_CA,pop_trt_CA)$ATT
  agg_bounds_lasso_CA_sub <- get_agg_ATT(Mu_trt_lasso_CA,Y1_CA,pop_trt_CA)$CI
  
  agg_ATT_gsc_CA_sub <- get_agg_ATT_gsc(Mu_trt_gsc_CA,Y1_CA,pop_trt_CA,n_trt_CA,fit_gsc)$ATT
  
  agg_ATT_svt_CA_sub <- get_agg_ATT_svt(Mu_trt_svt_CA,Y1_CA,pop_trt_CA,n_trt_CA)$ATT
  
  # Connecticut
  agg_ATT_ori_CT_sub <- get_agg_ATT(Mu_trt_ori_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_ori_CT_sub <-  get_agg_ATT(Mu_trt_ori_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_space_CT_sub <- get_agg_ATT(Mu_trt_space_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_space_CT_sub <-  get_agg_ATT(Mu_trt_space_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_ICAR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_ICAR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_AR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_AR_CT_sub <- get_agg_ATT(Mu_trt_spacetime_AR_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_shrink_CT_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_shrink_CT_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_lasso_CT_sub <- get_agg_ATT(Mu_trt_lasso_CT,Y1_CT,pop_trt_CT)$ATT
  agg_bounds_lasso_CT_sub <- get_agg_ATT(Mu_trt_lasso_CT,Y1_CT,pop_trt_CT)$CI
  
  agg_ATT_gsc_CT_sub <- get_agg_ATT_gsc(Mu_trt_gsc_CT,Y1_CT,pop_trt_CT,n_trt_CT,fit_gsc)$ATT
  
  agg_ATT_svt_CT_sub <- get_agg_ATT_svt(Mu_trt_svt_CT,Y1_CT,pop_trt_CT,n_trt_CT)$ATT
  
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI_ori_CA <- CI_fun(agg_bounds_ori_CA_sub)
  CI_space_CA <- CI_fun(agg_bounds_space_CA_sub)
  CI_ICAR_CA <- CI_fun(agg_bounds_ICAR_CA_sub)
  CI_AR_CA <- CI_fun(agg_bounds_AR_CA_sub)
  CI_shrink_CA <- CI_fun(agg_bounds_shrink_CA_sub)
  CI_lasso_CA <- CI_fun(agg_bounds_lasso_CA_sub)
  
  CI_ori_CT <- CI_fun(agg_bounds_ori_CT_sub)
  CI_space_CT <- CI_fun(agg_bounds_space_CT_sub)
  CI_ICAR_CT <- CI_fun(agg_bounds_ICAR_CT_sub)
  CI_AR_CT <- CI_fun(agg_bounds_AR_CT_sub)
  CI_shrink_CT <- CI_fun(agg_bounds_shrink_CT_sub)
  CI_lasso_CT <- CI_fun(agg_bounds_lasso_CT_sub)
  
  # Overall ATT
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq((m_trt-1),m) + (i*m))
  }
  
  Mu_trt_gsc <- as.vector(t(fit_gsc$Y.ct)[,(m_trt-1):m])
  
  Mu_trt_svt <- as.vector(t(fit_svt$X[1:(n_trt_CA + n_trt_CT),(m_trt-1):m]))
  
  agg_ATT_ori_sub <- get_agg_ATT(Mu_trt_ori[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_ori_sub <-  get_agg_ATT(Mu_trt_ori[,ind],Y1,pop_trt[ind])$CI
  
  agg_ATT_space_sub <- get_agg_ATT(Mu_trt_space[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_space_sub <-  get_agg_ATT(Mu_trt_space[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_space_sub <- NA
  #agg_bounds_space_sub <- NA
  
  agg_ATT_ICAR_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_ICAR_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_ICAR_sub <- NA
  #agg_bounds_ICAR_sub <- NA
  
  agg_ATT_AR_sub <- get_agg_ATT(Mu_trt_spacetime_AR[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_AR_sub <- get_agg_ATT(Mu_trt_spacetime_AR[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_AR_sub <- NA
  #agg_bounds_AR_sub <- NA
  
  agg_ATT_shrink_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_shrink_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_shrink_sub <- NA
  #agg_bounds_shrink_sub <- NA
  
  agg_ATT_lasso_sub <- get_agg_ATT(Mu_trt_lasso[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_lasso_sub <- get_agg_ATT(Mu_trt_lasso[,ind],Y1,pop_trt[ind])$CI
  
  agg_ATT_gsc_sub <- get_agg_ATT_gsc(Mu_trt_gsc,Y1,pop_trt[ind],n_trt,fit_gsc)$ATT
  agg_bounds_gsc_sub <- get_agg_ATT_gsc(Mu_trt_gsc,Y1,pop_trt[ind],n_trt,fit_gsc)$CI
  
  agg_ATT_svt_sub <- get_agg_ATT_svt(Mu_trt_svt,Y1,pop_trt[ind],n_trt)$ATT
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI_ori <- CI_fun(agg_bounds_ori_sub)
  CI_space <- CI_fun(agg_bounds_space_sub)
  CI_ICAR <- CI_fun(agg_bounds_ICAR_sub)
  CI_AR <- CI_fun(agg_bounds_AR_sub)
  CI_shrink <- CI_fun(agg_bounds_shrink_sub)
  CI_lasso <- CI_fun(agg_bounds_lasso_sub)
  CI_gsc <- CI_fun(agg_bounds_gsc_sub)
  
  
  
  
  # Make kable table
  Model <- c("Gsynth","SVT MC",rep("Vanilla",1),rep("Space",1),rep("Space-Time \n ICAR",1),rep("Space-Time \n AR",1), rep("Bayesian \n Lasso", 1),rep("Space-Time \n Shrinkage",1))
  df <- cbind(Model, 
              round(c(agg_ATT_gsc_CA_sub,agg_ATT_svt_CA_sub,agg_ATT_ori_CA_sub,agg_ATT_space_CA_sub,agg_ATT_ICAR_CA_sub,
                      agg_ATT_AR_CA_sub,agg_ATT_lasso_CA_sub,agg_ATT_shrink_CA_sub), 2),
              c("NA","NA",CI_ori_CA,CI_space_CA,CI_ICAR_CA,CI_AR_CA,CI_lasso_CA,CI_shrink_CA),
              round(c(agg_ATT_gsc_CT_sub,agg_ATT_svt_CT_sub,agg_ATT_ori_CT_sub,agg_ATT_space_CT_sub,agg_ATT_ICAR_CT_sub,
                      agg_ATT_AR_CT_sub,agg_ATT_lasso_CT_sub,agg_ATT_shrink_CT_sub), 2),
              c("NA","NA",CI_ori_CT,CI_space_CT,CI_ICAR_CT,CI_AR_CT,CI_lasso_CT,CI_shrink_CT), 
              round(c(agg_ATT_gsc_sub,agg_ATT_svt_sub,agg_ATT_ori_sub,agg_ATT_space_sub,agg_ATT_ICAR_sub,
                      agg_ATT_AR_sub,agg_ATT_lasso_sub,agg_ATT_shrink_sub), 2),
              c(CI_gsc,"NA",CI_ori,CI_space,CI_ICAR,CI_AR,CI_lasso,CI_shrink))
  
  rownames(df) <- NULL
  colnames(df) <- c("Model", "ATT", "95% CI", "ATT", "95% CI", "ATT", "95% CI")
  
  return(df)
}
#---------------------------------------------------------------------------------------------
# Function to create a table with the total ATT for CA and CT combined


ATT_table <- function(data_full_hisp, Mu_trt_ori, Mu_trt_space, Mu_trt_spacetime_ICAR, Mu_trt_spacetime_AR, 
         Mu_trt_space_comb_time_shrink, Mu_trt_lasso, fit_gsc){
  # Get indices for treated times
  ind <- c()
  for(i in 0:(n_trt-1)){
    ind <- append(ind,seq((m_trt-1),m) + (i*m))
  }
  
  Mu_trt_gsc <- as.vector(t(fit_gsc$Y.ct)[,(m_trt-1):m])
  
  agg_ATT_ori_sub <- get_agg_ATT(Mu_trt_ori[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_ori_sub <-  get_agg_ATT(Mu_trt_ori[,ind],Y1,pop_trt[ind])$CI
  
  agg_ATT_space_sub <- get_agg_ATT(Mu_trt_space[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_space_sub <-  get_agg_ATT(Mu_trt_space[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_space_sub <- NA
  #agg_bounds_space_sub <- NA
  
  agg_ATT_ICAR_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_ICAR_sub <- get_agg_ATT(Mu_trt_spacetime_ICAR[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_ICAR_sub <- NA
  #agg_bounds_ICAR_sub <- NA
  
  agg_ATT_AR_sub <- get_agg_ATT(Mu_trt_spacetime_AR[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_AR_sub <- get_agg_ATT(Mu_trt_spacetime_AR[,ind],Y1,pop_trt[ind])$CI
  
  agg_ATT_shrink_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_shrink_sub <- get_agg_ATT(Mu_trt_space_comb_time_shrink[,ind],Y1,pop_trt[ind])$CI
  #agg_ATT_shrink_sub <- NA
  #agg_bounds_shrink_sub <- NA
  
  agg_ATT_lasso_sub <- get_agg_ATT(Mu_trt_lasso[,ind],Y1,pop_trt[ind])$ATT
  agg_bounds_lasso_sub <- get_agg_ATT(Mu_trt_lasso[,ind],Y1,pop_trt[ind])$CI
  
  agg_ATT_gsc_sub <- get_agg_ATT_gsc(Mu_trt_gsc,Y1,pop_trt[ind],n_trt,fit_gsc)$ATT
  agg_bounds_gsc_sub <- get_agg_ATT_gsc(Mu_trt_gsc,Y1,pop_trt[ind],n_trt,fit_gsc)$CI
  
  CI_fun <- function(bounds){
    out <- paste0("(",round(bounds[1],2)," , ",round(bounds[2],2),")")
    return(out)
  }
  
  CI_ori <- CI_fun(agg_bounds_ori_sub)
  CI_space <- CI_fun(agg_bounds_space_sub)
  CI_ICAR <- CI_fun(agg_bounds_ICAR_sub)
  CI_AR <- CI_fun(agg_bounds_AR_sub)
  CI_shrink <- CI_fun(agg_bounds_shrink_sub)
  CI_lasso <- CI_fun(agg_bounds_lasso_sub)
  CI_gsc <- CI_fun(agg_bounds_gsc_sub)
  
  # Make kable table
  Model <- c("Gsynth",rep("Vanilla",1),rep("Space",1),rep("Space-Time \n ICAR",1),rep("Space-Time \n AR",1),rep("Space-Time \n Shrinkage",1), rep("Bayesian \n Lasso", 1))
  df <- cbind(Model, 
              round(c(agg_ATT_gsc_sub,agg_ATT_ori_sub,agg_ATT_space_sub,agg_ATT_ICAR_sub,
                      agg_ATT_AR_sub,agg_ATT_shrink_sub,agg_ATT_lasso_sub), 2),
              c(CI_gsc,CI_ori,CI_space,CI_ICAR,CI_AR,CI_shrink,CI_lasso))
  
  rownames(df) <- NULL
  colnames(df) <- c("Model", "ATT", "95% CI")
  
  table <- kbl(df) %>%
    kable_classic() 
  
  return(table)
}
#---------------------------------------------------------------------------------------------

# Function to plot Overall ATT over time
ATT_plot <- function(data_full_hisp, Mu_trt_ori, Mu_trt_space, Mu_trt_spacetime_ICAR, Mu_trt_spacetime_AR, 
                     Mu_trt_space_comb_time_shrink, Mu_trt_lasso, fit_gsc){
  # Calculate ATT
  # Vanilla Model
  ATT_ori <- get_ATT(Mu_trt_ori,Y1,pop_trt,n_trt,m_trt,m_trt)
  
  # Space Model
  ATT_space <- get_ATT(Mu_trt_space,Y1,pop_trt,n_trt,m_trt,m_trt)
  #ATT_space <- NA
  
  # Space-Time Model ICAR
  ATT_spacetime_ICAR <- get_ATT(Mu_trt_spacetime_ICAR,Y1,pop_trt,n_trt,m_trt,m_trt)
  #ATT_spacetime_ICAR <- NA
  
  # Space-Time Model AR
  ATT_spacetime_AR <- get_ATT(Mu_trt_spacetime_AR,Y1,pop_trt,n_trt,m_trt,m_trt)
  
  # Space-Time Shrink Model
  ATT_space_comb_time_shrink <- get_ATT(Mu_trt_space_comb_time_shrink,Y1,pop_trt,n_trt,m_trt,m_trt)
  #ATT_space_comb_time_shrink <- NA
  
  # Lasso Model
  ATT_lasso <- get_ATT(Mu_trt_lasso,Y1,pop_trt,n_trt,m_trt,m_trt)
  #ATT_lasso <- NA
  
  # Gsynth Model
  ATT_gsc <- fit_gsc$est.att[8:m,1]
  #ATT_gsc <- fit_gsc$att[(m-m_trt+1):m] # just treated years
  #gsc_bounds <- fit_gsc$est.att[(m-m_trt+1):m,c("CI.lower","CI.upper")] # inference for att
  gsc_bounds <- fit_gsc$est.att[8:m,3:4]
  
  # Summarize ATT
  # Vanilla Model
  med_ori <- rowMedians(ATT_ori,na.rm = TRUE)
  bounds_ori <- apply(ATT_ori, 1, function(x) quantile(x, probs = c(.025,.975),na.rm = TRUE))
  
  # Space Model
  med_space <- rowMedians(ATT_space,na.rm = TRUE)
  bounds_space <- apply(ATT_space, 1, function(x) quantile(x, probs = c(.025,.975),na.rm = TRUE))
  #med_space <- rep(0,m_trt)
  #bounds_space <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  # Space-Time ICAR Model 
  med_spacetime_ICAR <- rowMedians(ATT_spacetime_ICAR,na.rm = TRUE)
  bounds_spacetime_ICAR <- apply(ATT_spacetime_ICAR, 1, function(x) quantile(x, probs = c(.025,.975),na.rm = TRUE))
  #med_spacetime_ICAR <- rep(0,m_trt)
  #bounds_spacetime_ICAR <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  # Space-Time AR Model 
  med_spacetime_AR <- rowMedians(ATT_spacetime_AR,na.rm = TRUE)
  bounds_spacetime_AR <- apply(ATT_spacetime_AR, 1, function(x) quantile(x, probs = c(.025,.975),na.rm = TRUE))
  
  # Space-Time Shrinkage Model 
  med_space_comb_time_shrink <- rowMedians(ATT_space_comb_time_shrink)
  bounds_space_comb_time_shrink <- apply(ATT_space_comb_time_shrink, 1, function(x) quantile(x, probs = c(.025,.975)))
  # for now since this hasnt finished running
  #med_space_comb_time_shrink <- rep(0,m_trt)
  #bounds_space_comb_time_shrink <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  # Lasso Model
  med_lasso <- rowMedians(ATT_lasso, na.rm = TRUE)
  bounds_lasso <- apply(ATT_lasso, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  #med_lasso <- rep(0,m_trt)
  #bounds_lasso <- rbind(rep(0,m_trt),rep(0,m_trt)) 
  
  ## Plot
  model <- c(rep("Vanilla",m_trt),rep("Space",m_trt),rep("Space-Time \n ICAR",m_trt),rep("Space-Time \n AR",m_trt),rep("Space-Time \n Shrinkage",m_trt), rep("Bayesian \n Lasso", m_trt), rep("Gsynth",m_trt))
  year <- rep(c("1995","1996","1997","1998","1999","2000","2001","2002","2003"),7)
  ATT <- c(med_ori,med_space,med_spacetime_ICAR,med_spacetime_AR,med_space_comb_time_shrink,med_lasso,ATT_gsc)
  UB <- c(med_ori,med_space,med_spacetime_ICAR,med_spacetime_AR,med_space_comb_time_shrink,med_lasso,ATT_gsc)
  LB <-  c(bounds_ori[1,],bounds_space[1,],bounds_spacetime_ICAR[1,],bounds_spacetime_AR[1,],bounds_space_comb_time_shrink[1,], bounds_lasso[1,], gsc_bounds[,1])
  UB <-  c(bounds_ori[2,],bounds_space[2,],bounds_spacetime_ICAR[2,],bounds_spacetime_AR[2,],bounds_space_comb_time_shrink[2,],bounds_lasso[2,],gsc_bounds[,2])
  
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  ATT_df = ATT_df %>%
    mutate(Model = fct_relevel(Model, "Gsynth", "Vanilla", "Space", "Space-Time \n ICAR", "Space-Time \n AR", "Space-Time \n Shrinkage", "Bayesian \n Lasso"))
  
  plot <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = Model))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year, fill = Model), alpha = 0.1)+
    geom_hline(yintercept=0)+
    #scale_colour_manual("",values="blue")+
    #scale_fill_manual("",values="grey12")
    #theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+
    facet_wrap(~Model, nrow = 1) + 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Estimated ATT")
  # Observed Median Rate ") + geom_vline(xintercept = 1995, color = 'red')
  
  return(plot)
}

###################################################################################

ATT_full_plot <- function(data_full_hisp, Mu_trt, Model, trt_year, years){
  # Full ATT plot for CA
  
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")))
  
  Y1_CA <- treated1$CL_CASES
  
  n_trt_CA <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  ind <- c()
  for(i in 0:(n_trt_CA-1)){
    ind <- append(ind,seq((1),m) + (i*m))
  }
  
  Mu_trt_CA <- Mu_trt[,ind]
  
  # Subset Population
  pop_trt_CA <- treated1$POP
  
  
  # Calculate ATT for all MCMC simulations
  ATT_CA <- get_ATT_full(Mu_trt_CA,Y1_CA,pop_trt_CA,n_trt_CA,m_trt,m)
  
  # Summarize ATT
  med <- rowMedians(ATT_CA, na.rm = TRUE)
  bounds <- apply(ATT_CA, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  ## Plot
  model <- rep(Model,m)
  #year <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003")
  year <- years
  LB <-  c(bounds[1,])
  UB <-  c(bounds[2,])
  ATT <- c(med)
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  
  
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  plot_CA <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = "red"))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year), alpha = 0.1)+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=trt_year)+
    #scale_colour_manual(ATT,values="blue")+
    #scale_fill_manual("",values="grey12")
    theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+ 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Estimated ATT for California")
  
  plot_CA
  
  # Full ATT plot for CT
  
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "09")))
  
  Y1_CT <- treated1$CL_CASES
  
  n_trt_CT <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  ind_CT <- c()
  for(j in (i+1):(n_trt_CT+(i+1)-1)){
    ind_CT <- append(ind_CT,seq(1,m) + (j*m))
  }
  
  Mu_trt_CT <- Mu_trt[,ind_CT]
  
  # Subset Population
  pop_trt_CT <- treated1$POP
  
  
  # Calculate ATT for all MCMC simulations
  ATT_CT <- get_ATT_full(Mu_trt_CT,Y1_CT,pop_trt_CT,n_trt_CT,m_trt,m)
  
  # Summarize ATT
  med <- rowMedians(ATT_CT, na.rm = TRUE)
  bounds <- apply(ATT_CT, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  ## Plot
  model <- rep(Model,m)
  #year <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003")
  year <- years
  LB <-  c(bounds[1,])
  UB <-  c(bounds[2,])
  ATT <- c(med)
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  
  
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  plot_CT <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = "red"))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year), alpha = 0.1)+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=trt_year)+
    #scale_colour_manual(ATT,values="blue")+
    #scale_fill_manual("",values="grey12")
    theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+ 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Estimated ATT for Connecticut")
  
  plot_CT
  
  
  # Full ATT plot for both states
  
  treated1 <- data_full_hisp %>% filter((str_starts(FIPS, "06")) | (str_starts(FIPS, "09")))
  
  Y1 <- treated1$CL_CASES
  
  n_trt <- length(unique(treated1$FIPS))
  
  m <- length(unique(treated1$YEAR_DX))
  
  
  
  Mu_trt_CACT <- Mu_trt[,c(ind,ind_CT)]
  
  # Subset Population
  pop_trt <- treated1$POP
  
  
  # Calculate ATT for all MCMC simulations
  ATT_CACT <- get_ATT_full(Mu_trt_CACT,Y1,pop_trt,n_trt,m_trt,m)
  
  # Summarize ATT
  med <- rowMedians(ATT_CACT, na.rm = TRUE)
  bounds <- apply(ATT_CACT, 1, function(x) quantile(x, probs = c(.025,.975), na.rm = TRUE))
  
  ## Plot
  model <- rep(Model,m)
  #year <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001","2002","2003")
  year <- years
  LB <-  c(bounds[1,])
  UB <-  c(bounds[2,])
  ATT <- c(med)
  ATT_df <- tibble(model,as.numeric(year),as.numeric(round(ATT,3)),as.numeric(round(LB,3)),as.numeric(round(UB,3)))
  
  
  colnames(ATT_df) <- c("Model", "Year", "ATT", "LB","UB")
  
  plot_CACT <- ggplot(ATT_df) + geom_line(aes(y=ATT, x=Year, color = "red"))+
    geom_ribbon(aes(ymin=LB, ymax=UB, x=Year), alpha = 0.1)+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=trt_year)+
    #scale_colour_manual(ATT,values="blue")+
    #scale_fill_manual("",values="grey12")
    theme_bw()+
    theme(text=element_text(colour='black',size=14),
          axis.text.x = element_text(angle = 60, vjust = 0.9, hjust=1),
          legend.position = 'none')+ 
    xlab('Year')+ylab('ATT and 95% CI') +
    ggtitle("Overall Estimated ATT")
  
  plot_CACT
  
  library(cowplot)
  plot <-plot_grid(plot_CA, plot_CT, plot_CACT,nrow = 1)
  
  return(plot)
}
#########################################################################################
## CV Function for choosing lambda for SVT
## Inputs:
## - data_full_hisp: the dataset
## - k: number of folds for CV
## Outputs:
## - best_lambda

svt_cv <- function(data_full_hisp,k=10) {
  # Set up data matrix for fill.svt
  # Convert cases into rates 
  data_full_hisp$rate <- 100000 * (data_full_hisp$CL_CASES/data_full_hisp$POP)

  # Put NAs where treated units are
  Y0_miss <- data_full_hisp$rate
  Y0_miss[which(data_full_hisp$C==1)] <- NA
  Y0_miss <- matrix(Y0_miss, ncol = m, byrow=T)
  
  missing_data <- Y0_miss
  
  # Define the range of lambda values to consider
  lambda_values <- seq(0.1, 2, by = 0.1)
  
  # Create an empty vector to store cross-validation error
  cv_error <- numeric(length(lambda_values))
  
  # Get total elements
  n = nrow(missing_data) * ncol(missing_data)
  
  # Indicator matrix of missingness
  missing.matrix = is.na(missing_data)
  
  # Indices is all observed values
  obs.idx = which(as.vector(!missing.matrix))
  # Get observed set
  obs.data = as.vector(missing_data)[obs.idx] # returns vector
  
  ## Create the k folds of the data so there is no overlap ##
  library(caret)
  
  set.seed(501)
  folds.idx <- createFolds(obs.idx, k = 10, list = TRUE, returnTrain = FALSE) # returns a list of indices corresponding to obs.idx
  
  for (i in 1:k) {
    
    # Identify indices to be removed for validation
    remove.indices <- obs.idx[folds.idx[[i]]]
    
    # Select training data (put NAs where validation set is for testing)
    train_data = as.vector(missing_data)
    train_data[remove.indices] = NA
    
    # Put in matrix form
    train_data = matrix(train_data, nrow = nrow(missing_data), byrow = F)
    
    for (j in 1:length(lambda_values)) {
      # Fill missing values using SVT with current lambda value
      filled_data <- fill.SVT(train_data, lambda = lambda_values[j])$X
      
      # Calculate error on validation set (sum of squared errors)
      cv_error[j] <- cv_error[j] + sum(as.vector(filled_data)[remove.indices] - as.vector(missing_data)[remove.indices])^2
      
    }
    
  }
  
  # Average the errors across folds
  cv_error <- cv_error / k
  
  # Choose lambda with minimum cross-validation error
  best_lambda <- lambda_values[which.min(cv_error)]
  
  # Return best_lambda
  return(best_lambda)
}

