# This script takes the simulation results, calculates absolute percent biases
# of the ATT and counterfactual, and prepares them to inputted into a table.
# All absolute percent biases are calculated in the same format as ASC from
# Ben-Michael et al.

#############################################################################################################
## 0. Read command line arguments and load packages and functions                                          ##
##    Arguments should include:                                                                            ##
##    - rare: Rare or nonRare                                                                              ##
##    - smooth: smooth or nonSmooth                                                                        ##
##    - dir_in: directory input path for results from script 2                                             ##
##    - dir_out: directory output path for results                                                         ##
#############################################################################################################

# Read command line arguments 
args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }

# Load Packages
library(tidyverse)
library(matrixStats)

# load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Functions.R')
# load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/post_procesing_functions.R')

# load NM set up script
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/0_set_up.R') # We are basing the simulated data on the values from this script

#############################
## 1. Load and Set Up Data ##
#############################

# Load data
if(smooth == "smooth"){
  
  load(paste0(dir_in,"/sim_k7_",rare,"_smooth.RData"))
  res_7 <- data_big #100 sim
  
  load(paste0(dir_in,"/sim_k3_",rare,"_smooth.RData"))
  res_3 <- data_big #100 sim
  
  load(paste0(dir_in,"/sim_k1_",rare,"_smooth.RData"))
  res_1 <- data_big #100 sim
  
}else{
  load(paste0(dir_in,"/sim_k7_",rare,".RData"))
  res_7 <- data_big #100 sim
  
  load(paste0(dir_in,"/sim_k3_",rare,".RData"))
  res_3 <- data_big #100 sim
  
  load(paste0(dir_in,"/sim_k1_",rare,".RData"))
  res_1 <- data_big #100 sim
}


#remove data_big
rm(data_big)

# Set Up Data 
res_list_1 <- setup_sim_data(res_1)
list2env(setNames(res_list_1,c("Y0_matrix_1","Y1_matrix_1","Mu_trt_ori_1","Mu_trt_space_1",
                               "Mu_trt_spacetime_ICAR_1","Mu_trt_spacetime_AR_1", 
                               "Mu_trt_space_comb_time_shrink_1", "Mu_trt_lasso_1", 
                               "Y_pred_ori_1","Y_pred_space_1",
                               "Y_pred_spacetime_ICAR_1","Y_pred_spacetime_AR_1", 
                               "Y_pred_space_comb_time_shrink_1", "Y_pred_lasso_1",
                               "Y1_obs_rate_1","Y0_full_rate_1")), envir = environment()) 


res_list_3 <- setup_sim_data(res_3)
list2env(setNames(res_list_3,c("Y0_matrix_3","Y1_matrix_3","Mu_trt_ori_3","Mu_trt_space_3",
                               "Mu_trt_spacetime_ICAR_3","Mu_trt_spacetime_AR_3", 
                               "Mu_trt_space_comb_time_shrink_3", "Mu_trt_lasso_3", 
                               "Y_pred_ori_3","Y_pred_space_3",
                               "Y_pred_spacetime_ICAR_3","Y_pred_spacetime_AR_3", 
                               "Y_pred_space_comb_time_shrink_3", "Y_pred_lasso_3",
                               "Y1_obs_rate_3","Y0_full_rate_3")), envir = environment()) 

res_list_7 <- setup_sim_data(res_7)
list2env(setNames(res_list_7,c("Y0_matrix_7","Y1_matrix_7","Mu_trt_ori_7","Mu_trt_space_7",
                               "Mu_trt_spacetime_ICAR_7","Mu_trt_spacetime_AR_7", 
                               "Mu_trt_space_comb_time_shrink_7", "Mu_trt_lasso_7", 
                               "Y_pred_ori_7","Y_pred_space_7",
                               "Y_pred_spacetime_ICAR_7","Y_pred_spacetime_AR_7", 
                               "Y_pred_space_comb_time_shrink_7", "Y_pred_lasso_7",
                               "Y1_obs_rate_7","Y0_full_rate_7")), envir = environment()) 



# Get treated indices
ind <- c()
for(i in 0:5){
  ind <- append(ind,seq(9,15) + (i*15))
}

# Save true treatment effect
true_trt_effect_1 <- Y1_matrix_1 - Y0_matrix_1[ind,]
true_trt_effect_3 <- Y1_matrix_3 - Y0_matrix_3[ind,]
true_trt_effect_7 <- Y1_matrix_7 - Y0_matrix_7[ind,]
n_trt <- 6 # how we simulated the data

#get populations for treated counties
pop_trt <- as.vector(t(pop_mat[1:n_trt,]))

################################
## 3. Get convergence vectors ##
################################


# Check convergence 

conv <- conv_vec(Mu_trt_ori_1,Mu_trt_space_1,Mu_trt_spacetime_ICAR_1,
           Mu_trt_spacetime_AR_1,Mu_trt_lasso_1,
           Mu_trt_ori_3,Mu_trt_space_3,Mu_trt_spacetime_ICAR_3,
           Mu_trt_spacetime_AR_3,Mu_trt_lasso_3,
           Mu_trt_ori_7,Mu_trt_space_7,Mu_trt_spacetime_ICAR_7,
           Mu_trt_spacetime_AR_7,Mu_trt_lasso_7,
           Mu_trt_space_comb_time_shrink_7)


# Save convergence vector
if(smooth == "smooth"){
  filename_out = paste0(dir_out,"/conv_",rare,"_smooth.RData")
}else{
  filename_out = paste0(dir_out,"/conv_",rare,".RData")
}

print(filename_out)

save(conv, file = filename_out)





################################################################
## 4. Get Absolute Percent Bias of the Counterfactual vectors ##
################################################################

# Check if functions works
#medianAbsPercBiasATT_year(Mu_trt_space_comb_time_shrink_3,Y0_matrix_1,Y1_matrix_1, true_trt_effect_1, ind)

#medianAbsPercBiasY0_year(Mu_trt_lasso_1,Y0_matrix_1,ind)

# Mu_trt (not used in paper)

# Absolute Percent Bias of the Counterfactual
Y0_ori_1 <- medianAbsPercBiasY0_year(Mu_trt_ori_1,Y0_matrix_1,ind)$median
Y0_space_1 <- medianAbsPercBiasY0_year(Mu_trt_space_1,Y0_matrix_1,ind)$median
Y0_ICAR_1 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_ICAR_1,Y0_matrix_1,ind)$median
Y0_AR_1 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_AR_1,Y0_matrix_1,ind)$median
Y0_lasso_1 <- medianAbsPercBiasY0_year(Mu_trt_lasso_1,Y0_matrix_1,ind)$median

Y0_ori_3 <- medianAbsPercBiasY0_year(Mu_trt_ori_3,Y0_matrix_3,ind)$median
Y0_space_3 <- medianAbsPercBiasY0_year(Mu_trt_space_3,Y0_matrix_3,ind)$median
Y0_ICAR_3 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_ICAR_3,Y0_matrix_3,ind)$median
Y0_AR_3 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_AR_3,Y0_matrix_3,ind)$median
Y0_lasso_3 <- medianAbsPercBiasY0_year(Mu_trt_lasso_3,Y0_matrix_3,ind)$median

Y0_ori_7 <- medianAbsPercBiasY0_year(Mu_trt_ori_7,Y0_matrix_7,ind)$median
Y0_space_7 <- medianAbsPercBiasY0_year(Mu_trt_space_7,Y0_matrix_7,ind)$median
Y0_ICAR_7 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_ICAR_7,Y0_matrix_7,ind)$median
Y0_AR_7 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_AR_7,Y0_matrix_7,ind)$median
Y0_lasso_7 <- medianAbsPercBiasY0_year(Mu_trt_lasso_7,Y0_matrix_7,ind)$median
Y0_shrink_7 <- medianAbsPercBiasY0_year(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7,ind)$median

# Absolute Percent Bias of the ATT
ATT_ori_1 <- medianAbsPercBiasATT_year(Mu_trt_ori_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median
ATT_space_1 <- medianAbsPercBiasATT_year(Mu_trt_space_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median
ATT_ICAR_1 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_ICAR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median
ATT_AR_1 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_AR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median
ATT_lasso_1 <- medianAbsPercBiasATT_year(Mu_trt_lasso_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median

ATT_ori_3 <- medianAbsPercBiasATT_year(Mu_trt_ori_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median
ATT_space_3 <- medianAbsPercBiasATT_year(Mu_trt_space_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median
ATT_ICAR_3 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_ICAR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median
ATT_AR_3 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_AR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median
ATT_lasso_3 <- medianAbsPercBiasATT_year(Mu_trt_lasso_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median

ATT_ori_7 <- medianAbsPercBiasATT_year(Mu_trt_ori_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_space_7 <- medianAbsPercBiasATT_year(Mu_trt_space_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_ICAR_7 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_ICAR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_AR_7 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_AR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_lasso_7 <- medianAbsPercBiasATT_year(Mu_trt_lasso_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_shrink_7 <- medianAbsPercBiasATT_year(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median

# Upper bound
UB_Y0_ori_1 <- medianAbsPercBiasY0_year(Mu_trt_ori_1,Y0_matrix_1,ind)$UB
UB_Y0_space_1 <- medianAbsPercBiasY0_year(Mu_trt_space_1,Y0_matrix_1,ind)$UB
UB_Y0_ICAR_1 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_ICAR_1,Y0_matrix_1,ind)$UB
UB_Y0_AR_1 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_AR_1,Y0_matrix_1,ind)$UB
UB_Y0_lasso_1 <- medianAbsPercBiasY0_year(Mu_trt_lasso_1,Y0_matrix_1,ind)$UB

UB_Y0_ori_3 <- medianAbsPercBiasY0_year(Mu_trt_ori_3,Y0_matrix_3,ind)$UB
UB_Y0_space_3 <- medianAbsPercBiasY0_year(Mu_trt_space_3,Y0_matrix_3,ind)$UB
UB_Y0_ICAR_3 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_ICAR_3,Y0_matrix_3,ind)$UB
UB_Y0_AR_3 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_AR_3,Y0_matrix_3,ind)$UB
UB_Y0_lasso_3 <- medianAbsPercBiasY0_year(Mu_trt_lasso_3,Y0_matrix_3,ind)$UB

UB_Y0_ori_7 <- medianAbsPercBiasY0_year(Mu_trt_ori_7,Y0_matrix_7,ind)$UB
UB_Y0_space_7 <- medianAbsPercBiasY0_year(Mu_trt_space_7,Y0_matrix_7,ind)$UB
UB_Y0_ICAR_7 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_ICAR_7,Y0_matrix_7,ind)$UB
UB_Y0_AR_7 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_AR_7,Y0_matrix_7,ind)$UB
UB_Y0_lasso_7 <- medianAbsPercBiasY0_year(Mu_trt_lasso_7,Y0_matrix_7,ind)$UB
UB_Y0_shrink_7 <- medianAbsPercBiasY0_year(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7,ind)$UB

# Upper bound
UB_ATT_ori_1 <- medianAbsPercBiasATT_year(Mu_trt_ori_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB
UB_ATT_space_1 <- medianAbsPercBiasATT_year(Mu_trt_space_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB
UB_ATT_ICAR_1 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_ICAR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB
UB_ATT_AR_1 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_AR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB
UB_ATT_lasso_1 <- medianAbsPercBiasATT_year(Mu_trt_lasso_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB

UB_ATT_ori_3 <- medianAbsPercBiasATT_year(Mu_trt_ori_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB
UB_ATT_space_3 <- medianAbsPercBiasATT_year(Mu_trt_space_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB
UB_ATT_ICAR_3 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_ICAR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB
UB_ATT_AR_3 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_AR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB
UB_ATT_lasso_3 <- medianAbsPercBiasATT_year(Mu_trt_lasso_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB

UB_ATT_ori_7 <- medianAbsPercBiasATT_year(Mu_trt_ori_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_space_7 <- medianAbsPercBiasATT_year(Mu_trt_space_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_ICAR_7 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_ICAR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_AR_7 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_AR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_lasso_7 <- medianAbsPercBiasATT_year(Mu_trt_lasso_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_shrink_7 <- medianAbsPercBiasATT_year(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB

# Lower bound
LB_Y0_ori_1 <- medianAbsPercBiasY0_year(Mu_trt_ori_1,Y0_matrix_1,ind)$LB
LB_Y0_space_1 <- medianAbsPercBiasY0_year(Mu_trt_space_1,Y0_matrix_1,ind)$LB
LB_Y0_ICAR_1 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_ICAR_1,Y0_matrix_1,ind)$LB
LB_Y0_AR_1 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_AR_1,Y0_matrix_1,ind)$LB
LB_Y0_lasso_1 <- medianAbsPercBiasY0_year(Mu_trt_lasso_1,Y0_matrix_1,ind)$LB

LB_Y0_ori_3 <- medianAbsPercBiasY0_year(Mu_trt_ori_3,Y0_matrix_3,ind)$LB
LB_Y0_space_3 <- medianAbsPercBiasY0_year(Mu_trt_space_3,Y0_matrix_3,ind)$LB
LB_Y0_ICAR_3 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_ICAR_3,Y0_matrix_3,ind)$LB
LB_Y0_AR_3 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_AR_3,Y0_matrix_3,ind)$LB
LB_Y0_lasso_3 <- medianAbsPercBiasY0_year(Mu_trt_lasso_3,Y0_matrix_3,ind)$LB

LB_Y0_ori_7 <- medianAbsPercBiasY0_year(Mu_trt_ori_7,Y0_matrix_7,ind)$LB
LB_Y0_space_7 <- medianAbsPercBiasY0_year(Mu_trt_space_7,Y0_matrix_7,ind)$LB
LB_Y0_ICAR_7 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_ICAR_7,Y0_matrix_7,ind)$LB
LB_Y0_AR_7 <- medianAbsPercBiasY0_year(Mu_trt_spacetime_AR_7,Y0_matrix_7,ind)$LB
LB_Y0_lasso_7 <- medianAbsPercBiasY0_year(Mu_trt_lasso_7,Y0_matrix_7,ind)$LB
LB_Y0_shrink_7 <- medianAbsPercBiasY0_year(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7,ind)$LB

# Lower bound
LB_ATT_ori_1 <- medianAbsPercBiasATT_year(Mu_trt_ori_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB
LB_ATT_space_1 <- medianAbsPercBiasATT_year(Mu_trt_space_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB
LB_ATT_ICAR_1 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_ICAR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB
LB_ATT_AR_1 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_AR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB
LB_ATT_lasso_1 <- medianAbsPercBiasATT_year(Mu_trt_lasso_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB

LB_ATT_ori_3 <- medianAbsPercBiasATT_year(Mu_trt_ori_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB
LB_ATT_space_3 <- medianAbsPercBiasATT_year(Mu_trt_space_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB
LB_ATT_ICAR_3 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_ICAR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB
LB_ATT_AR_3 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_AR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB
LB_ATT_lasso_3 <- medianAbsPercBiasATT_year(Mu_trt_lasso_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB

LB_ATT_ori_7 <- medianAbsPercBiasATT_year(Mu_trt_ori_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_space_7 <- medianAbsPercBiasATT_year(Mu_trt_space_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_ICAR_7 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_ICAR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_AR_7 <- medianAbsPercBiasATT_year(Mu_trt_spacetime_AR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_lasso_7 <- medianAbsPercBiasATT_year(Mu_trt_lasso_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_shrink_7 <- medianAbsPercBiasATT_year(Mu_trt_space_comb_time_shrink_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB

# Preparing results for table 
Model <- c(rep(c("Vanilla", "Space", "Space-Time ICAR", "Space-Time AR", "Lasso"),3),"Space-Time Shrinkage")
k <- c(rep(1,5),rep(3,5),rep(7,6))
absPercBias_Y0 <- round(c(Y0_ori_1,Y0_space_1,Y0_ICAR_1,Y0_AR_1,Y0_lasso_1,
                          Y0_ori_3,Y0_space_3,Y0_ICAR_3,Y0_AR_3,Y0_lasso_3,
                          Y0_ori_7,Y0_space_7,Y0_ICAR_7,Y0_AR_7,Y0_lasso_7,Y0_shrink_7),2)
UB_absPercBias_Y0 <- round(c(UB_Y0_ori_1,UB_Y0_space_1,UB_Y0_ICAR_1,UB_Y0_AR_1,UB_Y0_lasso_1,
                             UB_Y0_ori_3,UB_Y0_space_3,UB_Y0_ICAR_3,UB_Y0_AR_3,UB_Y0_lasso_3,
                             UB_Y0_ori_7,UB_Y0_space_7,UB_Y0_ICAR_7,UB_Y0_AR_7,UB_Y0_lasso_7,UB_Y0_shrink_7),2)
LB_absPercBias_Y0 <- round(c(LB_Y0_ori_1,LB_Y0_space_1,LB_Y0_ICAR_1,LB_Y0_AR_1,LB_Y0_lasso_1,
                             LB_Y0_ori_3,LB_Y0_space_3,LB_Y0_ICAR_3,LB_Y0_AR_3,LB_Y0_lasso_3,
                             LB_Y0_ori_7,LB_Y0_space_7,LB_Y0_ICAR_7,LB_Y0_AR_7,LB_Y0_lasso_7,LB_Y0_shrink_7),2)
absPercBias_ATT <- round(c(ATT_ori_1,ATT_space_1,ATT_ICAR_1,ATT_AR_1,ATT_lasso_1,
                           ATT_ori_3,ATT_space_3,ATT_ICAR_3,ATT_AR_3,ATT_lasso_3,
                           ATT_ori_7,ATT_space_7,ATT_ICAR_7,ATT_AR_7,ATT_lasso_7,ATT_shrink_7),2)
UB_absPercBias_ATT <- round(c(UB_ATT_ori_1,UB_ATT_space_1,UB_ATT_ICAR_1,UB_ATT_AR_1,UB_ATT_lasso_1,
                              UB_ATT_ori_3,UB_ATT_space_3,UB_ATT_ICAR_3,UB_ATT_AR_3,UB_ATT_lasso_3,
                              UB_ATT_ori_7,UB_ATT_space_7,UB_ATT_ICAR_7,UB_ATT_AR_7,UB_ATT_lasso_7,UB_ATT_shrink_7),2)
LB_absPercBias_ATT <- round(c(LB_ATT_ori_1,LB_ATT_space_1,LB_ATT_ICAR_1,LB_ATT_AR_1,LB_ATT_lasso_1,
                              LB_ATT_ori_3,LB_ATT_space_3,LB_ATT_ICAR_3,LB_ATT_AR_3,LB_ATT_lasso_3,
                              LB_ATT_ori_7,LB_ATT_space_7,LB_ATT_ICAR_7,LB_ATT_AR_7,LB_ATT_lasso_7,LB_ATT_shrink_7),2)
df <- cbind(Model,k,absPercBias_Y0,LB_absPercBias_Y0,UB_absPercBias_Y0,absPercBias_ATT,LB_absPercBias_ATT,UB_absPercBias_ATT)

# Save results
if(smooth == "smooth"){
  filename_out = paste0(dir_out,"/results_",rare,"_smooth.RData")
}else{
  filename_out = paste0(dir_out,"/results_",rare,".RData")
}

save(df, file = filename_out)

#Y_pred: samples from the posterior predictive distribution (used in paper)
# Absolute Percent Bias of the Counterfactual
Y0_ori_1 <- medianAbsPercBiasY0_year(Y_pred_ori_1,Y0_matrix_1,ind)$median
Y0_space_1 <- medianAbsPercBiasY0_year(Y_pred_space_1,Y0_matrix_1,ind)$median
Y0_ICAR_1 <- medianAbsPercBiasY0_year(Y_pred_spacetime_ICAR_1,Y0_matrix_1,ind)$median
Y0_AR_1 <- medianAbsPercBiasY0_year(Y_pred_spacetime_AR_1,Y0_matrix_1,ind)$median
Y0_lasso_1 <- medianAbsPercBiasY0_year(Y_pred_lasso_1,Y0_matrix_1,ind)$median

Y0_ori_3 <- medianAbsPercBiasY0_year(Y_pred_ori_3,Y0_matrix_3,ind)$median
Y0_space_3 <- medianAbsPercBiasY0_year(Y_pred_space_3,Y0_matrix_3,ind)$median
Y0_ICAR_3 <- medianAbsPercBiasY0_year(Y_pred_spacetime_ICAR_3,Y0_matrix_3,ind)$median
Y0_AR_3 <- medianAbsPercBiasY0_year(Y_pred_spacetime_AR_3,Y0_matrix_3,ind)$median
Y0_lasso_3 <- medianAbsPercBiasY0_year(Y_pred_lasso_3,Y0_matrix_3,ind)$median

Y0_ori_7 <- medianAbsPercBiasY0_year(Y_pred_ori_7,Y0_matrix_7,ind)$median
Y0_space_7 <- medianAbsPercBiasY0_year(Y_pred_space_7,Y0_matrix_7,ind)$median
Y0_ICAR_7 <- medianAbsPercBiasY0_year(Y_pred_spacetime_ICAR_7,Y0_matrix_7,ind)$median
Y0_AR_7 <- medianAbsPercBiasY0_year(Y_pred_spacetime_AR_7,Y0_matrix_7,ind)$median
Y0_lasso_7 <- medianAbsPercBiasY0_year(Y_pred_lasso_7,Y0_matrix_7,ind)$median
Y0_shrink_7 <- medianAbsPercBiasY0_year(Y_pred_space_comb_time_shrink_7,Y0_matrix_7,ind)$median

# Absolute Percent Bias of the ATT
ATT_ori_1 <- medianAbsPercBiasATT_year(Y_pred_ori_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median
ATT_space_1 <- medianAbsPercBiasATT_year(Y_pred_space_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median
ATT_ICAR_1 <- medianAbsPercBiasATT_year(Y_pred_spacetime_ICAR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median
ATT_AR_1 <- medianAbsPercBiasATT_year(Y_pred_spacetime_AR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median
ATT_lasso_1 <- medianAbsPercBiasATT_year(Y_pred_lasso_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$median

ATT_ori_3 <- medianAbsPercBiasATT_year(Y_pred_ori_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median
ATT_space_3 <- medianAbsPercBiasATT_year(Y_pred_space_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median
ATT_ICAR_3 <- medianAbsPercBiasATT_year(Y_pred_spacetime_ICAR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median
ATT_AR_3 <- medianAbsPercBiasATT_year(Y_pred_spacetime_AR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median
ATT_lasso_3 <- medianAbsPercBiasATT_year(Y_pred_lasso_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$median

ATT_ori_7 <- medianAbsPercBiasATT_year(Y_pred_ori_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_space_7 <- medianAbsPercBiasATT_year(Y_pred_space_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_ICAR_7 <- medianAbsPercBiasATT_year(Y_pred_spacetime_ICAR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_AR_7 <- medianAbsPercBiasATT_year(Y_pred_spacetime_AR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_lasso_7 <- medianAbsPercBiasATT_year(Y_pred_lasso_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median
ATT_shrink_7 <- medianAbsPercBiasATT_year(Y_pred_space_comb_time_shrink_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$median

# Upper bound
UB_Y0_ori_1 <- medianAbsPercBiasY0_year(Y_pred_ori_1,Y0_matrix_1,ind)$UB
UB_Y0_space_1 <- medianAbsPercBiasY0_year(Y_pred_space_1,Y0_matrix_1,ind)$UB
UB_Y0_ICAR_1 <- medianAbsPercBiasY0_year(Y_pred_spacetime_ICAR_1,Y0_matrix_1,ind)$UB
UB_Y0_AR_1 <- medianAbsPercBiasY0_year(Y_pred_spacetime_AR_1,Y0_matrix_1,ind)$UB
UB_Y0_lasso_1 <- medianAbsPercBiasY0_year(Y_pred_lasso_1,Y0_matrix_1,ind)$UB

UB_Y0_ori_3 <- medianAbsPercBiasY0_year(Y_pred_ori_3,Y0_matrix_3,ind)$UB
UB_Y0_space_3 <- medianAbsPercBiasY0_year(Y_pred_space_3,Y0_matrix_3,ind)$UB
UB_Y0_ICAR_3 <- medianAbsPercBiasY0_year(Y_pred_spacetime_ICAR_3,Y0_matrix_3,ind)$UB
UB_Y0_AR_3 <- medianAbsPercBiasY0_year(Y_pred_spacetime_AR_3,Y0_matrix_3,ind)$UB
UB_Y0_lasso_3 <- medianAbsPercBiasY0_year(Y_pred_lasso_3,Y0_matrix_3,ind)$UB

UB_Y0_ori_7 <- medianAbsPercBiasY0_year(Y_pred_ori_7,Y0_matrix_7,ind)$UB
UB_Y0_space_7 <- medianAbsPercBiasY0_year(Y_pred_space_7,Y0_matrix_7,ind)$UB
UB_Y0_ICAR_7 <- medianAbsPercBiasY0_year(Y_pred_spacetime_ICAR_7,Y0_matrix_7,ind)$UB
UB_Y0_AR_7 <- medianAbsPercBiasY0_year(Y_pred_spacetime_AR_7,Y0_matrix_7,ind)$UB
UB_Y0_lasso_7 <- medianAbsPercBiasY0_year(Y_pred_lasso_7,Y0_matrix_7,ind)$UB
UB_Y0_shrink_7 <- medianAbsPercBiasY0_year(Y_pred_space_comb_time_shrink_7,Y0_matrix_7,ind)$UB

# Upper bound
UB_ATT_ori_1 <- medianAbsPercBiasATT_year(Y_pred_ori_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB
UB_ATT_space_1 <- medianAbsPercBiasATT_year(Y_pred_space_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB
UB_ATT_ICAR_1 <- medianAbsPercBiasATT_year(Y_pred_spacetime_ICAR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB
UB_ATT_AR_1 <- medianAbsPercBiasATT_year(Y_pred_spacetime_AR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB
UB_ATT_lasso_1 <- medianAbsPercBiasATT_year(Y_pred_lasso_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$UB

UB_ATT_ori_3 <- medianAbsPercBiasATT_year(Y_pred_ori_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB
UB_ATT_space_3 <- medianAbsPercBiasATT_year(Y_pred_space_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB
UB_ATT_ICAR_3 <- medianAbsPercBiasATT_year(Y_pred_spacetime_ICAR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB
UB_ATT_AR_3 <- medianAbsPercBiasATT_year(Y_pred_spacetime_AR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB
UB_ATT_lasso_3 <- medianAbsPercBiasATT_year(Y_pred_lasso_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$UB

UB_ATT_ori_7 <- medianAbsPercBiasATT_year(Y_pred_ori_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_space_7 <- medianAbsPercBiasATT_year(Y_pred_space_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_ICAR_7 <- medianAbsPercBiasATT_year(Y_pred_spacetime_ICAR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_AR_7 <- medianAbsPercBiasATT_year(Y_pred_spacetime_AR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_lasso_7 <- medianAbsPercBiasATT_year(Y_pred_lasso_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB
UB_ATT_shrink_7 <- medianAbsPercBiasATT_year(Y_pred_space_comb_time_shrink_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$UB

# Lower bound
LB_Y0_ori_1 <- medianAbsPercBiasY0_year(Y_pred_ori_1,Y0_matrix_1,ind)$LB
LB_Y0_space_1 <- medianAbsPercBiasY0_year(Y_pred_space_1,Y0_matrix_1,ind)$LB
LB_Y0_ICAR_1 <- medianAbsPercBiasY0_year(Y_pred_spacetime_ICAR_1,Y0_matrix_1,ind)$LB
LB_Y0_AR_1 <- medianAbsPercBiasY0_year(Y_pred_spacetime_AR_1,Y0_matrix_1,ind)$LB
LB_Y0_lasso_1 <- medianAbsPercBiasY0_year(Y_pred_lasso_1,Y0_matrix_1,ind)$LB

LB_Y0_ori_3 <- medianAbsPercBiasY0_year(Y_pred_ori_3,Y0_matrix_3,ind)$LB
LB_Y0_space_3 <- medianAbsPercBiasY0_year(Y_pred_space_3,Y0_matrix_3,ind)$LB
LB_Y0_ICAR_3 <- medianAbsPercBiasY0_year(Y_pred_spacetime_ICAR_3,Y0_matrix_3,ind)$LB
LB_Y0_AR_3 <- medianAbsPercBiasY0_year(Y_pred_spacetime_AR_3,Y0_matrix_3,ind)$LB
LB_Y0_lasso_3 <- medianAbsPercBiasY0_year(Y_pred_lasso_3,Y0_matrix_3,ind)$LB

LB_Y0_ori_7 <- medianAbsPercBiasY0_year(Y_pred_ori_7,Y0_matrix_7,ind)$LB
LB_Y0_space_7 <- medianAbsPercBiasY0_year(Y_pred_space_7,Y0_matrix_7,ind)$LB
LB_Y0_ICAR_7 <- medianAbsPercBiasY0_year(Y_pred_spacetime_ICAR_7,Y0_matrix_7,ind)$LB
LB_Y0_AR_7 <- medianAbsPercBiasY0_year(Y_pred_spacetime_AR_7,Y0_matrix_7,ind)$LB
LB_Y0_lasso_7 <- medianAbsPercBiasY0_year(Y_pred_lasso_7,Y0_matrix_7,ind)$LB
LB_Y0_shrink_7 <- medianAbsPercBiasY0_year(Y_pred_space_comb_time_shrink_7,Y0_matrix_7,ind)$LB

# Lower bound
LB_ATT_ori_1 <- medianAbsPercBiasATT_year(Y_pred_ori_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB
LB_ATT_space_1 <- medianAbsPercBiasATT_year(Y_pred_space_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB
LB_ATT_ICAR_1 <- medianAbsPercBiasATT_year(Y_pred_spacetime_ICAR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB
LB_ATT_AR_1 <- medianAbsPercBiasATT_year(Y_pred_spacetime_AR_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB
LB_ATT_lasso_1 <- medianAbsPercBiasATT_year(Y_pred_lasso_1,Y0_matrix_1,Y1_matrix_1,true_trt_effect_1,ind)$LB

LB_ATT_ori_3 <- medianAbsPercBiasATT_year(Y_pred_ori_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB
LB_ATT_space_3 <- medianAbsPercBiasATT_year(Y_pred_space_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB
LB_ATT_ICAR_3 <- medianAbsPercBiasATT_year(Y_pred_spacetime_ICAR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB
LB_ATT_AR_3 <- medianAbsPercBiasATT_year(Y_pred_spacetime_AR_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB
LB_ATT_lasso_3 <- medianAbsPercBiasATT_year(Y_pred_lasso_3,Y0_matrix_3,Y1_matrix_3,true_trt_effect_3,ind)$LB

LB_ATT_ori_7 <- medianAbsPercBiasATT_year(Y_pred_ori_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_space_7 <- medianAbsPercBiasATT_year(Y_pred_space_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_ICAR_7 <- medianAbsPercBiasATT_year(Y_pred_spacetime_ICAR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_AR_7 <- medianAbsPercBiasATT_year(Y_pred_spacetime_AR_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_lasso_7 <- medianAbsPercBiasATT_year(Y_pred_lasso_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB
LB_ATT_shrink_7 <- medianAbsPercBiasATT_year(Y_pred_space_comb_time_shrink_7,Y0_matrix_7,Y1_matrix_7,true_trt_effect_7,ind)$LB

# Preparing results for table 
Model <- c(rep(c("Vanilla", "Space", "Space-Time ICAR", "Space-Time AR", "Lasso"),3),"Space-Time Shrinkage")
k <- c(rep(1,5),rep(3,5),rep(7,6))
absPercBias_Y0 <- round(c(Y0_ori_1,Y0_space_1,Y0_ICAR_1,Y0_AR_1,Y0_lasso_1,
                          Y0_ori_3,Y0_space_3,Y0_ICAR_3,Y0_AR_3,Y0_lasso_3,
                          Y0_ori_7,Y0_space_7,Y0_ICAR_7,Y0_AR_7,Y0_lasso_7,Y0_shrink_7),2)
UB_absPercBias_Y0 <- round(c(UB_Y0_ori_1,UB_Y0_space_1,UB_Y0_ICAR_1,UB_Y0_AR_1,UB_Y0_lasso_1,
                             UB_Y0_ori_3,UB_Y0_space_3,UB_Y0_ICAR_3,UB_Y0_AR_3,UB_Y0_lasso_3,
                             UB_Y0_ori_7,UB_Y0_space_7,UB_Y0_ICAR_7,UB_Y0_AR_7,UB_Y0_lasso_7,UB_Y0_shrink_7),2)
LB_absPercBias_Y0 <- round(c(LB_Y0_ori_1,LB_Y0_space_1,LB_Y0_ICAR_1,LB_Y0_AR_1,LB_Y0_lasso_1,
                             LB_Y0_ori_3,LB_Y0_space_3,LB_Y0_ICAR_3,LB_Y0_AR_3,LB_Y0_lasso_3,
                             LB_Y0_ori_7,LB_Y0_space_7,LB_Y0_ICAR_7,LB_Y0_AR_7,LB_Y0_lasso_7,LB_Y0_shrink_7),2)
absPercBias_ATT <- round(c(ATT_ori_1,ATT_space_1,ATT_ICAR_1,ATT_AR_1,ATT_lasso_1,
                           ATT_ori_3,ATT_space_3,ATT_ICAR_3,ATT_AR_3,ATT_lasso_3,
                           ATT_ori_7,ATT_space_7,ATT_ICAR_7,ATT_AR_7,ATT_lasso_7,ATT_shrink_7),2)
UB_absPercBias_ATT <- round(c(UB_ATT_ori_1,UB_ATT_space_1,UB_ATT_ICAR_1,UB_ATT_AR_1,UB_ATT_lasso_1,
                              UB_ATT_ori_3,UB_ATT_space_3,UB_ATT_ICAR_3,UB_ATT_AR_3,UB_ATT_lasso_3,
                              UB_ATT_ori_7,UB_ATT_space_7,UB_ATT_ICAR_7,UB_ATT_AR_7,UB_ATT_lasso_7,UB_ATT_shrink_7),2)
LB_absPercBias_ATT <- round(c(LB_ATT_ori_1,LB_ATT_space_1,LB_ATT_ICAR_1,LB_ATT_AR_1,LB_ATT_lasso_1,
                              LB_ATT_ori_3,LB_ATT_space_3,LB_ATT_ICAR_3,LB_ATT_AR_3,LB_ATT_lasso_3,
                              LB_ATT_ori_7,LB_ATT_space_7,LB_ATT_ICAR_7,LB_ATT_AR_7,LB_ATT_lasso_7,LB_ATT_shrink_7),2)
df <- cbind(Model,k,absPercBias_Y0,LB_absPercBias_Y0,UB_absPercBias_Y0,absPercBias_ATT,LB_absPercBias_ATT,UB_absPercBias_ATT)


# Save results
if(smooth == "smooth"){
  filename_out = paste0(dir_out,"/results_",rare,"_smooth_Ypred.RData")
}else{
  filename_out = paste0(dir_out,"/results_",rare,"_Ypred.RData")
}


save(df, file = filename_out)