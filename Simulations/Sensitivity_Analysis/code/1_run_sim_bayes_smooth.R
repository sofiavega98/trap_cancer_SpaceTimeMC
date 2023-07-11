# This script generates simulated data and runs the data on models 1-6 from the paper
# Generates data with "true" census pop and NM spatial matrix 
# For this sensitivity analysis, we change the priors on tau^2. (i.e.tau^2~Ga(1,.01))
# See Supplementary Materials.

#############################################################################################################
## 0. Read command line arguments and load packages and functions                                          ##
##    Arguments should include:                                                                            ##
##    - est_k: estimated number of latent factors                                                          ##
##    - alpha: parameter for rarity                                                                        ##
##    - rho: degree of spatial corr                                                                        ##
##    - true_k: true k                                                                                     ##
##    - iter: number of HMC iterations                                                                     ##
##    - warmup: number of HMC iterations for warm-up.                                                      ##
##    - dir_out: directory output path for results                                                         ##
#############################################################################################################

# Read command line arguments 
args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }

# Load packages
library(tidyverse)
library(netdiffuseR)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)
library(mvtnorm)

# Load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Functions.R')

# Load prespecified NM spacial and temporal adjacency matrices
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/0_set_up.R') # We are basing the simulated data on the values from this script


########################
## 1. Load Models 1-6 ##
########################

original_mod <- stan_model("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/STAN Models/Vanilla_altPrior.stan",auto_write=F)
space_mod <- stan_model("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/STAN Models/Space_altPrior.stan",auto_write=F)
space_time_ICAR_mod <- stan_model("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/STAN Models/SpaceTimeICAR_altPrior.stan",auto_write=F)
space_time_AR_mod <- stan_model("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/STAN Models/SpaceTimeAR_altPrior.stan",auto_write=F)
space_comb_time_shrink_mod <- stan_model("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/STAN Models/SpaceTimeShrinkage_altPrior.stan",auto_write=F)
lasso_mod <- stan_model("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/STAN Models/SpaceTimeLasso_altPrior.stan",auto_write=F)

#######################
## 2. Set parameters ##
#######################
alpha = alpha
rho = rho
gamma_sd = 0.000015
psi_sd = 0.000015
rate = 0.5
tau2 = 0.01

true_k <- true_k
est_k <- est_k
k_max <- 7 # max k for shrinkage model
k_stan_vec <- c(est_k,k_max) #only true k for now


###################################
## 3. Set Up Storage for Results ##
###################################
#this will be a matrix with each col being a different simulated data set
#and each row being cl incidence for each of the treated counties at all time points
Y0_matrix <- matrix(NA,nrow = 90,ncol = 1) 
Y1_matrix <- matrix(NA, nrow = 42, ncol = 1)

########################
## 4. Run Simulations ##
########################

# Get simulated data
sim_data <- get_sim_spatial_rate_data_smooth(partial_W_nm,temp_adj_matrix, trt_fips,
                                      pop_mat, k=true_k, alpha = alpha , gamma_sd = gamma_sd,
                                      psi_sd = psi_sd, n_trt = 6, m_trt = 7, k_stan_in = k_stan_vec, 
                                      rho = rho, rate = rate, tau2 = tau2, simnum = simnum)

stan_in <- sim_data[[1]] #list of different stan_in with different k_max according to k_max_vec
Y0_full <- sim_data[[2]] # full counterfactual matrix
Y1_obs <- sim_data[[3]] # observed treated matrix
Y1_obs_rate <- sim_data[[4]] # full counterfactual matrix in rate
Y0_full_rate <- sim_data[[5]] # observed treated matrix in rate

# Store Y0_full
Y0_matrix[,1] <- as.vector(t(Y0_full[1:length(trt_fips),]))

# Store Y1_obs
Y1_matrix[,1] <- as.vector(t(Y1_obs))


#Fit Space-Time Shrinkage model
stan_in_1 <- stan_in[[2]] #fit is k = 7
fit_sim_space_comb_time_shrink <- sampling(object = space_comb_time_shrink_mod, 
                                           data = stan_in_1, chains = 1, iter = iter, warmup=warmup, 
                                           pars = c('Y_pred','Mu_trt'), init_r = .1, seed = 1)

# Fit on models 1-5
stan_in_1 <- stan_in[[1]] #k=est_k


# Fit on original model
fit_sim_ori <- sampling(object = original_mod, 
                        data = stan_in_1, chains = 1, iter = iter, warmup=warmup, 
                        pars = c('Y_pred','Mu_trt'), init_r = .1, seed = 1)


# Fit on space model
fit_sim_space <- sampling(object = space_mod, 
                          data = stan_in_1, chains = 1, iter = iter, warmup=warmup, 
                          pars = c('Y_pred','Mu_trt'), init_r = .1, seed = 1)


# Fit on space time ICAR model
fit_sim_spacetime_ICAR <- sampling(object = space_time_ICAR_mod, 
                                   data = stan_in_1, chains = 1, iter = iter, warmup=warmup, 
                                   pars = c('Y_pred','Mu_trt'), init_r = .1, seed = 1)


# Fit on space time AR model
fit_sim_spacetime_AR <- sampling(object = space_time_AR_mod, 
                                 data = stan_in_1, chains = 1, iter = iter, warmup=warmup, 
                                 pars = c('Y_pred','Mu_trt'), init_r = .1, seed = 1)

# Fit on lasso model
fit_sim_lasso <- sampling(object = lasso_mod, 
                          data = stan_in_1, chains = 1, iter = iter, warmup=warmup, 
                          pars = c('Y_pred','Mu_trt'), init_r = .1, seed = 1)

# Extract data
# original model
post_samp_ori <-rstan::extract(fit_sim_ori) # posterior samples 
Mu_trt_ori <- post_samp_ori$Mu_trt # Extract Mu_trt and store results
Y_pred_ori <- post_samp_ori$Y_pred # Extract Y_pred and store results

# space model
post_samp_space <-rstan::extract(fit_sim_space) # posterior samples 
Mu_trt_space <- post_samp_space$Mu_trt # Extract Mu_trt and store results
Y_pred_space <- post_samp_space$Y_pred # Extract Y_pred and store results

# ICAR model
post_samp_ICAR <-rstan::extract(fit_sim_spacetime_ICAR) # posterior samples 
Mu_trt_spacetime_ICAR <- post_samp_ICAR$Mu_trt # Extract Mu_trt and store results
Y_pred_spacetime_ICAR <- post_samp_ICAR$Y_pred # Extract Y_pred and store results

# AR model
post_samp_AR <-rstan::extract(fit_sim_spacetime_AR) # posterior samples 
Mu_trt_spacetime_AR <- post_samp_AR$Mu_trt # Extract Mu_trt and store results
Y_pred_spacetime_AR <- post_samp_AR$Y_pred # Extract Y_pred and store results

# shrinkage model
post_samp_space_comb_time_shrink<-rstan::extract(fit_sim_space_comb_time_shrink) # posterior samples 
Mu_trt_space_comb_time_shrink <- post_samp_space_comb_time_shrink$Mu_trt # Extract Mu_trt and store results
Y_pred_space_comb_time_shrink <- post_samp_space_comb_time_shrink$Y_pred # Extract Y_pred and store results

# lasso model 
post_samp_spacetime_lasso<-rstan::extract(fit_sim_lasso) # posterior samples 
Mu_trt_lasso <- post_samp_spacetime_lasso$Mu_trt # Extract Mu_trt and store results
Y_pred_lasso <- post_samp_spacetime_lasso$Y_pred # Extract Y_pred and store results

######################
## 5. Store Results ##
######################
results <- list(Y0_matrix=Y0_matrix,Y1_matrix = Y1_matrix, Mu_trt_ori = Mu_trt_ori, 
                Mu_trt_space=Mu_trt_space, Mu_trt_spacetime_ICAR= Mu_trt_spacetime_ICAR,
                Mu_trt_spacetime_AR= Mu_trt_spacetime_AR, Mu_trt_space_comb_time_shrink=Mu_trt_space_comb_time_shrink, 
                Mu_trt_lasso = Mu_trt_lasso,  
                Y_pred_ori = Y_pred_ori, 
                Y_pred_space=Y_pred_space, Y_pred_spacetime_ICAR= Y_pred_spacetime_ICAR,
                Y_pred_spacetime_AR= Y_pred_spacetime_AR, Y_pred_space_comb_time_shrink=Y_pred_space_comb_time_shrink, 
                Y_pred_lasso = Y_pred_lasso,
                Y1_obs_rate = Y1_obs_rate, Y0_full_rate = Y0_full_rate )
if(alpha==-5){
  rare="nonRare"
}else{
  rare="Rare"
}
save(results, file=paste0(dir_out,"/sim_k",est_k,"_",simnum,"_",rare,"_smooth.RData"))
