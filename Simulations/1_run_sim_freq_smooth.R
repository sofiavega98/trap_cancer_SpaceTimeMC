# Testing pseudo-spatial simulations using New Mexico 
# Generates data with "true" census pop and NM spatial matrix  
# Run on 3 models: MC in gsynth, GSC in gsynth, ASC in augsynth 

true_k <- 3
est_k <- c(1,3,7)

#############################################################################################################
## 0. Read command line arguments and load packages and functions                                          ##
##    Arguments should include:                                                                            ##
##    - alpha: parameter for rarity                                                                        ##
##    - rho: degree of spatial corr                                                                        ##
##    - true_k: true k                                                                                     ##
##    - dir_out: directory output path for results                                                         ##
#############################################################################################################

library(dplyr)
library(tidyr)
library(netdiffuseR)
library(gsynth)
library(augsynth)

# Read command line arguments 
args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }



#####################################################################
## 1. Load prespecified NM spacial and temporal adjacency matrices ##
#####################################################################

# Load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Functions.R')

# Load prespecified NM spacial and temporal adjacency matrices
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/0_set_up.R') # We are basing the simulated data on the values from this script


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
Y0_matrix <- matrix(NA,nrow = 90,ncol = 100) 
Y1_matrix <- matrix(NA, nrow = 42, ncol = 100)

fit_gsc_mc_1 <- list()
fit_gsc_mc_3 <- list()
fit_gsc_mc_7 <- list()
fit_gsc <- list()
fit_asc <- list()


########################
## 3. Run Simulations ##
########################

for(i in 1:100){
  simnum <- i
  
  # Get sim data
  sim_data <- get_sim_spatial_rate_data_smooth_freq(partial_W_nm,temp_adj_matrix, trt_fips,
                                             pop_mat, k=true_k, alpha = alpha , gamma_sd = gamma_sd,
                                             psi_sd = psi_sd, n_trt = 6, m_trt = 7, k_stan_in = k_stan_vec, 
                                             rho = rho, rate = rate, tau2 = tau2, simnum = simnum)
  Y0_obs <- sim_data[[1]]
  Y0_full <- sim_data[[2]]
  Y1_obs <- sim_data[[3]]
  
  # Store Y0_full
  Y0_matrix[,i] <- as.vector(t(Y0_full[1:length(trt_fips),]))
  
  # Store Y1_obs
  Y1_matrix[,i] <- as.vector(t(Y1_obs))
  
  
  # Get data for gsynth
  m_trt = 7
  n_trt = 6
  trt <- matrix(0,nrow=27,ncol=15)
  trt[1:(n_trt),9:15] <- 1
  treated = as.vector(t(trt))
  
  Y0_full_rate <- (Y0_obs/pop_mat) * 100000
  
  df_gsynth <- data.frame(rate = as.vector(t(Y0_full_rate)), treated = treated, 
                          FIPS = NM_new$FIPS,Year = NM_new$YEAR_DX)
  
  
  # Fit MC in gsynth
  set.seed(501)
  fit_gsc_mc_1[[i]] <- gsynth(rate ~ treated, data = df_gsynth, estimator = "mc", k=est_k[1], 
                            index = c("FIPS","Year"), force = "two-way", CV = TRUE, 
                            r = c(0, 5), se = FALSE,  nboots = 1000, parallel = FALSE)
  set.seed(501)
  fit_gsc_mc_3[[i]] <- gsynth(rate ~ treated, data = df_gsynth, estimator = "mc", k=est_k[2], 
                              index = c("FIPS","Year"), force = "two-way", CV = TRUE, 
                              r = c(0, 5), se = FALSE,  nboots = 1000, parallel = FALSE)
  set.seed(501)
  fit_gsc_mc_7[[i]] <- gsynth(rate ~ treated, data = df_gsynth, estimator = "mc", k=est_k[3], 
                              index = c("FIPS","Year"), force = "two-way", CV = TRUE, 
                              r = c(0, 5), se = FALSE,  nboots = 1000, parallel = FALSE)
  
  # Fit GSC in gsynth
  set.seed(501)
  fit_gsc[[i]] <- gsynth(rate ~ treated, data = df_gsynth, estimator = "ife", 
                         index = c("FIPS","Year"), force = "two-way", CV = TRUE, 
                         r = c(0, 5), se = FALSE,  nboots = 1000, parallel = FALSE)
  
  # Fit ASC in augsynth
  set.seed(501)
  fit_asc[[i]] <- augsynth(rate ~ treated, unit = FIPS, time = Year, data = df_gsynth,
                           progfunc = "Ridge", scm = T)
}

########################
## 4. Save Results    ##
########################

results <- list(Y0_matrix=Y0_matrix,Y1_matrix = Y1_matrix,fit_gsc_mc_1 = fit_gsc_mc_1,fit_gsc_mc_3 = fit_gsc_mc_3,
                fit_gsc_mc_7 = fit_gsc_mc_7,fit_gsc = fit_gsc, fit_asc = fit_asc)

if(alpha==-5){
  rare="nonRare"
}else{
  rare="Rare"
}
save(results, file=paste0(dir_out,"/sim_",rare,"_freq_smooth.RData"))
