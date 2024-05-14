# This script processes results from the application
# It also runs the application the MC through gsynth package

trt_year=1995 # treated year
data_path="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/Lymphoma_1988_2003_0-29.RData" # path to application data
adj_mat_path="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/lymphoma_adj_mat_1988_2003_0-29.RData" # path to adjacency matrix
dir_out="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Lymphoma/Results/" # directory output path for results


# Load packages
library(tidyverse)
library(gsynth)

# Load functions
source('~/cl_quasi_exp/Functions.R')
source('~/cl_quasi_exp/application_functions.R')

# Load data
load(data_path)
data_full_hisp <- data_full_hisp_noZero

######################################
## Load STAN Results for Models 1-6 ##
######################################

# Vanilla model
load(paste0(dir_out,"fit_vanilla_trt",trt_year,"_k3_smooth.RData"))
fit_ori <- fit
post_samp_ori <- rstan::extract(fit) # posterior samples 
Mu_trt_ori <- post_samp_ori$Y_pred # Extract Mu_trt and store results

# space model
load(paste0(dir_out,"fit_space_trt",trt_year,"_k3_smooth.RData"))
fit_space <- fit
post_samp_space <- rstan::extract(fit)  # posterior samples 
Mu_trt_space <- post_samp_space$Y_pred # Extract Mu_trt and store results

# ICAR model
load(paste0(dir_out,"fit_ICAR_trt",trt_year,"_k3_smooth.RData"))
fit_ICAR <- fit
post_samp_ICAR <- rstan::extract(fit) # posterior samples 
Mu_trt_ICAR <- post_samp_ICAR$Y_pred # Extract Mu_trt and store results

# AR model
load(paste0(dir_out,"fit_AR_trt",trt_year,"_k3_smooth.RData"))
fit_AR <- fit
post_samp_AR <- rstan::extract(fit) # posterior samples 
Mu_trt_AR <- post_samp_AR$Y_pred # Extract Mu_trt and store results

# lasso model 
load(paste0(dir_out,"fit_lasso_trt",trt_year,"_k3_smooth.RData"))
fit_lasso <- fit
post_samp_spacetime_lasso<- rstan::extract(fit)  # posterior samples 
Mu_trt_lasso <- post_samp_spacetime_lasso$Y_pred # Extract Mu_trt and store results

# shrinkage model
load(paste0(dir_out,"fit_shrink_trt",trt_year,"_k4_smooth.RData"))
fit_shrink<- fit
post_samp_shrink <- rstan::extract(fit)  # posterior samples 
Mu_trt_shrink <- post_samp_shrink$Y_pred # Extract Mu_trt and store results

############
## Set Up ##
############
# get treated counties
trt_fips <- data_full_hisp %>% filter(C == 1) %>%
  distinct(FIPS) %>% pull(FIPS)

# get Y0_obs
Y0_obs <- get_Y0_obs(data_full_hisp,trt_fips,trt_year = trt_year)[[1]]

# get Y(1) i.e. observed treated values
Y1 <- get_Y1_obs1(data_full_hisp) # 13*7=91 treated values; vector of length 91 (county1 time1, county1 time2,...)

# get population matrix
pop_mat <- get_pop_mat(data_full_hisp)

# get population matrix for treated counties
n_trt <- length(trt_fips)

#get populations for treated counties
pop_trt <- as.vector(t(pop_mat[1:n_trt,]))

# get number of treated time points
m_trt <- sum(is.na(Y0_obs)*1)/n_trt

# get number of time points
m <- length(unique(data_full_hisp$YEAR_DX))

# Get indices for treated times
ind <- c()
for(i in 0:(n_trt-1)){
  ind <- append(ind,seq((m_trt-1),m) + (i*m))
}

###################################
## Filter out -1 MCMC iterations ##
###################################
# Dimensions for model fits are number of MCMC iterations by number of treated units

# Vanilla model
if(length(which(Mu_trt_ori == -1)) > 0){
  
  # Identify the rows -1s (identify which MCMC iteration has the overflow problem)
  rows <- which(Mu_trt_ori == -1, arr.ind = TRUE)[,1] #[,1] since first column of the which statement identify rows
  
  # Filter out
  print("filtering")
  Mu_trt_ori <- Mu_trt_ori[-rows,]
  
}

# Space model
if(length(which(Mu_trt_space == -1)) > 0){
  
  # Identify the rows -1s (identify which MCMC iteration has the overflow problem)
  rows <- which(Mu_trt_space == -1, arr.ind = TRUE)[,1] #[,1] since first column of the which statement identify rows
  
  # Filter out
  print("filtering")
  Mu_trt_space <- Mu_trt_space[-rows,]
  
}

# ICAR model
if(length(which(Mu_trt_ICAR == -1)) > 0){
  
  # Identify the rows -1s (identify which MCMC iteration has the overflow problem)
  rows <- which(Mu_trt_ICAR == -1, arr.ind = TRUE)[,1] #[,1] since first column of the which statement identify rows
  
  # Filter out
  print("filtering")
  Mu_trt_ICAR <- Mu_trt_ICAR[-rows,]
  
}

# AR model
if(length(which(Mu_trt_AR == -1)) > 0){
  
  # Identify the rows -1s (identify which MCMC iteration has the overflow problem)
  rows <- which(Mu_trt_AR == -1, arr.ind = TRUE)[,1] #[,1] since first column of the which statement identify rows
  
  # Filter out
  print("filtering")
  Mu_trt_AR <- Mu_trt_AR[-rows,]
  
}

# Lasso model
if(length(which(Mu_trt_lasso == -1)) > 0){
  
  # Identify the rows -1s (identify which MCMC iteration has the overflow problem)
  rows <- which(Mu_trt_lasso == -1, arr.ind = TRUE)[,1] #[,1] since first column of the which statement identify rows
  
  # Filter out
  print("filtering")
  Mu_trt_lasso <- Mu_trt_lasso[-rows,]
  
}

# Shrink model
if(length(which(Mu_trt_shrink == -1)) > 0){
  
  # Identify the rows -1s (identify which MCMC iteration has the overflow problem)
  rows <- which(Mu_trt_shrink == -1, arr.ind = TRUE)[,1] #[,1] since first column of the which statement identify rows
  
  # Filter out
  print("filtering")
  Mu_trt_shrink <- Mu_trt_shrink[-rows,]
  
}


#################
## Diagnostics ##
#################
rhat_vec <- rhat_table_YPred_v2(fit_ori, fit_space, fit_ICAR, fit_AR, 
                                fit_shrink, fit_lasso, ind)
save(rhat_vec, file=paste0(dir_out,"rhat_vec_trt",trt_year,"_smooth_SensAnalysis.RData"))

##########################################
## Run Application on MC through gsynth ##
##########################################

# Get smoothed data using a GLM model
m = length(unique(data_full_hisp$YEAR_DX))
for(i in 1:length(unique(data_full_hisp$FIPS))){
  temp_county <- data_full_hisp %>%
    filter(FIPS == unique(data_full_hisp$FIPS)[i])
  fit <- glm(CL_CASES ~ ns(YEAR_DX,5) + offset(log(POP)), data = temp_county, family = "poisson")
  # Save smooth data as CL_CASES. This will keep in sync with processing functions
  data_full_hisp[seq((i-1)*m+1,m*i),]$CL_CASES <- round(predict.glm(fit, type = "response"))
}

# Convert cases into rates for gsynth
data_full_hisp$rate <- 100000 * (data_full_hisp$CL_CASES/data_full_hisp$POP)

# Run Gysnth
fit_gsc <- gsynth(rate ~ C + pct_hispanic_pop, data = data_full_hisp, estimator = "mc", 
                  index = c("FIPS","YEAR_DX"), force = "two-way", CV = TRUE, r = c(0, 5), 
                  se = TRUE, k=3, nboots = 1000, parallel = FALSE)


############################
## Save processed results ##
############################
res <- list(Mu_trt_ori_smooth = Mu_trt_ori, Mu_trt_space_smooth = Mu_trt_space, Mu_trt_ICAR_smooth = Mu_trt_ICAR, 
            Mu_trt_AR_smooth = Mu_trt_AR, Mu_trt_lasso_smooth = Mu_trt_lasso, Mu_trt_shrink_smooth = Mu_trt_shrink, 
            fit_gsc_smooth = fit_gsc)

save(res, file=paste0(dir_out,"processed_res_trt",trt_year,"_smooth_SensAnalysis.RData"))
