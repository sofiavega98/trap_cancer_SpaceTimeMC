# This script runs the models on the application data using the batch script

#############################################################################################################
## 0. Read command line arguments and load packages and functions                                          ##
##    Arguments should include:                                                                            ##
##    - k: number of latent factors                                                                        ##
##    - p: number of covariates                                                                            ##
##    - trt_year: treated year                                                                             ##
##    - model: which model needs to be run                                                                 ##
##     - In batch script, specify model=vanilla for Vanilla (Original) model, model=space for Space model  ##
##       model=ICAR for Space-Time ICAR model, model=AR for Space-Time AR model,                           ##
##       model=lasso for Space-Time Lasso model, model=shrink for Space-Time Shrinkage model               ##
##    - data_path: path to application data                                                                ##
##    - adj_mat_path: path to adjacency matrix for associated application data                             ##
##    - dir_out: directory output path for results                                                         ##
#############################################################################################################

args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }

#model = model
#trt_year = trt_year
#data_path = data_path
#adj_mat_path = adj_mat_path
#dir_out = dir_out

# For testing
#years <- c(1983:2003)

# Load packages
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)
library(tidyverse)
library("netdiffuseR")

# Load functions
source('~/cl_quasi_exp/Functions.R')

#####################################################################
## 1. Load prespecified spacial and temporal adjacency matrices    ##
#####################################################################

# Load data
#load(paste0('~/cl_quasi_exp/data/Lymphoma_',years[1],'_',years[length(years)],'_0-29.RData'))
load(data_path)
data_full_hisp <- data_full_hisp_noZero

# Load spacial adjacency matrix
load(adj_mat_path)
#load(paste0('~/cl_quasi_exp/data/lymphoma_adj_mat_',years[1],'_',years[length(years)],'_0-29.RData'))

########################
## 2. Load the models ##
########################

# Vanilla Model
if(model == 'vanilla'){
  stan_mod <- stan_model("~/cl_quasi_exp/STAN Models/SV_stan_cov_v2.stan",auto_write=F)
}

# Space Model
if(model == 'space'){
  stan_mod <- stan_model("~/cl_quasi_exp/STAN Models/SV_spatial_stan_cov_v2.stan",auto_write=F)
}

# Space-Time ICAR Model
if(model == 'ICAR'){
  stan_mod <- stan_model("~/cl_quasi_exp/STAN Models/SV_spatio-temporal_STAN_cov_v2.stan",auto_write=F)
}

# Space-Time AR Model
if(model == 'AR'){
  stan_mod <- stan_model("~/cl_quasi_exp/STAN Models/SV_spatio-temporal_AR_STAN_cov_v2.stan",auto_write=F)
}

# Space-Time Lasso Model
if(model == 'lasso'){
  stan_mod <- stan_model("~/cl_quasi_exp/STAN Models/SV_spatio-temporal_lasso_STAN_cov_v2.stan",auto_write=F)
}

# Space-Time Shrinkage Model
if(model == 'shrink'){
  stan_mod <- stan_model("~/cl_quasi_exp/STAN Models/SV_spatio-comb-temp-shrinkage_STAN_cov_v2.stan",auto_write=F)
}

#######################
## 3. Define stan_in ##
#######################

# define k
k = k

# Get treated counties
trt_fips <- data_full_hisp %>% filter(C == 1) %>%
  distinct(FIPS) %>% pull(FIPS)

# Get temporal adjacency matrix
m = length(unique(data_full_hisp$YEAR_DX))
n = length(unique(data_full_hisp$YEAR_DX))
temp_adj_matrix <- get_temp_adj_mat(m,n)

# Set Up covariates
p = p #number of covariates
x = data.frame(data_full_hisp$pct_hispanic_pop)
col_cov = rep(1:length(unique(data_full_hisp$YEAR_DX)), each = length(unique(data_full_hisp$FIPS)))
row_cov = rep(1:length(unique(data_full_hisp$FIPS)),m)

# Get smoothed data using a GLM model
for(i in 1:length(unique(data_full_hisp$FIPS))){
  temp_county <- data_full_hisp %>%
    filter(FIPS == unique(data_full_hisp$FIPS)[i])
  fit <- glm(CL_CASES ~ ns(YEAR_DX,5) + offset(log(POP)), data = temp_county, family = "poisson")
  # Save smooth data as CL_CASES. This will keep in sync with processing functions
  data_full_hisp[seq((i-1)*m+1,m*i),]$CL_CASES <- round(predict.glm(fit, type = "response"))
}


stan_in_1 <- get_space_time_cov_stan_in(k,data_full_hisp,trt_fips,partial_us_adj_mat,
                                        temp_adj_matrix,trt_year = trt_year,p=p,x=x, 
                                        row_cov=row_cov, col_cov=col_cov)
########################
## 4. Run Model       ##
########################

# Set seed
set.seed(1)

# Fit model
fit <- sampling(object = stan_mod, data = stan_in_1, chains = 1, iter = 2000, 
                warmup=1000, pars = c('Y_pred','Mu_trt'), init_r = .1)

# Save results
save(fit, file=paste0(dir_out,"/fit_",model,"_trt",trt_year,"_k",k,"_smooth.RData"))


