# This script takes the observed Y1 values and models the Y1s using CARBayes
# prior to assessing ATT. Then,
# This script takes the results from 2_process_results.R and creates
# tables and figures for the paper
# These results are displayed in the supplemental section of the paper

years=1988:2003
trt_year=1995 # treated year
data_path="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/Lymphoma_1988_2003_0-29.RData" # path to application data
adj_mat_path="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/lymphoma_adj_mat_1988_2003_0-29.RData" # path to adjacency matrix
dir_out="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Lymphoma/Results/" # directory output path for results
dir_out_paper="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Lymphoma/Manuscript/" # directory output path for results in manuscript


# Load packages
library(tidyverse)
library(xtable)
library(CARBayes)
library(CARBayesST)

# Load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/Functions.R')
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/application_functions.R')

# Load data
load(data_path)
data_full_hisp <- data_full_hisp_noZero


# Load spacial adjacency matrix
load(adj_mat_path)

############
## Set-Up ##
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

#######################
## Model observed Y1 ##
#######################
trt_df <- data_full_hisp %>% filter(C == 1)
Y <- trt_df$CL_CASES
HISP <- trt_df$pct_hispanic_pop
POP <- trt_df$POP
W <- partial_us_adj_mat[which(rownames(partial_us_adj_mat) %in% trt_fips),
                        which(rownames(partial_us_adj_mat) %in% trt_fips)]

# Re-order data 
# CARBayesST requires the data to be ordered so that the first K data points 
# relate to all the spatial units for time period 1 and so on
trt_df <- arrange(trt_df, YEAR_DX)
formula <- Y ~ HISP + offset(log(POP))

model <- ST.CARar(formula=formula, family="poisson",
         W=W, burnin=1000, n.sample=2000, AR = 1)

path <- "/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Lymphoma/Results"

# Save results
save(model, file=paste0(path,"/Y1_model.RData"))


##################################
## Get fitted values from model ##
##################################

fitted_values <- model$fitted.values
head(fitted_values) # not integers

trt_df$fitted_values <- fitted_values

# Save results
save(trt_df, file=paste0(path,"/Y1_fitted_values_df.RData"))

