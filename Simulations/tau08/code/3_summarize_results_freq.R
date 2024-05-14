# This script takes the simulation results, calculates absolute percent biases
# of the ATT and counterfactual, and prepares them to inputted into a table.
# All absolute percent biases are calculated in the same format as ASC from
# Ben-Michael et al.

# Load libraries
library(tidyverse)
library(augsynth)

# load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/Functions.R')
# load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/freq_post_processing_functions.R')

# load NM set up script
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/0_set_up.R') # We are basing the simulated data on the values from this script

# Directory to output results
dir_out <- "/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/tau08/Processed_Results"

#####################
## 1. Load results ##
#####################

load(paste0("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/tau08/Results/sim_Rare_freq.RData"))
results_Rare <- results
load(paste0("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/tau08/Results/sim_nonRare_freq.RData"))
results_nonRare <- results

load(paste0("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/tau08/Results/sim_Rare_freq_smooth.RData"))
results_Rare_smooth <- results
load(paste0("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/tau08/Results/sim_nonRare_freq_smooth.RData"))
results_nonRare_smooth <- results

rm(results)

###############
## 2. Set Up ##
###############

# Specify the indices of the treated units
ind <- c()
for(i in 0:5){
  ind <- append(ind,seq(9,15) + (i*15))
}

# Specify the indices of the treated units in full matrix
ind <- c()
for(i in 0:5){
  ind <- append(ind,seq(9,15) + (i*15))
}

# Get true treatment effects
true_trt_effect_Rare <- results_Rare$Y1_matrix - results_Rare$Y0_matrix[ind,]
true_trt_effect_nonRare <- results_nonRare$Y1_matrix - results_nonRare$Y0_matrix[ind,]

# Number of treated units
n_trt <- 6 # how we simulated the data

# Populations for treated counties
pop_trt <- as.vector(t(pop_mat[1:n_trt,]))

##########################
## 3. Summarize Results ##
##########################

# Rare and non-smooth
results <- results_Rare

# Matrix completion in gsynth
gsc_mc_1 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_1,results$Y0_matrix,ind)$median
gsc_mc_3 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_3,results$Y0_matrix,ind)$median
gsc_mc_7 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_7,results$Y0_matrix,ind)$median

# GSC in gsynth
gsc <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc,results$Y0_matrix,ind)$median

# ASC in augsynth
asc <- medianAbsPercBiasY0_year_asc(pop_trt,results$fit_asc,results$Y0_matrix,ind)$median

# ALS
als_1 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_1,results$Y0_matrix,ind)$median
als_3 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_3,results$Y0_matrix,ind)$median
als_7 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_7,results$Y0_matrix,ind)$median

# softImpute
softImpute <-medianAbsPercBiasY0_year_fill_soft(pop_trt,results$fit_softImpute,results$Y0_matrix,ind)$median

# Nuclear Norm minimization
nuclear <- medianAbsPercBiasY0_year_fill(pop_trt,results$fit_nuclear,results$Y0_matrix,ind)$median

# SVT
svt <- medianAbsPercBiasY0_year_fill(pop_trt,results$fit_svt,results$Y0_matrix,ind)$median

# Save results
df <- data.frame(gsc_mc_1, gsc_mc_3, gsc_mc_7, gsc, asc, 
                 als_1, als_3, als_7, softImpute,
                 nuclear, svt)
save(df, file = paste0(dir_out,"/results_Rare_freq.Rdata"))

#-------------------------------------------------------------------------

# nonRare and non-smooth
results <- results_nonRare

# Matrix completion in gsynth
gsc_mc_1 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_1,results$Y0_matrix,ind)$median
gsc_mc_3 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_3,results$Y0_matrix,ind)$median
gsc_mc_7 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_7,results$Y0_matrix,ind)$median

# GSC in gsynth
gsc <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc,results$Y0_matrix,ind)$median

# ASC in augsynth
asc <- medianAbsPercBiasY0_year_asc(pop_trt,results$fit_asc,results$Y0_matrix,ind)$median

# ALS
als_1 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_1,results$Y0_matrix,ind)$median
als_3 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_3,results$Y0_matrix,ind)$median
als_7 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_7,results$Y0_matrix,ind)$median

# softImpute
softImpute <-medianAbsPercBiasY0_year_fill_soft(pop_trt,results$fit_softImpute,results$Y0_matrix,ind)$median

# Nuclear Norm minimization
nuclear <- medianAbsPercBiasY0_year_fill(pop_trt,results$fit_nuclear,results$Y0_matrix,ind)$median

# SVT
svt <- medianAbsPercBiasY0_year_fill(pop_trt,results$fit_svt,results$Y0_matrix,ind)$median

# Save results
df <- data.frame(gsc_mc_1, gsc_mc_3, gsc_mc_7, gsc, asc, 
                 als_1, als_3, als_7, softImpute,
                 nuclear, svt)
save(df, file = paste0(dir_out,"/results_nonRare_freq.Rdata"))

#------------------------------------------------------------------------------
# Rare and smooth
results <- results_Rare_smooth

# Matrix completion in gsynth
gsc_mc_1 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_1,results$Y0_matrix,ind)$median
gsc_mc_3 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_3,results$Y0_matrix,ind)$median
gsc_mc_7 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_7,results$Y0_matrix,ind)$median

# GSC in gsynth
gsc <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc,results$Y0_matrix,ind)$median

# ASC in augsynth
asc <- medianAbsPercBiasY0_year_asc(pop_trt,results$fit_asc,results$Y0_matrix,ind)$median

# ALS
als_1 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_1,results$Y0_matrix,ind)$median
als_3 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_3,results$Y0_matrix,ind)$median
als_7 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_7,results$Y0_matrix,ind)$median

# softImpute
softImpute <-medianAbsPercBiasY0_year_fill_soft(pop_trt,results$fit_softImpute,results$Y0_matrix,ind)$median

# Nuclear Norm minimization
nuclear <- medianAbsPercBiasY0_year_fill(pop_trt,results$fit_nuclear,results$Y0_matrix,ind)$median

# SVT
svt <- medianAbsPercBiasY0_year_fill(pop_trt,results$fit_svt,results$Y0_matrix,ind)$median

# Save results
df <- data.frame(gsc_mc_1, gsc_mc_3, gsc_mc_7, gsc, asc, 
                 als_1, als_3, als_7, softImpute, nuclear, svt)
save(df, file = paste0(dir_out,"/results_Rare_smooth_freq.Rdata"))

#-------------------------------------------------------------------------

# nonRare and non-smooth
results <- results_nonRare_smooth

# Matrix completion in gsynth
gsc_mc_1 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_1,results$Y0_matrix,ind)$median
gsc_mc_3 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_3,results$Y0_matrix,ind)$median
gsc_mc_7 <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc_mc_7,results$Y0_matrix,ind)$median

# GSC in gsynth
gsc <- medianAbsPercBiasY0_year_gsc(pop_trt,results$fit_gsc,results$Y0_matrix,ind)$median

# ASC in augsynth
asc <- medianAbsPercBiasY0_year_asc(pop_trt,results$fit_asc,results$Y0_matrix,ind)$median

# ALS
als_1 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_1,results$Y0_matrix,ind)$median
als_3 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_3,results$Y0_matrix,ind)$median
als_7 <-medianAbsPercBiasY0_year_als(pop_trt,results$fit_als_7,results$Y0_matrix,ind)$median

# softImpute
softImpute <-medianAbsPercBiasY0_year_fill_soft(pop_trt,results$fit_softImpute,results$Y0_matrix,ind)$median

# Nuclear Norm minimization
nuclear <- medianAbsPercBiasY0_year_fill(pop_trt,results$fit_nuclear,results$Y0_matrix,ind)$median

# SVT
svt <- medianAbsPercBiasY0_year_fill(pop_trt,results$fit_svt,results$Y0_matrix,ind)$median

# Save results
df <- data.frame(gsc_mc_1, gsc_mc_3, gsc_mc_7, gsc, asc, 
                 als_1, als_3, als_7, softImpute,
                 nuclear, svt)
save(df, file = paste0(dir_out,"/results_nonRare_smooth_freq.Rdata"))
