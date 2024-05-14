# This script takes the results from 2_process_results.R and creates
# tables and figures for the paper

years=1988:2003
trt_year=1995 # treated year
data_path="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/Leukemia_1988_2003_0-9.RData" # path to application data
adj_mat_path="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/leukemia_adj_mat_1988_2003_0-9.RData" # path to adjacency matrix
dir_out="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Leukemia/Results/" # directory output path for results
dir_out_paper="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Leukemia/Manuscript/" # directory output path for results in manuscript

# Load packages
library(tidyverse)
library(xtable)
library(kableExtra)
library(matrixStats)
library(cowplot)

# Load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/Functions.R')
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/application_functions.R')

# Load data
load(data_path)
data_full_hisp <- data_full_hisp_noZero

# Load processed results
load(paste0(dir_out,"processed_res_trt",trt_year,".RData")) #non-smooth
res_nonSmooth <- res
load(paste0(dir_out,"processed_res_trt",trt_year,"_smooth.RData")) #smooth
res_Smooth <- res


############
## Set-Up ##
############
# get treated counties
trt_fips <- data_full_hisp %>% filter(C == 1) %>%
  distinct(FIPS) %>% pull(FIPS)

# get Y0_obs
Y0_obs <- get_Y0_obs(data_full_hisp,trt_fips,trt_year = trt_year)[[1]]

# Load modeled Y1 values
load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Leukemia/Results/Y1_fitted_values.RData")
Y1 <- model$fitted.values

# put modeled Y1 values into the observed Y1s
data_full_hisp$CL_CASES[which(data_full_hisp$C == 1)] <- Y1

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


################
## ATT Tables ##
################

# Non-Smooth
table <- ATT_CACT_overall_table(data_full_hisp, res_nonSmooth$Mu_trt_ori, res_nonSmooth$Mu_trt_space, 
               res_nonSmooth$Mu_trt_ICAR, res_nonSmooth$Mu_trt_AR, 
               res_nonSmooth$Mu_trt_shrink, res_nonSmooth$Mu_trt_lasso, 
               res_nonSmooth$fit_gsc,trt_year)
ATT_table <- save_kable(table, file =  paste0(dir_out_paper,"ATT_table_trt",trt_year,"_nonSmooth_modeled.tex"),float = FALSE)

# Smooth
table <- ATT_CACT_overall_df(data_full_hisp, res_Smooth$Mu_trt_ori, res_Smooth$Mu_trt_space, 
                             res_Smooth$Mu_trt_ICAR, res_Smooth$Mu_trt_AR, 
                             res_Smooth$Mu_trt_shrink, res_Smooth$Mu_trt_lasso, 
                             res_Smooth$fit_gsc,trt_year)
#ATT_table <- save_kable(table, file =  paste0(dir_out_paper,"ATT_table_trt",trt_year,"_Smooth_SensAnalysis.tex"),float = FALSE)
# Save to combine in a different script
save(table, file = paste0(dir_out_paper,"ATT_table",trt_year,"_Smooth_modeled.RData"))

#########################
## ATT Plots Over-Time ##
#########################

# Non-smooth; Model = AR

ATT_full_plot(data_full_hisp, res_nonSmooth$Mu_trt_AR, Model = "AR",trt_year,years) %>%
  ggsave(filename = paste0(dir_out_paper,"ATT_plot_AR_trt",trt_year,"_nonSmooth_modeled.png"),
         height = 5,
         width = 12,
         units = "in")

# Smooth; Model = AR

ATT_full_plot(data_full_hisp, res_Smooth$Mu_trt_AR, Model = "AR",trt_year,years) %>%
  ggsave(filename = paste0(dir_out_paper,"ATT_plot_AR_trt",trt_year,"_Smooth_modeled.png"),
         height = 5,
         width = 12,
         units = "in")


#############################
## Raw vs. Estimated Plots ##
#############################
# get raw values
raw <- data_full_hisp$CL_CASES[1:(n_trt*m)]

raw_v_est_plot(data_full_hisp,pop_trt,res_nonSmooth$Mu_trt_ori,treated,n_trt,m,raw)

raw_v_est_plot(data_full_hisp,pop_trt,res_Smooth$Mu_trt_ori_smooth,treated,n_trt,m,raw)


raw_v_est_plot(data_full_hisp,pop_trt,res_nonSmooth$Mu_trt_space,treated,n_trt,m,raw)

raw_v_est_plot(data_full_hisp,pop_trt,res_Smooth$Mu_trt_space_smooth,treated,n_trt,m,raw)


raw_v_est_plot(data_full_hisp,pop_trt,res_nonSmooth$Mu_trt_ICAR,treated,n_trt,m,raw)

raw_v_est_plot(data_full_hisp,pop_trt,res_Smooth$Mu_trt_ICAR_smooth,treated,n_trt,m,raw)


raw_v_est_plot(data_full_hisp,pop_trt,res_nonSmooth$Mu_trt_AR,treated,n_trt,m,raw)

raw_v_est_plot(data_full_hisp,pop_trt,res_Smooth$Mu_trt_AR_smooth,treated,n_trt,m,raw)


raw_v_est_plot(data_full_hisp,pop_trt,res_nonSmooth$Mu_trt_lasso,treated,n_trt,m,raw)

raw_v_est_plot(data_full_hisp,pop_trt,res_Smooth$Mu_trt_lasso_smooth,treated,n_trt,m,raw)


raw_v_est_plot(data_full_hisp,pop_trt,res_nonSmooth$Mu_trt_shrink,treated,n_trt,m,raw)

raw_v_est_plot(data_full_hisp,pop_trt,res_Smooth$Mu_trt_shrink_smooth,treated,n_trt,m,raw)
