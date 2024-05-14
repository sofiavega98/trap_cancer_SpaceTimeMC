# This script takes the results from 2_process_results.R and creates
# tables and figures for the paper

years=1983:2003
trt_year=1990 # treated year
data_path="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/Lymphoma_1983_2003_0-29.RData" # path to application data
adj_mat_path="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/data/lymphoma_adj_mat_1983_2003_0-29.RData" # path to adjacency matrix
dir_out="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Sensitivity_Analysis_1983_2003/Lymphoma/Results/" # directory output path for results
dir_out_paper="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Sensitivity_Analysis_1983_2003/Lymphoma/Manuscript/" # directory output path for results in manuscript

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
Y0_obs <- get_Y0_obs(data_full_hisp,trt_fips,trt_year = "1995")[[1]]

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


#################
## Diagnostics ##
#################

# R-hat table

# Load R-hat vectors
load(paste0(dir_out,"rhat_vec_trt",trt_year,".RData"))
rhat_vec_nonSmooth <- rhat_vec
load(paste0(dir_out,"rhat_vec_trt",trt_year,"_smooth.RData"))
rhat_vec_Smooth <- rhat_vec

# Combine columns
df <- data.frame(K = as.integer(c(rep(2,5),3)), "Non-Smoothed" = rhat_vec_nonSmooth,"Smoothed" = rhat_vec_Smooth)
colnames(df) <- c("K","Non-Smoothed", "Smoothed")

# Save to combine with Leukemia in a different script
save(df, file = paste0(dir_out,"rhat_df_trt",trt_year,".RData"))

################
## ATT Tables ##
################

# Non-Smooth
table <- ATT_CACT_overall_table(data_full_hisp, res_nonSmooth$Mu_trt_ori, res_nonSmooth$Mu_trt_space, 
               res_nonSmooth$Mu_trt_ICAR, res_nonSmooth$Mu_trt_AR, 
               res_nonSmooth$Mu_trt_shrink, res_nonSmooth$Mu_trt_lasso, 
               res_nonSmooth$fit_gsc,trt_year)
ATT_table <- save_kable(table, file =  paste0(dir_out_paper,"ATT_table_trt",trt_year,"_nonSmooth.tex"),float = FALSE)

# Smooth
table <- ATT_CACT_overall_table(data_full_hisp, res_Smooth$Mu_trt_ori, res_Smooth$Mu_trt_space, 
                                res_Smooth$Mu_trt_ICAR, res_Smooth$Mu_trt_AR, 
                                res_Smooth$Mu_trt_shrink, res_Smooth$Mu_trt_lasso, 
                                res_Smooth$fit_gsc,trt_year)
ATT_table <- save_kable(table, file =  paste0(dir_out_paper,"ATT_table_trt",trt_year,"_Smooth.tex"),float = FALSE)

#########################
## ATT Plots Over-Time ##
#########################

# Non-smooth; Model = AR

ATT_full_plot(data_full_hisp, res_nonSmooth$Mu_trt_AR, Model = "AR",trt_year,years) %>%
  ggsave(filename = paste0(dir_out_paper,"ATT_plot_AR_trt",trt_year,"_nonSmooth.png"),
         height = 5,
         width = 12,
         units = "in")

# Smooth; Model = AR

ATT_full_plot(data_full_hisp, res_Smooth$Mu_trt_AR, Model = "AR",trt_year,years) %>%
  ggsave(filename = paste0(dir_out_paper,"ATT_plot_AR_trt",trt_year,"_Smooth.png"),
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

raw_v_est_plot_svt(data_full_hisp,pop_trt,res_nonSmooth$fit_svt$X,treated,n_trt,m,raw)

raw_v_est_plot_gsc(data_full_hisp,pop_trt,res_nonSmooth$fit_gsc,treated,n_trt,m,raw)

####################################################
## Make Scree-Plots to Determine k For Each Model ##
####################################################

# Perform PCA
data_cases <- matrix(data_full_hisp$CL_CASES, nrow=length(unique(data_full_hisp$FIPS)), byrow = TRUE)
pca_results <- prcomp(data_cases, scale = TRUE)

# Calculate total variance explained by each principal component
var_explained = pca_results$sdev^2 / sum(pca_results$sdev^2)

# Create scree plot (save 544 347)
qplot(c(1:length(var_explained)), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  theme_bw() +
  #ggtitle("Scree Plot") +
  ylim(0, 1) 