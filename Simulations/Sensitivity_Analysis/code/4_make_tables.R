# This script makes and saves tables used in the paper

#############################################################################################################
## 0. Read command line arguments and load packages and functions                                          ##
##    Arguments should include:                                                                            ##
##    - rare: Rare or nonRare                                                                              ##
##    - smooth: smooth or nonSmooth                                                                        ##
##    - dir_in: directory input path for results from script 2                                             ##
##    - dir_out: directory output path for results                                                         ##
#############################################################################################################

dir_in <- "/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/Sensitivity_Analysis/Processed_Results/"
dir_in_freq <- "/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/Main/Processed_Results/"

library(tidyverse)
library(kableExtra)

#############################
## 1. Load and Set Up Data ##
#############################

# Convergence Vectors
load(paste0(dir_in, "conv_Rare.RData"))
conv_Rare <- data.frame(conv)
load(paste0(dir_in, "conv_nonRare.RData"))
conv_nonRare <- data.frame(conv)
load(paste0(dir_in, "conv_Rare_smooth.RData"))
conv_Rare_smooth <- data.frame(conv)
load(paste0(dir_in, "conv_nonRare_smooth.RData"))
conv_nonRare_smooth <- data.frame(conv)

# Results from the posterior predictive
load(paste0(dir_in, "results_Rare_Ypred.RData"))
results_Rare <- data.frame(df)
load(paste0(dir_in, "results_nonRare_Ypred.RData"))
results_nonRare <- data.frame(df)
load(paste0(dir_in, "results_Rare_smooth_Ypred.RData"))
results_Rare_smooth <- data.frame(df)
load(paste0(dir_in, "results_nonRare_smooth_Ypred.RData"))
results_nonRare_smooth <- data.frame(df)


# Load frequensist results
load(paste0(dir_in_freq,"results_Rare_freq.Rdata"))
results_Rare_freq <- df
load(paste0(dir_in_freq,"results_nonRare_freq.Rdata"))
results_nonRare_freq <- df
load(paste0(dir_in_freq,"results_Rare_smooth_freq.Rdata"))
results_Rare_freq_smooth <- df
load(paste0(dir_in_freq,"results_nonRare_smooth_freq.Rdata"))
results_nonRare_smooth_freq <- df

####################################
## 2. Make LaTeX tables for paper ##
####################################

## Convergence Table

conv_df <- data.frame(conv_nonRare,conv_nonRare_smooth$conv,conv_Rare$conv,conv_Rare_smooth$conv) %>%
  mutate(Model = as_factor(Model),
         Model = fct_relevel(Model, "Original", "Space", "Space-Time \n ICAR", "Space-Time \n AR(1)", "Shrinkage \n Lasso", "Space-Time \n Shrinkage"))%>%
  arrange(Model) 

K <- c(rep(c(1,3,7),5),7)
conv_df <- data.frame(conv_df$Model, K,conv_df[,c(3:6)])
colnames(conv_df) <- c("Model","K","Non-Rare & Non-Smoothed", "Non-Rare & Smoothed",
                       "Rare & Non-Smoothed", "Rare & Smoothed")

conv_df <- conv_df %>% kbl(format = "latex") %>%
  kable_classic()


conv_table <- save_kable(conv_df, file =  paste0(dir_in,"conv_table.tex"),float = FALSE)


## Absolute percent bias of counterfactual table
Model <- c(results_nonRare$Model, rep("Athey MC",3), "GSC", "ASC")
K <- c(results_nonRare$k, 1,3,7, "NA", "NA")
results_df <- data.frame( Model, K,
                         c(results_nonRare$absPercBias_Y0, t(round(results_nonRare_freq,2))),
                         c(results_nonRare_smooth$absPercBias_Y0, t(round(results_nonRare_smooth_freq,2))),
                         c(results_Rare$absPercBias_Y0, t(round(results_Rare_freq,2))),
                         c(results_Rare_smooth$absPercBias_Y0, t(round(results_Rare_freq_smooth,2)))) %>% 
  mutate(Model = as_factor(Model),
         Model = fct_relevel(Model, "Vanilla", "Space", "Space-Time ICAR", "Space-Time AR", "Lasso", "Space-Time Shrinkage")
         ) %>%
              arrange(Model)

view(results_df)
colnames(results_df) <- c("Model","K","Non-Rare & Non-Smoothed", "Non-Rare & Smoothed",
                          "Rare & Non-Smoothed", "Rare & Smoothed")

results_df <- results_df %>%
  kbl(format = "latex") %>%
  kable_classic()


results_table <- save_kable(results_df, file =  paste0(dir_in,"absPercBiasY0_table.tex"),float = FALSE)
