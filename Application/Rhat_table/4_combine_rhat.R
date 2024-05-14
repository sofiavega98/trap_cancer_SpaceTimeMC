# This script takes the Rhat values from both data sets and combines them

# Load packages
library(tidyverse)
library(kableExtra)


# Load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/Functions.R')
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/application_functions.R')

# Load data
load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Leukemia/Results/rhat_df_trt1995.RData")
df_Leukemia <- df

load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Lymphoma/Results/rhat_df_trt1995.RData")
df_Lymphoma <- df

df <- data.frame(df_Lymphoma, df_Leukemia$`Non-Smoothed`,df_Leukemia$Smoothed)

df_table <- kbl(round(df,2), format = "latex") %>%
  kable_classic()

save_kable(df_table, file =  "/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Rhat_table/rhat_table.tex",float = FALSE)

########################################
## Sensitivity Analysis: Increasing K ##
########################################

# Here we combine Rhat tables to make a table with all Rhat values

# Load data
load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Leukemia/Results/rhat_df_trt1995_SensAnalysis.RData")
df_Leukemia_SensAnalysis <- df

load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Lymphoma/Results/rhat_df_trt1995_SensAnalysis.RData")
df_Lymphoma_SensAnalysis <- df


Model <- rep(c("Vanilla","Space","Space-Time ICAR","Space-Time AR","Space-Time Lasso","Space-Time Shrinkage"),2)
results_df <- data.frame( Model, 
                          rbind(round(df_Lymphoma,2),round(df_Lymphoma_SensAnalysis,2)), 
                                rbind(round(df_Leukemia[,2:3],2),round(df_Leukemia_SensAnalysis[,2:3],2))) %>% 
  mutate(Model = as_factor(Model),
         Model = fct_relevel(Model, "Vanilla", "Space", "Space-Time ICAR", "Space-Time AR", "Space-Time Lasso", "Space-Time Shrinkage")
  ) %>%
  arrange(Model)

rownames(results_df) <- NULL
colnames(results_df)[3:6] <- c("Lymphoma & Non-Smoothed","Lymphoma & Smoothed",
                               "Leukemia & Non-Smoothed","Leukemia & Smoothed")

results_df <- kbl(results_df, format = "latex") %>%
  kable_classic()

save_kable(results_df, file =  "/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Application/Rhat_table/rhat_table_combined.tex",float = FALSE)

