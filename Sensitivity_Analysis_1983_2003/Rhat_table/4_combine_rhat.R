# This script takes the Rhat values from both data sets and combines them

# Load packages
library(tidyverse)
library(kableExtra)


# Load functions
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/Functions.R')
source('/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/functions/application_functions.R')

# Load data
load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Sensitivity_Analysis_1983_2003/Leukemia/Results/rhat_df_trt1990.RData")
df_Leukemia <- df

load("/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Sensitivity_Analysis_1983_2003/Lymphoma/Results/rhat_df_trt1990.RData")
df_Lymphoma <- df

df <- data.frame(df_Lymphoma, df_Leukemia$`Non-Smoothed`,df_Leukemia$Smoothed)

df_table <- kbl(round(df,2), format = "latex") %>%
  kable_classic()

save_kable(df_table, file =  "/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Sensitivity_Analysis_1983_2003/Rhat_table/rhat_table.tex",float = FALSE)
