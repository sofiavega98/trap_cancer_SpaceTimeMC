# This script combines the results from script 3 to make a combined table
# with the results from multiple latent factors for the smoothed results

# Load packages
library(tidyverse)
library(kableExtra)

# Load data
load("~/cl_quasi_exp/Application/Leukemia/Manuscript/ATT_table1995_Smooth.RData")
df <- table

load("~/cl_quasi_exp/Application/Leukemia/Manuscript/ATT_table1995_Smooth_SensAnalysis.RData")
df_SensAnalysis <- table


Model <- rep(c("Athey MC","Vanilla","Space","Space-Time ICAR","Space-Time AR","Space-Time Lasso","Space-Time Shrinkage"),2)
K <- c(rep(2,6),3,rep(3,6),4)
results_df <- data.frame( Model, K,
                          rbind(df[,2:7], 
                                df_SensAnalysis[,2:7])) %>% 
  mutate(Model = as_factor(Model),
         Model = fct_relevel(Model, "Athey MC", "Vanilla", "Space", "Space-Time ICAR", "Space-Time AR", "Space-Time Lasso", "Space-Time Shrinkage")
  ) %>%
  arrange(Model)

rownames(results_df) <- NULL
colnames(results_df) <- c("Model", "ATT", "95% CI", "ATT", "95% CI", "ATT", "95% CI")
table <- kbl(results_df, format = "latex") %>%
  kable_classic() %>%
  add_header_above(c(" " , "California" = 2, "Connecticut" = 2, "Overall" = 2))

save_kable(table, file =  "~/cl_quasi_exp/Application/Leukemia/Manuscript/ATT_trt1995_smooth_combined.tex",float = FALSE)
