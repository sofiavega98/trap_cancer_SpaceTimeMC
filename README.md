# trap_cancer_SpaceTimeMC

Code for running analysis in: https://arxiv.org/abs/2307.09546


This respository contains code used to conduct simulation studies comparing 6 Bayesian methods for matrix completion (Vanilla, Space, Space-Time ICAR, Space-Time AR, Space-Time Lasso, and Space-Time Shrinkage) with 3 exisiting methods (Athey's matrix completion (MC) found in the $\texttt{gysnth}$ R package, generalized synthetic control (GSC) found in the $\texttt{gysnth}$ R package, and augmented synthetic control (ASC) found in the $\texttt{augysnth}$ R package) on rare and non-rare simulated data. 

This repository also contains an application to childhood and young adult hematologic cancers.

Scripts containing functions used for the Simulations and Applications are found in `\functions\`.

## Simulations
Code related to running simulations is found in `\Simulations\`
* `0_set_up.R` is script to set up simulations. This script identifies the NM counties in the SEER data, finds the partial adjacency matrix, and identifies 6 counties as treated
* `1_sim_bayes.R` is a script to generate simulated data and run the data on models 1-6 from the paper.
* `1_sim_bayes_smooth.R` is a script to generate simulated data and run the data on models 1-6 from the paper, but includes a temporal smoothing pre-processing step.
* `1_sim_freq.R` is a script to generate simulated data and run the data on the frequentist models Athey MC, GSC, and ASC.
* `1_sim_freq_smooth.R` is a script to generate simulated data and run the data on the frequentist models Athey MC, GSC, and ASC. but includes a temporal smoothing pre-processing step.
* `2_process_results.R` is a script to take the individual Bayesian simulation results and saves them as one
* `3_summarize_results_filter.R` is a script to take the simulation results, calculate absolute percent biases of the counterfactual, and prepare results to inputted into a table. All absolute percent biases are calculated in the same format as ASC from Ben-Michael et al.
* `\code\3_summarize_results_freq.R` is a script to take the frequentist simulation results, calculate absolute percent biases of the counterfactual, and prepare results to inputted into a table. All absolute percent biases are calculated in the same format as ASC from Ben-Michael et al.
* `\code\4_make_tables.R` is a script to create LaTeX tables used in the paper.
* `\tau1\` contains shell scripts ending in ".sh" to launch jobs on the FASRC cluster for the results in the main paper.
* `\tau08\` contains shell scripts ending in ".sh" to launch jobs on the FASRC cluster for the results in the supplemental section of the paper.

## Application
Code related to running the models on cancer data is found in `\Application\` and `\Sensitivity_Analysis_1983_2003\`. Within each folder, are folders for each outcome, `/Leukemia/` and `/Lymphoma/` with associated code and shell scripts inside. 
* `/code/0_cleaning_data.R` is a script to clean the dat sets.
* `/code/1_run_model.R` is script to run the models on the application data.
* `/code/1_run_models_smooth.R` is script to run the models on the application data with a temporal pre-processing step.
* `/code/2_process_results.R` is a script that processes in from the application. It also runs the application the MC through $\texttt{gsynth}$ package.
* `/code/2_process_results_smooth.R` is a script that processes results in the application with the temporal smoothing pre-processing step. It also runs the application the MC through $\texttt{gsynth}$ package.
* `/code/2_process_results_SensAnalysis.R` is a script that processes results in the application with increased K. It also runs the application the MC through $\texttt{gsynth}$ package.
* `/code/2_process_results_smooth_SensAnalysis.R` is a script that processes results in the application with increased K and with the temporal smoothing pre-processing step. It also runs the application the MC through $\texttt{gsynth}$ package.
* `/code/3_tables_figures.R` is a script that takes the results from `2_process_results.R` and `2_process_results_smooth.R` and create tables and figures for the paper.
* `/code/3_tables_figures_SensAnalysis.R` is a script that takes the results from `2_process_results_SensAnalysis.R` and `2_process_results_smooth_SensAnalysis.R` and create tables and figures for the paper.
* `/Rhat_table/4_combine_rhat.R` is a script that takes Rhat values from the models applied to both data sets and combines results in a table.
* `/code/5_combine_tables_smooth.R` is script that combines the results from `3_tables_figures.R` and `3_tables_figures_SensAnalysis.R` to make a combined table with the results from multiple choices of latent factors for the smoothed results.

## STAN Models
STAN code related to Models 1-6 is found in `\STAN Models\`.
* `Vanilla.stan` is a script to run Model 1 (Vanilla) for the simulations.
* `Space.stan` is a script to run Model 2 (Space) for the simulations.
* `SpaceTimeICAR.stan` is a script to run Model 3 (Space-Time ICAR) for the simulations.
* `SpaceTimeAR.stan` is a script to run Model 4 (Space-Time AR) for the simulations.
* `SpaceTimeLasso.stan` is a script to run Model 5 (Space-Time Lasso) for the simulations.
* `SpaceTimeShrinkage.stan` is a script to run Model 6 (Space-Time Shrinkage) for the simulations.


* `Vanilla_cov.stan` is a script to run Model 1 (Vanilla) for the application. This script allows for the inclusion of covariates.
* `Space_cov.stan` is a script to run Model 2 (Space) for the application. This script allows for the inclusion of covariates.
* `SpaceTimeICAR_cov.stan` is a script to run Model 3 (Space-Time ICAR) for the application. This script allows for the inclusion of covariates.
* `SpaceTimeAR_cov.stan` is a script to run Model 4 (Space-Time AR) for the application. This script allows for the inclusion of covariates.
* `SpaceTimeLasso_cov.stan` is a script to run Model 5 (Space-Time Lasso) for the application. This script allows for the inclusion of covariates.
* `SpaceTimeShrinkage_cov.stan` is a script to run Model 6 (Space-Time Shrinkage) for the application. This script allows for the inclusion of covariates.


  

