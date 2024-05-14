# This script takes the inidividual simulation results and saves them as one

#############################################################################################################
## 0. Read command line arguments and load packages and functions                                          ##
##    Arguments should include:                                                                            ##
##    - est_k: estimated number of latent factors                                                          ##
##    - rare: Rare or nonRare                                                                              ##
##    - smooth: smooth or nonSmooth                                                                        ##
##    - dir_out: directory output path for results                                                         ##
#############################################################################################################
# Read command line arguments 

args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }

library(tidyverse)

#est_k=1
#rare="Rare"
#smooth="smooth"
#dir_out="/n/holylfs05/LABS/nethery_lab/Users/svega/trap_cancer_mc/Simulations/tau08/Results" # directory output path for results


if(smooth == "smooth"){
  filenames <- mapply(c,
                      expand_grid(y=c(seq(1,100))) %>% unname()) %>% 
    as_tibble(.name_repair = 'unique') %>%
    rename(
           simnum2 = 1) %>%
    mutate(filename = paste0(dir_out,"/sim_k",est_k,"_",simnum2,"_",rare,"_smooth.RData")) %>%
    select(filename) %>% as_vector() 
  
}else{
  filenames <- mapply(c,
                      expand_grid(y=c(seq(1,100))) %>% unname()) %>% 
    as_tibble(.name_repair = 'unique') %>%
    rename(
      simnum2 = 1) %>%
    mutate(filename = paste0(dir_out,"/sim_k",est_k,"_",simnum2,"_",rare,".RData")) %>%
    select(filename) %>% as_vector() 
}


data_big <- lapply(filenames, function(x) get(load(x)))

names(data_big) <- filenames %>% str_extract(., "job_\\d_\\d")

if(smooth == "smooth"){
  filename_out = paste0(dir_out,"/sim_k",est_k,"_",rare,"_smooth.RData")
}else{
  filename_out = paste0(dir_out,"/sim_k",est_k,"_",rare,".RData")
}

save(data_big, file=filename_out)
