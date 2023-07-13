# Script of Functions primarily used in pre-processing

## Get Y1_obs1:
# This functions pulls the observed treated values using the indicator variable
get_Y1_obs1 <- function(data) {
  #Input:
  ## data: the data set with columns FIPS,CL_CASES, YEAR_DX, C
  ## C is an indicator column indicating treatment
  #Output:
  ## Y1_obs: a vector of treated outcomes
  
  Y1_obs1 <- data[which(data$C==1),]$CL_CASES
  
  return(Y1_obs1)
}


#----------------------------------------------------------------------
## Get Y0_obs1 matrix: 
### This obtains a matrix of observed non treated values
### Columns will be years, rows will be counties. 
### NAs will be placed where counties are treated

get_Y0_obs <- function(data,trt_fips, 
                        trt_year = "1996") {
  # Input:
  ## data: the data set with columns FIPS, CL_CASES, YEAR_DX
  ## trt_fips: a vector of treated counties
  ## trt_year: the year where treatment starts (default "1996")
  # Output:
  ## Y0_obs1 matrix of observed nontreated values with NAs for trt values
  ## Y0_obs2 matrix of observed nontreated values with 99999 for trt values
  
  #observed non treated values
  Y0_obs1 <- pivot_wider(data,id_cols = FIPS, 
                         names_from = YEAR_DX, values_from = CL_CASES) 
  
  #put NAs where counties are "treated" (Note: the +1 accounts for the FIPS column)
  Y0_obs1[which(Y0_obs1$FIPS %in% trt_fips),
          (which(unique(data$YEAR_DX)==trt_year)+1):(length(unique(data$YEAR_DX))+1)] <-NA
  #remove FIPS column
  Y0_obs1 <- Y0_obs1[,-1]
  
  ## stan doesn't allow missing values so we put a numeric "placeholder" in the missings
  Y0_obs2 <- Y0_obs1
  Y0_obs2[is.na(Y0_obs1)]<-(99999)
  
  return(list(Y0_obs1,Y0_obs2))
}

#---------------------------------------------------------------------
## Get Population matrix:
# This functions pivots the population column in the data to be a matrix
# with counties in the rows and years in the columns

get_pop_mat <- function(data) {

  pop_mat <- matrix(data$POP,nrow=length(unique(data$FIPS)),
                    ncol=length(unique(data$YEAR_DX)),byrow=T)
  
  return(pop_mat)
}



#---------------------------------------------------------------------
## Get y_miss:

get_y_miss <- function(Y0_obs1) {
  #Input: 
  ## Y0_obs1: a matrix of observed  values with NAs for trt values
  #Output: 
  ## y_miss: a matrix with 1s for trted values and 0s for non-treated values
  y_miss <- as.matrix(is.na(Y0_obs1),
                      nrow=nrow(Y0_obs1),ncol=ncol(Y0_obs1),byrow=T)*1
  return(y_miss)
}


#--------------------------------------------------------------------

#Transform adjacency matrix to an edge list

get_edg <- function(adj_mat) {
  edg <- as.data.frame(adjmat_to_edgelist(adj_mat,undirected=getOption("diffnet.undirected")))
  #get the lower traingular part
  edg <- edg[which(edg$alter < edg$ego, arr.ind = TRUE),] 
  N = max(edg$ego)
  node1 <- edg$alter
  node2 <- edg$ego
  N_edges <- nrow(edg)
  
  return(list(edg,N,node1,node2,N_edges))
}

#------------------------------------------------------------------------------
# This functions creates temporal adjacency matrix

get_temp_adj_mat <- function(m,n) {
  temp_adj_matrix <- matrix(0,nrow=m,ncol=m)
  for(i in 1:(m-1)) {
    temp_adj_matrix[i,i+1] <- 1
    temp_adj_matrix[i+1,i] <- 1
  }
  return(temp_adj_matrix)
}


#------------------------------------------------------------------------------
## Get the space-time stan_in
# This function put the data in the format to work with our STAN models
# This function is specifically for for non-simulated data with no covariates

get_space_time_stan_in <- function(k,data,trt_fips,adj_mat,temp_adj_matrix,trt_year = "1996") {
  
  Y0_obs1 <- get_Y0_obs(data,trt_fips,trt_year)[[1]]
  y_miss <- get_y_miss(Y0_obs1)
  Y0_obs2 <- get_Y0_obs(data,trt_fips,trt_year)[[2]]
  pop_mat <- get_pop_mat(data)
  
  ## put adjacency matrices in node form ##
  N <- get_edg(adj_mat)[[2]]
  s_node1 <- get_edg(adj_mat)[[3]]
  s_node2 <- get_edg(adj_mat)[[4]]
  N_edges <- get_edg(adj_mat)[[5]]
  
  M <- get_edg(temp_adj_matrix)[[2]]
  t_node1 <- get_edg(temp_adj_matrix)[[3]]
  t_node2 <- get_edg(temp_adj_matrix)[[4]]
  M_edges <- get_edg(temp_adj_matrix)[[5]]
  
  ## set up the list of data for stan call ##
  stan_in <- list(k=k, ## k is user-specified number of factors ##
                  m = ncol(Y0_obs1), # m corresponds to T in our notation (number of time points) ##
                  n = nrow(Y0_obs1), # n is number of counties ##
                  nmiss = length(which(is.na(Y0_obs1))), #total number of missings
                  n_exp =  length(trt_fips), # number of exposed counties
                  y_miss = y_miss, ## y_miss is n-length indicator of exposure      (missingness in Y(0)) at time m ##
                  y=Y0_obs2, ## y is the observed Y(0) matrix ##
                  offs=pop_mat, ## offs is n x m matrix of population size offsets for each county ##
                  #offs=pop
                  N = N,
                  s_node1 = s_node1,
                  s_node2 = s_node2,
                  N_edges = N_edges,
                  M = M,
                  t_node1 = t_node1,
                  t_node2 = t_node2,
                  M_edges = M_edges
  )
  
  return(stan_in)
}

#---------------------------------------------------------------
## Get the space-time stan_in
# This function put the data in the format to work with our STAN models
# This function is specifically for for non-simulated data with covariates

get_space_time_cov_stan_in <- function(k,data,trt_fips,adj_mat,temp_adj_matrix,
                                       trt_year = "1996",p, x, row_cov, col_cov) {
  
  Y0_obs1 <- get_Y0_obs(data,trt_fips,trt_year)[[1]]
  y_miss <- get_y_miss(Y0_obs1)
  Y0_obs2 <- get_Y0_obs(data,trt_fips,trt_year)[[2]]
  pop_mat <- get_pop_mat(data)
  
  ## put adjacency matrices in node form ##
  N <- get_edg(adj_mat)[[2]]
  s_node1 <- get_edg(adj_mat)[[3]]
  s_node2 <- get_edg(adj_mat)[[4]]
  N_edges <- get_edg(adj_mat)[[5]]
  
  M <- get_edg(temp_adj_matrix)[[2]]
  t_node1 <- get_edg(temp_adj_matrix)[[3]]
  t_node2 <- get_edg(temp_adj_matrix)[[4]]
  M_edges <- get_edg(temp_adj_matrix)[[5]]
  
  ## set up the list of data for stan call ##
  stan_in <- list(k=k, ## k is user-specified number of factors ##
                  m = ncol(Y0_obs1), # m corresponds to T in our notation (number of time points) ##
                  n = nrow(Y0_obs1), # n is number of counties ##
                  nmiss = length(which(is.na(Y0_obs1))), #total number of missings
                  n_exp =  length(trt_fips), # number of exposed counties
                  y_miss = y_miss, ## y_miss is n-length indicator of exposure      (missingness in Y(0)) at time m ##
                  y=Y0_obs2, ## y is the observed Y(0) matrix ##
                  offs=pop_mat, ## offs is n x m matrix of population size offsets for each county ##
                  N = N,
                  s_node1 = s_node1,
                  s_node2 = s_node2,
                  N_edges = N_edges,
                  M = M,
                  t_node1 = t_node1,
                  t_node2 = t_node2,
                  M_edges = M_edges,
                  p=p,
                  x=x, 
                  row_cov=row_cov, 
                  col_cov=col_cov
                  )
  
  return(stan_in)
}




#######--------------------------------------------------

# Get spatially simulated  data
# Changing DGF to be in terms of rates per 100,000
get_sim_spatial_rate_data <- function(partial_W,temp_adj_matrix, trt_fips,pop_mat,k,alpha,
                                      gamma_sd=.15,psi_sd=.015, n_trt,m_trt, k_stan_in, rho = .6, 
                                      rate = .5, tau2 = .01, simnum = simnum) {
  #####################################################################################
  #### Input:                                                                      ####
  ####  partial_W: spatial adjacency matrix (n x n)                                ####
  ####  temp_adj_matrix: temporal adjacency matrix (m x m)                         ####
  ####  trt_fips: vector of fips codes for treated counties                        ####
  ####  pop_mat : population matrix (n x m)                                        ####
  ####  k: user specified number of latent factors                                 ####
  ####  alpha: user specified global intercept                                     ####
  ####  gamma_sd: county intercept sd                                              ####
  ####  psi_sd: time intercept sd                                                  ####
  ####  n_trt: number of treated counties                                          ####
  ####  m_trt: number of treated time points                                       ####
  ####  k_stan_in: a vector of k_max-values for the STAN model                     ####
  ####  rho: degree of space/time correlation                                      ####
  ####  rate: treatment effect in terms of rate                                    ####
  ####  simnum: seed according to simnum                                           ####
  #### Output:                                                                     ####
  ####  stan_in_list: list of data for stan call                                   ####
  ####  Y0_full: the simulated "true" Y(0) matrix with no missingness (n x m)      ####     
  ####  Y1_obs: the generated oserved Y(1) values after treatment (n_trt x m_trt)  ####
  #####################################################################################
  
  #####################################################################
  ## Generate Mean Structure: Should be the same for each simulation ## 
  #####################################################################
  
  library("mvtnorm")
  
  ## use a common seed to generate the parameters fixed across sims ##
  set.seed(4)
  #set.seed(123)
  #set.seed(10)
  
  #### Define variables  
  n <- nrow(partial_W) # number of counties
  m <- nrow(temp_adj_matrix) # number of time points
  k <- true_k # user specified number of latent factors
  n_trt <- n_trt # number of treated counties
  m_trt <- m_trt # number of treated years
  W <- partial_W # spatial adjacency matrix
  D <- temp_adj_matrix # temporal adjacency matrix
  tau2 <- tau2 #precision parameter
  
  #### Create the spatial precision matrix based on the CAR prior proposed by Leroux et al. (2000)
  Q.W <- rho * (diag(apply(W, 2, sum)) - W) + (1-rho) * diag(rep(1, n))
  #### Create the spatial covariance matrix and simulate a set of spatial random effects
  Q.W.inv <- solve(Q.W)
  
  #### Create the temporal covariance matrix and simulate a set of temporal random effects
  Q.D <- rho * (diag(apply(D, 2, sum)) - D) + (1-rho) * diag(rep(1, m))
  Q.D.inv <- solve(Q.D)
  
  #### Simulate factor loadings 
  L <- t(rmvnorm(n = k,mean = rep(0, m), sigma = tau2 * Q.D.inv )) #(m x k)
  
  #### Simulate latent factors 
  FS <- t(rmvnorm(n = k,mean = rep(0, n), sigma = tau2 * Q.W.inv )) #(n x k)
  
  #### Simulate Intercepts
  ## global intercept ##
  alpha<-alpha #(-4)
  ## county-specific intercepts ##
  gamma_i<- rnorm(n,sd=gamma_sd)
  ## time-specific intercepts ##
  psi_t<- rnorm(m,sd=psi_sd)
  
  #### Get Mean Structure
  mean_struct<-exp(alpha+
                     matrix(rep(gamma_i,m),nrow=n)+
                     matrix(rep(psi_t,n),nrow=n,byrow = T)+
                     FS%*%t(L)+
                     log(pop_mat))
  
  ###################
  ## Generate Data ##
  ###################
  ## use a job-specific seed to generate the random components ##
  set.seed(simnum)
  
  #### Generate the "true" Y(0) matrix with no missingness 
  Y0_full<-Y0_obs<-matrix(rpois(m*n,lambda=c(mean_struct)),nrow=n,ncol=m)
  
  #### Insert missingness into the Y(0) matrix to correspond to the observed Y(0) ##
  Y0_obs[1:n_trt,(m-m_trt+1):m]<-NA
  
  #### Convert full counterfactual matrix into rate per 100000
  Y0_full_rate <- (Y0_full/pop_mat) * 100000
  
  #### Generate the Y(1) in terms of rate per 100000
  Y1_obs_rate <- Y0_full_rate
  Y1_obs_rate[1:n_trt,(m-m_trt+1):m] <- Y0_full_rate[1:n_trt,(m-m_trt+1):m] + rate
  
  #### Convert back into counts
  Y1_obs_mat <- (Y1_obs_rate/100000) * pop_mat
  Y1_obs <- Y1_obs_mat[1:n_trt,(m-m_trt+1):m]
  
  ##########################
  ## Set Up Data for STAN ##
  ##########################
  #### Get variables for STAN 
  y_miss <- get_y_miss(Y0_obs)
  Y0_full2 <- Y0_full ## STAN doesn't allow missingness
  Y0_full2[is.na(Y0_obs)]<-(99999)
  
  N <- get_edg(partial_W)[[2]]
  s_node1 <- get_edg(partial_W)[[3]]
  s_node2 <- get_edg(partial_W)[[4]]
  N_edges <- get_edg(partial_W)[[5]]
  
  M <- get_edg(temp_adj_matrix)[[2]]
  t_node1 <- get_edg(temp_adj_matrix)[[3]]
  t_node2 <- get_edg(temp_adj_matrix)[[4]]
  M_edges <- get_edg(temp_adj_matrix)[[5]]
  
  stan_in_list <- list()
  for(i in 1:length(k_stan_in)){
    k_stan <- k_stan_in[i]
    stan_in_list[[i]] <- list(k=k_stan, ## k is user-specified number of factors ##
                              m = m, # m corresponds to T in our notation (number of time points) ##
                              n = n, # n is number of counties ##
                              nmiss = length(which(is.na(Y0_obs))), #total number of missings
                              n_exp =  length(trt_fips), # number of exposed counties
                              y_miss = y_miss, ## y_miss is n-length indicator of exposure     (missingness in Y(0)) at time m ##
                              y=Y0_full2, ## y is the observed Y(0) matrix ##
                              offs=pop_mat, ## offs is n x m matrix of population size offsets for each county ##
                              N = N,
                              s_node1 = s_node1,
                              s_node2 = s_node2,
                              N_edges = N_edges,
                              M = M,
                              t_node1 = t_node1,
                              t_node2 = t_node2,
                              M_edges = M_edges
    )
  }
  
  return(list(stan_in_list,Y0_full,Y1_obs,Y1_obs_rate,Y0_full_rate))
}

#############################################################################################################

# Get spatially simulated  data
# Changing DGF to be in terms of rates per 100,000
library(splines)
get_sim_spatial_rate_data_smooth <- function(partial_W,temp_adj_matrix, trt_fips,pop_mat,k,alpha,
                                      gamma_sd=.15,psi_sd=.015, n_trt,m_trt, k_stan_in, rho = .6, 
                                      rate = .5, tau2 = .01, simnum = simnum) {
  #####################################################################################
  #### Input:                                                                      ####
  ####  partial_W: spatial adjacency matrix (n x n)                                ####
  ####  temp_adj_matrix: temporal adjacency matrix (m x m)                         ####
  ####  trt_fips: vector of fips codes for treated counties                        ####
  ####  pop_mat : population matrix (n x m)                                        ####
  ####  k: user specified number of latent factors                                 ####
  ####  alpha: user specified global intercept                                     ####
  ####  gamma_sd: county intercept sd                                              ####
  ####  psi_sd: time intercept sd                                                  ####
  ####  n_trt: number of treated counties                                          ####
  ####  m_trt: number of treated time points                                       ####
  ####  k_stan_in: a vector of k_max-values for the STAN model                     ####
  ####  rho: degree of space/time correlation                                      ####
  ####  rate: treatment effect in terms of rate                                    ####
  ####  simnum: seed according to simnum                                           ####
  #### Output:                                                                     ####
  ####  stan_in_list: list of data for stan call                                   ####
  ####  Y0_full: the simulated "true" Y(0) matrix with no missingness (n x m)      ####     
  ####  Y1_obs: the generated oserved Y(1) values after treatment (n_trt x m_trt)  ####
  #####################################################################################
  
  #####################################################################
  ## Generate Mean Structure: Should be the same for each simulation ## 
  #####################################################################
  
  library("mvtnorm")
  
  ## use a common seed to generate the parameters fixed across sims ##
  set.seed(4)
  #set.seed(123)
  #set.seed(10)
  
  #### Define variables  
  n <- nrow(partial_W) # number of counties
  m <- nrow(temp_adj_matrix) # number of time points
  k <- true_k # user specified number of latent factors
  n_trt <- n_trt # number of treated counties
  m_trt <- m_trt # number of treated years
  W <- partial_W # spatial adjacency matrix
  D <- temp_adj_matrix # temporal adjacency matrix
  tau2 <- tau2 #precision parameter
  
  #### Create the spatial precision matrix based on the CAR prior proposed by Leroux et al. (2000)
  Q.W <- rho * (diag(apply(W, 2, sum)) - W) + (1-rho) * diag(rep(1, n))
  #### Create the spatial covariance matrix and simulate a set of spatial random effects
  Q.W.inv <- solve(Q.W)
  
  #### Create the temporal covariance matrix and simulate a set of temporal random effects
  Q.D <- rho * (diag(apply(D, 2, sum)) - D) + (1-rho) * diag(rep(1, m))
  Q.D.inv <- solve(Q.D)
  
  #### Simulate factor loadings 
  L <- t(rmvnorm(n = k,mean = rep(0, m), sigma = tau2 * Q.D.inv )) #(m x k)
  
  #### Simulate latent factors 
  FS <- t(rmvnorm(n = k,mean = rep(0, n), sigma = tau2 * Q.W.inv )) #(n x k)
  
  #### Simulate Intercepts
  ## global intercept ##
  alpha<-alpha #(-4)
  ## county-specific intercepts ##
  gamma_i<- rnorm(n,sd=gamma_sd)
  ## time-specific intercepts ##
  psi_t<- rnorm(m,sd=psi_sd)
  
  #### Get Mean Structure
  mean_struct<-exp(alpha+
                     matrix(rep(gamma_i,m),nrow=n)+
                     matrix(rep(psi_t,n),nrow=n,byrow = T)+
                     FS%*%t(L)+
                     log(pop_mat))
  
  ###################
  ## Generate Data ##
  ###################
  ## use a job-specific seed to generate the random components ##
  set.seed(simnum)
  
  #### Generate the "true" Y(0) matrix with no missingness 
  Y0_smooth<-Y0_full<-Y0_obs<-matrix(rpois(m*n,lambda=c(mean_struct)),nrow=n,ncol=m)
  
  # Get smoothed data using a GLM model
  for(i in 1:length(unique(trt_fips))){
    CL_CASES <- Y0_full[i,] #current county is ith row
    YEAR_DX <- seq(1988,2002)
    POP <- pop_mat[i,] #current county is ith row
    fit <- glm(CL_CASES ~ ns(YEAR_DX,5) + offset(log(POP)), family = "poisson")
    # Save smooth data as CL_CASES. This will keep in sync with processing functions
    Y0_smooth[i,] <- round(predict.glm(fit, type = "response"))
  }
  
  
  #### Insert missingness into the Y(0) matrix to correspond to the observed Y(0) ##
  Y0_obs[1:n_trt,(m-m_trt+1):m]<-NA
  
  #### Convert full counterfactual matrix into rate per 100000
  Y0_full_rate <- (Y0_full/pop_mat) * 100000
  
  #### Generate the Y(1) in terms of rate per 100000
  Y1_obs_rate <- Y0_full_rate
  Y1_obs_rate[1:n_trt,(m-m_trt+1):m] <- Y0_full_rate[1:n_trt,(m-m_trt+1):m] + rate
  
  #### Convert back into counts
  Y1_obs_mat <- (Y1_obs_rate/100000) * pop_mat
  Y1_obs <- Y1_obs_mat[1:n_trt,(m-m_trt+1):m]
  
  ##########################
  ## Set Up Data for STAN ##
  ##########################
  #### Get variables for STAN 
  y_miss <- get_y_miss(Y0_obs)
  Y0_full2 <- Y0_smooth ## STAN doesn't allow missingness
  Y0_full2[is.na(Y0_obs)]<-(99999)
  
  N <- get_edg(partial_W)[[2]]
  s_node1 <- get_edg(partial_W)[[3]]
  s_node2 <- get_edg(partial_W)[[4]]
  N_edges <- get_edg(partial_W)[[5]]
  
  M <- get_edg(temp_adj_matrix)[[2]]
  t_node1 <- get_edg(temp_adj_matrix)[[3]]
  t_node2 <- get_edg(temp_adj_matrix)[[4]]
  M_edges <- get_edg(temp_adj_matrix)[[5]]
  
  stan_in_list <- list()
  for(i in 1:length(k_stan_in)){
    k_stan <- k_stan_in[i]
    stan_in_list[[i]] <- list(k=k_stan, ## k is user-specified number of factors ##
                              m = m, # m corresponds to T in our notation (number of time points) ##
                              n = n, # n is number of counties ##
                              nmiss = length(which(is.na(Y0_obs))), #total number of missings
                              n_exp =  length(trt_fips), # number of exposed counties
                              y_miss = y_miss, ## y_miss is n-length indicator of exposure     (missingness in Y(0)) at time m ##
                              y=Y0_full2, ## y is the observed Y(0) matrix ##
                              offs=pop_mat, ## offs is n x m matrix of population size offsets for each county ##
                              N = N,
                              s_node1 = s_node1,
                              s_node2 = s_node2,
                              N_edges = N_edges,
                              M = M,
                              t_node1 = t_node1,
                              t_node2 = t_node2,
                              M_edges = M_edges
    )
  }
  
  return(list(stan_in_list,Y0_full,Y1_obs,Y1_obs_rate,Y0_full_rate))
}


#############################################################################################################

# Get spatially simulated  data
# Changing DGF to be in terms of rates per 100,000
library(splines)
get_sim_spatial_rate_data_smooth_freq <- function(partial_W,temp_adj_matrix, trt_fips,pop_mat,k,alpha,
                                             gamma_sd=.15,psi_sd=.015, n_trt,m_trt, k_stan_in, rho = .6, 
                                             rate = .5, tau2 = .01, simnum = simnum) {
  #####################################################################################
  #### Input:                                                                      ####
  ####  partial_W: spatial adjacency matrix (n x n)                                ####
  ####  temp_adj_matrix: temporal adjacency matrix (m x m)                         ####
  ####  trt_fips: vector of fips codes for treated counties                        ####
  ####  pop_mat : population matrix (n x m)                                        ####
  ####  k: user specified number of latent factors                                 ####
  ####  alpha: user specified global intercept                                     ####
  ####  gamma_sd: county intercept sd                                              ####
  ####  psi_sd: time intercept sd                                                  ####
  ####  n_trt: number of treated counties                                          ####
  ####  m_trt: number of treated time points                                       ####
  ####  k_stan_in: a vector of k_max-values for the STAN model                     ####
  ####  rho: degree of space/time correlation                                      ####
  ####  rate: treatment effect in terms of rate                                    ####
  ####  simnum: seed according to simnum                                           ####
  #### Output:                                                                     ####
  ####  stan_in_list: list of data for stan call                                   ####
  ####  Y0_full: the simulated "true" Y(0) matrix with no missingness (n x m)      ####     
  ####  Y1_obs: the generated oserved Y(1) values after treatment (n_trt x m_trt)  ####
  #####################################################################################
  
  #####################################################################
  ## Generate Mean Structure: Should be the same for each simulation ## 
  #####################################################################
  
  library("mvtnorm")
  
  ## use a common seed to generate the parameters fixed across sims ##
  set.seed(4)
  #set.seed(123)
  #set.seed(10)
  
  #### Define variables  
  n <- nrow(partial_W) # number of counties
  m <- nrow(temp_adj_matrix) # number of time points
  k <- true_k # user specified number of latent factors
  n_trt <- n_trt # number of treated counties
  m_trt <- m_trt # number of treated years
  W <- partial_W # spatial adjacency matrix
  D <- temp_adj_matrix # temporal adjacency matrix
  tau2 <- tau2 #precision parameter
  
  #### Create the spatial precision matrix based on the CAR prior proposed by Leroux et al. (2000)
  Q.W <- rho * (diag(apply(W, 2, sum)) - W) + (1-rho) * diag(rep(1, n))
  #### Create the spatial covariance matrix and simulate a set of spatial random effects
  Q.W.inv <- solve(Q.W)
  
  #### Create the temporal covariance matrix and simulate a set of temporal random effects
  Q.D <- rho * (diag(apply(D, 2, sum)) - D) + (1-rho) * diag(rep(1, m))
  Q.D.inv <- solve(Q.D)
  
  #### Simulate factor loadings 
  L <- t(rmvnorm(n = k,mean = rep(0, m), sigma = tau2 * Q.D.inv )) #(m x k)
  
  #### Simulate latent factors 
  FS <- t(rmvnorm(n = k,mean = rep(0, n), sigma = tau2 * Q.W.inv )) #(n x k)
  
  #### Simulate Intercepts
  ## global intercept ##
  alpha<-alpha #(-4)
  ## county-specific intercepts ##
  gamma_i<- rnorm(n,sd=gamma_sd)
  ## time-specific intercepts ##
  psi_t<- rnorm(m,sd=psi_sd)
  
  #### Get Mean Structure
  mean_struct<-exp(alpha+
                     matrix(rep(gamma_i,m),nrow=n)+
                     matrix(rep(psi_t,n),nrow=n,byrow = T)+
                     FS%*%t(L)+
                     log(pop_mat))
  
  ###################
  ## Generate Data ##
  ###################
  ## use a job-specific seed to generate the random components ##
  set.seed(simnum)
  
  #### Generate the "true" Y(0) matrix with no missingness 
  Y0_smooth<-Y0_full<-Y0_obs<-matrix(rpois(m*n,lambda=c(mean_struct)),nrow=n,ncol=m)
  
  # Get smoothed data using a GLM model
  for(i in 1:length(unique(trt_fips))){
    CL_CASES <- Y0_full[i,] #current county is ith row
    YEAR_DX <- seq(1988,2002)
    POP <- pop_mat[i,] #current county is ith row
    fit <- glm(CL_CASES ~ ns(YEAR_DX,5) + offset(log(POP)), family = "poisson")
    # Save smooth data as CL_CASES. This will keep in sync with processing functions
    Y0_smooth[i,] <- round(predict.glm(fit, type = "response"))
  }
  
  
  #### Insert missingness into the Y(0) matrix to correspond to the observed Y(0) ##
  Y0_obs[1:n_trt,(m-m_trt+1):m]<-NA
  
  #### Convert full counterfactual matrix into rate per 100000
  Y0_full_rate <- (Y0_full/pop_mat) * 100000
  
  #### Generate the Y(1) in terms of rate per 100000
  Y1_obs_rate <- Y0_full_rate
  Y1_obs_rate[1:n_trt,(m-m_trt+1):m] <- Y0_full_rate[1:n_trt,(m-m_trt+1):m] + rate
  
  #### Convert back into counts
  Y1_obs_mat <- (Y1_obs_rate/100000) * pop_mat
  Y1_obs <- Y1_obs_mat[1:n_trt,(m-m_trt+1):m]
  
  ##########################
  ## Set Up Data for STAN ##
  ##########################
  #### Get variables for STAN 
  y_miss <- get_y_miss(Y0_obs)
  Y0_full2 <- Y0_smooth ## STAN doesn't allow missingness
  Y0_full2[is.na(Y0_obs)]<-Y1_obs_mat[1:n_trt,(m-m_trt+1):m]
 
  
  Y0_obs <- Y0_full2
  
  return(list(Y0_obs = Y0_obs, Y0_full = Y0_full, Y1_obs  = Y1_obs))
}


# Get spatially simulated  data
# Changing DGF to be in terms of rates per 100,000
get_sim_spatial_rate_data_freq <- function(partial_W,temp_adj_matrix, trt_fips,pop_mat,k,alpha,
                                      gamma_sd=.15,psi_sd=.015, n_trt,m_trt, k_stan_in, rho = .6, 
                                      rate = .5, tau2 = .01, simnum = simnum) {
  #####################################################################################
  #### Input:                                                                      ####
  ####  partial_W: spatial adjacency matrix (n x n)                                ####
  ####  temp_adj_matrix: temporal adjacency matrix (m x m)                         ####
  ####  trt_fips: vector of fips codes for treated counties                        ####
  ####  pop_mat : population matrix (n x m)                                        ####
  ####  k: user specified number of latent factors                                 ####
  ####  alpha: user specified global intercept                                     ####
  ####  gamma_sd: county intercept sd                                              ####
  ####  psi_sd: time intercept sd                                                  ####
  ####  n_trt: number of treated counties                                          ####
  ####  m_trt: number of treated time points                                       ####
  ####  k_stan_in: a vector of k_max-values for the STAN model                     ####
  ####  rho: degree of space/time correlation                                      ####
  ####  rate: treatment effect in terms of rate                                    ####
  ####  simnum: seed according to simnum                                           ####
  #### Output:                                                                     ####
  ####  stan_in_list: list of data for stan call                                   ####
  ####  Y0_full: the simulated "true" Y(0) matrix with no missingness (n x m)      ####     
  ####  Y1_obs: the generated oserved Y(1) values after treatment (n_trt x m_trt)  ####
  #####################################################################################
  
  #####################################################################
  ## Generate Mean Structure: Should be the same for each simulation ## 
  #####################################################################
  
  library("mvtnorm")
  
  ## use a common seed to generate the parameters fixed across sims ##
  set.seed(4)
  #set.seed(123)
  #set.seed(10)
  
  #### Define variables  
  n <- nrow(partial_W) # number of counties
  m <- nrow(temp_adj_matrix) # number of time points
  k <- true_k # user specified number of latent factors
  n_trt <- n_trt # number of treated counties
  m_trt <- m_trt # number of treated years
  W <- partial_W # spatial adjacency matrix
  D <- temp_adj_matrix # temporal adjacency matrix
  tau2 <- tau2 #precision parameter
  
  #### Create the spatial precision matrix based on the CAR prior proposed by Leroux et al. (2000)
  Q.W <- rho * (diag(apply(W, 2, sum)) - W) + (1-rho) * diag(rep(1, n))
  #### Create the spatial covariance matrix and simulate a set of spatial random effects
  Q.W.inv <- solve(Q.W)
  
  #### Create the temporal covariance matrix and simulate a set of temporal random effects
  Q.D <- rho * (diag(apply(D, 2, sum)) - D) + (1-rho) * diag(rep(1, m))
  Q.D.inv <- solve(Q.D)
  
  #### Simulate factor loadings 
  L <- t(rmvnorm(n = k,mean = rep(0, m), sigma = tau2 * Q.D.inv )) #(m x k)
  
  #### Simulate latent factors 
  FS <- t(rmvnorm(n = k,mean = rep(0, n), sigma = tau2 * Q.W.inv )) #(n x k)
  
  #### Simulate Intercepts
  ## global intercept ##
  alpha<-alpha #(-4)
  ## county-specific intercepts ##
  gamma_i<- rnorm(n,sd=gamma_sd)
  ## time-specific intercepts ##
  psi_t<- rnorm(m,sd=psi_sd)
  
  #### Get Mean Structure
  mean_struct<-exp(alpha+
                     matrix(rep(gamma_i,m),nrow=n)+
                     matrix(rep(psi_t,n),nrow=n,byrow = T)+
                     FS%*%t(L)+
                     log(pop_mat))
  
  ###################
  ## Generate Data ##
  ###################
  ## use a job-specific seed to generate the random components ##
  set.seed(simnum)
  
  #### Generate the "true" Y(0) matrix with no missingness 
  Y0_full<-Y0_obs<-matrix(rpois(m*n,lambda=c(mean_struct)),nrow=n,ncol=m)
  
  #### Insert missingness into the Y(0) matrix to correspond to the observed Y(0) ##
  Y0_obs[1:n_trt,(m-m_trt+1):m]<-NA
  
  #### Convert full counterfactual matrix into rate per 100000
  Y0_full_rate <- (Y0_full/pop_mat) * 100000
  
  #### Generate the Y(1) in terms of rate per 100000
  Y1_obs_rate <- Y0_full_rate
  Y1_obs_rate[1:n_trt,(m-m_trt+1):m] <- Y0_full_rate[1:n_trt,(m-m_trt+1):m] + rate
  
  #### Convert back into counts
  Y1_obs_mat <- (Y1_obs_rate/100000) * pop_mat
  Y1_obs <- Y1_obs_mat[1:n_trt,(m-m_trt+1):m]
  
  ##########################
  ## Set Up Data for STAN ##
  ##########################
  #### Get variables for STAN 
  y_miss <- get_y_miss(Y0_obs)
  Y0_full2 <- Y0_full ## STAN doesn't allow missingness
  Y0_full2[is.na(Y0_obs)]<-Y1_obs_mat[1:n_trt,(m-m_trt+1):m]

  
  Y0_obs <- Y0_full2
  
  return(list(Y0_obs = Y0_obs, Y0_full = Y0_full, Y1_obs  = Y1_obs))
}
