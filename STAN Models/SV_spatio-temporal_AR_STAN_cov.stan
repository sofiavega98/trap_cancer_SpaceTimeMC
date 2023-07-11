// spatio-temporal model with AR(1) prior on L 
data {
  int<lower = 0> m; // number of time points (T in the grant notation)
  int<lower = 0> k; // number of latent factors
  int<lower = 0> n; // number of counties
  int<lower = 0> nmiss; // total number of missings
  int<lower = 0> n_exp; // number of exposed counties
  int<lower=0, upper=1> y_miss[n,m]; //indicator of missingness in each element of Y
  int<lower=0> y[n,m]; //Y(0) matrix with counties in rows and times in columns
  int<lower=0> offs[n,m]; //offset (population size) nxm matrix with counties in rows and times in columns
  //Spatial part
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1,upper=N> s_node1[N_edges]; //// s_node1[i] adjacent to s_node2[i]
  int<lower=1, upper=N> s_node2[N_edges];  // and s_node1[i] < s_node2[i]
  //Temporal part
  int<lower=0> M; // number of time points after first time point is removed
  int<lower=0> M_edges;
  int<lower=1,upper=N> t_node1[M_edges]; //// t_node1[i] adjacent to t_node2[i]
  int<lower=1, upper=N> t_node2[M_edges];  // and t_node1[i] < t_node2[i]
  //accounting for covariates
  int<lower=0> p; // observed number of covariates
  matrix[n*m,p] x; // observed covariates (in long-form, ordered by row and column indices)
  int<lower = 1, upper = n> row_cov[n*m]; // row indices for covariates
  int<lower = 1, upper = m> col_cov[n*m]; // column indices for covariates

}

parameters {
  real alpha; //Global intercept
  vector[n] d0; //County deviations from global intercept
  vector[m] c0; //time-specific deviations from global intercept
  matrix[n,k] FS; //Factor scores matrix (U in the grant notation)
  matrix[m, k] L; //factor loadings matrix (V in the grant notation)
  real<lower=0> phi; //nb scale parameter
  // adding parameters for the AR(1) prior
  real alpha_AR;
  real beta_AR;
  real<lower=0> sigma_AR; //AR(1) scale parameter-- is this right??
  //accounting for covariates
  vector[p] beta_x; // coefficients for observed covariates
}

transformed parameters {
  matrix[n,m] Ups; //intermediate predictor
  matrix<lower=0>[n,m] Mu; //poisson mean
  
  //accounting for covariates
  vector[m] eta[n]; //means of fixed effect terms
  
  //accounting for covariates
  // build the matrix of fixed-effect means
  for (j in 1:n){
    matrix[m,p] x_j;
    for (i in 1:(n*m)){
      if (row_cov[i]==j) x_j[col_cov[i]]=x[i];
    }
    eta[j]=x_j*beta_x;
  }
  
  // Predictor
  Ups = FS * L';
  // Poisson means
  for(j in 1:m){
    for (i in 1:n){
      Mu[i,j] = exp(alpha + d0[i] + c0[j] + Ups[i,j] + eta[i,j] + log(offs)[i,j]);
    }
  } 
  
}

model {
  //target += -0.5 * dot_self(d0[node1]-d0[node2]); //puts a prior on d0
  // put priors on the other variables to see if there's a change
  alpha ~ normal(0,1);
  d0 ~ normal(0,1);
  c0 ~ normal(0,1);
  
  // the following computes the prior on FS on the unit scale with sd = 1
  for (i in 1:k) {
    target += -0.5 * dot_self((FS[s_node1,i])' - (FS[s_node2,i])');
  }
  //AR prior on L
  for (i in 1:k) {
    L[1,i] ~ normal(0,10);
    L[2:m,i] ~ normal(alpha_AR + beta_AR' * L[1:(m - 1),i],sigma_AR); //inv returns 1/x for each element
  }
  
  for(i in 1:n){
    for (j in 1:m){
      if (1-y_miss[i,j]) y[i,j] ~ neg_binomial_2(Mu[i,j],phi); //Likelihood contribution when y isn't missing
    }
  }
}

generated quantities{
  int<lower=0> Y_pred[nmiss]; //Compute the predictions for treated units at treated times
  real<lower=0> Mu_trt[n_exp*m]; //Extract the expected value for all treated units (all time periods)
  {
    int idy=0;
    int idz=0;
    for (i in 1:n){
      if (sum(y_miss[i,])>0){
        Mu_trt[(idz*m+1):(idz*m+m)]=to_array_1d(Mu[i,]);
        idz=idz+1;
      }
      for (j in 1:m){
        if (y_miss[i,j]){
          Y_pred[idy+1]=neg_binomial_2_rng(Mu[i,j],phi);
          idy=idy+1;
        }
      }
    }
  }
}



