// spatio-temporal model with ICAR prior on L 
// Includes a soft-sum constraint for ICAR prior and tau^2~gamma(1,.01)
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
  
}

parameters {
  real alpha; //Global intercept
  vector[n] d0; //County deviations from global intercept
  vector[m] c0; //time-specific deviations from global intercept
  matrix[n,k] FS; //Factor scores matrix (U in the grant notation)
  matrix[m, k] L; //factor loadings matrix (V in the grant notation)
  real<lower=0> phi; //nb scale parameter
  
  vector<lower=0>[k] lambdaFS1; //parameter for laplace prior
  vector<lower=0>[k] lambdaFS2; //parameter for laplace prior
  vector<lower=0>[k] lambdaL1; //parameter for laplace prior
  vector<lower=0>[k] lambdaL2; //parameter for laplace prior
  
  //real lambdaFS1; //parameter for laplace prior
  //real lambdaFS2; //parameter for laplace prior
  //real lambdaL1; //parameter for laplace prior
  //real lambdaL2; //parameter for laplace prior
}

transformed parameters {
  matrix[n,m] Ups; //intermediate predictor
  matrix<lower=0>[n,m] Mu; //poisson mean
  
  // Predictor
  Ups = FS * L';
  // Poisson means
  for(j in 1:m){
    for (i in 1:n){
      Mu[i,j] = exp(alpha + d0[i] + c0[j] + Ups[i,j] + log(offs)[i,j]);
    }
  } 
  
}

model {
  //target += -0.5 * dot_self(d0[node1]-d0[node2]); //puts a prior on d0
  
  // gamma distribution on hyperparameters lambda from Casella et al 2010 paper
  lambdaFS1 ~ gamma(1,0.1);
  lambdaFS2 ~ gamma(1,0.1);
  lambdaL1 ~ gamma(1,0.1);
  lambdaL1 ~ gamma(1,0.1);
  
  
  
  // the following computes the prior on FS 
  //for (i in 1:k) {
  //  target += (sum((-lambdaFS1*(fabs(((FS[s_node1,i])' - (FS[s_node2,i])')))))) + sum((-lambdaFS2*fabs((FS[,i])')));
  //  print("target = ", target())
  //}

  
  
  // prior on L
  //for (i in 1:k) {
  //  target += (sum((-lambdaL1*(fabs(((L[t_node1,i])' - (L[t_node2,i])')))))) + sum((-lambdaL2*fabs((L[,i])')));
  //  print("target = ", target())
  //}
  
  
  
  
  // the following computes the prior on FS 
  for (i in 1:k) {
    target += log((prod(exp(-lambdaFS1[i]*(fabs(((FS[s_node1,i])' - (FS[s_node2,i])')))))) * prod(exp(-lambdaFS2[i]*fabs((FS[,i])'))));
  }

  // prior on L
  for (i in 1:k) {
   target += log((prod(exp(-lambdaL1[i]*(fabs(((L[t_node1,i])' - (L[t_node2,i])')))))) * prod(exp(-lambdaL2[i]*fabs((L[,i])'))));
  }
  
  
  for(i in 1:n){
    for (j in 1:m){
      if (1-y_miss[i,j]) y[i,j] ~ neg_binomial_2(Mu[i,j],phi); //Likelihood contribution when y isn't missing
    }
  }
}

generated quantities{
  int<lower=0> Y_pred[n_exp*m]; //Compute the predictions for treated units at treated times
  real<lower=0> Mu_trt[n_exp*m]; //Extract the expected value for all treated units (all time periods)
  {
    int idy=0;
    int idz=0;
    for (i in 1:n){
      if (sum(y_miss[i,])>0){
        Mu_trt[(idz*m+1):(idz*m+m)]=to_array_1d(Mu[i,]); // save mean values
        idz=idz+1;
        
        for (j in 1:m){
          Y_pred[idy+1]=neg_binomial_2_rng(Mu[i,j],phi); // Calculate and save posterior prediction
          idy=idy+1;
        }
      }
      
    }
  }
}




