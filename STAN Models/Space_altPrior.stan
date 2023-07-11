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
  int<lower=1,upper=N> s_node1[N_edges]; //// node1[i] adjacent to node2[i]
  int<lower=1, upper=N> s_node2[N_edges];  // and node1[i] < node2[i]
}

parameters {
  real alpha; //Global intercept
  vector[n] d0; //County deviations from global intercept
  vector[m] c0; //time-specific deviations from global intercept
  matrix[n,k] FS; //Factor scores matrix (U in the grant notation)
  matrix[m, k] L; //factor loadings matrix (V in the grant notation)
  //vector[n] FS;
  //vector[m] L;
  real<lower=0> phi; //nb scale parameter
  real <lower=0> tau_ICAR; //precision parameter for ICAR prior
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
  tau_ICAR ~ gamma(1,.01);
  // the following computes the prior on FS on the unit scale with sd = 1
  for (i in 1:k) {
    target += -(0.5 * sqrt(tau_ICAR)) * dot_self((FS[s_node1,i])' - (FS[s_node2,i])');
    // soft sum-to-zero constraint on phi)
    sum(FS[,i]) ~ normal(0, 0.001 * n);  // equivalent to mean(phi) ~ normal(0,0.001)
    
    //print(target());
    //target += target() - target();// didn't work
  }
  //target += -0.5 * columns_dot_self(FS[node1,1:k] - FS[node2,1:k]);
  //real log_p;
  
  for (i in 1:k){
    L[,i] ~ normal(0,1);
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




