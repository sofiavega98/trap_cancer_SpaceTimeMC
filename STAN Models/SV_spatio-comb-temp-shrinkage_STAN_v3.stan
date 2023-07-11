// Same as SV_spatio-temporal_STAN but with multiplicative gamma shrinkage prior combined with AR1
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
  
  // adding parameters for the multiplicative gamma shrinkage prior
  matrix<lower = 0>[m-1,k] p; //hyperprior for factor scores (phi in Bhattacharya paper)
  real<lower = 0> delta1;
  vector<lower = 0>[k-1] delta2; 
  real<lower = 0> a1;
  real<lower = 1> a2;
  // adding parameters for the AR(1) prior
  real alpha_AR;
  real beta_AR; 
  
  real <lower=0> tau_ICAR; //precision parameter for ICAR prior
}

transformed parameters {
  matrix[n,m] Ups; //intermediate predictor
  matrix<lower=0>[n,m] Mu; //poisson mean
  
  //added
  vector<lower = 0>[k] tau;
  vector[k] delta = append_row(delta1, delta2);
  
  for ( i in 1:k ) tau[i] = prod(delta[1:i]);
  
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
  // added combined gamma shrinkage prior with AR(1)
  
  for (i in 1:k) {
    L[1,i] ~ normal(0,10);
    L[2:m,i] ~ normal(alpha_AR + beta_AR' * L[1:(m - 1),i], inv(sqrt((p[,i] * tau[i])))); //inv returns 1/x for each element
  }
  to_vector(p) ~ gamma(1.5, 1.5);
  delta1 ~ gamma(a1, 1);
  delta2 ~ gamma(a2,1);
  tau_ICAR ~ gamma(1,.01);
  //target += -0.5 * dot_self(d0[node1]-d0[node2]); //puts a prior on d0
  
  
  // the following computes the prior on FS on the unit scale with sd = 1
  for (i in 1:k) {
    target += -(0.5 * sqrt(tau_ICAR)) * dot_self((FS[s_node1,i])' - (FS[s_node2,i])');
    // soft sum-to-zero constraint on phi)
    sum(FS[,i]) ~ normal(0, 0.001 * n);  // equivalent to mean(phi) ~ normal(0,0.001)
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





