data {
  int<lower = 0> m; // number of time points (T in the grant notation)
  int<lower = 0> k; // number of latent factors
  int<lower = 0> n; // number of counties
  int<lower = 0> nmiss; // total number of missings
  int<lower = 0> n_exp; // number of exposed counties
  int<lower=0, upper=1> y_miss[n,m]; //indicator of missingness in each element of Y
  int<lower=0> y[n,m]; //Y(0) matrix with counties in rows and times in columns
  int<lower=0> offs[n,m]; //offset (population size) nxm matrix with counties in rows and times in columns
}

parameters {
  real alpha; //Global intercept
  vector[n] d0; //County deviations from global intercept
  vector[m] c0; //time-specific deviations from global intercept
  matrix[n,k] FS; //Factor scores matrix (U in the grant notation)
  matrix[m, k] L; //factor loadings matrix (V in the grant notation)
  real<lower=0> log_phi; //nb scale parameter
}

transformed parameters {
  matrix[n,m] Ups; //intermediate predictor
  matrix<lower=0>[n,m] Mu; //poisson mean
  
  // Predictor
  Ups = FS * L';
  // Poisson means
  for(i in 1:n){
    for (j in 1:m){
      Mu[i,j] = exp(alpha + d0[i] + c0[j] + Ups[i,j] + log(offs)[i,j]);
    }
  } 
  
}

model {
  log_phi ~ normal(0, 1);
  
  for (i in 1:k){
    FS[,i] ~ normal(0,1);
  }
  for (i in 1:k){
    L[,i] ~ normal(0,1);
  }
  
  for(i in 1:n){
    for (j in 1:m){
      if (1-y_miss[i,j]) y[i,j] ~ neg_binomial_2(Mu[i,j],exp(log_phi)); //Likelihood contribution when y isn't missing
    }
  }
}

generated quantities{
  int Y_pred[n_exp*m]; //Compute the predictions for treated units at treated times
  real<lower=0> Mu_trt[n_exp*m]; //Extract the expected value for all treated units (all time periods)
  {
    int idy=0;
    int idz=0;
    for (i in 1:n){
      if (sum(y_miss[i,])>0){
        Mu_trt[(idz*m+1):(idz*m+m)]=to_array_1d(Mu[i,]); // save mean values
        idz=idz+1;
        
        for (j in 1:m){
          
          if (log(Mu[i,j])>20){ // set some reasonable threshold here
            Y_pred[idy+1]=-1; // set to some value that you can "filter" out
          } else {
              Y_pred[idy+1]=neg_binomial_2_rng(Mu[i,j],exp(log_phi)); // Calculate and save posterior prediction
          }
          idy=idy+1;
        }
      }
      
    }
  }
}

