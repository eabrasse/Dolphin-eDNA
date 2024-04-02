data {
  int<lower=0> N; //distinguish between total samples N (which will include replicates, perhaps, of qPCR) and Nt (unique time points)
  //int Nt; //unique time points
  vector[N] y; //observations
  //int<lower=1, upper= Nt> time_idx[N];
} 

parameters {
  real<lower=0, upper=1> alpha;  //otherwise unexplained autoregression
  real<lower=0> x0; //initial condition
  real<lower=0> beta; // dolphin impact on DNA
  real<lower=0> theta; // scale by which variance increases with mean
  
}

transformed parameters{
  vector<lower=0>[N] x; //mean concentration present. if not on log scale, x cannot be less than 0
  vector<lower=0>[N] v; // observation variance

  //begin time series w initial condition
  x[1] = x0; 
  v[1] = theta*x0;
  
  for (t in 2:N){

    x[t] = x[t-1]*alpha + beta;    // dolphin DNA sources and sinks  
    
    v[t] = theta*x[t]; // variance

  }
  

}

model {
  for (i in 1:N){
    
    y[i] ~ normal(x[i], v[i]); //normal dist with variance tied to mean
    
  }


  // PRIORS
  theta ~ gamma(1,3); // scale of variance to mean; mean of 0.33 based on matrix samples
  
  //DNA autoregression
  x0 ~ normal(100,50);  //gamma(2,1); //based on scale of first observed data points
  alpha ~ beta(1,5); //otherwise unexplained autoregression in eDNA observations. Essentially, 1 - (loss term)
  beta ~ lognormal(5.5,.5); //lognormal(2,2); //normal(3,1); //gamma(1,1); //background dolphin impact on DNA concentration

}

generated quantities{
  vector[N] ysim; //simulated observations, given model

  for (t in 1:N){
    ysim[t] = normal_rng((x[t]), v[t]);
  }

}
