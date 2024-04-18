data {
  int<lower=0> N; //distinguish between total samples N (which will include replicates, perhaps, of qPCR) and Nt (unique time points)
  vector[N] y; //observations
  vector[N] tide_height_frac; // fractional change due to tides
  int N_dolphins[N];
  vector[N] over_height;
} 

parameters {
  real<lower=0, upper=1> alpha;  //otherwise unexplained autoregression
  real gamma;
  real<lower=0> x0; //initial condition
  real<lower=0> beta; // dolphin impact on DNA
  real<lower=0> theta; // scale by which variance increases with mean
  vector<lower=0,upper=1>[N] p; // number of dolphins in high shedding state
  real<lower=0> psi; // magnitude of high shedding
  vector[N] w; // process variance
  real<lower=0> sig_w; // deviation of process variance
  
}

transformed parameters{
  vector<lower=0>[N] x; //mean concentration present. if not on log scale, x cannot be less than 0
  vector<lower=0>[N] v; // observation variance

  //begin time series w initial condition
  x[1] = x0; 
  v[1] = theta*x0;
  
  for (t in 2:N){
    
    x[t] = x[t-1]*(alpha+gamma*tide_height_frac[t-1]) + over_height[t-1]*N_dolphins[t-1]*(beta+p[t-1]*psi) + w[t-1];
    
    v[t] = theta*x[t]; // variance

  }
  

}

model {
  for (i in 1:N){
    
    y[i] ~ normal(x[i], v[i]); //normal dist with variance tied to mean
    
    w[i] ~ normal(0, sig_w); //normal dist with unknown variance
    
  }

  // PRIORS
  theta ~ gamma(1,3); // scale of variance to mean; mean of 0.33 based on matrix samples
  sig_w ~ lognormal(5.5,.5); // same prior as beta, just to start
  
  //DNA autoregression
  x0 ~ normal(100,50);  //gamma(2,1); //based on scale of first observed data points
  alpha ~ beta(1,5); //otherwise unexplained autoregression in eDNA observations. Essentially, 1 - (loss term)
  gamma ~ normal(0,10); // scale of tidal height effect on loss term
  beta ~ lognormal(5.5,.5); //lognormal(2,2); //normal(3,1); //gamma(1,1); //background dolphin impact on DNA concentration
  psi ~ normal(max(y),theta*max(y));
  p ~ beta(1,10);


}

generated quantities{
  vector[N] ysim; //simulated observations, given model
  vector[N] log_likelihood;

  for (t in 1:N){
    ysim[t] = normal_rng((x[t]), v[t]);
    
    log_likelihood[t] = normal_lpdf(y[t] | x[t], v[t]);
  }
  

}
