data{
  int<lower=1> N1; //number of detections
  int<lower=0, upper=1> z[N1]; //occupancy
  
  int<lower=1> N2; //number of plots with a detection
  int<lower=1> K; //number of cover classes 
  int<lower=1, upper=K> y_int[N2]; //non-zero cover scores
  ordered[K-1] cuts; //upper threshold for cover class
  int<lower=1> J; //number of groups, here species
  int<lower=1, upper=J> grp_id[N1]; // group for spp detection
  int<lower=1, upper=J> grp_id2[N2]; // group for non-zero cover 
}

parameters{
  real logit_mu_raw[J];
  real logit_psi_raw[J];
  
  real mean_mu;
  real<lower=0> sigma_mu;
  real mean_psi;
  real<lower=0> sigma_psi;
  
  real<lower=0, upper=1> y[N1];
  real<lower=0> phi[J];
}

transformed parameters{
  real<lower=0, upper=1> mu[J];
  real<lower=0, upper=1> psi[J];
  
  real<lower=0> alpha[J];
  real<lower=0> beta[J];
  //real logit_p[N1];
  matrix[J, K] theta;
  
  for(j in 1:J){
    mu[j] = inv_logit(mean_mu + logit_mu_raw[j] * sigma_mu);
    psi[j] = inv_logit(mean_psi + logit_psi_raw[j] * sigma_psi);
  }
  
  for(j in 1:J){
    alpha[j] = mu[j] * phi[j];  // for beta cdf - first arg
    beta[j] = (1-mu[j]) * phi[j]; // for beta cdf - 2nd arg
    
    theta[j, 1] = beta_cdf(cuts[1], alpha[j], beta[j]); // theta is the beta cdf for 
    for(k in 2:(K-1)){
      theta[j, k] = beta_cdf(cuts[k], alpha[j], beta[j]) - 
        beta_cdf(cuts[(k-1)], alpha[j], beta[j]);
    }
    theta[j, K] = 1-beta_cdf(cuts[(K-1)], alpha[j], beta[j]);
  }

}

model{
  for(i in 1:N1){
    z[i] ~ bernoulli(psi[grp_id[i]]);
  }
  
  for(i in 1:N2){
    y_int[i] ~ categorical(to_vector(theta[grp_id2[i], ])); //Increment target log probability density with categorical_lpmf(y | theta) dropping constant additive terms.
  }
  
  for(j in 1:J){
    logit_mu_raw[j] ~ normal(0, 1);
    logit_psi_raw[j] ~ normal(0, 1);
    phi[j] ~ cauchy(0, 2.5);
  }
  
  mean_mu ~ cauchy(0, 2.5);
  mean_psi ~ cauchy(0, 2.5);
  
  sigma_mu ~ cauchy(0, 2.5);
  sigma_psi ~ cauchy(0, 2.5);
  
 
}
