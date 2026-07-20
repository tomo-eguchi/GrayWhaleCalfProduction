data {
  int<lower=1> n_obs;
  int<lower=1> n_years;
  int<lower=1> n_weeks;
  
  array[n_obs] int<lower=0> count_obs;
  vector<lower=0>[n_obs] effort; // Added raw effort vector
  vector[n_obs] log_offset; 
  array[n_obs] int<lower=1, upper=n_years> year_idx;
  array[n_obs] int<lower=1, upper=n_weeks> week_idx;
}

parameters {
  real<lower=0, upper=1> p_obs; 
  real mu_year;
  
  real<lower=0> sigma_year;
  real<lower=0> sigma_week;
  real<lower=0.01> phi; // Bounded away from zero for stability
  
  vector[n_years] year_eff_raw;
  vector[n_weeks] week_eff_raw;
}

transformed parameters {
  vector[n_years] year_eff;
  vector[n_weeks] week_eff;
  vector[n_obs] log_alpha;

  year_eff = mu_year + sigma_year * year_eff_raw;
  week_eff = sigma_week * week_eff_raw; 
  
  for (i in 1:n_obs) {
    log_alpha[i] = year_eff[year_idx[i]] + week_eff[week_idx[i]] + log(p_obs) + log_offset[i];
  }
}

model {
  year_eff_raw ~ normal(0, 1);  // this is eta_{year}
  week_eff_raw ~ normal(0, 1);  // this is eta_{week}
  p_obs ~ normal(0.889, 0.06375);
  mu_year ~ normal(0, 5);
  sigma_year ~ normal(0, 2); 
  sigma_week ~ normal(0, 2);
  phi ~ gamma(2, 0.1); 

  // Likelihood loop with zero-effort safety switch
  for (i in 1:n_obs) {
    if (effort[i] > 0) {
      count_obs[i] ~ neg_binomial_2_log(log_alpha[i], phi);
    }
    // If effort == 0, it skips completely! Perfect gradient safety.
  }
}

generated quantities {
  vector[n_obs] log_lik;
  vector[n_obs] mu_true;
  vector[n_obs] p_obs_corr;
  array[n_obs] int count_rep;
  
  for (i in 1:n_obs) {
    mu_true[i] = exp(year_eff[year_idx[i]] + week_eff[week_idx[i]]);
    
    if (effort[i] > 0) {
      log_lik[i] = neg_binomial_2_log_lpmf(count_obs[i] | log_alpha[i], phi);
      p_obs_corr[i] = p_obs * exp(log_offset[i]); 
	  
	  // Generate a simulated observation from the posterior distribution
      count_rep[i] = neg_binomial_2_log_rng(log_alpha[i], phi);
    } else {
      log_lik[i] = 0; // Zero-effort days have no predictive likelihood
      p_obs_corr[i] = 0; // 0% chance of observing a calf because nobody looked
	  
	  // Replicated observed count is structurally zero because nobody looked
      count_rep[i] = 0;
    }
  }
}
