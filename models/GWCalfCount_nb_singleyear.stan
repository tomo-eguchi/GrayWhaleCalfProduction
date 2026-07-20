// =====================================================================
// Single-year Negative Binomial model  (NO cross-year pooling)
// =====================================================================
// Purpose: isolate the YEAR-LEVEL shrinkage. Fit one season at a time with
// the SAME NB likelihood as the hierarchical model, so the Poisson-vs-NB
// confound is removed. Compare this fit's estimate for a given year against
// the hierarchical model's estimate for that year:
//     gap  =  year-level pull (shrinkage), cleanly isolated.
//
// How to use: pass the rows for a SINGLE season only. Re-index that season's
// weeks to 1..n_weeks (a contiguous range) before passing them in.
//
// Expectation: data-rich years should converge and roughly agree with the
// hierarchical fit; for 2025/2026 (~26 sightings, ~60% zeros) phi, sigma_week
// and the weekly deviations are weakly identified, so divergences / low ESS /
// R-hat > 1.01 are likely. That non-convergence is itself a reportable result:
// it quantifies how little a single sparse season can say without borrowing.
// =====================================================================

data {
  int<lower=1> n_obs;
  int<lower=1> n_weeks;

  array[n_obs] int<lower=0> count_obs;
  vector<lower=0>[n_obs] effort;
  vector[n_obs] log_offset;
  array[n_obs] int<lower=1, upper=n_weeks> week_idx;
}

parameters {
  real<lower=0, upper=1> p_obs;
  real beta0;                       // single-year log-level (replaces year_eff)
  real<lower=0> sigma_week;         // within-year weekly SD (no sharing across years)
  real<lower=0.01> phi;
  vector[n_weeks] week_eff_raw;
}

transformed parameters {
  vector[n_weeks] week_eff;
  vector[n_obs]   log_mu_true;
  vector[n_obs]   log_alpha;

  week_eff = sigma_week * week_eff_raw;     // zero-centred: beta0 carries the level

  for (i in 1:n_obs) {
    log_mu_true[i] = beta0 + week_eff[week_idx[i]];
    log_alpha[i]   = log_mu_true[i] + log(p_obs) + log_offset[i];
  }
}

model {
  week_eff_raw ~ normal(0, 1);
  p_obs        ~ normal(0.889, 0.06375);
  beta0        ~ normal(0, 5);
  sigma_week   ~ normal(0, 2);             // half-normal via lower=0
  phi          ~ gamma(2, 0.1);

  for (i in 1:n_obs)
    if (effort[i] > 0)
      count_obs[i] ~ neg_binomial_2_log(log_alpha[i], phi);
}

generated quantities {
  vector[n_obs] log_lik;
  vector[n_obs] mu_true;
  vector[n_obs] p_obs_corr;
  array[n_obs] int count_rep;

  for (i in 1:n_obs) {
    mu_true[i] = exp(log_mu_true[i]);
    if (effort[i] > 0) {
      log_lik[i]    = neg_binomial_2_log_lpmf(count_obs[i] | log_alpha[i], phi);
      p_obs_corr[i] = p_obs * exp(log_offset[i]);
      count_rep[i]  = neg_binomial_2_log_rng(log_alpha[i], phi);
    } else {
      log_lik[i]    = 0;
      p_obs_corr[i] = 0;
      count_rep[i]  = 0;
    }
  }
}
