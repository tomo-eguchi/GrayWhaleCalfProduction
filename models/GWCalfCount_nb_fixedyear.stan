// =====================================================================
// Decoupled model:  shared seasonal shape  +  UNPOOLED year levels
// =====================================================================
// Identical to the hierarchical NB model EXCEPT year_eff is given a wide,
// independent prior instead of being shrunk toward a common mean. The terms
// mu_year and sigma_year are removed.
//
//   - week_eff is STILL shared across years -> sparse years keep borrowing the
//     SEASONAL SHAPE, which is what makes them estimable in the first place.
//   - each year's LEVEL is now data-driven with NO pull toward the multi-decade
//     average.
//
// This is the production-candidate fix for the shrinkage concern: it keeps the
// estimability you need without the upward pull on recent low years. Compare
// year_eff[2025], year_eff[2026] here against the hierarchical model to read off
// the year-level shrinkage directly.
//
// Note on identifiability: log_mu_true = year_eff + week_eff has a soft additive
// ridge (add c to every year_eff, subtract c from every week_eff). It is broken
// the same way as your original model -- week_eff is zero-centred by its prior,
// so year_eff absorbs the level. If you see the level and shape trading off
// (correlated posteriors, low ESS on year_eff/sigma_week), add a hard
// sum-to-zero constraint on week_eff (sum_to_zero_vector in Stan >= 2.36, or a
// soft sum(week_eff) ~ normal(0, 0.001 * n_weeks) penalty).
// =====================================================================

data {
  int<lower=1> n_obs;
  int<lower=1> n_years;
  int<lower=1> n_weeks;

  array[n_obs] int<lower=0> count_obs;
  vector<lower=0>[n_obs] effort;
  vector[n_obs] log_offset;
  array[n_obs] int<lower=1, upper=n_years> year_idx;
  array[n_obs] int<lower=1, upper=n_weeks> week_idx;
}

parameters {
  real<lower=0, upper=1> p_obs;
  vector[n_years] year_eff;        // FREE, unpooled (was mu_year + sigma_year * raw)
  real<lower=0> sigma_week;
  real<lower=0.01> phi;
  // vector[n_weeks] week_eff_raw;
  sum_to_zero_vector[n_weeks] week_eff_raw;   // replaces vector[n_weeks]
}

transformed parameters {
  vector[n_weeks] week_eff;
  vector[n_obs]   log_mu_true;
  vector[n_obs]   log_alpha;

  week_eff = sigma_week * week_eff_raw;     // shared seasonal shape, zero-centred

  for (i in 1:n_obs) {
    log_mu_true[i] = year_eff[year_idx[i]] + week_eff[week_idx[i]];
    log_alpha[i]   = log_mu_true[i] + log(p_obs) + log_offset[i];
  }
}

model {
  year_eff     ~ normal(0, 5);            // wide, INDEPENDENT: no cross-year shrinkage
  week_eff_raw ~ normal(0, 1);
  p_obs        ~ normal(0.889, 0.06375);
  sigma_week   ~ normal(0, 2);
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
