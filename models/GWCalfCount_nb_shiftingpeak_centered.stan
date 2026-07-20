// =====================================================================
// Smooth seasonal pulse -- MEAN-CENTERED -- year-varying peak timing
// =====================================================================
// Fixes the amp <-> year_eff ridge that caused the tree-depth blowups and the
// marginal mu_peak convergence in the baseline-anchored version.
//
// Change: the Gaussian pulse is centered to have mean ZERO across the week grid
// (per year), so year_eff becomes the MEAN log-level and no longer trades off
// against amp. Direct analog of the sum_to_zero fix that cleaned up fixed-year.
//
// Year levels are UNPOOLED (flat), matching the fixed-year model, so this
// isolates the effect of freeing the timing. Abundance estimates should be
// invariant to the centering; only the amp/year_eff coordinates change.
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

transformed data {
  real week_mid = (n_weeks + 1) / 2.0;
}

parameters {
  real<lower=0, upper=1> p_obs;

  vector[n_years] year_eff;          // UNPOOLED mean log-level (flat prior)

  real<lower=0> amp;                 // log peak-to-trough amplitude (shared)
  real<lower=0.5> width;             // pulse width in weeks (shared)
  real mu_peak;                      // mean peak week across years
  real<lower=0> sigma_peak;          // between-year SD of peak timing (weeks)
  vector[n_years] peak_raw;          // non-centerd year peak offsets

  real<lower=0.01> phi;
}

transformed parameters {
  vector[n_years] peak_week = mu_peak + sigma_peak * peak_raw;
  vector[n_obs]   log_mu_true;
  vector[n_obs]   log_alpha;

  {
    // per-year mean of the raw pulse over the full week grid (for centering)
    vector[n_years] gbar;
    for (y in 1:n_years) {
      real acc = 0;
      for (k in 1:n_weeks)
        acc += exp(-0.5 * square((k - peak_week[y]) / width));
      gbar[y] = acc / n_weeks;
    }
    for (i in 1:n_obs) {
      real g    = exp(-0.5 * square((week_idx[i] - peak_week[year_idx[i]]) / width));
      real seas = amp * (g - gbar[year_idx[i]]);     // mean-zero across weeks
      log_mu_true[i] = year_eff[year_idx[i]] + seas;
      log_alpha[i]   = log_mu_true[i] + log(p_obs) + log_offset[i];
    }
  }
}

model {
  year_eff   ~ normal(0, 5);
  amp        ~ normal(0, 5);           // widened from N(0,3): data wants amp ~ 8
  width      ~ gamma(2, 0.5);
  mu_peak    ~ normal(week_mid, n_weeks / 4.0);
  sigma_peak ~ normal(0, 2);
  peak_raw   ~ normal(0, 1);
  p_obs      ~ normal(0.889, 0.06375);
  phi        ~ gamma(2, 0.1);

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
