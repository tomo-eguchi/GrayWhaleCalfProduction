// =====================================================================
// Smooth seasonal pulse with YEAR-VARYING peak timing
// =====================================================================
// Replaces the n_weeks free weekly deviations (week_eff) with a parametric
// migration curve: a Gaussian bump in week-space with SHARED amplitude and
// width and a YEAR-SPECIFIC peak week. This:
//   (a) costs far fewer effective parameters, so the seasonal shape stays
//       estimable even in sparse years; and
//   (b) directly tests whether migration timing has shifted -- inspect
//       peak_week[y] across years, and check whether sigma_peak is bounded
//       away from 0.
//
// The seasonal multiplier is exp(seas), seas in [0, amp]:
//   peak week      -> seas = amp           (rate = exp(year_eff + amp))
//   far off-peak   -> seas ~ 0             (rate = exp(year_eff))
// so year_eff is now the OFF-PEAK baseline log-level, and amp is the log
// peak-to-baseline ratio (shared across years).
//
// Year LEVEL keeps the original hierarchical structure (mu_year, sigma_year).
// To combine with unpooled year levels, swap the year block for the one in
// GWCalfCount_nb_fixedyear.stan.
//
// The Gaussian kernel is symmetric. If posterior predictive checks show tail
// misfit (migration curves often have a longer right tail), replace it with an
// asymmetric kernel (separate left/right widths, or a skew-normal shape) -- a
// small extension that does not change the rest of the model.
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
  real week_mid = (n_weeks + 1) / 2.0;     // season midpoint, centres the peak prior
}

parameters {
  real<lower=0, upper=1> p_obs;

  // --- year level (hierarchical, as in the original) ---
  real mu_year;
  real<lower=0> sigma_year;
  vector[n_years] year_eff_raw;

  // --- smooth seasonal pulse ---
  real<lower=0> amp;                 // log peak-to-baseline ratio (shared)
  real<lower=0.5> width;             // pulse width in weeks (shared); away from 0
  real mu_peak;                      // mean peak week across years
  real<lower=0> sigma_peak;          // between-year SD of peak timing (weeks)
  vector[n_years] peak_raw;          // non-centred year peak offsets

  real<lower=0.01> phi;
}

transformed parameters {
  vector[n_years] year_eff;
  vector[n_years] peak_week;
  vector[n_obs]   log_mu_true;
  vector[n_obs]   log_alpha;

  year_eff  = mu_year + sigma_year * year_eff_raw;
  peak_week = mu_peak + sigma_peak * peak_raw;     // year-specific migration peak

  for (i in 1:n_obs) {
    real w    = week_idx[i];
    real seas = amp * exp(-0.5 * square((w - peak_week[year_idx[i]]) / width));
    log_mu_true[i] = year_eff[year_idx[i]] + seas;
    log_alpha[i]   = log_mu_true[i] + log(p_obs) + log_offset[i];
  }
}

model {
  // year level
  year_eff_raw ~ normal(0, 1);
  mu_year      ~ normal(0, 5);
  sigma_year   ~ normal(0, 2);

  // seasonal pulse
  amp        ~ normal(0, 3);              // half-normal; peak/baseline ~ exp(amp)
  width      ~ gamma(2, 0.5);             // mean ~4 weeks, stays away from 0
  mu_peak    ~ normal(week_mid, n_weeks / 4.0);
  sigma_peak ~ normal(0, 2);              // half-normal; allows multi-week timing shifts
  peak_raw   ~ normal(0, 1);

  p_obs ~ normal(0.889, 0.06375);
  phi   ~ gamma(2, 0.1);

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
