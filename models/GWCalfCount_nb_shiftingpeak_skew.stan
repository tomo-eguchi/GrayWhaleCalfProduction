// =====================================================================
// Smooth seasonal pulse -- MEAN-CENTERED, ASYMMETRIC (split-normal)
// =====================================================================
// Same as GWCalfCount_nb_smoothpeak_centered.stan, but the single `width` is
// split into a pre-peak (width_l) and post-peak (width_r) width. This lets the
// migration curve have a longer tail on one side -- the long RIGHT tail typical
// of cow-calf passage would show up as width_r > width_l.
//
// Only fit this AFTER the by-week PPC on the symmetric model shows a systematic
// tail misfit. The generated quantity `width_skew = width_r - width_l` is the
// asymmetry: a 95% interval excluding 0 is direct evidence the skew is real.
// Compare against the symmetric model with LOO (same likelihood, same estimand).
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

  vector[n_years] year_eff;

  real<lower=0> amp;
  real<lower=0.5> width_l;           // pre-peak (left) width
  real<lower=0.5> width_r;           // post-peak (right) width; expect > width_l
  real mu_peak;
  real<lower=0> sigma_peak;
  vector[n_years] peak_raw;

  real<lower=0.01> phi;
}

transformed parameters {
  vector[n_years] peak_week = mu_peak + sigma_peak * peak_raw;
  vector[n_obs]   log_mu_true;
  vector[n_obs]   log_alpha;

  {
    vector[n_years] gbar;
    for (y in 1:n_years) {
      real acc = 0;
      for (k in 1:n_weeks) {
        real dwk = k - peak_week[y];
        real ww  = dwk <= 0 ? width_l : width_r;
        acc += exp(-0.5 * square(dwk / ww));
      }
      gbar[y] = acc / n_weeks;
    }
    for (i in 1:n_obs) {
      real dw   = week_idx[i] - peak_week[year_idx[i]];
      real ww   = dw <= 0 ? width_l : width_r;
      real g    = exp(-0.5 * square(dw / ww));
      real seas = amp * (g - gbar[year_idx[i]]);
      log_mu_true[i] = year_eff[year_idx[i]] + seas;
      log_alpha[i]   = log_mu_true[i] + log(p_obs) + log_offset[i];
    }
  }
}

model {
  year_eff   ~ normal(0, 5);
  amp        ~ normal(0, 5);
  width_l    ~ gamma(2, 0.5);
  width_r    ~ gamma(2, 0.5);
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
  real width_skew = width_r - width_l;   // >0 (CI excluding 0) => longer right tail

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
