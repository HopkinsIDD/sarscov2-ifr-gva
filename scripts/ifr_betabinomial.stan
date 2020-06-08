// Beta-Binomial model
// https://mc-stan.org/docs/2_18/stan-users-guide/reparameterizations.html

data {
  int<lower=1> N;         // number of serosurveys
  int<lower=1> M;         // number of posterior draws of seroprevalence
  int<lower=0> deaths[N]; // cumulative number of deaths at serosurvey dates
  int<lower=0> I[N, M];   // number of infected at risk of dying
}
parameters {
  real<lower=0,upper=1> phi;
  real<lower=0.1> lambda;
  real<lower=0, upper=1> IFR;
}
transformed parameters {
  real<lower=0> alpha = lambda * phi;
  real<lower=0> beta = lambda * (1 - phi);
}
model {
  real lp[N, M];
  
  for (n in 1:N) {
    for (m in 1:M) {
      lp[n, m] = binomial_lpmf(deaths[n] | I[n, m], IFR);
    }
  }
  
  for (n in 1:N) {
    target += -log(M) + log_sum_exp(to_vector(lp[n, ]));
  }
  
  phi ~ beta(1, 6.5); // median of prior ~ 0.1
  lambda ~ pareto(0.1, 1.5);
  IFR ~ beta(alpha, beta);
}