data {
  int<lower = 1> n;
  int<lower = 1> nS;
  int<lower = 1> nT;
  int<lower = 1> nST;
  int<lower = 1> nrep;
  int<lower = 1> p;
  real X[nS,nT,nrep,p];
  int<lower = 0> y[nS,nT,nrep];
  int<lower = 1> T2;
  int<lower = 1> S1;
  vector[n] log_offset;
  matrix[nST,nST] Qst;
}
parameters {
  vector[p] beta;
  vector[nST] phi;
  real<lower = 0> tau;
}

model {
  phi ~ multi_normal_prec(0, tau * Qst);
  beta ~ normal(0, 1);
  tau ~ gamma(2, 2);
  for(i in 1:nS)
  for(j in 1:nT)
    y[i,j,] ~ poisson_log(X[i,j,,] * beta + phi[i*(nT-1)+j] + log_offset[i,j,]);
}
