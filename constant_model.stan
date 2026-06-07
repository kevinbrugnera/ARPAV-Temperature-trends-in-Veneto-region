data {
int<lower=1> N;
vector[N] T; 
vector[N] time; 

int<lower=1> N_future;
vector[N_future] time_future;
}



parameters {
real T0;
real<lower=0> sigma;
}


model {
T0 ~ normal(0, 20); 
sigma ~ normal(0,10);

T ~ normal(T0, sigma);
}



generated quantities {
array[N_future] real T_rep_future;
  for (i in 1:N_future)
    T_rep_future[i] = normal_rng(T0, sigma);

vector[N] log_lik;
for (i in 1:N)
log_lik[i] = normal_lpdf(T[i] | T0, sigma);
}

