data {
int<lower=1> N;
vector[N] T; # use standardized 
vector[N] time; # use series: from -1 to 1 at equal steps

int<lower=1> N_future;
vector[N_future] time_future;
}



parameters {
real T0;
real T_dot;
real<lower=0> sigma;

}


model {
T0 ~ normal(0, 20); # check if correct variance??
T_dot ~ normal(0, 20);
sigma ~ normal(0,10);

T ~ normal(T0 + T_dot * time, sigma);
}



generated quantities {
real T_rep_future[N_future];
  for (i in 1:N_future)
    T_rep_future[i] = normal_rng(T0 + T_dot * time_future[i], sigma);

vector[N] log_lik;
for (i in 1:N)
log_lik[i] = normal_lpdf(T[i] | T0 + T_dot * time[i], sigma);
}

