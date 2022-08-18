data {
  int<lower = 1> N_trials;
  int<lower = 0,upper = N_trials> ans[4];

}
parameters {
    simplex[4] theta;
}

model {
    target += dirichlet_lpdf(theta | rep_vector(0.1, 4));
    target += multinomial_lpmf(ans |theta  );
}

