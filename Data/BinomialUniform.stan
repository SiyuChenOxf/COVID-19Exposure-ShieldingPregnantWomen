data {
  // Number of data points
  int N;
  // Number of successes
  int n;
}

parameters {
  real<lower=0, upper=1> theta;
}

model {  
  theta ~ beta(1, 1);
  n ~ binomial(N, theta);
}

