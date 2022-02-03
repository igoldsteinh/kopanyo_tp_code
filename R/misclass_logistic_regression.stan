//This model is for logistic regression with misclassified responses on a single covariate
//here z is the observed response, possibly misclassified, y is the true response
//we assume there is no dependence between misclassification and the covariate
//we assume that sensitivity and specificity of misclassification are known
data {
  int<lower=0> N;
  //misclassified response
  int<lower=0, upper=1> z[N];
  //covariate of interese
  real x[N];
  //prob(z=0 | y=0)
  real<lower=0> specificity;
  //prob(z = 1 | y=1)
  real<lower=0> sensitivity;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real beta0;
  real beta1;
}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  real  prob_y1[N];
  real  prob_y0[N];
  real  p[N];
  

  beta0 ~ normal(0,1);
  beta1 ~ normal(0,1);
  
  for (n in 1:N) {
    prob_y1[n] = inv_logit(beta0 + beta1*x[n]);
    prob_y0[n] = 1-prob_y1[n];
    p[n] = sensitivity*prob_y1[n] + (1-specificity)*prob_y0[n];

  }

  for (n in 1:N) {
      z[n] ~ bernoulli(p[n]);
  }
}

