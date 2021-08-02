// An Introduction to Stan
// Michael Betancourt
// March 2020

// 1. Data block

data {
  int N; // specify that N is an integer
  real y[N]; // specify that the obserational space includes integer, defines an array
}

transformed data {
  real sigma = 1; // internal data block, like constant
}

// 2. Parameters block

parameters {
  real theta1; // the model configuration space is theta which defined by the real, defines an array
  real theta2;
}

transformed parameters { // function to define expectation values of interest
  real x;
  { // this scope enable the definition of x but to forget the local value of y (need to defined in the data block)
    real y = 5;
    x = exp(y);
  }
  real mu = baseline(theta1, theta2, 0.5); # new variable defined by the two variables define above
}

// 3. model block defines the target log probability density function

model {
  target += normal_lpdf(theta1 | 0, 1); // each element is added to the target function 
  target += normal_lpdf(theta2 | 0, 1); // the order of each line does not matter, addition is cumulative 
  target += normal_lpdf(y | mu, sigma); 
}
// or more complex

model {
  target += normal_lpdf(theta1 | 0, 1); // defines the independent observational model
  for (n in 1:N) //  add up each contributing log probability density function
    target += normal_lpdf(y[n] | theta, 1); // into a global target accumulator variable
}

generated quantities { // compute posterior predictive samples
  int y_predict = normal_rng(mu, sigma);
}
