// An Introduction to Stan
// Michael Betancourt
// March 2020

// 1. Data block

data {
  int N; // specify that N is an integer
  real y[N]; // specify that the obserational space includes integer
}

// 2. Parameters block

parameters {
  real theta; // the model configuration space is theta which defined by the real 
}

// 3. model block defines the target log probability density function

model {
  target += normal_lpdf(theta | 0, 1); // defines the independent observational model
  for (n in 1:N) //  add up each contributing log probability density function
    target += normal_lpdf(y[n] | theta, 1); // into a global target accumulator variable
}
