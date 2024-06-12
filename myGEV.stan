// inspired by brms 2.16.1
functions {
  /* generalized extreme value log-PDF for a single response
   * Args: 
   *   y: the response value 
   *   mu: location parameter
   *   sigma: scale parameter
   *   xi: shape parameter
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
   real gen_extreme_value_lpdf(real y, real mu, real sigma, real xi) { 
     real x = (y - mu) / sigma;
     if (xi == 0) {
       return  -log(sigma) - x - exp(-x);
     } else {
       real t = 1 + xi * x;
       real inv_xi = 1 / xi;
       if (t <= 0){
         return log(0); 
         } else{ return -log(sigma) - (1 + inv_xi) * log(t) - pow(t, -inv_xi); 
       }
      }
   }
  /* generalized extreme value log-CDF for a single response
   * Args: 
   *   y: a quantile 
   *   mu: location parameter
   *   sigma: scale parameter
   *   xi: shape parameter
   * Returns:  
   *   log(P(Y <= y))
   */ 
   real gen_extreme_value_lcdf(real y, real mu, real sigma, real xi) { 
     real x = (y - mu) / sigma;
     if (xi == 0) {
       return - exp(-x);
     } else {
       return - pow(1 + xi * x, - 1 / xi);
     }
   }
  /* generalized extreme value log-CCDF for a single response
   * Args: 
   *   y: a quantile
   *   mu: location parameter
   *   sigma: scale parameter
   *   xi: shape parameter
   * Returns:  
   *   log(P(Y > y))
   */ 
   real gen_extreme_value_lccdf(real y, real mu, real sigma, real xi) { 
     return log1m_exp(gen_extreme_value_lcdf(y | mu, sigma, xi));
   }
     /* generalized extreme value log-CCDF for a single response
   * Args: 
   *   y: a quantile
   *   mu: location parameter
   *   sigma: scale parameter
   *   xi: shape parameter
   * Returns:  
   *   q such that P(Y > q) = p 
   */ 
   real gen_extreme_value_qfun(real p, real mu, real sigma, real xi) { 
     if (xi == 0) {
       return mu - sigma * log(-log(p));
     } else {
       return mu + sigma*(pow(-log(p), -xi) - 1)/xi;
     }
   }
   real gen_extreme_value_rng(real mu, real sigma, real xi) { 
       real randomp = uniform_rng(0,1); 
       return gen_extreme_value_qfun(randomp,mu,sigma,xi);
   }
  // /* scale auxiliary parameter xi to a suitable region
  //  * expecting sigma to be a scalar
  //  * Args:
  //  *   xi: unscaled shape parameter
  //  *   y: response values
  //  *   mu: location parameter
  //  *   sigma: scale parameter
  //  * Returns:
  //  *   scaled shape parameter xi
  //  */
  // real scale_xi(real xi, vector y, vector mu, real sigma) {
  //   vector[rows(y)] x = (y - mu) / sigma;
  //   vector[2] bounds = [-inv(min(x)), -inv(max(x))]';
  //   real lb = min(bounds);
  //   real ub = max(bounds);
  //   return inv_logit(xi) * (ub - lb) + lb;
  // }
  // /* scale auxiliary parameter xi to a suitable region
  //  * expecting sigma to be a vector
  //  * Args:
  //  *   xi: unscaled shape parameter
  //  *   y: response values
  //  *   mu: location parameter
  //  *   sigma: scale parameter
  //  * Returns:
  //  *   scaled shape parameter xi
  //  */
  // real scale_xi_vector(real xi, vector y, vector mu, vector sigma) {
  //   vector[rows(y)] x = (y - mu) ./ sigma;
  //   vector[2] bounds = [-inv(min(x)), -inv(max(x))]';
  //   real lb = min(bounds);
  //   real ub = max(bounds);
  //   return inv_logit(xi) * (ub - lb) + lb;
  // }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  real mean_mu;
  real sd_mu; 
  real mean_lsigma;
  real sd_lsigma; 
  real mean_xi;
  real sd_xi; 
}
// transformed data {
// }
parameters {
  real mu;  // location
  real lsigma; // log-scale
  real xi;  // shape
}
transformed parameters {
  real sigma; // scale
  sigma = exp(lsigma);
}
model {
  // likelihood including constants
  //vector[N] mu = mu + rep_vector(0.0, N);
  //vector[N] sigma = sigma + rep_vector(0.0, N);
  //vector[N] xi = xi + rep_vector(0.0, N);
  //for (n in 1:N) {
      // apply the inverse link function
  //    sigma[n] = exp(sigma[n]);
  //}
  // for (n in 1:N) {
  //   // apply the inverse link function
  //   mu[n] = exp(mu[n]);
  // }
  for (n in 1:N) {
      target += gen_extreme_value_lpdf(Y[n] | mu, sigma, xi);
  }

  // priors including constants
  target += student_t_lpdf(mu | 3, mean_mu, sd_mu);
  target += student_t_lpdf(lsigma | 3, mean_lsigma, sd_lsigma);
  target += normal_lpdf(xi | mean_xi, sd_xi);
}
// generated quantities {
//  real theq;
//  theq = gen_extreme_value_qfun(0.01, mu, sigma, xi);
// }

// for prior predictive checks  
// generated quantities {
//   real mu = student_t_rng(3, mean_mu, sd_mu);
//   real xi = normal_rng(mean_xi, sd_xi);
//   real sigma = exp(student_t_rng(3,mean_lsigma, sd_lsigma));
//   vector[N] y_sim;
//   for (i in 1:N) {
//     y_sim[i] = gen_extreme_value_rng(mu, sigma, xi);
//   }
// }
