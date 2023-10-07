// Stan model for simple linear regression

data {
    int <lower = 1> n; // Sample size
    int <lower = 1> p; // regression coefficient size, p
    int <lower = 1> psi; // splines coefficient size, psi
    real <lower = 1> u; // large threshold value, u
    vector[n] zero; // zero vector for multivariate normal
    matrix[n,p] bs.linear; // x.origin
    matrix[n, p*psi] bs.nonlinear // thin plate splines basis
    matrix[n] y; // extreme response
}

parameters {
    vector[n] alpha; // tail index
    vector[p] theta; // linear predictor
    // matrix[psi, p] gamma;
    real lambda1; // lasso penalty
    // real lambda2; // group lasso penalty
    // vector[p] <lower=0> tau; //
}
transformed parameters {
}

model {
    lambda1 ~ dgamma(1, 10);
    lambda2 ~ dgamma(0.1, 0.1);
    for (j in 1:p)
        theta[j] ~ double_exponential(0, lambda1);
    for (i in 1:n)
        alpha[i] <- bs.linear[i] * theta;

    y ~ pareto(u, alpha[i]);
}
//generated quantities {} // The posterior predictive distribution
