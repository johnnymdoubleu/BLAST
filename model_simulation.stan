// Stan model for simple linear regression
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp; 
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector[n] y; // extreme response
    real <lower=0> atau;
}

parameters {
    vector[newp] theta; // linear predictor
    vector[psi] gamma[p]; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real sigma; //
    vector[p] tau;
}

transformed parameters {
    vector[n] alpha; // tail index
    vector[n] newalpha; // tail index
    matrix[n, p] gsmooth; // nonlinear component
    matrix[n, p] newgsmooth; // nonlinear component
    for (j in 1:p){
        gsmooth[,j] <- bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgsmooth[,j] <- xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    for (i in 1:n){
        alpha[i] <- exp(theta[1] + dot_product(bsLinear[i], theta[2:newp]) + (gsmooth[i,] *rep_vector(1, p)));
        newalpha[i] <- exp(theta[1] + dot_product(xholderLinear[i], theta[2:newp]) + (newgsmooth[i,] *rep_vector(1, p)));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 0.0001, 1);
    target += gamma_lpdf(lambda2 | 0.000001, 1);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    target += double_exponential_lpdf(theta[1] | 0, lambda1); // target += normal_lpdf(theta[1] | 0, 0.01);
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, (square(lambda2)/2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), (diag_matrix(rep_vector(1, psi)) * tau[j] * sigma));
    }
}
generated quantities {
    // Used in Posterior predictive check
    vector[n] log_lik;
    real y_rep[n] = pareto_rng(u, alpha);
    for (i in 1:n) {
        log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}

