// Stan model for simple linear regression
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp; //parameter and the intercept
    int <lower=1> psi; // splines coefficient size
    int <lower=1> sc; // soft constrained parameter size
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
    array[p] vector[psi] gamma; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real <lower=0> sigma; //
    array[p] real <lower=0> tau;
}

transformed parameters {
    array[n] real <lower=0> alpha; // tail index
    matrix[n, p] gsmooth; // nonlinear component
    array[n] real <lower=0> newalpha; // tail index
    matrix[n, p] newgsmooth; // nonlinear component
    array[2] vector[p] gammasc; // simplex scaled
    for (i in 1:p){
        gammasc[1, i] = sum(bsNonlinear[,(((i-1)*psi)+1)] * gamma[i, 1]);
        gammasc[2, i] = sum(bsNonlinear[,(((i-1)*psi)+psi)] * gamma[i, psi]);
    };

    for (j in 1:p){
        gsmooth[,j] = bsNonlinear[,(((j-1)*psi)+2):(((j-1)*psi)+psi-1)] * gamma[j, 2:(psi-1)] + (bsNonlinear[,(((j-1)*psi)+1)] * gamma[1, j]) + (bsNonlinear[,(((j-1)*psi)+1)] * gamma[2, j]);
        newgsmooth[,j] = xholderNonlinear[,(((j-1)*psi)+2):(((j-1)*psi)+psi-1)] * gamma[j, 2:(psi-1)] + (xholderNonlinear[,(((j-1)*psi)+1)] * gamma[1, j]) + (xholderNonlinear[,(((j-1)*psi)+1)] * gamma[2, j]);
    };
    for (i in 1:n){
        alpha[i] = exp(theta[1] + dot_product(bsLinear[i], theta[2:newp]) + (gsmooth[i,] * rep_vector(1, p)));
        newalpha[i] = exp(theta[1] + dot_product(xholderLinear[i], theta[2:newp]) + (newgsmooth[i,] * rep_vector(1, p)));        
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 0.1, 1);
    target += gamma_lpdf(lambda2 | 0.1, 1);
    target += normal_lpdf(theta[1] | 0, 10);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01); // target += double_exponential_lpdf(theta[1] | 0, lambda1)
    target += (p * log(lambda1) + (p * psi * log(lambda2)));
    for (i in 1:2){
        target += normal_lpdf(sum(gammasc[i]) | 0, 0.001*p);
    }
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, (lambda2/sqrt(2)));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, sc), diag_matrix(rep_vector(1, sc)) * sqrt(tau[j]) * sqrt(sigma)); //if (j < 2 && j > 3) {targpsiet += normal_lpdf(sum(gamma[j]) | 0, 0.001*psi)}
    }
}
generated quantities {
    // Used in Posterior predictive check
    vector[n] log_lik;
    array[n] real y_rep = pareto_rng(rep_vector(u, n), alpha);
    for (i in 1:n) {
        log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}

