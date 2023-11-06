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
    vector[p] theta; // linear predictor
    vector[psi] gamma[p]; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real <lower=0> sigma; //
    vector<lower=0>[p] tau;
}

transformed parameters {
    vector<lower=0>[n] alpha; // tail index
    matrix[n, p] gsmooth; // nonlinear component
    vector<lower=0>[n] newalpha; // tail index
    matrix[n, p] newgsmooth; // nonlinear component
    for (j in 1:p){
        gsmooth[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgsmooth[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    for (i in 1:n){
        alpha[i] = exp(dot_product(bsLinear[i], theta[1:p]) + (gsmooth[i,] * rep_vector(1, p)));
        newalpha[i] = exp(dot_product(xholderLinear[i], theta[1:p]) + (newgsmooth[i,] * rep_vector(1, p)));        
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    };
    target += gamma_lpdf(lambda1 | 1, 10);
    target += gamma_lpdf(lambda2 | 0.1, 1);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    for (j in 1:p){
        target += double_exponential_lpdf(theta[j] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, (square(lambda2)/2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * tau[j] * sigma);
    };
}
generated quantities {
    // Used in Posterior predictive check
    vector[n] log_lik;
    real y_rep[n] = pareto_rng(rep_vector(u, n), alpha);
    for (i in 1:n) {
        log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}

