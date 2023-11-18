// Stan model for simple linear regression
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp; 
    real <lower=0> u; // large threshold value
    matrix[n, p] x; // train dataset
    matrix[n,p] xholder; // test dataset
    vector[n] y; // extreme response
}

parameters {
    vector[newp] beta; // linear predictor
    real <lower=0> lambda1; // lasso penalty
}

transformed parameters {
    array[n] real <lower=0> alpha; // tail index
    array[n] real <lower=0> newalpha; // new tail index
    
    for (i in 1:n){
        alpha[i] = exp(beta[1] + dot_product(x[i], beta[2:newp]));
        newalpha[i] = exp(beta[1] + dot_product(xholder[i], beta[2:newp]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 0.1, 0.1);
    target += normal_lpdf(beta[1] | 0, 10);
    target += newp * log(lambda1);
    for (j in 1:p){
        target += double_exponential_lpdf(beta[(j+1)] | 0, lambda1);
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

