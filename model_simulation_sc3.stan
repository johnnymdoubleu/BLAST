// Stan model for simple linear regression
functions{
    real halft_lpdf(real y, real c){
        // Half-t distribution log pdf
        return ((c+1)/2) * log(1+((y^2)/c));
    }

    real burr_rng(real c){
        return ((1-uniform_rng(0,1))^(-1)-1)^(1/c);
    }
}

data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp;
    real u; // <lower=0> 
    int <lower=1> psi; // splines coefficient size
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    array[n] real y; // extreme responses
    real <lower=0> atau;
}

parameters {
    real theta; // intercept term
    array[p] vector[psi] gamma; // splines coefficient
    real <lower=0> lambda; // group lasso penalty
    real <lower=0> sigma; //
    array[p] real <lower=0> tau;
}

transformed parameters {
    array[n] real <lower=0> alpha; // tail index
    matrix[n, p] gsmooth; // nonlinear component
    array[n] real <lower=0> newalpha; // tail index
    matrix[n, p] newgsmooth; // nonlinear component
    for (j in 1:p){
        gsmooth[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgsmooth[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    for (i in 1:n){
        alpha[i] = exp(theta + sum(gsmooth[i,]));
        newalpha[i] = exp(theta + sum(newgsmooth[i,]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += student_t_lpdf(y[i] | alpha[i], 0, 1); // student_t_lpdf(y[i] | alpha[i], 0, 1) halft_lpdf(y[i] | alpha[i]) pareto_lpdf(y[i]|u, alpha[i])
        target += -1*log(1-student_t_cdf(u, alpha[i], 0, 1));
    };
    target += gamma_lpdf(lambda | 0.01, 1000);
    target += normal_lpdf(theta | 0, 10);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    target += (p * psi * log(lambda));
    for (j in 1:p){
        target += gamma_lpdf(tau[j] | atau, lambda/sqrt(2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * sqrt(tau[j]) * sqrt(sigma));
    };
}
generated quantities {
    // Used in Posterior predictive check
    vector[n] log_lik;
    array[n] real y_rep = student_t_rng(alpha, rep_vector(0, n),rep_vector(1, n));
    for (i in 1:n) {
        log_lik[i] = student_t_lpdf(y[i] | alpha[i], 0, 1);
    }
}

