// Stan model for simple linear regression
functions{
    real burr_lpdf(real y, real c){
        // Burr distribution log pdf
        return log(c)+((c-1)*log(y)) - ((1+1)*log1p(y^c));
    }

    real burr_rng(real c){
        return ((1-uniform_rng(0,1))^(-1)-1)^(1/c);
    }
}

data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp; 
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector[n] y; // extreme response
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
        target += burr_lpdf(y[i] | alpha[i]);
    };
    target += gamma_lpdf(lambda | 0.1, 100);
    target += normal_lpdf(theta | 0, 1);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    target += (p * psi * log(lambda));
    for (j in 1:p){
        target += gamma_lpdf(tau[j] | atau, lambda/sqrt(2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * sqrt(tau[j]) * sqrt(sigma));
    };
}
