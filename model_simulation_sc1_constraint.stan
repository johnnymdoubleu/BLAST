// Stan model for BRSTIR Pareto Uncorrelated Samples
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    array[n] real <lower=1> y; // extreme response
    real <lower=0> atau;
    matrix[2, (2*p)] basisFL;
    array[(p*2)] int indexFL;
}
parameters {
    real theta; // linear predictor
    vector[psi] gamma[p]; // splines coefficient 
    real <lower=0> lambda; // group lasso penalty
    real sigma;
    array[p] real <lower=0> tau;          
}
transformed parameters {
    array[n] real <lower=0> alpha; // covariate-adjusted tail index
    matrix[n, p] gnl; // nonlinear component
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgnl; // nonlinear component

    for (j in 1:p){
        gnl[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };

    for (i in 1:n){
        alpha[i] = exp(theta + sum(gnl[i,])); 
        newalpha[i] = exp(theta + sum(newgnl[i,]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += normal_lpdf(theta | 0, 100);
    target += gamma_lpdf(lambda | 0.01, 0.01);
    target += (p * psi * log(lambda)/2);
    for (j in 1:p){
        target += inv_gamma_lpdf(sigma | 0.01, 0.01); 
        target += gamma_lpdf(tau[j] | atau, sqrt(lambda/2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * tau[j] * sigma);
    }
}

