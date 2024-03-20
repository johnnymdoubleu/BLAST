// Stan model for BRSTIR Student-t Uncorrelated Samples

data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; //
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // 
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector[n] y; // extreme response
    real <lower=0> atau;
    array[n] real <lower=0> trueAlpha; // True Alpha
}

parameters {
    vector[(p+1)] theta; // linear predictor
    vector[psi] gamma[p]; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real <lower=0> sigma; //
    array[p] real <lower=0> tau;
}

transformed parameters {
    array[n] real <lower=0> alpha; // tail index
    array[n] real <lower=0> se; // squared error
    real <lower=0> mse; // mean squared error
    matrix[n, p] gnl; // nonlinear component
    matrix[n, p] gl; // linear component
    matrix[n, p] gsmooth; // linear component
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgnl; // nonlinear component
    matrix[n, p] newgl; // linear component
    matrix[n, p] newgsmooth; // linear component
    
    for (j in 1:p){
        gnl[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        gl[,j] = bsLinear[,j] * theta[j+1];
        newgl[,j] = xholderLinear[,j] * theta[j+1];
        gsmooth[,j] = gl[,j] + gnl[,j];
        newgsmooth[,j] = newgl[,j] + newgnl[,j];
    };
    for (i in 1:n){
        alpha[i] = exp(theta[1] + sum(gnl[i,]) + sum(gl[i,]));
        newalpha[i] = exp(theta[1] + sum(newgnl[i,]) + sum(newgl[i,]));
        se[i] = pow((newalpha[i]-trueAlpha[i]), 2);
    };
    mse = mean(se);
}


model {
    // likelihood
    for (i in 1:n){
        target += student_t_lpdf(y[i] | alpha[i], 0, 1);
        target += -1*log(1-student_t_cdf(u, alpha[i], 0, 1));
    }
    target += gamma_lpdf(lambda1 | 1, 10);
    target += gamma_lpdf(lambda2 | 0.1, 10);
    target += normal_lpdf(theta[1] | 0, 100);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01); // target += double_exponential_lpdf(theta[1] | 0, lambda1)
    target += (p * log(lambda1) + (p * psi * log(lambda2)));
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, lambda2/sqrt(2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * sqrt(tau[j]) * sqrt(sigma));
    }
}

