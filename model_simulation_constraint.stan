// Stan model for BRSTIR Pareto Uncorrelated Samples
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    array[n] real <lower=1> y; // extreme response
    real <lower=0> atau;
    matrix[2, (2*p)] basisFL;
    array[(p*2)] int indexFL;
}
parameters {
    real theta; // intercept term
    vector[p] beta; // linear coefficient
    vector[(psi-2)] gammaTemp[p]; // constraint splines coefficient from 2 to psi-1
    real <lower=0> lambda; // group lasso penalty
    real sigma;
    array[p] real <lower=0> tau;
    real <lower=0, upper=1> pie;
}
transformed parameters {
    array[n] real <lower=0> alpha; // covariate-adjusted tail index
    vector[(psi+1)] gamma[p]; // splines coefficient 
    vector[2] gammaFL[p]; 
    matrix[2, p] subgnl;
    matrix[n, p] gnl; // nonlinear component
    matrix[n, p] gl; // linear component
    matrix[n, p] gsmooth; // linear component
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgnl; // nonlinear component
    matrix[n, p] newgl; // linear component
    matrix[n, p] newgsmooth; // linear component

    for(j in 1:p){
        gamma[j][3:(psi)] = gammaTemp[j][1:(psi-2)];
        subgnl[,j] = bsNonlinear[indexFL[(((j-1)*2)+1):(((j-1)*2)+2)], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j];
        gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * subgnl[,j] * -1;
        gamma[j][2] = gammaFL[j][1];
        gamma[j][(psi+1)] = gammaFL[j][2];
        gamma[j][1] = beta[j];
    };
    for (j in 1:p){
        gnl[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j][2:(psi+1)];
        newgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j][2:(psi+1)];
        gl[,j] = bsLinear[,j] * gamma[j][1];
        newgl[,j] = xholderLinear[,j] * gamma[j][1];
        gsmooth[,j] = gl[,j] + gnl[,j];
        newgsmooth[,j] = newgl[,j] + newgnl[,j];
    };

    for (i in 1:n){
        alpha[i] = exp(theta + sum(gsmooth[i,])); 
        newalpha[i] = exp(theta + sum(newgsmooth[i,]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda | 1, 1e-6);
    target += normal_lpdf(theta | 0, 100);
    target += inv_gamma_lpdf(sigma | 1, 1e-6);
    target += (-p * (psi+1) * log(lambda)/2);
    target += beta_lpdf(pie | 1, 1);
    for (j in 1:p){
        target += gamma_lpdf(tau[j] | atau, (sqrt(lambda)/2));
        target += log_mix(pie, multi_normal_lpdf(gamma[j] | rep_vector(0, (psi+1)), diag_matrix(rep_vector(1, (psi+1))) * tau[j] * sigma), multi_normal_lpdf(gamma[j] | rep_vector(0, (psi+1)), diag_matrix(rep_vector(1, (psi+1))) * 0.1));
    }
}

