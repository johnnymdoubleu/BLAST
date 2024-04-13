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
    vector[(p+1)] theta; // linear predictor
    vector[(psi-2)] gammaTemp[p]; // constraint splines coefficient from 2 to psi-1
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    array[p] real <lower=0> tau1;
    array[p] real <lower=0> tau2;
    real <lower=0, upper=1> pie;
}
transformed parameters {
    array[n] real <lower=0> alpha; // covariate-adjusted tail index
    vector[psi] gamma[p]; // splines coefficient 
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
        gamma[j][2:(psi-1)] = gammaTemp[j][1:(psi-2)];
        subgnl[,j] = bsNonlinear[indexFL[(((j-1)*2)+1):(((j-1)*2)+2)], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j];
        gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * subgnl[,j] * -1;
        gamma[j][1] = gammaFL[j][1];
        gamma[j][psi] = gammaFL[j][2];  
    };
    for (j in 1:p){
        gnl[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        gl[,j] = bsLinear[,j] * theta[j+1];
        newgl[,j] = xholderLinear[,j] * theta[j+1];
        gsmooth[,j] = gl[,j] + gnl[,j];
        newgsmooth[,j] = newgl[,j] + newgnl[,j];
    };

    for (i in 1:n){
        alpha[i] = exp(theta[1] + sum(gsmooth[i,])); 
        newalpha[i] = exp(theta[1] + sum(newgsmooth[i,]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += normal_lpdf(theta[1] | 0, 100);
    target += gamma_lpdf(lambda1 | 0.01, 0.01);
    target += gamma_lpdf(lambda2 | 1, 1e-6);
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, sqrt(lambda1));
    }
    for (j in 1:p){
        target += gamma_lpdf(tau1[j] | atau, sqrt(lambda2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * tau1[j]);
        target += -(psi * log(lambda2)/2);
    }
}

