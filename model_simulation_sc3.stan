// Stan model for simple linear regression
functions{
    real halft_lpdf(real y, real c){
        // halft distribution log pdf
        return ((c+1)/2) * log(1+((y^2)/c));
    }

    real burr_rng(real c){
        return ((1-uniform_rng(0,1))^(-1)-1)^(1/c);
    }
}

data {
    int <lower=1> n; // Sample size
    real <lower=0> u;
}

parameters {
    array[n] real <lower=0> newy;
}

model {
    // likelihood
    for (i in 1:n){
        target += student_t_lpdf(newy[i] | 2, 0, 1); // student_t_lpdf(y[i] | alpha[i], 0, 1) halft_lpdf(y[i] | alpha[i]) pareto_lpdf(y[i]|u, alpha[i])
        target += -1*log(1-student_t_cdf(u, 2, 0, 1));
    };
}

