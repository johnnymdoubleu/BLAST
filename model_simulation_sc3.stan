// Stan model for simple linear regression
data {
    int <lower=1> n; // Sample size
    real <lower=0> u;
}

parameters {
    array[n] real <lower=u> newy;
}

model {
    // likelihood
    for (i in 1:n){
        target += student_t_lpdf(newy[i] | 2, 0, 1); // student_t_lpdf(y[i] | alpha[i], 0, 1) halft_lpdf(y[i] | alpha[i]) pareto_lpdf(y[i]|u, alpha[i])
        target += -1*log(1-student_t_cdf(u, 2, 0, 1));
    };
}

