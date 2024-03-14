library(npreg)
suppressMessages(library(tidyverse))
library(ReIns)
library(rstan)

# Scenario D
total.iter <- 500

n <- 15000
psi <- 10
threshold <- 0.95
p <- 5
newp <- p+1
no.theta <- 1

# Function to generate Gaussian copula
C <- diag(p)
## Generate sample
gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
    for (ps in 1:psi){
        if(j %in% c(1,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
        else if(j==7){
            if(ps <= (psi/2)){gamma.origin[ps, j] <- -0.1}
            else{gamma.origin[ps, j] <- -0.1}
        }
        else {
            if(ps <= (psi/2)){gamma.origin[ps, j] <- -0.1}
            else{gamma.origin[ps, j] <- -0.1}
        }
    }
}

theta.origin <- c(0.5, 0, -0.2, -0.2, 0, 0)

write("// Stan model for BRSTIR Burr Uncorrelated Samples
functions{
    real burr_lpdf(real y, real c){
        // Burr distribution log pdf
        return log(c)+((c-1)*log(y)) - ((1+1)*log1p(y^c));
    }

    real burr_cdf(real y, real c){
        // Bur distribution cdf
        return 1 - (1 + y^c)^(-1);
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
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector[n] y; // extreme response
    real <lower=0> atau;
}

parameters {
    vector[newp] theta; // linear predictor
    array[p] vector[psi] gamma; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real <lower=0> sigma; //
    array[p] real <lower=0> tau;
}

transformed parameters {
    array[n] real <lower=0> alpha; // tail index
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
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += burr_lpdf(y[i] | alpha[i]); // pareto_lpdf(y[i] | u, alpha[i]) burr_lpdf(y[i] | alpha[i], 1) student_t_lpdf(y[i] | alpha[i], 0, 1)
        target += -1*log(1-burr_cdf(u, alpha[i]));
    }
    target += gamma_lpdf(lambda1 | 0.1, 0.1);
    target += gamma_lpdf(lambda2 | 0.1, 0.1);
    target += normal_lpdf(theta[1] | 0, 100);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01); // target += double_exponential_lpdf(theta[1] | 0, lambda1)
    target += (newp * log(lambda1) + (p * psi * log(lambda2)));
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, lambda2/sqrt(2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * sqrt(tau[j]) * sqrt(sigma));
    }
}
"
, "model_simulation_scD.stan")

newgsmooth.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
alpha.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
alpha.lower.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
alpha.upper.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))

for(iter in 1:total.iter){
    n <- 15000
    x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
    xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
    for(i in 1:p){
        knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
        tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
        bs.linear <- cbind(bs.linear, tps[,1:no.theta])
        bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
    }
   
    g.nonlinear.origin <- g.linear.origin <- g.origin <- matrix(, nrow = n, ncol = p)
    for(j in 1:p){
        g.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
        g.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
        g.origin[, j] <- g.linear.origin[,j] + g.nonlinear.origin[,j]
    }

    alp.origin <- y.origin <- NULL
    for(i in 1:n){
        alp.origin[i] <- exp(theta.origin[1] + sum(g.origin[i,]))
        y.origin[i] <- rburr(1, m=1, s=alp.origin[i], f=1)
    }

    u <- quantile(y.origin, threshold)
    x.origin <- x.origin[which(y.origin>u),]
    y.origin <- y.origin[y.origin > u]
    n <- length(y.origin)

    xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
    newx <- seq(0, 1, length.out=n)
    xholder <- bs.x <- matrix(, nrow = n, ncol = p)
    for(i in 1:p){
        xholder[,i] <- seq(0, 1, length.out = n)  
        test.knot <- seq(0, 1, length.out = psi)
        splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
        xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
        xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
        knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
        tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
        bs.linear <- cbind(bs.linear, tps[,1:no.theta])
        bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)]) 
    }

    g.nonlinear.new <- g.linear.new <- g.new <- g.nonlinear.origin <- g.linear.origin <- g.origin <- matrix(, nrow = n, ncol = p)
    for(j in 1:p){
        g.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
        g.nonlinear.origin[,j] <- bs.nonlinear[1:n, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
        g.origin[, j] <- g.linear.origin[,j] + g.nonlinear.origin[,j]
        g.linear.new[,j] <- xholder.linear[, j] * theta.origin[j+1]
        g.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
        g.new[,j] <- g.linear.new[,j] + g.nonlinear.new[,j]
    }

    true.alpha <- alp.new <- alp.origin <- NULL
    for(i in 1:n){
        alp.origin[i] <- exp(theta.origin[1] + sum(g.origin[i,]))
        alp.new[i] <- exp(theta.origin[1] + sum(g.new[i,]))
    }

    data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                        atau = ((psi+1)/2),
                        bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                        xholderLinear = xholder.linear, 
                        xholderNonlinear = xholder.nonlinear)

    init.alpha <- list(list(gamma = array(rep(0, (psi*p)), dim=c(psi, p)),
                            theta = rep(0, (p+1)), 
                            tau = rep(0.1, p), sigma = 0.1, 
                            lambda1 = 0.1, lambda2 = 0.1),
                      list(gamma = array(rep(0.02, (psi*p)), dim=c(psi, p)),
                            theta = rep(0.01, (p+1)), 
                            tau = rep(0.01, p), sigma = 0.001,
                            lambda1 = 0.01, lambda2 = 0.1),
                      list(gamma = array(rep(0.01, (psi*p)), dim=c(psi, p)),
                            theta = rep(0.05, (p+1)), 
                            tau = rep(0.01, p), sigma = 0.01,
                            lambda1 = 0.1, lambda2 = 0.01))

    fit1 <- stan(
        file = "model_simulation_scD.stan",  # Stan program
        data = data.stan,    # named list of data
        init = init.alpha,      # initial value
        chains = 3,             # number of Markov chains
        iter = 2000,            # total number of iterations per chain
        cores = parallel::detectCores(), # number of cores (could use one per chain)
        refresh = 500             # no progress shown
    )
    newgsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
    newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary

    alpha.lower.container[,iter] <- (newalpha.samples[,4])
    alpha.container[,iter] <- (newalpha.samples[,5])
    alpha.upper.container[,iter] <- (newalpha.samples[,6])
    newgsmooth.container[,iter] <- as.vector(matrix(newgsmooth.samples[,5], nrow = n, byrow=TRUE))
}

alpha.container$x <- seq(0,1, length.out = n)
alpha.container$true <- alp.new
alpha.container <- cbind(alpha.container, t(apply(alpha.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
colnames(alpha.container)[(dim(alpha.container)[2]-2):(dim(alpha.container)[2])] <- c("q1","q2","q3")
alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
alpha.container$q1 <- apply(alpha.lower.container[,1:total.iter], 1, quantile, c(.5))
alpha.container$q3 <- apply(alpha.upper.container[,1:total.iter], 1, quantile, c(.5))
alpha.container <- as.data.frame(alpha.container)

equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

newgsmooth.container$x <- seq(0,1, length.out = n)
newgsmooth.container$true <- as.vector(g.new)
newgsmooth.container <- cbind(newgsmooth.container, t(apply(newgsmooth.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
colnames(newgsmooth.container)[(dim(newgsmooth.container)[2]-2):(dim(newgsmooth.container)[2])] <- c("q1","q2","q3")
newgsmooth.container$mean <- rowMeans(newgsmooth.container[,1:total.iter])
newgsmooth.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]"))
newgsmooth.container <- as.data.frame(newgsmooth.container)
         
save(alpha.container, newgsmooth.container, file = (paste0("./",Sys.Date(),"_",total.iter,"_MC_scD.Rdata")))