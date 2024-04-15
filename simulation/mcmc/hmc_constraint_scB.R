library(npreg)
library(tmvnsim)
library(reshape2)
library(splines2)
library(scales)
library(MASS)
library(npreg)
library(Pareto)
suppressMessages(library(tidyverse))
library(JOPS)
library(readxl)
library(gridExtra)
library(colorspace)
library(corrplot)
library(ReIns)
library(evir)
library(rstan)
library(ggmcmc)
library(MCMCvis)
#library(cmdstanr)
library(ggh4x)

#Scenario 1
set.seed(10)
# set.seed(6)

n <- 15000
psi <- 10
threshold <- 0.95
p <- 5
no.theta <- 1
simul.no <- 50

# Function to generate Gaussian copula
C <- matrix(c(1, 0.3, 0.5, 0.3, 0.3,
            0.3, 1, 0.95, 0.4, 0.4,
            0.5, 0.95, 1, 0.5, 0.1,
            0.3, 0.4, 0.5 , 1, 0.5,
            0.3, 0.4, 0.5, 0.5, 1), nrow = p)
xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))

# x.origin <- apply(x.origin, 2, sort, decreasing=F)
end.holder <- basis.holder <- matrix(, nrow = 2, ncol =0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                        matrix(c(which.min(x.origin[,i]),
                                 which.max(x.origin[,i])), ncol=2))
}
for(i in 1:p){
  knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
  tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  bs.linear <- cbind(bs.linear, tps[,1:no.theta])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
  basis.holder <- cbind(basis.holder, 
                        solve(t(matrix(c(tps[index.holder[i,1], no.theta+1],
                                         tps[index.holder[i,1], no.theta+psi],
                                         tps[index.holder[i,2], no.theta+1],
                                         tps[index.holder[i,2], no.theta+psi]), 
                                       nrow = 2, ncol = 2))))
  end.holder <- cbind(end.holder, 
                      matrix(c(tps[index.holder[i,1], no.theta+1],
                               tps[index.holder[i,1], no.theta+psi],
                               tps[index.holder[i,2], no.theta+1],
                               tps[index.holder[i,2], no.theta+psi]), 
                             nrow = 2, ncol = 2))
}


## Generate sample
gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
  for (ps in 1:psi){
    if(j %in% c(1,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
    else {
      if(ps == 1 || ps == psi){gamma.origin[ps, j] <- 0}
      else{gamma.origin[ps, j] <- -200}
    }
  }
}
theta.origin <- c(-0.5, 0, -0.2, -0.2, 0, 0)

f.sub.origin <- matrix(, nrow = 2, ncol = p)
for(j in 1:p){
  f.sub.origin[,j] <- as.matrix(bs.nonlinear[index.holder[j,], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))], nrow = 2) %*% gamma.origin[(2:(psi-1)), j]
  gamma.origin[c(1,psi),j] <- -1 * basis.holder[,(((j-1)*2)+1):(((j-1)*2)+2)] %*% as.matrix(f.sub.origin[,j], nrow=2)
}

f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
  f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
  f.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
  f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
  alp.origin[i] <- exp(theta.origin[1] + sum(f.origin[i,]))
  y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i]) 
}

u <- quantile(y.origin, threshold)
excess.index <- which(y.origin>u)
x.origin <- as.matrix(x.origin[excess.index,])
# bs.nonlinear <- bs.nonlinear[excess.index,]
# bs.linear <- bs.linear[excess.index,]

y.origin <- y.origin[y.origin > u]
n <- length(y.origin)

xholder.nonlinear <- xholder.linear <-  matrix(,nrow=n, ncol=0)
bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
newx <- seq(0, 1, length.out=n)
xholder <- bs.x <- matrix(, nrow = n, ncol = p)
end.holder <- basis.holder <- matrix(, nrow = 2, ncol =0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                        matrix(c(which.min(x.origin[,i]),
                                 which.max(x.origin[,i])), ncol=2))
}
# # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# # x.origin <- as.data.frame(sapply(as.data.frame(x.origin), FUN = range01))
for(i in 1:p){
  # xholder[,i] <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = n)  
  # test.knot <- seq(min(xholder[,i]), max(xholder[,i]), length.out = psi)
  xholder[,i] <- seq(0, 1, length.out = n)  
  test.knot <- seq(0, 1, length.out = psi)
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
  tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  basis.holder <- cbind(basis.holder, 
                        solve(t(matrix(c(tps[index.holder[i,1], no.theta+1],
                                         tps[index.holder[i,1], no.theta+psi],
                                         tps[index.holder[i,2], no.theta+1],
                                         tps[index.holder[i,2], no.theta+psi]), 
                                       nrow = 2, ncol = 2))))
  # end.holder <- cbind(end.holder, 
  #             matrix(c(tps[index.holder[i,1], no.theta+1],
  #                 tps[index.holder[i,1], no.theta+psi],
  #                 tps[index.holder[i,2], no.theta+1],
  #                 tps[index.holder[i,2], no.theta+psi]), 
  #                nrow = 2, ncol = 2)) 
  bs.linear <- cbind(bs.linear, tps[,1:no.theta])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)]) 
}


f.nonlinear.new <- f.linear.new <- f.new <- f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)

for(j in 1:p){
  f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
  f.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
  f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
  f.linear.new[,j] <- xholder.linear[, j] * theta.origin[j+1]
  f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
  f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}

true.alpha <- alp.new <- alp.origin <- NULL
for(i in 1:n){
  alp.origin[i] <- exp(theta.origin[1] + sum(f.origin[i,]))
  alp.new[i] <- exp(theta.origin[1] + sum(f.new[i,]))
}

write("// Stan model for BRSTIR Pareto Uncorrelated Samples
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
    //vector[(psi)] gammaTemp[p];
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    array[p] real <lower=0> tau;
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
        gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * subgnl[,j] * (-1);
        gamma[j][1] = gammaFL[j][1];
        gamma[j][psi] = gammaFL[j][2];  
    };
    //gamma=gammaTemp;
    
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
    target += gamma_lpdf(lambda1 | 1, 1e-3);
    target += gamma_lpdf(lambda2 | 1, 1e-3);
    target += (2*p*log(lambda2));
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, lambda2^2*0.5);
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * (1/tau[j]));
    }
}
"
, "model_simulation_sc1_constraint.stan")
# for ( i in 1:psi){
# target += inv_gamma_lpdf(sigma | 0.1, 0.1); 
#     target += double_exponential_lpdf(gamma[j][i] | 0, sqrt(lambda2));
#  target += double_exponential_lpdf(theta[(j+1)] | 0, sqrt(lambda1));
# }
# target += beta_lpdf(pie | 3, 1);
# target += ((p * log(lambda1)/2) + (p * psi * log(lambda2)/2));
# target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * (1/tau1[j]));
# target += normal_lpdf(tau2[j] | 0, 0.01);
# target += log_mix(pie, multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * (1/tau1[j])), multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * (1/tau2[j])));

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                  atau = ((psi+1)/2), basisFL = basis.holder,
                  indexFL = as.vector(t(index.holder)),
                  bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                  xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear)

init.alpha <- list(list(gammaTemp = array(rep(2, ((psi-2)*p)), dim=c((psi-2),p)),
                        theta = rep(0, (p+1)), tau = rep(0.1, p),
                        # pie = 0.05,tau2 = rep(0.1, p),
                        lambda1 = 0.1, lambda2 = 1),
                   list(gammaTemp = array(rep(-1, ((psi-2)*p)), dim=c((psi-2),p)),
                        theta = rep(0, (p+1)), tau = rep(0.001, p),
                        #tau2 = rep(0.5, p),pie = 0.75,
                        lambda1 = 100, lambda2 = 100),
                   list(gammaTemp = array(rep(-3, ((psi-2)*p)), dim=c((psi-2),p)),
                        theta = rep(0.1, (p+1)), tau = rep(0.5, p),
                        # pie = 0.5,  tau2 = rep(0.01, p),
                        lambda1 = 5, lambda2 = 55))

# init.alpha <- list(list(gammaTemp = array(rep(2, ((psi)*p)), dim=c(p,(psi))),
#                         theta = rep(0, (p+1)), tau1 = rep(0.1, p), 
#                         # pie = 0.05,tau2 = rep(0.1, p),
#                         lambda1 = 0.1, lambda2 = 1),
#                    list(gammaTemp = array(rep(-1, ((psi)*p)), dim=c(p,(psi))),
#                         theta = rep(0, (p+1)), tau1 = rep(0.001, p), 
#                         #tau2 = rep(0.5, p),pie = 0.75,
#                         lambda1 = 100, lambda2 = 100),
#                    list(gammaTemp = array(rep(-3, ((psi)*p)), dim=c(p,(psi))),
#                         theta = rep(0.1, (p+1)), tau1 = rep(0.5, p),
#                         # pie = 0.5,  tau2 = rep(0.01, p),
#                         lambda1 = 5, lambda2 = 55))
                   
# setwd("C:/Users/Johnny Lee/Documents/GitHub")
fit1 <- stan(
  file = "model_simulation_sc1_constraint.stan",  # Stan program
  # file = "model_BRSTIR.stan",  # Stan program
  data = data.stan,    # named list of data
  init = init.alpha,      # initial value
  chains = 3,             # number of Markov chains
  # warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = parallel::detectCores(), # number of cores (could use one per chain)
  refresh = 500             # no progress shown
  # control=list(adapt_delta=0.99, max_treedepth=15)
)

posterior <- extract(fit1)

plot(fit1, plotfun = "trace", pars = c("theta"), nrow = 3)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_theta_trace_sc1-wi.pdf"), width=10, height = 7.78)
# plot(fit1, plotfun = "trace", pars = c("lambda1", "lambda2"), nrow = 2)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_lambda_sc1-wi.pdf"), width=10, height = 7.78)

tau.samples <- summary(fit1, par=c("tau"), probs = c(0.05,0.5, 0.95))$summary
# pie.samples <- summary(fit1, par=c("pie"), probs = c(0.05,0.5, 0.95))$summary
theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
newgl.samples <- summary(fit1, par=c("newgl"), probs = c(0.05, 0.5, 0.95))$summary
newgnl.samples <- summary(fit1, par=c("newgnl"), probs = c(0.05, 0.5, 0.95))$summary
newgsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary
subgnl.samples <- summary(fit1, par=c("subgnl"), probs = c(0.05,0.5, 0.95))$summary
gammafl.samples <- summary(fit1, par=c("gammaFL"), probs = c(0.05,0.5, 0.95))$summary

gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q2 <- gamma.samples[,1]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q2 <- theta.samples[,5]
theta.q3 <- theta.samples[,6]

# array(gamma.q2, dim=c(psi,p))
# sampled <- end.holder[,1:2] %*% gammafl.samples[1:2, 5]
# trued <- as.matrix(c(bs.nonlinear[index.holder[1,1], 2:(psi-1)] %*% gamma.samples[2:(psi-1), 5], bs.nonlinear[index.holder[1,2], 2:(psi-1)] %*% gamma.samples[2:(psi-1), 5]), nrow = 2)
# sampled - trued
# sampled
# trued
df.theta <- data.frame("seq" = seq(1, (p+1)),
                       "true" = theta.origin,
                       "m" = theta.q2,
                       "l" = theta.q1,
                       "u" = theta.q3)
df.theta$covariate <- factor(0:p)
df.theta$labels <- factor(0:p)
ggplot(df.theta, aes(x = covariate, y=m, color = covariate)) + ylab("") + xlab('') +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_point(size = 5) + geom_point(aes(y = true), color="red", size = 4) +
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, linewidth =1.2) + 
  scale_x_discrete(labels = c(expression(bold(theta[0])),
                              expression(bold(theta[1])),
                              expression(bold(theta[2])),
                              expression(bold(theta[3])),
                              expression(bold(theta[4])),
                              expression(bold(theta[5])))) + 
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.text = element_text(size=25),
        legend.margin=margin(0,0,0,-10),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        axis.text.x = element_text(hjust=0.35),
        axis.text = element_text(size = 28))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_theta_sc1-wi.pdf"), width=10, height = 7.78)

df.gamma <- data.frame("seq" = seq(1, (psi*p)), 
                       "true" = as.vector(gamma.origin),
                       "m" = as.vector(gamma.q2),
                       "l" = as.vector(gamma.q1),
                       "u" = as.vector(gamma.q3))
df.gamma$covariate <- factor(rep(seq(1, 1 + nrow(df.gamma) %/% psi), each = psi, length.out = nrow(df.gamma)))
df.gamma$labels <- factor(1:(psi*p))
ggplot(df.gamma, aes(x =labels, y = m, color = covariate)) + 
  geom_errorbar(aes(ymin = l, ymax = u),alpha = 0.4, width = 4, linewidth = 1.2) +
  geom_point(aes(y=true), size =4, color ="red")+
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_point(size = 4) + ylab("") + xlab("" ) + #ylim(-15,15) +
  # geom_ribbon(aes(ymin = l, ymax = u)) +
  scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+7), 
                   label = c(expression(bold(gamma[1])), 
                             expression(bold(gamma[2])), 
                             expression(bold(gamma[3])), 
                             expression(bold(gamma[4])), 
                             expression(bold(gamma[5]))),
                   expand=c(0,3)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=25),
        legend.margin=margin(0,0,0,-10),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        axis.text.x = element_text(hjust=0.5),
        axis.text = element_text(size = 28))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_gamma_sc1-wi.pdf"), width=10, height = 7.78)

g.linear.mean <- as.vector(matrix(newgl.samples[,1], nrow = n, byrow=TRUE))
g.linear.q1 <- as.vector(matrix(newgl.samples[,4], nrow = n, byrow=TRUE))
g.linear.q2 <- as.vector(matrix(newgl.samples[,5], nrow = n, byrow=TRUE))
g.linear.q3 <- as.vector(matrix(newgl.samples[,6], nrow = n, byrow=TRUE))
g.nonlinear.mean <- as.vector(matrix(newgnl.samples[,1], nrow = n, byrow=TRUE))
g.nonlinear.q1 <- as.vector(matrix(newgnl.samples[,4], nrow = n, byrow=TRUE))
g.nonlinear.q2 <- as.vector(matrix(newgnl.samples[,5], nrow = n, byrow=TRUE))
g.nonlinear.q3 <- as.vector(matrix(newgnl.samples[,6], nrow = n, byrow=TRUE))
g.smooth.mean <- as.vector(matrix(newgsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(newgsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(newgsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(newgsmooth.samples[,6], nrow = n, byrow=TRUE))

equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

data.smooth <- data.frame("x"=newx,
                          "true" = as.vector(f.new),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]")),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=true, colour = "True"), linewidth=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab(expression(c)) +
  facet_wrap(covariates ~ ., scales = "free_x", nrow = 5,
             labeller = label_parsed, strip.position = "left") + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
         fill = guide_legend(order = 1)) + #ylim(-0.55, 0.5) +
  #   scale_y_continuous(breaks=c(-0.6,-0.3,-0,0.3)) + 
  # scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.margin=margin(t = 1, unit='cm'),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text.y = element_text(size = 25, colour = "black", angle = 0, face = "bold.italic"),
        strip.placement = "outside",
        axis.title.x = element_text(size = 35),
        axis.text = element_text(size=18))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_smooth_sc1-wi.pdf"), width=12.5, height = 15)

data.linear <- data.frame("x"=newx,
                          "true" = as.vector(f.linear.new),
                          "post.mean" = as.vector(g.linear.mean),
                          "q1" = as.vector(g.linear.q1),
                          "q2" = as.vector(g.linear.q2),
                          "q3" = as.vector(g.linear.q3),
                          "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]")),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

ggplot(data.linear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=true, colour = "True"), linewidth=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab(expression(c)) +
  facet_wrap(covariates ~ ., scales = "free_x", nrow = 5,
             labeller = label_parsed, strip.position = "left") + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
         fill = guide_legend(order = 1)) + ylim(-0.55, 0.5) +
  # scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.margin=margin(t = 1, unit='cm'),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text.y = element_text(size = 25, colour = "black", angle = 0, face = "bold.italic"),
        strip.placement = "outside",
        axis.title.x = element_text(size = 35),
        axis.text = element_text(size=18))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_linear_sc1-wi.pdf"), width=12.5, height = 15)


data.nonlinear <- data.frame("x"=newx,
                             "true" = as.vector(f.nonlinear.new),
                             "post.mean" = as.vector(g.nonlinear.mean),
                             "q1" = as.vector(g.nonlinear.q1),
                             "q2" = as.vector(g.nonlinear.q2),
                             "q3" = as.vector(g.nonlinear.q3),
                             "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]")),
                             "fakelab" = rep(1, (p*n)),
                             "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=true, colour = "True"), linewidth=2) + 
  geom_line(aes(y=post.mean, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab(expression(c)) +
  facet_wrap(covariates ~ ., scales = "free_x", nrow = 5,
             labeller = label_parsed, strip.position = "left") + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
         fill = guide_legend(order = 1)) + #ylim(-1, 1) +
  # scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.margin=margin(t = 1, unit='cm'),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text.y = element_text(size = 25, colour = "black", angle = 0, face = "bold.italic"),
        strip.placement = "outside",
        axis.title.x = element_text(size = 35),
        axis.text = element_text(size=18))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_nonlinear_sc1-wi.pdf"), width=12.5, height = 15)

data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "true" = (alp.new),
                            "post.mean" = (newalpha.samples[,1]),
                            "post.median" = (newalpha.samples[,5]),
                            "q1" = (newalpha.samples[,4]),
                            "q3" = (newalpha.samples[,6]))
# "post.mean" = sort(alpha.smooth.new),
# "post.median" = sort(newalpha.samples[,5]),
# "q1" = sort(alpha.smooth.q1),
# "q3" = sort(alpha.smooth.q3))

ggplot(data.scenario, aes(x=newx)) + 
  ylab(expression(alpha(c*bold("1")))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2) + 
  geom_line(aes(y=post.mean, col = "Posterior Median"), linewidth=1.5) +
  scale_color_manual(values=c("steelblue", "red")) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  theme_minimal(base_size = 30) + #ylim(0.5,2.5)+
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_test_sc1-wi.pdf"), width=10, height = 7.78)


mcmc.alpha <- posterior$alpha
len <- dim(mcmc.alpha)[1]
r <- matrix(, nrow = n, ncol = 30)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 30
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y.origin[i], u, alpha = mcmc.alpha[round(runif(1,1,len)),i]))
  }
}
lgrid <- n
grid <- qnorm(ppoints(lgrid))
# qqnorm(r[, 1])
# points(grid, quantile(r[, 1], ppoints(lgrid), type = 2), 
#     xlim = c(-3, 3), col = "red")
traj <- matrix(NA, nrow = T, ncol = lgrid)
for (t in 1:T){
  traj[t, ] <- quantile(r[, t], ppoints(lgrid), type = 2)
}
l.band <- apply(traj, 2, quantile, prob = 0.025)
trajhat <- apply(traj, 2, quantile, prob = 0.5)
u.band <- apply(traj, 2, quantile, prob = 0.975)

ggplot(data = data.frame(grid = grid, l.band = l.band, trajhat = trajhat, 
                         u.band = u.band)) + 
  geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
              fill = "steelblue",
              alpha = 0.4, linetype = "dashed") + 
  geom_line(aes(x = grid, y = trajhat), colour = "steelblue", linetype = "dashed", linewidth = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.2) + 
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
  theme_minimal(base_size = 20) +
  theme(text = element_text(size = 20)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_qqplot_sc1-wi.pdf"), width=10, height = 7.78)

# lambda.container <- data.frame("x" = seq(0, max(posterior$lambda2), length.out = 1000),
#                                "GamDist" = dgamma(seq(0, max(posterior$lambda2), length.out = 1000), 1, 1e-3),
#                                "lambda.post" = posterior$lambda2)


# ggplot(data = lambda.container, aes(x = x)) + ylab("density") + xlab("lambdas") + labs(col = "") +
#   geom_line(aes(x=x, y=GamDist), color = "red", linewidth = 0.7) +
#   geom_density(aes(x=lambda.post), color = "steelblue", linewidth = 1) +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#         axis.text = element_text(size = 35))
