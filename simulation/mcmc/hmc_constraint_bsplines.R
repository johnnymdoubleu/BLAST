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
library(cmdstanr)
library(ggh4x)

#Scenario 1
# set.seed(10)
# set.seed(6)


n <- 10000
psi <- 2
threshold <- 0.95
p <- 2
no.theta <- 1
simul.no <- 50

# Function to generate Gaussian copula
C <- diag(p)
xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))

# x.origin <- apply(x.origin, 2, sort, decreasing=F)
# end.holder <- basis.holder <- matrix(, nrow = 2, ncol =0)
# index.holder <- matrix(, nrow = 0, ncol = 2)
# for(i in 1:p){
#   index.holder <- rbind(index.holder, 
#                       matrix(c(which.min(x.origin[,i]),
#                               which.max(x.origin[,i])), ncol=2))
# }
for(i in 1:p){
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = (psi+2))[2:(psi+1)]  
    bsp <- bSpline(x.origin[,i], knots= knots, degree = 3)
    bs.nonlinear <- cbind(bs.nonlinear, bsp)
    # basis.holder <- cbind(basis.holder, 
    #         solve(t(matrix(c(bsp[index.holder[i,1], 1],
    #                 bsp[index.holder[i,1], psi],
    #                 bsp[index.holder[i,2], 1],
    #                 bsp[index.holder[i,2], psi]), 
    #                 nrow = 2, ncol = 2))))
}


## Generate sample
gamma.origin <- matrix(, nrow = (psi+3), ncol = p)
for(j in 1:p){
    for (ps in 1:(psi+3)){
        if(j %in% c(1,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
        else {
            gamma.origin[ps, j] <- -2.5
        }
    }
}
theta.origin <- -0.01

f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    # f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*(psi+3))+1):(((j-1)*(psi+3))+(psi+3))] %*% gamma.origin[,j]
    # f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp(theta.origin + sum(f.nonlinear.origin[i,]))
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
    test.knot <- seq(min(xholder[,i]), max(xholder[,i]), length.out = (psi+2))[2:(psi+1)]
    splines <- bSpline(xholder[,i], knots = test.knot, degree = 3)    
    xholder.nonlinear <- cbind(xholder.nonlinear, splines)
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = (psi+2))[2:(psi+1)]    
    bsp <- bSpline(x.origin[,i], knots = knots, degree = 3)
    # basis.holder <- cbind(basis.holder, 
    #         solve(t(matrix(c(bsp[index.holder[i,1], 1],
    #                 bsp[index.holder[i,1], psi],
    #                 bsp[index.holder[i,2], 1],
    #                 bsp[index.holder[i,2], psi]), 
    #                 nrow = 2, ncol = 2))))
    # end.holder <- cbind(end.holder, 
    #             matrix(c(tps[index.holder[i,1], no.theta+1],
    #                 tps[index.holder[i,1], no.theta+psi],
    #                 tps[index.holder[i,2], no.theta+1],
    #                 tps[index.holder[i,2], no.theta+psi]), 
    #                nrow = 2, ncol = 2))
    bs.nonlinear <- cbind(bs.nonlinear, bsp) 
}


f.nonlinear.new <- f.linear.new <- f.new <- f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
f.sub.origin <- matrix(, nrow = 2, ncol = p)
# for(j in 1:p){
#     f.sub.origin[,j] <- as.matrix(bs.nonlinear[index.holder[j,], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))], nrow = 2) %*% gamma.origin[(2:(psi-1)), j]
#     gamma.origin[c(1,psi),j] <- -1 * basis.holder[,(((j-1)*2)+1):(((j-1)*2)+2)] %*% as.matrix(f.sub.origin[,j], nrow=2)
# }
for(j in 1:p){
    # f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*(psi+3))+1):(((j-1)*(psi+3))+(psi+3))] %*% gamma.origin[,j]
    # f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
    # f.linear.new[,j] <- xholder.linear[, j] * theta.origin[j+1]
    f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*(psi+3))+1):(((j-1)*(psi+3))+(psi+3))] %*% gamma.origin[,j]
    # f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}

true.alpha <- alp.new <- alp.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp(theta.origin + sum(f.nonlinear.origin[i,]))
    alp.new[i] <- exp(theta.origin + sum(f.nonlinear.new[i,]))
}

write("// Stan model for BRSTIR Pareto Uncorrelated Samples
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n, ((psi+3)*p)] bsNonlinear; // thin plate splines basis
    matrix[n, ((psi+3)*p)] xholderNonlinear; // thin plate splines basis    
    array[n] real <lower=1> y; // extreme response
    real <lower=0> atau;
}
parameters {
    real theta; // linear predictor
    vector[(psi+3)] gamma[p]; // splines coefficient 
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
        gnl[,j] = bsNonlinear[,(((j-1)*(psi+3))+1):(((j-1)*(psi+3))+(psi+3))] * gamma[j];
        newgnl[,j] = xholderNonlinear[,(((j-1)*(psi+3))+1):(((j-1)*(psi+3))+(psi+3))] * gamma[j];
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
    target += gamma_lpdf(lambda | 1, 0.0001);
    target += (p * (psi+3) * log(lambda)/2);
    for (j in 1:p){
        target += inv_gamma_lpdf(sigma | 0.01, 0.01); 
        target += gamma_lpdf(tau[j] | atau, sqrt(lambda/2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, (psi+3)), diag_matrix(rep_vector(1, (psi+3))) * tau[j] * sigma);
    }
}
"
, "model_simulation_sc1_constraint.stan")
    # matrix[2, (2*p)] basisFL;
    # array[(p*2)] int indexFL;
    # vector[2] gammaFL[p]; 
    # matrix[2, p] subgnl;
    # for(j in 1:p){
    #     gamma[j][2:(psi-1)] = gammaTemp[j][1:(psi-2)];
    #     subgnl[,j] = bsNonlinear[indexFL[(((j-1)*2)+1):(((j-1)*2)+2)], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j];
    #     gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * subgnl[,j] * -1;
    #     gamma[j][1] = gammaFL[j][1];
    #     gamma[j][psi] = gammaFL[j][2];
    # };

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                    atau = (((psi+3)+1)/2), 
                    # basisFL = basis.holder, indexFL = as.vector(t(index.holder)),
                    bsNonlinear = bs.nonlinear, xholderNonlinear = xholder.nonlinear)

init.alpha <- list(list(gamma = array(rep(0, ((psi+3)*p)), dim=c((psi+3),p)),
                        theta = 0, tau = rep(0.1, p), sigma = 0.1, lambda = 0.1)
                #   list(gammaTemp = array(rep(-0.2, ((psi-2)*p)), dim=c((psi-2),p)),
                #         theta = rep(0.01, (p+1)),
                #         tau = rep(0.01, p), sigma = 0.001,
                #         # lambda1 = 0.1, 
                #         lambda2 = 0.001),
                #   list(gammaTemp = array(rep(-0.5, ((psi-2)*p)), dim=c((psi-2),p)),
                #         theta = rep(-0.05, (p+1)),
                #         tau = rep(0.5, p), sigma = 0.01,
                #         #lambda1 = 0.01, 
                #         lambda2 = 0.05)
                        )
# setwd("C:/Users/Johnny Lee/Documents/GitHub")
fit1 <- stan(
    file = "model_simulation_sc1_constraint.stan",  # Stan program
    # file = "model_BRSTIR.stan",  # Stan program
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    chains = 1,             # number of Markov chains
    # warmup = 1000,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 500             # no progress shown
)

posterior <- extract(fit1)

plot(fit1, plotfun = "trace", pars = c("theta"), nrow = 3)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_theta_trace_sc1-wi.pdf"), width=10, height = 7.78)
# plot(fit1, plotfun = "trace", pars = c("lambda1", "lambda2"), nrow = 2)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_lambda_sc1-wi.pdf"), width=10, height = 7.78)

# tau.samples <- summary(fit1, par=c("tau"), probs = c(0.05,0.5, 0.95))$summary
theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
newgnl.samples <- summary(fit1, par=c("newgnl"), probs = c(0.05, 0.5, 0.95))$summary
newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary

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

df.gamma <- data.frame("seq" = seq(1, ((psi+3)*p)), 
                  "true" = as.vector(gamma.origin),
                  "m" = as.vector(gamma.q2),
                  "l" = as.vector(gamma.q1),
                  "u" = as.vector(gamma.q3))
df.gamma$covariate <- factor(rep(seq(1, 1 + nrow(df.gamma) %/% (psi+3)), each = (psi+3), length.out = nrow(df.gamma)))
df.gamma$labels <- factor(1:((psi+3)*p))
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

g.nonlinear.mean <- as.vector(matrix(newgnl.samples[,1], nrow = n, byrow=TRUE))
g.nonlinear.q1 <- as.vector(matrix(newgnl.samples[,4], nrow = n, byrow=TRUE))
g.nonlinear.q2 <- as.vector(matrix(newgnl.samples[,5], nrow = n, byrow=TRUE))
g.nonlinear.q3 <- as.vector(matrix(newgnl.samples[,6], nrow = n, byrow=TRUE))


equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}


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
          fill = guide_legend(order = 1)) + 
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
  theme_minimal(base_size = 30) +
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

lambda.container <- data.frame("x" = seq(0, max(posterior$lambda), length.out = 1000),
                        "GamDist" = dgamma(seq(0, max(posterior$lambda), length.out = 1000), 0.01, 0.01),
                        "lambda.post" = posterior$lambda)

                        
ggplot(data = lambda.container, aes(x = x)) + ylab("density") + xlab("lambdas") + labs(col = "") +
    geom_line(aes(x=x, y=GamDist), color = "red", linewidth = 0.7) +
    geom_density(aes(x=lambda.post), color = "steelblue", linewidth = 1) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                axis.text = element_text(size = 35))
