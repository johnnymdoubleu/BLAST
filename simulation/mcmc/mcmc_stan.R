library(npreg)
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
# library(ggplotify)

#Scenario 1
# set.seed(2)

n <- 5000
psi <- 20
threshold <- 0.90
p <- 5
no.theta <- 1
simul.no <- 50

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)

sample_meanvector <- runif(p,0,1)
sample_covariance_matrix <- matrix(NA, nrow = p, ncol = p)
diag(sample_covariance_matrix) <- 1
# set.seed(666)

mat_Sim <- matrix(data = NA, nrow = p, ncol = p)
U <- runif(n = p) * 0.5

for(i in 1 : p)
{
  if(i <= (p*1.3/2))
  {
    U_Star <- pmin(U + 0.25 * runif(n = p), 0.99999)
    
  }else
  {
    U_Star <- pmin(pmax(U + sample(c(0, 1), size = p, replace = TRUE) * runif(n = p), 0.00001), 0.99999)
  }
  
  mat_Sim[, i] <- qnorm(U_Star)  
}

cor_Mat <- cor(mat_Sim)
sample_covariance_matrix <- cor_Mat * (p/2)

## create multivariate normal distribution
x.origin <- mvrnorm(n = n, mu = sample_meanvector, Sigma = sample_covariance_matrix)


# x.origin <- cbind(replicate(p, runif(n, 0, 1)))

corrplot.mixed(cor(x.origin),
                upper = "circle",
                lower = "number",
                addgrid.col = "black")
for(i in 1:p){
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    # tps <- bSpline(x = x.origin[,i], Boundary.knots = c(min(x.origin[,i]),max(x.origin[,i])), df = psi, intercept = TRUE)
    # tps <- mSpline(x.origin[,i], df=psi, Boundary.knots = range(x.origin[,i]), degree = 3, intercept=TRUE)
    #   bs.x <- cbind(bs.x, tps)
    bs.linear <- cbind(bs.linear, tps[,1:no.theta])
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
}

gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
    for (ps in 1:psi){
        if(j %in% c(2,3,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
        else if(j==7){
            if(ps <= (psi/2)){gamma.origin[ps, j] <- 0.01}
            else{gamma.origin[ps, j] <- 0.01}
        }
        else {
            if(ps <= (psi/2)){gamma.origin[ps, j] <- 0.01}
            else{gamma.origin[ps, j] <- 0.01}
        }
    }
}

# theta.origin <- matrix(, nrow = 2, ncol = p)
# for(j in 1:p){
#     for (k in 1:2){
#         if(j %in% c(2,4,5,6,9,10)){theta.origin[k, j] <- 0}
#         else if(j==7){
#             if(k==1){theta.origin[k, j] <- 0.5}
#             else{theta.origin[k, j] <- -0.3}
#         }
#         else {
#             if(k==1){theta.origin[k,j] <- -0.2}
#             else{theta.origin[k,j] <- 0.8}
#         }
#     }
# }
theta.origin <- c(0, 0.04, 0, 0, 0, 0)

f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
    f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp((theta.origin[1]*p) + sum(f.origin[i,]))
    y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i])
}

u <- quantile(y.origin, threshold)
x.origin <- x.origin[which(y.origin>u),]
x.origin <- scale(x.origin)
y.origin <- y.origin[y.origin > u]
n <- length(y.origin)

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
newx <- seq(0, 1, length.out=n)
xholder <- bs.x <- matrix(, nrow = n, ncol = p)
for(i in 1:p){
    xholder[,i] <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = n)  
    test.knot <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
    splines <- basis.tps(newx, test.knot, m=2, rk=FALSE, intercept = FALSE)
    xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
    xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    bs.linear <- cbind(bs.linear, tps[,1:no.theta])
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
}

f.nonlinear.new <- f.linear.new <- f.new <- f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
    f.linear.new[,j] <- xholder.linear[, j] * theta.origin[j+1]
    f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}

true.alpha <- alp.new <- alp.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp((theta.origin[1]) + sum(f.origin[i,]))
    alp.new[i] <- exp((theta.origin[1]) + sum(f.new[i,]))
}

# lambda1 ~ gamma(1, 1.78);
# intercept ~ double_exponential(0, lambda1);
# lambda2 ~ gamma(0.1, 0.1);
# sigma ~ inv_gamma(0.01, 0.01);
# for (j in 1:p){
#     theta[j] ~ double_exponential(0, lambda1);
#     tau[j] ~ gamma(atau, square(lambda2));
#     gamma[j] ~ multi_normal(rep_vector(0, psi), diag_matrix(rep_vector(1,psi)) * tau[j] * sigma);
# }
# // likelihood
# for (i in 1:n){
#     target += pareto_lpdf(y[i] | u, alpha[i]);
# }

write("// Stan model for simple linear regression
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
    vector[psi] gamma[p]; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real <lower=-1, upper = 1> sigma; //
    vector[p] tau;
}

transformed parameters {
    vector[n] alpha; // tail index
    vector[n] newalpha; // tail index
    matrix[n, p] gsmooth; // nonlinear component
    matrix[n, p] newgsmooth; // nonlinear component
    cov_matrix[psi] covmat[p]; // covariance
    for(j in 1:p){
        covmat[j] <- diag_matrix(rep_vector(1, psi)) * tau[j] * sigma;
    };
    for (j in 1:p){
        gsmooth[,j] <- bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgsmooth[,j] <- xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    for (i in 1:n){
        alpha[i] <- exp(theta[1] + dot_product(bsLinear[i], theta[2:newp]) + (gsmooth[i,] *rep_vector(1, p)));
        newalpha[i] <- exp(theta[1] + dot_product(xholderLinear[i], theta[2:newp]) + (newgsmooth[i,] *rep_vector(1, p)));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 0.1, 2);
    target += gamma_lpdf(lambda2 | 1, 1);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    target += normal_lpdf(theta[1] | 0, 0.01);
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, (square(lambda2)/2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), covmat[j]);
    }
}
//generated quantities {} // Used in Posterior predictive check"
, "model_simulation.stan")

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), newp = (p+1),
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear)

set_cmdstan_path(path = NULL)
#> CmdStan path set to: /Users/jgabry/.cmdstan/cmdstan-2.32.2

# Create a CmdStanModel object from a Stan program,
# here using the example model that comes with CmdStan
file <- file.path(cmdstan_path(), "model_simulation.stan")

init.alpha <- list(list(gamma = array(rep(0,(psi*p)), dim=c(psi, p)),
                        theta = rep(0, (p+1)), 
                        tau = rep(0.01, p), sigma = 0.001, 
                        lambda1 = 0.1, lambda2 = 0.1),
                  list(gamma = array(rep(0.02,(psi*p)), dim=c(psi, p)),
                        theta = rep(0.001, (p+1)), 
                        tau = rep(0.01, p), sigma = 0.001,
                        lambda1 = 0.1, lambda2 = 0.1),
                  list(gamma = array(rep(-0.02, (psi*p)), dim=c(psi, p)),
                        theta = rep(-0.02, (p+1)), 
                        tau = rep(0.01, p), sigma = 0.001,
                        lambda1 = 0.1, lambda2 = 0.1))

# stanc("C:/Users/Johnny Lee/Documents/GitHub/BRSTIR/application/model1.stan")
fit1 <- stan(
    file = "model_simulation.stan",  # Stan program
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    # init_r = 1,
    chains = 3,             # number of Markov chains
    warmup = 2000,          # number of warmup iterations per chain
    iter = 6000,            # total number of iterations per chain
    cores = 4,              # number of cores (could use one per chain)
    refresh = 1             # no progress shown
)

# saveRDS(fit1, file=paste0("./BRSTIR/application/",Sys.Date(),"_stanfit.rds"))
posterior <- extract(fit1)
str(posterior)

# print(as.mcmc(fit1), pars=c("alpha", "gamma", "intercept", "theta", "lambda1", "lambda2","lp__"), probs=c(.05,.5,.95))
plot(fit1, plotfun = "trace", pars = c("theta"), nrow = 3)
# ggsave(paste0("./BRSTIR/simulation/results/",Sys.Date(),"_mcmc_theta_trace.pdf"), width=10, height = 7.78)
plot(fit1, plotfun = "trace", pars = c("lambda1", "lambda2"), nrow = 2)
# ggsave(paste0("./BRSTIR/simulation/results/",Sys.Date(),"_mcmc_lambda.pdf"), width=10, height = 7.78)


theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary

gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q3 <- theta.samples[,6]


# theta.samples <- data.frame(apply(posterior$theta, 2, summary))


df.theta <- data.frame("seq" = seq(1, (p+1)),
                        "true" = theta.origin,
                        "m" = theta.post.mean,
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
#   scale_color_discrete(labels = c(expression(theta[0]),colnames(fwi.scaled))) + 
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

# ggsave(paste0("./BRSTIR/simulation/results/",Sys.Date(),"_mcmc_theta.pdf"), width=10, height = 7.78)

# ggplot(data.frame(group = factor(1:(p+1)), m=theta.post.mean, l = theta.q1, u = theta.q3), 
#        aes(group)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1.4) +
#   geom_point(aes(x = group, y = m), size = 4.5) + 
#   #geom_point(aes(x = group, y = beta), shape=8, size = 4.5, col="red")+
#   geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, size = 1.2) + 
#   labs(x = "Regression coefficients", y = "") + 
#   ylim(-5,5) + 
#   scale_x_discrete(labels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC"))+
#                               #,expression(beta[8]),
#                               #expression(beta[9]))) + 
#   theme_minimal(base_size = 30) + 
#   theme(text = element_text(size = 30), 
#         axis.text.x = element_text(angle = 0, hjust = 0.5))

df.gamma <- data.frame("seq" = seq(1, (psi*p)), 
                  "true" = as.vector(gamma.origin),
                  "m" = as.vector(gamma.post.mean),
                  "l" = as.vector(gamma.q1),
                  "u" = as.vector(gamma.q3))
df.gamma$covariate <- factor(rep(seq(1, 1 + nrow(df.gamma) %/% psi), each = psi, length.out = nrow(df.gamma)))
df.gamma$labels <- factor(1:(psi*p))
ggplot(df.gamma, aes(x =labels, y = m, color = covariate)) + 
  geom_point(size = 4) + ylab("") + xlab("" ) + #ylim(-15,15) +
  # geom_ribbon(aes(ymin = l, ymax = u)) +
  geom_errorbar(aes(ymin = l, ymax = u), width = 4, linewidth = 1.2) + 
  geom_point(size = 4, color = "black") + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+7), 
                    label = c(expression(bold(gamma[1])), 
                              expression(bold(gamma[2])), expression(bold(gamma[3])), expression(bold(gamma[4])), expression(bold(gamma[5]))),
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
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_gamma.pdf"), width=10, height = 7.78)

g.nonlinear.q1 <- g.linear.q1 <- g.q1 <- g.nonlinear.q3 <- g.linear.q3 <- g.q3 <- g.nonlinear.new <- g.linear.new <- g.new <- matrix(, nrow = n, ncol=p)
g.smooth.q1 <- g.smooth.q3 <- g.smooth.new <- alpha.new <- NULL
for (j in 1:p){
  g.linear.new[,j] <- xholder.linear[,j] * theta.post.mean[(j+1)]
  g.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.post.mean, nrow=psi)[,j] 
  g.new[1:n, j] <- g.linear.new[,j] + g.nonlinear.new[,j]
  g.linear.q1[,j] <- xholder.linear[,j] * theta.q1[(j+1)]
  g.nonlinear.q1[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q1, nrow=psi)[,j] 
  g.q1[1:n, j] <- g.linear.q1[,j] + g.nonlinear.q1[,j]
  g.linear.q3[,j] <- xholder.linear[,j] * theta.q3[(j+1)]
  g.nonlinear.q3[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q3, nrow=psi)[,j] 
  g.q3[1:n, j] <- g.linear.q3[,j] + g.nonlinear.q3[,j]
}

for(i in 1:n){
  g.smooth.new[i] <- theta.post.mean[1] + sum(g.new[i,])
  g.smooth.q1[i] <- theta.q1[1] + sum(g.q1[i,])
  g.smooth.q3[i] <- theta.q3[1] + sum(g.q3[i,])
}

### Plotting linear and nonlinear components
# post.mean <- as.vector(apply(as.data.frame(matrix(alpha.summary[(n+1):(n+(n*p)),1], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q1 <- as.vector(apply(as.data.frame(matrix(alpha.summary[(n+1):(n+(n*p)),4], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q3 <- as.vector(apply(as.data.frame(matrix(alpha.summary[(n+1):(n+(n*p)),5], nrow = n, ncol = p)), 2, sort, decreasing=F))
equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

data.smooth <- data.frame("x"=c(1:n),
                          "true" = as.vector(f.new),
                          "post.mean" = as.vector(g.new),
                          "q1" = as.vector(g.q1),
                          "q3" = as.vector(g.q3),
                          "covariates" = gl(p, n, (p*n)),
                          "replicate" = gl(2, n, (p*n)))
ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=true, colour = covariates, linetype = "True"), linewidth=2) + 
  geom_line(aes(y=post.mean, colour = covariates, linetype = "Posterior Mean"), linewidth=2) + 
  ylab("") + xlab ("Smooth Functions") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  facet_grid(covariates ~ .) +
  scale_linetype_manual("functions",values=c("Posterior Mean"=3,"True"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))

# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_smooth.pdf"), width=10.5, height = 15)
data.linear <- data.frame("x"=c(1:n),
                          "true" = as.vector(f.linear.new),
                          "post.mean" = as.vector(g.linear.new),
                          "q1" = as.vector(g.linear.q1),
                          "q3" = as.vector(g.linear.q3),
                          "covariates" = gl(p, n, (p*n)),
                          "replicate" = gl(2, n, (p*n)))
ggplot(data.linear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=true, colour = covariates, linetype = "True"), linewidth=2) + 
  geom_line(aes(y=post.mean, colour = covariates, linetype = "Posterior Mean"), linewidth=2) + xlab("Linear Components") + ylab("") +
  facet_grid(covariates ~ ., scales = "free_y") + 
  scale_linetype_manual("functions",values=c("Posterior Mean"=3,"True"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_linear.pdf"), width=10.5, height = 15)
# post.mean <- as.vector(apply(as.data.frame(matrix(alpha.summary[((n+(n*p))+1):(n+(2*n*p)),1], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q1 <- as.vector(apply(as.data.frame(matrix(alpha.summary[((n+(n*p))+1):(n+(2*n*p)),4], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q3 <- as.vector(apply(as.data.frame(matrix(alpha.summary[((n+(n*p))+1):(n+(2*n*p)),5], nrow = n, ncol = p)), 2, sort, decreasing=F))

data.nonlinear <- data.frame("x"=c(1:n),
                          "true" = as.vector(f.nonlinear.new),
                          "post.mean" = as.vector(g.nonlinear.new),
                          "q1" = as.vector(g.nonlinear.q1),
                          "q3" = as.vector(g.nonlinear.q3),
                          "covariates" = gl(p, n, (p*n)),
                          "replicate" = gl(2, n, (p*n)))
ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=true, colour = covariates, linetype = "True"), linewidth=2) + 
  geom_line(aes(y=post.mean, colour = covariates, linetype = "MCMC"), linewidth=2) + xlab("Nonlinear Components") + ylab("") +
  facet_grid(covariates ~ ., scales = "free_y") + 
  scale_linetype_manual("functions",values=c("MCMC"=3,"True"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=45),
        legend.margin=margin(0,0,0,-10),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_nonlinear.pdf"), width=12.5, height = 15)

data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "true" = sort(alp.new),
                            "post.mean" = sort(newalpha.samples[,1]),
                            "post.median" = sort(newalpha.samples[,5]),
                            "q1" = sort(newalpha.samples[,4]),
                            "q3" = sort(newalpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(x))) + xlab(expression(x)) + labs(col = "") + 
  geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=2, linetype = 2) + 
  ylim(0, (max(data.scenario$post.mean)+10)) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linetype=2, linewidth=2) +
  # geom_line(aes(y=chain2, col = "Chian 2"), linetype=3) +
  # facet_grid(covariates ~ .) + 
  # scale_y_continuous(breaks=c(0)) + 
  scale_color_manual(values = c("blue","#e0b430","red"))+
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position="top", 
        legend.key.size = unit(1, 'cm'),
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_alpha.pdf"), width=10, height = 7.78)

mcmc.gamma <- posterior$gamma
gamma.container <- as.data.frame(matrix(NA, nrow = 20, ,ncol = 0))
for(i in 1:(dim(mcmc.gamma)[1]/3)){
    gamma.container <- cbind(gamma.container, mcmc.gamma[i,1,])    
}
# corrplot.mixed(cor(gamma.container),
#                 upper = "circle",
#                 lower = "number",
#                 addgrid.col = "black")

mcmc.alpha <- posterior$alpha
len <- dim(mcmc.alpha)[1]
r <- matrix(, nrow = n, ncol = 30)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 30
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y.origin[i], u, alpha = mcmc.alpha[round(runif(1,1,len)),i]))
    # r[i, t] <- qnorm(pPareto(y[i], u, alpha = alpha.new[i]))
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
  #geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
  #           color = "lightgrey", fill = "lightgrey",
  #          alpha = 0.4, linetype = "dashed") + 
  geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
              # color = "darkgrey", fill = "darkgrey",
              alpha = 0.4, linetype = "dashed") + 
  geom_line(aes(x = grid, y = trajhat), linetype = "dashed", linewidth = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.2) + 
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
  theme_minimal(base_size = 20) +
  theme(text = element_text(size = 20)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_qqplot.pdf"), width=10, height = 7.78)              
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))

cat("sc1_Alp Done")

