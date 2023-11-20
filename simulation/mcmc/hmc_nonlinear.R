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
# library(simstudy)
# library(ggplotify)

#Scenario 1
set.seed(36)
# set.seed(50)

# n <- 5000
# psi <- 20
# threshold <- 0.90
# p <- 10
# no.theta <- 1
# simul.no <- 50

# xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
# x.origin <- cbind(replicate(p, runif(n, 0, 1)))
# for(i in 1:p){
#     knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
#     tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
#     # tps <- mSpline(x.origin[,i], df=psi, Boundary.knots = range(x.origin[,i]), degree = 3, intercept=TRUE)
#     #   bs.x <- cbind(bs.x, tps)
#     bs.linear <- cbind(bs.linear, tps[,1:no.theta])
#     bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
# }

# gamma.origin <- matrix(, nrow = psi, ncol = p)
# for(j in 1:p){
#     for (ps in 1:psi){
#         if(j %in% c(2,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
#         else if(j==7){
#             if(ps <= (psi/2)){gamma.origin[ps, j] <- 1}
#             else{gamma.origin[ps, j] <- 1}
#         }
#         else {
#             if(ps <= (psi/2)){gamma.origin[ps, j] <- 1}
#             else{gamma.origin[ps, j] <- 1}
#         }
#     }
# }

# theta.origin <- c(-0.1, 0.8, 0, 0.8, 0, 0, 0, -0.3, 0.8, 0, 0)

n <- 5000
psi <- 20
threshold <- 0.90
p <- 5
no.theta <- 1
simul.no <- 50

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)


C <- matrix(c(1, 0.3, 0.5, 0.3, 0.3,
              0.3, 1, 0.95, 0.4, 0.4,
              0.5, 0.95, 1, 0.5, 0.1,
              0.3, 0.4, 0.5 , 1, 0.5,
              0.3, 0.4, 0.5, 0.5, 1), nrow = p)
x.origin <- tmvnsim(n = n, k = p, lower = rep(0, p), means = rep(0, p), sigma = C)$samp
# x.origin <- scale(x.origin)

# corrplot.mixed(cor(x.origin),
#                 upper = "circle",
#                 lower = "number",
#                 addgrid.col = "black")

for(i in 1:p){
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    bs.linear <- cbind(bs.linear, tps[,1:no.theta])
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
}

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

theta.origin <- c(-0.5, 0, -0.2, -0.2, 0, 0)

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
x.origin <- x.origin[which(y.origin>u),]
# x.bs <- x.origin
y.origin <- y.origin[y.origin > u]
n <- length(y.origin)

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
newx <- seq(0, 1, length.out=n)
xholder <- bs.x <- matrix(, nrow = n, ncol = p)
for(i in 1:p){
    # xholder[,i] <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = n)  
    # test.knot <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
    xholder[,i] <- seq(0, 0.1, length.out = n)
    test.knot <- seq(0, 0.1, length.out = psi)
    splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
    xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
    xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    bs.linear <- cbind(bs.linear, tps[,1:no.theta])
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
}

f.nonlinear.new <- f.linear.new <- f.new <- f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    # f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    # f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
    # f.linear.new[,j] <- xholder.linear[, j] * theta.origin[j+1]
    f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    # f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}

true.alpha <- alp.new <- alp.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp(theta.origin[1] + sum(f.nonlinear.origin[i,]))
    alp.new[i] <- exp(theta.origin[1] + sum(f.nonlinear.new[i,]))
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
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector[n] y; // extreme response
    real <lower=0> atau;
}

parameters {
    array[p] vector[psi] gamma; // splines coefficient
    real theta0; // intercept term
    real <lower=0> lambda; // group lasso penalty
    real <lower=0> sigma; //
    array[p] real <lower=0> tau;
}

transformed parameters {
    array[n] real <lower=0> alpha; // tail index
    matrix[n, p] gsmooth; // nonlinear component
    array[n] real <lower=0> newalpha; // tail index
    matrix[n, p] newgsmooth; // nonlinear component

    for (j in 1:p){
        gsmooth[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgsmooth[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    for (i in 1:n){
        alpha[i] = exp(theta0 + sum(gsmooth[i,])); //gsmooth[i,] * rep_vector(1, p)
        newalpha[i] = exp(theta0 + sum(newgsmooth[i,])); // newalpha[i] = exp(newgsmooth[i,] * rep_vector(1, p))
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda | 0.1, 0.1);
    target += normal_lpdf(theta0 | 0, 10);
    target += inv_gamma_lpdf(sigma | 0.1, 0.1);
    target += (p * psi * log(lambda));
    for (j in 1:p){
        target += gamma_lpdf(tau[j] | atau, lambda/sqrt(2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * sqrt(tau[j]) * sqrt(sigma));
    }
}
generated quantities {
    // Used in Posterior predictive check
    vector[n] log_lik;
    array[n] real y_rep = pareto_rng(rep_vector(u, n), alpha);
    for (i in 1:n) {
        log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}
"
, "model_simulation_sc3.stan")

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), newp = (p+1),
                    bsNonlinear = bs.nonlinear,
                    xholderNonlinear = xholder.nonlinear)

# set_cmdstan_path(path = NULL)
#> CmdStan path set to: /Users/jgabry/.cmdstan/cmdstan-2.32.2

# Create a CmdStanModel object from a Stan program,
# here using the example model that comes with CmdStan
# file <- file.path(cmdstan_path(), "model_simulation.stan")

init.alpha <- list(list(gamma = array(rep(0, (psi*p)), dim=c(psi, p)),
                        tau = rep(0.1, p), sigma = 0.1, lambda = 0.1),
                  list(gamma = array(rep(0.02, (psi*p)), dim=c(psi, p)),
                        tau = rep(0.1, p), sigma = 0.1, lambda = 0.1),
                  list(gamma = array(rep(0.05, (psi*p)), dim=c(psi, p)),
                        tau = rep(0.1, p), sigma = 0.1, lambda = 0.1))

fit1 <- stan(
    file = "model_simulation_sc3.stan",  # Stan program
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    chains = 3,             # number of Markov chains
    warmup = 1000,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = 4,              # number of cores (could use one per chain)
    refresh = 500             # no progress shown
)

# saveRDS(fit1, file=paste0("./BRSTIR/application/",Sys.Date(),"_stanfit.rds"))
posterior <- extract(fit1)
str(posterior)

# plot(fit1, plotfun = "trace", pars = c("theta"), nrow = 3)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_theta_trace_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_theta_trace_sc3-nl.pdf"), width=10, height = 7.78)
plot(fit1, plotfun = "trace", pars = c("theta0", "lambda"), nrow=2)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_theta_lambda_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_theta_lambda_sc3-nl.pdf"), width=10, height = 7.78)

theta.samples <- summary(fit1, par=c("theta0"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary
summary(fit1, par=c("sigma"), probs = c(0.05,0.5, 0.95))$summary
gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q3 <- theta.samples[,6]


df.gamma <- data.frame("seq" = seq(1, (psi*p)), 
                  "true" = as.vector(gamma.origin),
                  "m" = as.vector(gamma.post.mean),
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
                              expression(bold(gamma[5])), 
                              expression(bold(gamma[6])), 
                              expression(bold(gamma[7])), 
                              expression(bold(gamma[8])), 
                              expression(bold(gamma[9])), 
                              expression(bold(gamma[10]))),
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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_gamma_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_gamma_sc3-nl.pdf"), width=10, height = 7.78)

g.nonlinear.q1 <- g.linear.q1 <- g.q1 <- g.nonlinear.q3 <- g.linear.q3 <- g.q3 <- g.nonlinear.new <- g.linear.new <- g.new <- matrix(, nrow = n, ncol=p)
alpha.smooth.q1 <- alpha.smooth.q3 <- alpha.smooth.new <- alpha.new <- NULL
for (j in 1:p){
#   g.linear.new[,j] <- xholder.linear[,j] * theta.post.mean[(j+1)]
  g.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.post.mean, nrow=psi)[,j] 
#   g.new[1:n, j] <- g.linear.new[,j] + g.nonlinear.new[,j]
#   g.linear.q1[,j] <- xholder.linear[,j] * theta.q1[(j+1)]
  g.nonlinear.q1[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q1, nrow=psi)[,j] 
#   g.q1[1:n, j] <- g.linear.q1[,j] + g.nonlinear.q1[,j]
#   g.linear.q3[,j] <- xholder.linear[,j] * theta.q3[(j+1)]
  g.nonlinear.q3[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q3, nrow=psi)[,j] 
#   g.q3[1:n, j] <- g.linear.q3[,j] + g.nonlinear.q3[,j]
}

for(i in 1:n){
  alpha.smooth.new[i] <- exp(theta.post.mean + sum(g.nonlinear.new[i,]))
  alpha.smooth.q1[i] <- exp(theta.q1 + sum(g.nonlinear.q1[i,]))
  alpha.smooth.q3[i] <- exp(theta.q3 + sum(g.nonlinear.q3[i,]))
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

# data.smooth <- data.frame("x"=c(1:n),
#                           "true" = as.vector(f.new),
#                           "post.mean" = as.vector(g.nonlinear.new),
#                           "q1" = as.vector(g.nonlinear.q1),
#                           "q3" = as.vector(g.nonlinear.q3),
#                           "covariates" = gl(p, n, (p*n)),
#                           "replicate" = gl(2, n, (p*n)))
# ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
#   geom_line(aes(y=true, colour = covariates, linetype = "True"), linewidth=2) + 
#   geom_line(aes(y=post.mean, colour = covariates, linetype = "Posterior Mean"), linewidth=2) + 
#   ylab("") + xlab ("Smooth Functions") +
#   # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
#   facet_grid(covariates ~ .) +
#   scale_linetype_manual("functions",values=c("Posterior Mean"=3,"True"=1)) +
#   scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
#   theme(plot.title = element_text(hjust = 0.5, size = 30),
#         legend.position = "none",
#         plot.margin = margin(0,0,0,-10),
#         strip.text = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size=33),
#         axis.title.x = element_text(size = 35))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_smooth_sc2-nl.pdf"), width=10.5, height = 15)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_smooth_sc3-nl.pdf"), width=10.5, height = 15)

data.nonlinear <- data.frame("x"=c(1:n),
                          "true" = as.vector(f.nonlinear.new),
                          "post.mean" = as.vector(g.nonlinear.new),
                          "q1" = as.vector(g.nonlinear.q1),
                          "q3" = as.vector(g.nonlinear.q3),
                          "covariates" = gl(p, n, (p*n)),
                          "replicate" = gl(2, n, (p*n)))

                        
plot.nonlinear <- ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  # geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=true, colour = covariates, linetype = "True"), linewidth=2) + 
  geom_line(aes(y=post.mean, colour = covariates, linetype = "MCMC"), linewidth=2) + xlab("Nonlinear Components") + ylab("") +
  facet_wrap(covariates ~ ., scale="free_y", nrow = p)  +
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
plot.nonlinear
# plot.nonlinear + facetted_pos_scales(y = list(
#     covariates == "1" ~ ylim(-0.01, 0.01),
#     covariates == "2" ~ ylim(-0.06, 0),
#     covariates == "3" ~ ylim(-0.09, 0),
#     covariates == "4" ~ ylim(-0.001, 0.05),
#     covariates == "5" ~ ylim(-0.01, 0.01)))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_nonlinear_sc2-nl.pdf"), width=12.5, height = 15)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_nonlinear_sc3-nl.pdf"), width=12.5, height = 15)
data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "true" = sort(alp.origin),
                            "post.mean" = sort(alpha.samples[,1]),
                            "post.median" = sort(alpha.samples[,5]),
                            "q1" = sort(alpha.samples[,4]),
                            "q3" = sort(alpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(x))) + xlab(expression(x)) + labs(col = "") + 
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
  geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=2, linetype = 2) +
  # ylim(0, (max(data.scenario$post.mean)+10)) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linetype=2, linewidth=2) +
  # geom_line(aes(y=chain2, col = "Chian 2"), linetype=3) +
  # facet_grid(covariates ~ .) + 
  # scale_y_continuous(breaks=c(0)) + 
  scale_color_manual(values = c("blue","#e0b430","red"))+
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position="top", 
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=18),
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_alpha_train_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_alpha_train_sc3-nl.pdf"), width=10, height = 7.78)

data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "true" = sort(alp.new),
                            "post.mean" = sort(newalpha.samples[,1]),
                            "post.median" = sort(newalpha.samples[,5]),
                            "q1" = sort(newalpha.samples[,4]),
                            "q3" = sort(newalpha.samples[,6]))
                            # "post.mean" = sort(alpha.smooth.new),
                            # "post.median" = sort(newalpha.samples[,5]),
                            # "q1" = sort(alpha.smooth.q1),
                            # "q3" = sort(alpha.smooth.q3))                            

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(x))) + xlab(expression(x)) + labs(col = "") + 
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
  geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=2, linetype = 2) + 
  # ylim(0, (max(data.scenario$q3)+10)) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linetype=2, linewidth=2) +
  # geom_line(aes(y=chain2, col = "Chian 2"), linetype=3) +
  # facet_grid(covariates ~ .) + 
  # scale_y_continuous(breaks=c(0)) + 
  scale_color_manual(values = c("blue","#e0b430","red"))+
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position="top", 
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=18),
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_alpha_test_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_alpha_test_sc3-nl.pdf"), width=10, height = 7.78)
# mcmc.gamma <- posterior$gamma
# gamma.container <- as.data.frame(matrix(NA, nrow = 20, ,ncol = 0))
# for(i in 1:(dim(mcmc.gamma)[1]/3)){
#     gamma.container <- cbind(gamma.container, mcmc.gamma[i,1,])    
# }
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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_qqplot_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_qqplot_sc3-nl.pdf"), width=10, height = 7.78)           
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))

# for (i in 1:n){
#     alpha[i] <- exp(theta[1] + dot_product(bsLinear[i], theta[2:newp]) + (gsmooth[i,] * rep_vector(1, p)));    
# };

# target += gamma_lpdf(lambda1 | 1, 10);
# target += gamma_lpdf(lambda2 | 1, 100);
cat("sc1_Alp Done")
# set_cmdstan_path(path = NULL)
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
    array[p] vector[psi] gamma; // splines coefficient array[p] real <lower=0> tau
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real <lower=0> sigma; //
    array[p] real <lower=0> tau;
}

transformed parameters {
    array[n] real <lower=0> alpha; // tail index
    matrix[n, p] gsmooth; // nonlinear component
    array[n] real <lower=0> newalpha; // tail index
    matrix[n, p] newgsmooth; // nonlinear component
    for (j in 1:p){
        gsmooth[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgsmooth[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    for (i in 1:n){
        alpha[i] = exp(theta[1] + dot_product(bsLinear[i], theta[2:newp]) + (gsmooth[i,] * rep_vector(1, p)));
        newalpha[i] = exp(theta[1] + dot_product(xholderLinear[i], theta[2:newp]) + (newgsmooth[i,] * rep_vector(1, p)));        
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 0.1, 0.01);
    target += gamma_lpdf(lambda2 | 0.1, 0.01);
    target += normal_lpdf(theta[1] | 0, square(100));
    target += inv_gamma_lpdf(sigma | 0.01, 0.01); // target += double_exponential_lpdf(theta[1] | 0, lambda1)
    target += (p * log(lambda1) + (p*psi) * log(lambda2));
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, (square(lambda2)/2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * tau[j] * sigma);
    }
}
generated quantities {
    // Used in Posterior predictive check
    vector[n] log_lik;
    array[n] real y_rep = pareto_rng(rep_vector(u, n), alpha);
    for (i in 1:n) {
        log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}
"
, "C:/Users/Johnny Lee/Documents/.cmdstan/cmdstan-2.33.1/examples/model_simulation_sc3.stan")

file <- file.path(cmdstan_path(), "examples/model_simulation_sc3.stan")
mod <- cmdstan_model(file)

op <- mod$optimize(
  data = list(y = as.vector(y.origin), u = u, p = p, 
                      n= n, psi = psi, atau = ((psi+1)/2), newp = (p+1),
                      bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                      xholderLinear = xholder.linear, 
                      xholderNonlinear = xholder.nonlinear),
  init = 0,
  # init = list(list(gamma = t(gamma.origin),
  #                   theta = theta.origin)),
  # init = list(list(gamma = array(rep(0.01, (psi*p)), dim=c(p, psi)),
  #                 theta = rep(0.1, (p+1)), 
  #                 tau = rep(0.1, p), sigma = 0.1, 
  #                 lambda1 = 0.01, lambda2 = 0.01)),
  iter = 3500,
  algorithm = "newton",
  refresh = 50
)


# sm <- stan_model(model_code = stan.code)
# op <- optimizing(sm, data = list(y = as.vector(y.origin), u = u, p = p, 
#                     n= n, psi = psi, atau = ((psi+1)/2), newp = (p+1),
#                     bsLinear = bs.linear, bsNonlinear = bs.nonlinear),
#                     # xholderLinear = xholder.linear, 
#                     # xholderNonlinear = xholder.nonlinear),
#               init = list(gamma = array(rep(0.01, (psi*p)), dim=c(psi, p)),
#                         theta = rep(0.1, (p+1)), 
#                         tau = rep(0.1, p), sigma = 0.1, 
#                         lambda1 = 1, lambda2 = 1),
#               iter = 3500,
#               algorithm = "LBFGS",
#               verbose = TRUE)


# theta.map <- op$par[1:(p+1)]
# gamma.map <- as.vector(matrix(op$par[(p+1+1):(p+1+(psi*p))]))
# gamma.map <- as.vector(t(matrix(gamma.map, nrow=p)))
# lambda.map <- op$par[(p+2+(psi*p)):(p+3+(psi*p))]

theta.map <- as.vector(op$draws(variables = "theta"))
gamma.map <- as.vector(t(matrix(op$draws(variables = "gamma"), nrow = p)))
lambda.map <- as.vector(op$draws(variables = c("lambda1", "lambda2")))
alpha.map <- as.vector(op$draws(variables = "alpha"))
newalpha.map <- as.vector(op$draws(variables = "newalpha"))
# alpha.map <- op$par[(p+p+5+(psi*p)):(p+p+4+n+(psi*p))]
# newalpha.map <- op$par[(p+p+5+n+(psi*p)):(p+p+4+n+n+(psi*p))]

systime <- Sys.time()
Sys.time()
systime <- chartr(":","-",systime)
date <- gsub("-","", substr(systime, 1, 10))
time <- substr(systime, 12, 20)


df.theta <- data.frame("seq" = seq(1, p+1),
                  theta.map,
                  "theta.true" = theta.origin)
# df.theta$covariate <- factor(rep(seq(1, 1 + nrow(df.theta) %/% no.theta), each = no.theta, length.out = nrow(df.theta)))
df.theta$covariate <- factor(0:p)
df.theta$labels <- factor(0:p)
ggplot(df.theta, aes(x = labels)) + ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_point(aes(y = theta.map, col = covariate), size = 6) + 
  geom_point(aes(y = theta.true), color="red", size = 4) +
  # labs(title=expression("MAP vs True for"~theta)) + 
  scale_x_discrete(labels = c(expression(bold(theta[0])),
                              expression(bold(theta[1])),
                              expression(bold(theta[2])),
                              expression(bold(theta[3])),
                              expression(bold(theta[4])),
                              expression(bold(theta[5])),
                              expression(bold(theta[6])),
                              expression(bold(theta[7])),
                              expression(bold(theta[8])),
                              expression(bold(theta[9])),
                              expression(bold(theta[10]))))+
  theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size=30),
          # axis.ticks.x = element_blank(),
          axis.text.x = element_text(hjust=0.3),
          axis.text = element_text(size = 30),
          panel.grid.minor.x = element_blank())
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_theta_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_theta_sc3-nl.pdf"), width=10, height = 7.78)

df <- data.frame("seq" = seq(1, (psi*p)), 
                  gamma.map, 
                  "gamma.true" = as.vector(gamma.origin))
df$covariate <- factor(rep(seq(1, 1 + nrow(df) %/% psi), each = psi, length.out = nrow(df)))
df$labels <- factor(1:(psi*p))

ggplot(df, aes(x =labels , y = gamma.map, col = covariate)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("")+
  geom_point(aes(y = gamma.true), color = "red", size=3) + 
  geom_point(size = 4) + ylab("") +
  scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+15), 
                    label = c(expression(bold(gamma[1])), 
                              expression(bold(gamma[2])), 
                              expression(bold(gamma[3])), 
                              expression(bold(gamma[4])), 
                              expression(bold(gamma[5])), 
                              expression(bold(gamma[6])),
                              expression(bold(gamma[7])),
                              expression(bold(gamma[8])),
                              expression(bold(gamma[9])),
                              expression(bold(gamma[10]))), 
                    expand=c(0,5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size=30),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(hjust=1),
          axis.text = element_text(size = 30),
          panel.grid.major.x = element_blank())
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_gamma_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_gamma_sc3-nl.pdf"), width=10, height = 7.78)


f.nonlinear.origin <- f.linear.origin <- f.origin <- f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
for(j in 1:p){
    f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
    f.linear.new[,j] <- bs.linear[, j] * theta.map[j+1]
    f.nonlinear.new[,j] <- bs.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map,nrow =20)[,j]
    f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}

func.linear.new <- func.nonlinear.new <- func.linear.origin <- func.nonlinear.origin <- func.new <- func.origin <- matrix(, nrow=n, ncol=0)
for (j in 1:p){
  func.origin <- cbind(func.origin, sort(f.origin[,j]))
  func.linear.origin <- cbind(func.linear.origin, sort(f.linear.origin[,j]))
  func.nonlinear.origin <- cbind(func.nonlinear.origin, sort(f.nonlinear.origin[,j]))
  func.new <- cbind(func.new, sort(f.new[,j]))
  func.linear.new <- cbind(func.linear.new, sort(f.linear.new[,j]))
  func.nonlinear.new <- cbind(func.nonlinear.new, sort(f.nonlinear.new[,j]))  
}


covariates <- gl(p, n, (p*n))
replicate <- gl(2, n, (p*n))
func.df <- data.frame(x = as.vector(apply(x.origin, 2, sort, method = "quick")),
                        origin=as.vector(func.origin),
                        origin.linear=as.vector(func.linear.origin),
                        origin.nonlinear=as.vector(func.nonlinear.origin), 
                        new=as.vector(func.new),
                        new.linear=as.vector(func.linear.new),
                        new.nonlinear=as.vector(func.nonlinear.new),
                        covariates=covariates, 
                        replicate=replicate)

equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=origin, colour = covariates, linetype = "true"), linewidth=2) + 
  geom_line(aes(y=new, colour = covariates, linetype = "MAP"), linewidth=2) + 
  ylab("") + xlab ("Smooth Functions") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  facet_grid(covariates ~ ., scale = "free_y") + 
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_linear_train_sc2-nl.pdf"), width=10.5, height = 15)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_linear_train_sc3-nl.pdf"), width=10.5, height = 15)

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=origin.linear, colour = covariates, linetype = "true"), linewidth=2) + 
  geom_line(aes(y=new.linear, colour = covariates, linetype = "MAP"), linewidth=2) + 
  ylab ("") + facet_grid(covariates ~ ., scale = "free_y") + xlab("Linear Components") + 
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,-20,0,-30),
        strip.text = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_linear_train_sc2-nl.pdf"), width=10.5, height = 15)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_linear_train_sc3-nl.pdf"), width=10.5, height = 15)


ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=origin.nonlinear, colour = covariates, linetype = "true"), linewidth=2) + 
  geom_line(aes(y=new.nonlinear, colour = covariates, linetype = "MAP"), linewidth=2) +
  ylab ("") + facet_grid(covariates ~ ., scale = "free_y") + xlab("Nonlinear Components") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=45),
        legend.margin=margin(0,0,0,-10),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_nonlinear_train_sc2-nl.pdf"), width=12.5, height = 15)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_nonlinear_train_sc3-nl.pdf"), width=12.5, height = 15)

old.alpha <- NULL
for(i in 1:n){
  old.alpha[i] <- exp(theta.map[1] + sum(f.new))
}

data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "trueAlp" = sort(alp.origin),
                            "mapAlp" = sort(alpha.map),
                            "optimAlp" = sort(old.alpha))
ggplot(data = data.scenario, aes(x = constant)) + 
  ylab(expression(alpha(x))) + xlab("") +
  geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
  geom_line(aes(y = mapAlp, col = "MAP Alpha"), linewidth = 2.5, linetype = 2) +
  labs(col = "") +
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
    scale_color_manual(values = c("#e0b430", "red"))+
    theme(text = element_text(size = 15),
            legend.position="bottom", legend.key.size = unit(1, 'cm'),
            axis.text = element_text(size = 20),
            legend.margin=margin(-15,-15,-15,-15),
            legend.box.margin=margin(-25,0,20,0))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_alpha_train_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_alpha_train_sc3-nl.pdf"), width=10, height = 7.78)

#### Test Set
f.nonlinear.origin <- f.linear.origin <- f.origin <- f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
true.alpha <- new.alpha <- NULL
for (j in 1:p){
  f.linear.origin[,j] <- xholder.linear[,j] * theta.origin[j+1]
  f.nonlinear.origin[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[, j]
  f.origin[,j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
  f.linear.new[,j] <- xholder.linear[,j] * theta.map[j+1]
  f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}
for(i in 1:n){
  new.alpha[i] <- exp(theta.map[1] + sum(f.new[i,]))
}

func.linear.new <- func.nonlinear.new <- func.linear.origin <- func.nonlinear.origin <- func.new <- func.origin <- matrix(, nrow=n, ncol=0)
for (j in 1:p){
  func.origin <- cbind(func.origin, f.origin[,j])
  func.linear.origin <- cbind(func.linear.origin, f.linear.origin[,j])
  func.nonlinear.origin <- cbind(func.nonlinear.origin, f.nonlinear.origin[,j])
  func.new <- cbind(func.new, f.new[,j])
  func.linear.new <- cbind(func.linear.new, f.linear.new[,j])
  func.nonlinear.new <- cbind(func.nonlinear.new, f.nonlinear.new[,j])
}

covariates <- gl(p, n, (p*n))
replicate <- gl(2, n, (p*n))
func.df <- data.frame(x = as.vector(apply(x.origin, 2, sort, method = "quick")),
                        origin=as.vector(func.origin),
                        origin.linear=as.vector(func.linear.origin),
                        origin.nonlinear=as.vector(func.nonlinear.origin), 
                        new=as.vector(func.new),
                        new.linear=as.vector(func.linear.new),
                        new.nonlinear=as.vector(func.nonlinear.new),
                        covariates=covariates, 
                        replicate=replicate)

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=origin, colour = covariates, linetype = "true"), linewidth = 2) + 
  geom_line(aes(y=new, colour = covariates, linetype = "MAP"), linewidth = 2) + 
  ylab ("") + facet_grid(covariates ~ .) + xlab("Smooth Functions") +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_smooth_test_sc2-nl.pdf"), width=10.5, height = 15)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_smooth_test_sc3-nl.pdf"), width=10.5, height = 15)

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) +  
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) +  
  geom_line(aes(y=origin.linear, colour = covariates, linetype = "true"), linewidth = 2) + 
  geom_line(aes(y=new.linear, colour = covariates, linetype = "MAP"), linewidth = 2) + 
  facet_grid(covariates ~ .) + xlab("Linear Component") + ylab ("") +
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,-20,0,-30),
        strip.text = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_linear_test_sc2-nl.pdf"), width=10, height = 15)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_linear_test_sc3-nl.pdf"), width=10, height = 15)

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=origin.nonlinear, colour = covariates, linetype = "true"), linewidth = 2) + 
  geom_line(aes(y=new.nonlinear, colour = covariates, linetype = "MAP"), linewidth=2) + ylab ("") + xlab("Nonlinear Component") + facet_grid(covariates ~ .) + 
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=45),
        legend.margin=margin(0,0,0,-10),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_nonlinear_test_sc2-nl.pdf"), width=12, height = 15)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_nonlinear_test_sc3-nl.pdf"), width=12, height = 15)


data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "trueAlp" = sort(alp.new),
                            "mapAlp" = sort(newalpha.map))

ggplot(data = data.scenario, aes(x = constant)) + 
  ylab(expression(alpha(x))) + xlab("") +
  geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
  geom_line(aes(y = mapAlp, col = "MAP Alpha"), linewidth = 2.5, linetype = 2) +
  labs(col = "") + #ylim(0, 5) +
  theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
  theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
  scale_color_manual(values = c("#e0b430", "red"))+
  theme(text = element_text(size = 15),
          legend.position="bottom", legend.key.size = unit(1, 'cm'),
          axis.text = element_text(size = 20),
          legend.margin=margin(-15,-15,-15,-15),
          legend.box.margin=margin(-25,0,20,0))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_alpha_test_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_alpha_test_sc3-nl.pdf"), width=10, height = 7.78)

# Randomized quantile residuals
r <- matrix(NA, nrow = n, ncol = 20)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 20
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y.origin[i], u, alpha = alpha.map[i]))
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
              color = "lightgrey", fill = "lightgrey",
              alpha = 0.4, linetype = "dashed") + 
  geom_line(aes(x = grid, y = trajhat), linetype = "dashed", linewidth = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.2) + 
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
  theme(text = element_text(size = 30)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_qqplot_sc2-nl.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_map_qqplot_sc3-nl.pdf"), width=10, height = 7.78)
