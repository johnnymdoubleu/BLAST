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


n <- 5000

beta.origin <- c(0.2, 0, 0.8, 0, 0,-0.1, 0, 0, 0,-0.4)
newp <- length(beta.origin)
p <- length(beta.origin)-1
X <- replicate(p, runif(n))
# X <- scale(X)

y <- alpha.origin <- alpha.new <- NULL
for(i in 1:n){
  y[i] <- rPareto(1, 1, alpha = exp(beta.origin[1] + sum(X[i, ] * beta.origin[2:newp])))
}

x.origin <- X[which(y>quantile(y,0.9)),]
u <- quantile(y,0.9)

y.origin <- as.matrix(y[y>quantile(y,0.9)])
n <- length(y.origin)
xholder <- replicate(p, seq(0,1, length.out = n))
for(i in 1:n){
  alpha.origin[i] <- exp(beta.origin[1] + sum(x.origin[i, ] * beta.origin[2:newp]))
  alpha.new[i] <- exp(beta.origin[1] + sum(xholder[i,] * beta.origin[2:newp]))
}

write("// Stan model for simple linear regression
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp; 
    real <lower=0> u; // large threshold value
    matrix[n, p] x; // train dataset
    matrix[n,p] xholder; // test dataset
    vector[n] y; // extreme response
}

parameters {
    vector[newp] beta; // linear predictor
    real <lower=0> lambda1; // lasso penalty
}

transformed parameters {
    array[n] real <lower=0> alpha; // tail index
    array[n] real <lower=0> newalpha; // new tail index
    
    for (i in 1:n){
        alpha[i] = exp(beta[1] + dot_product(xholder[i], beta[2:newp]));
        newalpha[i] = exp(beta[1] + dot_product(xholder[i], beta[2:newp]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 0.1, 0.1);
    target += normal_lpdf(beta[1] | 0, 10);
    target += newp * log(lambda1);
    for (j in 1:p){
        target += double_exponential_lpdf(beta[(j+1)] | 0, lambda1);
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
, "model_simulation_sc1.stan")

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, 
                    newp = newp, x = x.origin, xholder = xholder)

init.alpha <- list(list(beta = rep(0, (p+1)), lambda1 = 0.1),
                  list(beta = rep(0.01, (p+1)), lambda1 = 0.01),
                  list(beta = rep(0.05, (p+1)), lambda1 = 0.1))

fit1 <- stan(
    file = "model_simulation_sc1.stan",  # Stan program
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    chains = 3,             # number of Markov chains
    warmup = 1000,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = 4,              # number of cores (could use one per chain)
    refresh = 500             # no progress shown
)

posterior <- extract(fit1)
str(posterior)

plot(fit1, plotfun = "trace", pars = c("beta"), nrow = 3)
# ggsave(paste0("./simulation/results/",Sys.Date(),n,"_mcmc_theta_trace_sc1-wi.pdf"), width=10, height = 7.78)

plot(fit1, plotfun = "trace", pars = c("lambda1"), nrow = 2)

# ggsave(paste0("./simulation/results/",Sys.Date(),n,"_mcmc_lambda_sc1-wi.pdf"), width=10, height = 7.78)


beta.samples <- summary(fit1, par=c("beta"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary

beta.post.mean <- beta.samples[,1]
beta.q1 <- beta.samples[,4]
beta.q3 <- beta.samples[,6]

df.beta <- data.frame("seq" = seq(1, (p+1)),
                        "true" = beta.origin,
                        "m" = beta.post.mean,
                        "l" = beta.q1,
                        "u" = beta.q3)
df.beta$covariate <- factor(0:p)
df.beta$labels <- factor(0:p)
ggplot(df.beta, aes(x = covariate, y=m, color = covariate)) + ylab("") + xlab('') +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_point(size = 5) + geom_point(aes(y = true), color="red", size = 4) +
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, linewidth =1.2) + 
  scale_x_discrete(labels = c(expression(bold(beta[0])),
                              expression(bold(beta[1])),
                              expression(bold(beta[2])),
                              expression(bold(beta[3])),
                              expression(bold(beta[4])),
                              expression(bold(beta[5])),
                              expression(bold(beta[6])),
                              expression(bold(beta[7])),
                              expression(bold(beta[8])),
                              expression(bold(beta[9])),
                              expression(bold(beta[10])))) + 
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

# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_theta_sc2-wi.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_theta_sc3-wi.pdf"), width=10, height = 7.78)

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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_gamma_sc2-wi.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_gamma_sc3-wi.pdf"), width=10, height = 7.78)

g.nonlinear.q1 <- g.linear.q1 <- g.q1 <- g.nonlinear.q3 <- g.linear.q3 <- g.q3 <- g.nonlinear.new <- g.linear.new <- g.new <- matrix(, nrow = n, ncol=p)
alpha.smooth.q1 <- alpha.smooth.q3 <- alpha.smooth.new <- alpha.new <- NULL
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
  alpha.smooth.new[i] <- exp(theta.post.mean[1] + sum(g.new[i,]))
  alpha.smooth.q1[i] <- exp(theta.q1[1] + sum(g.q1[i,]))
  alpha.smooth.q3[i] <- exp(theta.q3[1] + sum(g.q3[i,]))
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
plot.smooth <- ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=true, colour = covariates, linetype = "True"), linewidth=2) + 
  geom_line(aes(y=post.mean, colour = covariates, linetype = "Posterior Mean"), linewidth=2) + 
  ylab("") + xlab ("Smooth Functions") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  facet_grid(covariates ~ ., scales = "free_y") +
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

plot.smooth + facetted_pos_scales(y = list(
    covariates == "1" ~ ylim(0, 1.4), 
    covariates == "2" ~ ylim(0, 1.3),
    covariates == "3" ~ ylim(0, 1.5),
    covariates == "4" ~ ylim(-0.68, 0.05),
    covariates == "5" ~ ylim(-0.2, 0.15),
    covariates == "6" ~ ylim(-0.1, 0.8),
    covariates == "7" ~ ylim(-0.5, 1),
    covariates == "8" ~ ylim(0, 1.4),
    covariates == "9" ~ ylim(-0.81, 0.1),
    covariates == "10" ~ ylim(-0.7, 0.1)
    ))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_smooth_sc2-wi.pdf"), width=10.5, height = 15)

# plot.smooth + facetted_pos_scales(y = list(
#     covariates == "1" ~ ylim(-0.01, 0.38),
#     covariates == "2" ~ ylim(-0.35, 0),
#     covariates == "3" ~ ylim(-0.35, 0),
#     covariates == "4" ~ ylim(-0.1, 0.05),
#     covariates == "5" ~ ylim(-0.35, 0)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_smooth_sc3-wi.pdf"), width=10.5, height = 15)
data.linear <- data.frame("x"=c(1:n),
                          "true" = as.vector(f.linear.new),
                          "post.mean" = as.vector(g.linear.new),
                          "q1" = as.vector(g.linear.q1),
                          "q3" = as.vector(g.linear.q3),
                          "covariates" = gl(p, n, (p*n)),
                          "replicate" = gl(2, n, (p*n)))
plot.linear <- ggplot(data.linear, aes(x=x, group=interaction(covariates, replicate))) + 
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

plot.linear + facetted_pos_scales(y = list(
    covariates == "1" ~ ylim(0, 1.4), 
    covariates == "2" ~ ylim(0, 0.28),
    covariates == "3" ~ ylim(0, 1.5),
    covariates == "4" ~ ylim(-0.27, 0.15),
    covariates == "5" ~ ylim(-0.2, 0.15),
    covariates == "6" ~ ylim(-0.1, 0.28),
    covariates == "7" ~ ylim(-0.33, 0.33),
    covariates == "8" ~ ylim(0, 1.4),
    covariates == "9" ~ ylim(-0.28, 0.1),
    covariates == "10" ~ ylim(-0.7, 0.1)
    ))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_linear_sc2-wi.pdf"), width=10.5, height = 15)

# plot.linear + facetted_pos_scales(y = list(
#     covariates == "1" ~ ylim(-0.01, 0.38),
#     covariates == "2" ~ ylim(-0.35, 0),
#     covariates == "3" ~ ylim(-0.35, 0),
#     covariates == "4" ~ ylim(-0.1, 0.05),
#     covariates == "5" ~ ylim(-0.35, 0)))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_linear_sc3-wi.pdf"), width=10.5, height = 15)

data.nonlinear <- data.frame("x"=c(1:n),
                          "true" = as.vector(f.nonlinear.new),
                          "post.mean" = as.vector(g.nonlinear.new),
                          "q1" = as.vector(g.nonlinear.q1),
                          "q3" = as.vector(g.nonlinear.q3),
                          "covariates" = gl(p, n, (p*n)),
                          "replicate" = gl(2, n, (p*n)))
plot.nonlinear <- ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
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

plot.nonlinear + facetted_pos_scales(y = list(
    covariates == "1" ~ ylim(0, 1.4), 
    covariates == "2" ~ ylim(-0.02, 0.6),
    covariates == "3" ~ ylim(0, 1.5),
    covariates == "4" ~ ylim(-0.33, 0.33),
    covariates == "5" ~ ylim(-0.33, 0.33),
    covariates == "6" ~ ylim(-0.4, 0.01),
    covariates == "7" ~ ylim(-0.2, 1.4),
    covariates == "8" ~ ylim(0, 1.4),
    covariates == "9" ~ ylim(-0.33, 0.33),
    covariates == "10" ~ ylim(-0.7, 0.1)
    ))        
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_nonlinear_sc2-wi.pdf"), width=12.5, height = 15)

# plot.nonlinear + facetted_pos_scales(y = list(
#     covariates == "1" ~ ylim(-0.03, 0.05), 
#     covariates == "2" ~ ylim(-0.09, 0),
#     covariates == "3" ~ ylim(-0.09, 0),
#     covariates == "4" ~ ylim(-0.03, 0.03),
#     covariates == "5" ~ ylim(-0.03, 0.03)))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_nonlinear_sc3-wi.pdf"), width=12.5, height = 15)
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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_alpha_train_sc2-wi.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_alpha_train_sc3-wi.pdf"), width=10, height = 7.78)

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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_alpha_test_sc2-wi.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_alpha_test_sc3-wi.pdf"), width=10, height = 7.78)
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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_qqplot_sc2-wi.pdf"), width=10, height = 7.78)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_mcmc_qqplot_sc3-wi.pdf"), width=10, height = 7.78)