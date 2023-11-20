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
# set.seed(9)

total.iter <- 10

n <- 5000
threshold <- 0.9
beta.origin <- c(0.2, 0, 0.8, 0, 0,-0.1, 0, 0, 0,-0.4)
newp <- length(beta.origin)
p <- length(beta.origin)-1

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
        alpha[i] = exp(beta[1] + dot_product(x[i], beta[2:newp]));
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

beta.container <- as.data.frame(matrix(, nrow = newp, ncol= total.iter))
linear.container <- nonlinear.container <- f.container <- lapply(1:simul.no, data.frame)
alpha.container <- as.data.frame(matrix(, nrow=n, ncol = total.iter))
for(iter in 1:total.iter){
    X <- replicate(p, runif(n))
    # X <- scale(X)

    y <- alpha.origin <- alpha.new <- NULL
    for(i in 1:n){
        y[i] <- rPareto(1, 1, alpha = exp(beta.origin[1] + sum(X[i, ] * beta.origin[2:newp])))
    }

    x.origin <- X[which(y>quantile(y, threshold)),]
    u <- quantile(y, threshold)

    y.origin <- as.matrix(y[y>quantile(y, threshold)])
    n <- length(y.origin)
    xholder <- replicate(p, seq(0, 1, length.out = n))
    for(i in 1:n){
    alpha.origin[i] <- exp(beta.origin[1] + sum(x.origin[i, ] * beta.origin[2:newp]))
    alpha.new[i] <- exp(beta.origin[1] + sum(xholder[i,] * beta.origin[2:newp]))
    }
    data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, 
                    newp = newp, x = x.origin, xholder = xholder)

    init.alpha <- list(list(beta = rep(0, (p+1)), lambda1 = 0.01),
                    list(beta = rep(0.01, (p+1)), lambda1 = 0.01),
                    list(beta = rep(-0.05, (p+1)), lambda1 = 0.01))

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

    beta.samples <- summary(fit1, par=c("beta"), probs = c(0.05,0.5, 0.95))$summary
    lambda.samples <- summary(fit1, par=c("lambda1"), probs = c(0.05,0.5, 0.95))$summary
    alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
    newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary

    alpha.container[,iter] <- sort(newalpha.samples[,1])
    beta.container[,iter] <- beta.samples[,1]
}

alpha.container$x <- seq(0,1, length.out = n)
alpha.container$true <- alpha.new
alpha.container <- cbind(alpha.container, t(apply(alpha.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
colnames(alpha.container)[(dim(alpha.container)[2]-2):(dim(alpha.container)[2])] <- c("q1","q2","q3")
alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
alpha.container <- as.data.frame(alpha.container)
plt <- ggplot(data = alpha.container, aes(x = x)) + ylab(expression(alpha(x))) + xlab("")
# for(i in 1:total.iter){
#   # print(.data[[names(data.scenario)[i]]])
#   plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.4,linewidth = 0.7)
# }

print(plt + geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) + 
        geom_line(aes(y=true, col = "True"), linewidth = 1.8) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.8, linetype = 2) +
        # geom_line(aes(y = post.check, col=paste0("Simulated Alpha: ",n,"/",psi,"/",threshold)), linewidth = 1.5) +
        theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        labs(col = "") + #ylim(0, 150) +
        scale_color_manual(values = c("#e0b430", "red"))+
        theme(text = element_text(size = 15),
                legend.position="bottom", legend.key.size = unit(1, 'cm'),
                axis.text = element_text(size = 20),
                legend.margin=margin(-15,-15,-15,-15),
                legend.box.margin=margin(-25,0,20,0)))


resg <- gather(beta.container,
               key = "group",
               names(beta.container),
               value = "values")
resg$group1 <- factor(rep(1:newp, total.iter))
somelines <- data.frame(value=c(as.vector(beta.origin)),boxplot.nr=c(1:(newp)))
ggplot(resg, aes(group=group1, x = group1, y = values)) + ylim(-1,1) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) +
  geom_boxplot() + #coord_cartesian(ylim=c(-1,1))+
  geom_segment(data=somelines,aes(x=boxplot.nr-0.5, xend=boxplot.nr+0.5, 
                                  y=value,yend=value),inherit.aes=FALSE,color="red",linewidth=1.5)+
  # facet_wrap( ~ group2, labeller = label_parsed) +
  labs(x = "", y = "") + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust=0.8),
        axis.text = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size=30))

# X <- replicate(p, runif(n))
# # X <- scale(X)

# y <- alpha.origin <- alpha.new <- NULL
# for(i in 1:n){
#   y[i] <- rPareto(1, 1, alpha = exp(beta.origin[1] + sum(X[i, ] * beta.origin[2:newp])))
# }

# x.origin <- X[which(y>quantile(y, threshold)),]
# u <- quantile(y, threshold)

# y.origin <- as.matrix(y[y>quantile(y, threshold)])
# n <- length(y.origin)
# xholder <- replicate(p, seq(0, 1, length.out = n))
# for(i in 1:n){
#   alpha.origin[i] <- exp(beta.origin[1] + sum(x.origin[i, ] * beta.origin[2:newp]))
#   alpha.new[i] <- exp(beta.origin[1] + sum(xholder[i,] * beta.origin[2:newp]))
# }

# write("// Stan model for simple linear regression
# data {
#     int <lower=1> n; // Sample size
#     int <lower=1> p; // regression coefficient size
#     int <lower=1> newp; 
#     real <lower=0> u; // large threshold value
#     matrix[n, p] x; // train dataset
#     matrix[n,p] xholder; // test dataset
#     vector[n] y; // extreme response
# }

# parameters {
#     vector[newp] beta; // linear predictor
#     real <lower=0> lambda1; // lasso penalty
# }

# transformed parameters {
#     array[n] real <lower=0> alpha; // tail index
#     array[n] real <lower=0> newalpha; // new tail index
    
#     for (i in 1:n){
#         alpha[i] = exp(beta[1] + dot_product(x[i], beta[2:newp]));
#         newalpha[i] = exp(beta[1] + dot_product(xholder[i], beta[2:newp]));
#     };
# }

# model {
#     // likelihood
#     for (i in 1:n){
#         target += pareto_lpdf(y[i] | u, alpha[i]);
#     }
#     target += gamma_lpdf(lambda1 | 0.1, 0.1);
#     target += normal_lpdf(beta[1] | 0, 100);
#     target += newp * log(lambda1);
#     for (j in 1:p){
#         target += double_exponential_lpdf(beta[(j+1)] | 0, lambda1);
#     }
# }
# generated quantities {
#     // Used in Posterior predictive check
#     vector[n] log_lik;
#     array[n] real y_rep = pareto_rng(rep_vector(u, n), alpha);
#     for (i in 1:n) {
#         log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
#     }
# }
# "
# , "model_simulation_sc1.stan")

# data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, 
#                     newp = newp, x = x.origin, xholder = xholder)

# init.alpha <- list(list(beta = rep(0, (p+1)), lambda1 = 0.01),
#                   list(beta = rep(0.01, (p+1)), lambda1 = 0.01),
#                   list(beta = rep(-0.05, (p+1)), lambda1 = 0.01))

# fit1 <- stan(
#     file = "model_simulation_sc1.stan",  # Stan program
#     data = data.stan,    # named list of data
#     init = init.alpha,      # initial value
#     chains = 3,             # number of Markov chains
#     warmup = 1000,          # number of warmup iterations per chain
#     iter = 2000,            # total number of iterations per chain
#     cores = 4,              # number of cores (could use one per chain)
#     refresh = 500             # no progress shown
# )

# posterior <- extract(fit1)
# str(posterior)

# plot(fit1, plotfun = "trace", pars = c("beta"), nrow = 3)
# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_theta_trace_sc1-wi.pdf"), width=10, height = 7.78)

# plot(fit1, plotfun = "trace", pars = c("lambda1"), nrow = 2)
# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_lambda_sc1-wi.pdf"), width=10, height = 7.78)


# beta.samples <- summary(fit1, par=c("beta"), probs = c(0.05,0.5, 0.95))$summary
# lambda.samples <- summary(fit1, par=c("lambda1"), probs = c(0.05,0.5, 0.95))$summary
# alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
# newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary

# beta.post.mean <- beta.samples[,1]
# beta.q1 <- beta.samples[,4]
# beta.q3 <- beta.samples[,6]

# df.beta <- data.frame("seq" = seq(1, (p+1)),
#                         "true" = beta.origin,
#                         "m" = beta.post.mean,
#                         "l" = beta.q1,
#                         "u" = beta.q3)
# df.beta$covariate <- factor(0:p)
# df.beta$labels <- factor(0:p)
# ggplot(df.beta, aes(x = covariate, y=m, color = covariate)) + ylab("") + xlab('') +
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_point(size = 5) + geom_point(aes(y = true), color="red", size = 4) +
#   geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, linewidth =1.2) + 
#   scale_x_discrete(labels = c(expression(bold(beta[0])),
#                               expression(bold(beta[1])),
#                               expression(bold(beta[2])),
#                               expression(bold(beta[3])),
#                               expression(bold(beta[4])),
#                               expression(bold(beta[5])),
#                               expression(bold(beta[6])),
#                               expression(bold(beta[7])),
#                               expression(bold(beta[8])),
#                               expression(bold(beta[9])),
#                               expression(bold(beta[10])))) + 
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#           legend.text.align = 0,
#           legend.title = element_blank(),
#           legend.text = element_text(size=25),
#           legend.margin=margin(0,0,0,-10),
#           legend.box.margin=margin(-10,0,-10,0),
#           plot.margin = margin(0,0,0,-20),
#           axis.text.x = element_text(hjust=0.35),
#           axis.text = element_text(size = 28))

# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_beta_sc1-wi.pdf"), width=10, height = 7.78)

# data.scenario <- data.frame("x" = c(1:n),
#                             "constant" = newx,
#                             "true" = sort(alpha.origin),
#                             "post.mean" = sort(alpha.samples[,1]),
#                             "post.median" = sort(alpha.samples[,5]),
#                             "q1" = sort(alpha.samples[,4]),
#                             "q3" = sort(alpha.samples[,6]))

# ggplot(data.scenario, aes(x=x)) + 
#   ylab(expression(alpha(x))) + xlab(expression(x)) + labs(col = "") + 
#   geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
#   geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
#   geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=2, linetype = 2) +
#   # ylim(0, (max(data.scenario$post.mean)+10)) +
#   geom_line(aes(y=post.median, col = "Posterior Median"), linetype=2, linewidth=2) +
#   # geom_line(aes(y=chain2, col = "Chian 2"), linetype=3) +
#   # facet_grid(covariates ~ .) + 
#   # scale_y_continuous(breaks=c(0)) + 
#   scale_color_manual(values = c("blue","#e0b430","red"))+
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 30),
#         legend.position="top", 
#         legend.key.size = unit(1, 'cm'),
#         legend.text = element_text(size=18),
#         plot.margin = margin(0,0,0,-10),
#         strip.text = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size=33),
#         axis.title.x = element_text(size = 35))
# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_train_sc1-wi.pdf"), width=10, height = 7.78)

# data.scenario <- data.frame("x" = c(1:n),
#                             "constant" = newx,
#                             "true" = sort(alpha.new),
#                             "post.mean" = sort(newalpha.samples[,1]),
#                             "post.median" = sort(newalpha.samples[,5]),
#                             "q1" = sort(newalpha.samples[,4]),
#                             "q3" = sort(newalpha.samples[,6]))
#                             # "post.mean" = sort(alpha.smooth.new),
#                             # "post.median" = sort(newalpha.samples[,5]),
#                             # "q1" = sort(alpha.smooth.q1),
#                             # "q3" = sort(alpha.smooth.q3))                            

# ggplot(data.scenario, aes(x=x)) + 
#   ylab(expression(alpha(x))) + xlab(expression(x)) + labs(col = "") + 
#   geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
#   geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
#   geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=2, linetype = 2) + 
#   # ylim(0, (max(data.scenario$q3)+10)) +
#   geom_line(aes(y=post.median, col = "Posterior Median"), linetype=2, linewidth=2) +
#   # geom_line(aes(y=chain2, col = "Chian 2"), linetype=3) +
#   # facet_grid(covariates ~ .) + 
#   # scale_y_continuous(breaks=c(0)) + 
#   scale_color_manual(values = c("blue","#e0b430","red"))+
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 30),
#         legend.position="top", 
#         legend.key.size = unit(1, 'cm'),
#         legend.text = element_text(size=18),
#         plot.margin = margin(0,0,0,-10),
#         strip.text = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size=33),
#         axis.title.x = element_text(size = 35))
# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_test_sc1-wi.pdf"), width=10, height = 7.78)

# mcmc.alpha <- posterior$alpha
# len <- dim(mcmc.alpha)[1]
# r <- matrix(, nrow = n, ncol = 30)
# # beta <- as.matrix(mcmc[[1]])[, 1:7] 
# T <- 30
# for(i in 1:n){
#   for(t in 1:T){
#     r[i, t] <- qnorm(pPareto(y.origin[i], u, alpha = mcmc.alpha[round(runif(1,1,len)),i]))
#   }
# }
# lgrid <- n
# grid <- qnorm(ppoints(lgrid))
# traj <- matrix(NA, nrow = T, ncol = lgrid)
# for (t in 1:T){
#   traj[t, ] <- quantile(r[, t], ppoints(lgrid), type = 2)
# }
# l.band <- apply(traj, 2, quantile, prob = 0.025)
# trajhat <- apply(traj, 2, quantile, prob = 0.5)
# u.band <- apply(traj, 2, quantile, prob = 0.975)

# ggplot(data = data.frame(grid = grid, l.band = l.band, trajhat = trajhat, 
#                          u.band = u.band)) + 
#   geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
#               alpha = 0.4, linetype = "dashed") + 
#   geom_line(aes(x = grid, y = trajhat), linetype = "dashed", linewidth = 1.2) + 
#   geom_abline(intercept = 0, slope = 1, linewidth = 1.2) + 
#   labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
#   theme_minimal(base_size = 20) +
#   theme(text = element_text(size = 20)) + 
#   coord_fixed(xlim = c(-3, 3),  
#               ylim = c(-3, 3))
# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_qqplot_sc1-wi.pdf"), width=10, height = 7.78)
