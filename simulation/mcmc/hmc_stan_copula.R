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
# set.seed(22)
set.seed(3)

n <- 2500
psi <- 20
threshold <- 0.9
p <- 5
no.theta <- 1
simul.no <- 50

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)

## Function to generate Gaussian copula
C <- matrix(c(1, 0.3, 0.5, 0.3, 0.3,
            0.3, 1, 0.95, 0.4, 0.4,
            0.5, 0.95, 1, 0.5, 0.1,
            0.3, 0.4, 0.5 , 1, 0.5,
            0.3, 0.4, 0.5, 0.5, 1), nrow = p)    
## Generate sample
x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
# x.origin <- scale(x.origin)
# pairs(x.origin, diag.panel = function(x){
#           h <- hist(x, plot = FALSE)
#           rect(head(h$breaks, -1), 0, tail(h$breaks, -1), h$counts/max(h$counts))})



for(i in 1:p){
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    # bs.linear <- cbind(bs.linear, tps[,1:no.theta])
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
}

gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
    for (ps in 1:psi){
        if(j %in% c(1,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
        else if(j==2){
            if(ps <= (psi/2)){gamma.origin[ps, j] <- -0.5}
            else{gamma.origin[ps, j] <- -0.5}
        }
        else {
            if(ps <= (psi/2)){gamma.origin[ps, j] <- 1.7}
            else{gamma.origin[ps, j] <- 1.7}
        }
    }
}
theta.origin <- -1.6

f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp(theta.origin + sum(f.nonlinear.origin[i,]))
    y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i])
}

u <- quantile(y.origin, threshold)
x.origin <- x.origin[which(y.origin>u),]
# x.bs <- x.origin
# x.origin <- scale(x.origin)
y.origin <- y.origin[y.origin > u]
n <- length(y.origin)

# corrplot.mixed(cor(x.origin),
#                 upper = "circle",
#                 lower = "number",
#                 addgrid.col = "black")

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
newx <- seq(0, 1, length.out=n)
xholder <- bs.x <- matrix(, nrow = n, ncol = p)
for(i in 1:p){
    # xholder[,i] <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = n)  
    # test.knot <- seq(min(xholder[,i]), max(xholder[,i]), length.out = psi)
    xholder[,i] <- seq(0, 1, length.out = n)  
    test.knot <- seq(0, 1, length.out = psi)
    splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
    xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)]) 
}

f.nonlinear.new <- f.nonlinear.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
}

true.alpha <- alp.new <- alp.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp(theta.origin + sum(f.nonlinear.origin[i,]))
    alp.new[i] <- exp(theta.origin + sum(f.nonlinear.new[i,]))
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
    real theta; // intercept term
    array[p] vector[psi] gamma; // splines coefficient
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
        alpha[i] = exp(theta + sum(gsmooth[i,]));
        newalpha[i] = exp(theta + sum(newgsmooth[i,]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    };
    target += gamma_lpdf(lambda | 0.1, 0.1);
    target += normal_lpdf(theta | 0, 10);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    target += (p * psi * log(lambda));
    for (j in 1:p){
        target += gamma_lpdf(tau[j] | atau, lambda/sqrt(2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * sqrt(tau[j]) * sqrt(sigma));
    };
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
                        theta = 0, tau = rep(0.1, p), sigma = 0.1, 
                        lambda = 0.1),
                  list(gamma = array(rep(0.02, (psi*p)), dim=c(psi, p)),
                        theta = -1, tau = rep(0.01, p), sigma = 0.001,
                        lambda = 0.1),
                  list(gamma = array(rep(0.01, (psi*p)), dim=c(psi, p)),
                        theta = -2, tau = rep(0.01, p), sigma = 0.01,
                        lambda = 1))

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
# str(posterior)

plot(fit1, plotfun = "trace", pars = c("theta", "lambda"), nrow = 2)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_lambda_sc3-wi.pdf"), width=10, height = 7.78)


theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda"), probs = c(0.05,0.5, 0.95))$summary
# alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary

newgsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary


theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q2 <- theta.samples[,5]
theta.q3 <- theta.samples[,6]
gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q2 <- gamma.samples[,5]
gamma.q3 <- gamma.samples[,6]

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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_gamma_sc3-wi.pdf"), width=10, height = 7.78)

g.nonlinear.q1 <- g.nonlinear.q2  <- g.nonlinear.q3 <- g.nonlinear.new <- matrix(, nrow = n, ncol=p)
alpha.smooth.q1 <- alpha.smooth.q2 <- alpha.smooth.q3 <- alpha.smooth.new <- alpha.new <- NULL
for (j in 1:p){
  g.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.post.mean, nrow=psi)[,j]
  g.nonlinear.q1[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q1, nrow=psi)[,j]
  g.nonlinear.q2[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q2, nrow=psi)[,j] 
  g.nonlinear.q3[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q3, nrow=psi)[,j] 
}

for(i in 1:n){
  alpha.smooth.new[i] <- exp(theta.post.mean + sum(g.nonlinear.new[i,]))
  alpha.smooth.q1[i] <- exp(theta.q1 + sum(g.nonlinear.q1[i,]))
  alpha.smooth.q2[i] <- exp(theta.q2 + sum(g.nonlinear.q2[i,]))
  alpha.smooth.q3[i] <- exp(theta.q3 + sum(g.nonlinear.q3[i,]))
}

equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

g.nonlinear.new <- as.vector(matrix(newgsmooth.samples[,1], nrow = n, byrow=TRUE))
g.nonlinear.q1 <- as.vector(matrix(newgsmooth.samples[,4], nrow = n, byrow=TRUE))
g.nonlinear.q2 <- as.vector(matrix(newgsmooth.samples[,5], nrow = n, byrow=TRUE))
g.nonlinear.q3 <- as.vector(matrix(newgsmooth.samples[,6], nrow = n, byrow=TRUE))


data.nonlinear <- data.frame("x"=newx,
                          "true" = as.vector(f.nonlinear.new),
                          "post.mean" = as.vector(g.nonlinear.new),
                          "q1" = as.vector(g.nonlinear.q1),
                          "q2" = as.vector(g.nonlinear.q2),
                          "q3" = as.vector(g.nonlinear.q3),
                          "covariates" = gl(p, n, (p*n), labels = c("g[1](x[1])", "g[2](x[2])", "g[3](x[3])", "g[4](x[4])", "g[5](x[5])")),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))
# plot.nonlinear <- ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_line(aes(y=true, colour = "true"), linewidth=2) + 
#   geom_line(aes(y=q2, colour = covariates), linewidth=1.5) + 
#   xlab("Nonlinear Components") + ylab("") +
#   facet_grid(covariates ~ ., scales = "free_y") + 
#   scale_color_manual(values=c("pink", "#F5C710", "#61D04F", "#2297E6", "purple", "red")) +
#   scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 15),
#         axis.ticks = element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(size=45),
#         legend.margin=margin(0,0,0,-10),
#         legend.box.margin=margin(-10,0,-10,0),
#         plot.margin = margin(0,0,0,-20),
#         strip.text = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size=33),
#         axis.title.x = element_text(size = 35))

# plot.nonlinear
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_nonlinear_sc3-nl.pdf"), width=12.5, height = 15)
ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=true, colour = "True"), linewidth=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab(expression(x[j])) +
  facet_wrap(covariates ~ ., scales = "free_x", nrow=5,
              labeller = label_parsed, strip.position = "left") + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
#   scale_linetype_manual("functions",values=c("MCMC"=1,"True"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
  # coord_capped_cart(bottom='both', left='both', xlim=c(0, 1)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        # axis.ticks = element_blank(),
        legend.position="top",
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.margin=margin(t = 1, unit='cm'),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text.y = element_text(size = 25, colour = "black", angle = 0, face = "bold.italic"),
        # strip.text.x.bottom = element_text(size=18, face = "bold.italic"),
        strip.placement = "outside",
        # strip.text.x = ,
        # strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"),
        # axis.text.x = element_blank(),
        axis.title.x = element_text(size = 35),
        axis.text = element_text(size=18))
#   facetted_pos_scales(y = list(
#         covariates == "1" ~ ylim(-0.03, 0.05), 
#         covariates == "2" ~ ylim(-0.09, 0.01),
#         covariates == "3" ~ ylim(-0.09, 0.01),
#         covariates == "4" ~ ylim(-0.03, 0.03),
#         covariates == "5" ~ ylim(-0.03, 0.03)))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_nonlinear(CI)_sc3-nl.pdf"), width=12.5, height = 16)

data.nonlinear <- data.frame("x"=as.vector(xholder),
                          "true" = as.vector(f.nonlinear.new),
                          "post.mean" = g.nonlinear.new,
                          "q1" = g.nonlinear.q1,
                          "q2" = g.nonlinear.q2,
                          "q3" = g.nonlinear.q3,
                          "covariates" = gl(p, n, (p*n), labels = c("f[1](x[1])", "g[2](x[2])", "g[3](x[3])", "g[4](x[4])", "g[5](x[5])")),
                          "replicate" = gl(2, n, (p*n)))
data.charge <- data.nonlinear[(1*n):(2*n),]
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_nonlinear_sc3-nl.pdf"), width=12.5, height = 15)
ggplot(data.charge, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=true, colour = "True"), linewidth=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  xlab(expression(x[1])) + ylab(expression(g[1](x[1]))) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) + ylim(-0.5,1) +
#   scale_linetype_manual("functions",values=c("MCMC"=1,"True"=1)) +
  # scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        # axis.ticks = element_blank(),
        legend.position="none",
        # legend.position="top",
        # legend.title = element_blank(),
        # legend.text = element_text(size=20),
        # legend.margin=margin(t = 1, unit='cm'),
        # legend.box.margin=margin(-10,0,-10,0),
        # plot.margin = margin(0,0,0,-20),
        strip.text.y = element_text(size = 18, colour = "black", angle = 0, face = "bold.italic"),
        # strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"),
        # axis.text.x = element_blank(),
        axis.title.x = element_text(size = 35),
        axis.text.y = element_text(size=25))


data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "true" = sort(alp.new),
                            "post.mean" = sort(newalpha.samples[,1]),
                            "post.median" = sort(newalpha.samples[,5]),
                            "q1" = sort(newalpha.samples[,4]),
                            "q3" = sort(newalpha.samples[,6]))
                            # "post.mean" = sort(alpha.smooth.new),
                            # "post.median" = sort(alpha.smooth.q2),
                            # "q1" = sort(alpha.smooth.q1),
                            # "q3" = sort(alpha.smooth.q3))                            
# data.scenario <- data.frame("x" = c(1:n),
#                             "constant" = newx,
#                             "true" = alp.new,
#                             "post.mean" = newalpha.samples[,1],
#                             "post.median" = newalpha.samples[,5],
#                             "q1" = newalpha.samples[,4],
#                             "q3" = newalpha.samples[,6])

ggplot(data.scenario, aes(x=constant)) + 
  ylab(expression(alpha(x))) + xlab(expression(c)) + labs(col = "") + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2) + 
  geom_line(aes(y = true, col = "True"), linewidth = 2) +
  ylim(0, 1.5) +
#   geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=2, linetype = 2) + 
#   ylim(0, (max(data.scenario$q3)+10)) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  # geom_line(aes(y=chain2, col = "Chian 2"), linetype=3) +
  # facet_grid(covariates ~ .) + 
  # scale_y_continuous(breaks=c(0)) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue","red"))+
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position="top", 
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=20),
        plot.margin = margin(0,0,0,-1),
        strip.text = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        # axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))

ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_test_sc3-nl.pdf"), width=10, height = 7.78)


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

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_qqplot_sc3-wi.pdf"), width=10, height = 7.78)

# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))
