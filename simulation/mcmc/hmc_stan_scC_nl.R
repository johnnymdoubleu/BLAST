library(npreg)
library(scales)
library(MASS)
suppressMessages(library(tidyverse))
library(readxl)
library(corrplot)
library(rstan)
library(ggmcmc)
library(MCMCvis)
library(cmdstanr)
library(ggh4x)
library(LaplacesDemon)
library(extraDistr)

# library(simstudy)
# library(ggplotify)

#Scenario 1
# set.seed(22)
set.seed(36)

n <- 5000
psi <- 20
threshold <- 0.9
p <- 1
no.theta <- 1
simul.no <- 50

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
# Function to generate Gaussian copula
C <- diag(p)                
## Generate sample
x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
# x.origin <- scale(x.origin)

for(i in 1:p){
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
}

gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
    for (ps in 1:psi){
        if(j %in% c(2,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
        else if(j==7){
            if(ps <= (psi/2)){gamma.origin[ps, j] <- 1}
            else{gamma.origin[ps, j] <- 1}
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
    y.origin[i] <- rt(1, df = alp.origin[i])
    # y.origin[i] <- rht(1, nu = alp.origin[i])
    # y.origin[i] <- rhalft(1, nu = alp.origin[i], scale=1)
    # y.origin[i] <- rst(1, 0, 1, nu = alp.origin[i])
}

u <- quantile(y.origin, threshold)
x.origin <- x.origin[which(y.origin>u),]
# x.origin <- scale(x.origin)
y.origin <- y.origin[y.origin > u]
n <- length(y.origin)

alp.origin <- alp.origin[which(y.origin>u)]
bs.nonlinear <- bs.nonlinear[which(y.origin>u),]
# pairs(x.origin, diag.panel = function(x){
#           h <- hist(x, plot = FALSE)
#           rect(head(h$breaks, -1), 0, tail(h$breaks, -1), h$counts/max(h$counts))})

# corrplot.mixed(cor(x.origin),
#                 upper = "circle",
#                 lower = "number",
#                 addgrid.col = "black")

xholder.nonlinear <- xholder.linear <- matrix(,nrow=n, ncol=0)
# bs.nonlinear <- matrix(,nrow=n, ncol=0)
newx <- seq(0, 1, length.out=n)
xholder <- bs.x <- matrix(, nrow = n, ncol = p)
for(i in 1:p){
    # xholder[,i] <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = n)  
    # test.knot <- seq(min(xholder[,i]), max(xholder[,i]), length.out = psi)
    xholder[,i] <- seq(0, 1, length.out = n)  
    test.knot <- seq(0, 1, length.out = psi)
    splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
    xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
    # knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
    # tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    # bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
}

f.nonlinear.new <- f.nonlinear.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
}

alp.new <- NULL
# alp.origin <- NULL
for(i in 1:n){
    # alp.origin[i] <- exp(theta.origin + sum(f.nonlinear.origin[i,]))
    alp.new[i] <- exp(theta.origin + sum(f.nonlinear.new[i,]))
}


write("// Stan model for simple linear regression
functions{
    real halft_lpdf(real y, real c){
        // Burr distribution log pdf
        return ((c+1)/2) * log(1+((y^2)/c));
    }

    real burr_rng(real c){
        return ((1-uniform_rng(0,1))^(-1)-1)^(1/c);
    }
}

data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp;
    real <lower=0> u;
    int <lower=1> psi; // splines coefficient size
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    array[n] real <lower=0> y; // extreme responses
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
        target += student_t_lpdf(y[i] | alpha[i], 0, 1); // student_t_lpdf(y[i] | alpha[i], 0, 1) halft_lpdf(y[i] | alpha[i]) pareto_lpdf(y[i]|u, alpha[i])
        target += -1*log(1-student_t_cdf(u, alpha[i], 0, 1));
    };
    target += gamma_lpdf(lambda | 0.1, 1000);
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
    array[n] real y_rep = student_t_rng(alpha, rep_vector(0, n),rep_vector(1, n));
    for (i in 1:n) {
        log_lik[i] = student_t_lpdf(y[i] | alpha[i], 0, 1);
    }
}
"
, "model_simulation_sc3.stan")

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), newp = (p+1),
                    bsNonlinear = bs.nonlinear,
                    xholderNonlinear = xholder.nonlinear)

# init.alpha <- list(list(gamma = gamma.map, theta = theta.map, tau = tau.map, sigma = sigma.map, lambda = lambda.map),
# list(gamma = gamma.map, theta = theta.map, tau = tau.map, sigma = sigma.map, lambda = lambda.map),
# list(gamma = gamma.map, theta = theta.map, tau = tau.map, sigma = sigma.map, lambda = lambda.map))

init.alpha <- list(list(gamma = array(rep(0, (psi*p)), dim=c(psi, p)),
                        theta = 0, tau = rep(0.1, p), sigma = 0.1, 
                        lambda = 0.1),
                  list(gamma = array(rep(0.5, (psi*p)), dim=c(psi, p)),
                        theta = 0.5, tau = rep(0.1, p), sigma = 0.1,
                        lambda = 0.1),
                  list(gamma = array(rep(1, (psi*p)), dim=c(psi, p)),
                        theta = -1, tau = rep(0.1, p), sigma = 0.1,
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

plot(fit1, plotfun = "trace", pars = c("theta", "lambda"), nrow = 2)
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_lambda_sc1-nl.pdf"), width=10, height = 7.78)


theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda"), probs = c(0.05,0.5, 0.95))$summary

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
                  "pm" = as.vector(gamma.post.mean),
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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_gamma_sc1-nk.pdf"), width=10, height = 7.78)

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
                          "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]", "g[6]", "g[7]", "g[8]", "g[9]", "g[10]")),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]", "x[6]", "x[7]", "x[8]", "x[9]", "x[10]")))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_nonlinear_sc1-nl.pdf"), width=12.5, height = 15)

ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=true, colour = "True"), linewidth=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab("") +
  facet_wrap(covariates ~ ., scales = "free_x", nrow = 5,
              labeller = label_parsed, strip.position = "left") + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) + ylim(-0.5, 1) +
#   scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) +
#   scale_y_continuous(breaks=c(-0.1, 0.2, 0.6)) + 
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

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_nonlinear(CI)_sc1-nl.pdf"), width=20, height = 16)

data.nonlinear <- data.frame("x"=as.vector(xholder),
                          "true" = as.vector(f.nonlinear.new),
                          "post.mean" = g.nonlinear.new,
                          "q1" = g.nonlinear.q1,
                          "q2" = g.nonlinear.q2,
                          "q3" = g.nonlinear.q3,
                          "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]", "g[6]", "g[7]", "g[8]", "g[9]", "g[10]")),
                          "replicate" = gl(2, n, (p*n)))
data.charge <- data.nonlinear[(1*n):(2*n),]

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
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position="none",
        strip.text.y = element_text(size = 18, colour = "black", angle = 0, face = "bold.italic"),
        axis.title.x = element_text(size = 35),
        axis.text.y = element_text(size=25))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_nonlinear_X_sc1-nl.pdf"), width=12.5, height = 15)

data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "true" = sort(alp.new),
                            "post.mean" = sort(newalpha.samples[,1]),
                            "post.median" = sort(newalpha.samples[,5]),
                            "q1" = sort(newalpha.samples[,4]),
                            "q3" = sort(newalpha.samples[,6]))

ggplot(data.scenario, aes(x=constant)) + 
  ylab(expression(alpha(bold(c1)))) + xlab(expression(c)) + labs(col = "") + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  geom_line(aes(y = true, col = "True"), linewidth = 2) +
  ylim(0, 1.5) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue","red")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position="none", 
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=20),
        plot.margin = margin(0,0,0,-1),
        strip.text = element_blank(),
        axis.title.x = element_text(size = 35))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_test_sc1-nl.pdf"), width=10, height = 7.78)


mcmc.alpha <- posterior$alpha

len <- dim(mcmc.alpha)[1]
r <- matrix(, nrow = n, ncol = 30)
T <- 30
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pst(y.origin[i], 0, 1, mcmc.alpha[round(runif(1,1,len)),i], lower.tail = FALSE))
  }
}
lgrid <- n
grid <- qnorm(ppoints(lgrid))

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
        target += student_t_lpdf(y[i] | alpha[i], 0, 1);
    };
    target += gamma_lpdf(lambda | 0.1, 100);
    target += normal_lpdf(theta | 0, 100);
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
    array[n] real y_rep = student_t_rng(alpha, rep_vector(u, n),rep_vector(1, n));
    for (i in 1:n) {
        log_lik[i] = student_t_lpdf(y[i] | alpha[i], 0, 1);
    }
}
"
, "C:/Users/Johnny Lee/Documents/.cmdstan/cmdstan-2.33.1/examples/model_simulation_sc3.stan")

file <- file.path(cmdstan_path(), "examples/model_simulation_sc3.stan")
mod <- cmdstan_model(file)

op <- mod$optimize(
  data = list(y = as.vector(y.origin), u = u, p = p, 
                      n= n, psi = psi, atau = ((psi+1)/2), newp = (p+1), bsNonlinear = bs.nonlinear,
                      xholderNonlinear = xholder.nonlinear),
#   init = 0,
  # init = list(list(gamma = t(gamma.origin),
  #                   theta = theta.origin)),
  init = list(list(gamma = array(rep(0.01, (psi*p)), dim=c(p, psi)),
                  theta = 0.1, lambda = 0.1,
                  tau = rep(0.1, p), sigma = 0.1)),
  iter = 3500,
  algorithm = "lbfgs",
  refresh = 50
)

theta.map <- as.vector(op$draws(variables = "theta"))
gamma.map <- as.vector(t(matrix(op$draws(variables = "gamma"), nrow = p)))
lambda.map <- as.vector(op$draws(variables = c("lambda")))
# alpha.map <- as.vector(op$draws(variables = "alpha"))
newalpha.map <- as.vector(op$draws(variables = "newalpha"))
sigma.map <- op$draws(variables = "sigma")
tau.map <- as.vector(op$draws(variables = "tau"))
