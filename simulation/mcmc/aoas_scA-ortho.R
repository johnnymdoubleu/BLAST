# library(npreg)
library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(qqboxplot)

#Scenario A-2
# set.seed(10)

n <- 15000
psi <- 10
threshold <- 0.95
p <- 5

# Function to generate Gaussian copula
C <- diag(p)
xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
x.origin <- data.frame(x.origin)


# ## Generate sample
psi <- psi -2
gamma.origin <- matrix(, nrow=(psi), ncol=p)
for(j in 1:p){
  if(j %in% c(1,4,5)){gamma.origin[,j] <- rep(0, psi)}
  else{gamma.origin[,j] <- rep(0, psi)} #runif(psi, 0.25, 0.75)}
}

theta.origin <- c(0.5, 0, 0.25, 0.25, 0, 0)
f <- function(x){
  0.25*x[,2] + x[,2]^2+ 
  0.25*x[,3] - 1.5*x[,3]^3
}
alp.origin <- exp(rep(theta.origin[1], n) + f(x.origin))
y.origin <- rPareto(n, rep(1,n), alpha = alp.origin)

u <- quantile(y.origin, threshold)
excess.index <- which(y.origin>u)
x.origin <- x.origin[excess.index,]
y.origin <- y.origin[y.origin > u]
n <- length(y.origin)
newx <- seq(0,1,length.out = n)
f.l <- function(x) {0.25*x}
f2.nl <- function(x) {x^2}
f2 <- f.l(newx) + f2.nl(newx)
f3.nl <- function(x) {-1.5*x^3}
f3 <- f.l(newx) + f3.nl(newx)
grid.zero <- rep(0, n)
g.new <- c(grid.zero, f2, f3, grid.zero, grid.zero)
l.new <- c(grid.zero, f.l(newx), f.l(newx), grid.zero, grid.zero)
nl.new <- c(grid.zero, f2.nl(newx), f3.nl(newx), grid.zero, grid.zero)

covariates <- colnames(x.origin) 
p <- length(covariates)
n_train <- nrow(x.origin)

bs.linear <- model.matrix(~ ., data = x.origin)
qr.x <- base::qr(bs.linear)
q.x <- qr.Q(qr.x)
group.map <- c()
Z.list <- list()
Gamma_list <- list()      
keep_cols_list <- list()  
sm_spec_list <- list() 
scale_stats_list <- list()

# --- 2. TRAINING LOOP ---
for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  
  sm_expr <- call("s", as.symbol(var_name), bs = "tp", k = (psi+2))
  sm_spec <- eval(sm_expr)
  sm <- smoothCon(sm_spec, data = x.origin, knots = NULL)[[1]]
  Z_raw <- as.matrix(sm$X)
  Gamma <- t(q.x) %*% Z_raw 
  Z_orth <- Z_raw - q.x %*% Gamma
  keep_cols <- colSums(abs(Z_orth)) > 1e-9
  Z_clean <- Z_orth[, keep_cols, drop = FALSE]
  
  col_sds <- apply(Z_clean, 2, sd)
  col_sds[col_sds < 1e-9] <- 1 # Safety
  Z_clean <- scale(Z_clean, center = FALSE, scale = col_sds)
  Z.list[[i]] <- Z_clean
  group.map <- c(group.map, rep(i, ncol(Z_clean)))
  Gamma_linear <- backsolve(qr.R(qr.x), Gamma)
  
  Gamma_list[[i]] <- Gamma_linear
  keep_cols_list[[i]] <- keep_cols
  sm_spec_list[[i]] <- sm
  scale_stats_list[[i]] <- col_sds
}

bs.nonlinear <- do.call(cbind, Z.list)


grid_seq <- seq(0, 1, length.out = n)
grid_Z_list <- list() 

xholder <- matrix(0, nrow = n, ncol = p)
colnames(xholder) <- covariates

for (i in 1:p) {
  var_name <- covariates[i] 
  xholder[,i] <- grid_seq
  grid_df <- data.frame(matrix(0, nrow = n, ncol = p))
  colnames(grid_df) <- covariates # Must match x.origin
  grid_df[[var_name]] <- grid_seq 
  X_grid <- model.matrix(~ ., data = grid_df)
  Z_grid_raw <- PredictMat(sm_spec_list[[i]], grid_df)
  Gamma_curr <- Gamma_list[[i]]
  Z_grid_orth <- Z_grid_raw - X_grid %*% Gamma_curr
  keep_cols <- keep_cols_list[[i]]
  Z_grid_final <- Z_grid_orth[, keep_cols, drop = FALSE]
  Z_grid_final <- scale(Z_grid_final, center = FALSE, scale = scale_stats_list[[i]])
  
  grid_Z_list[[i]] <- Z_grid_final
}

xholder.linear <- model.matrix(~ ., data = data.frame(xholder))
xholder.nonlinear <- do.call(cbind, grid_Z_list)

sum(t(bs.linear) %*% bs.nonlinear)

true.alpha <- alp.new <- alp.origin <- NULL
# alp.origin <- exp(rep(theta.origin[1], n) + f(bs.linear[,c(2:(p+1))]))
alp.new <- exp(rep(theta.origin[1], n) + f(xholder.linear[,c(2:(p+1))]))

X_means <- colMeans(bs.linear[,c(2:(p+1))])
bs.linear[,c(2:(p+1))] <- scale(bs.linear[,c(2:(p+1))], center = X_means, scale = FALSE)

model.stan <- "// Stan model for BLAST Pareto Samples
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n, p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n, p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    array[n] real <lower=1> y; // extreme response
    real <lower=0> atau;
    vector[p] X_means;
}

parameters {
    vector[(p+1)] theta; // linear predictor
    array[p] vector[psi] gamma_raw;
    array[p] real <lower=0> lambda1; // lasso penalty //array[p] real <lower=0> 
    array[p] real <lower=0> lambda2; // lambda2 group lasso penalty
    array[p] real <lower=0> tau;
    real <lower=0> rho;
}

transformed parameters {
    array[n] real <lower=0> alpha; // covariate-adjusted tail index
    real theta0;
    array[p] vector[psi] gamma;
    {
      matrix[n, p] gsmooth; // linear component
      for (j in 1:p){
        gamma[j] = gamma_raw[j] * sqrt(tau[j]);
        gsmooth[,j] = bsLinear[,j] * theta[j+1] + bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
      };
      theta0 = theta[1] - dot_product(X_means, theta[2:(p+1)]);
      for (i in 1:n){
        alpha[i] = exp(theta[1] + sum(gsmooth[i,]));
      };
    }
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += normal_lpdf(theta[1] | 0, 10);
    target += gamma_lpdf(rho | 2, 2);
    for (j in 1:p){
        target += gamma_lpdf(lambda1[j] | 0.1, 0.1); 
        target += gamma_lpdf(lambda2[j] | 0.1, 0.1);
        target += double_exponential_lpdf(theta[(j+1)] | 0, 1/(lambda1[j]*rho));
        target += gamma_lpdf(tau[j] | atau, square(lambda2[j]*rho)*0.5);
        target += std_normal_lpdf(gamma_raw[j]); // target += multi_normal_lpdf(gamma_raw[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)))
    }
}

generated quantities {
    // Used in Posterior predictive check    
    vector[n] log_lik;
    array[n] real <lower=0> gridalpha; // new tail index
    matrix[n, p] gridgnl; // nonlinear component
    matrix[n, p] gridgl; // linear component
    matrix[n, p] gridgsmooth; // linear component

    for(i in 1:n){
      log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
    for (j in 1:p){
        gridgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        gridgl[,j] = xholderLinear[,j] * theta[j+1];
        gridgsmooth[,j] = gridgl[,j] + gridgnl[,j];
    };

    for (i in 1:n){
        gridalpha[i] = exp(theta[1] + sum(gridgsmooth[i,]));
    };
}
"

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                  atau = ((psi+1)/2), X_means = X_means,
                  bsLinear = bs.linear[,c(2:(p+1))], bsNonlinear = bs.nonlinear,
                  xholderLinear = xholder.linear[,c(2:(p+1))], xholderNonlinear = xholder.nonlinear)

init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p,psi)),
                        theta = rep(-0.1, (p+1)), tau = rep(0.1, p), rho = 1, 
                        lambda1 = rep(0.1, p), lambda2 = rep(1, p)),
                   list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p,psi)),
                        theta = rep(0.05, (p+1)), tau = rep(2, p), rho = 1,
                        lambda1 = rep(2, p), lambda2 = rep(5, p)),
                   list(gamma_raw = array(rep(-0.75, (psi*p)), dim=c(p,psi)),
                        theta = rep(0.01, (p+1)), tau = rep(1.5, p), rho = 1,
                        lambda1 = rep(0.5,p), lambda2= rep(5, p)))

system.time(fit1 <- stan(
  model_code = model.stan,  # Stan program
  data = data.stan,    # named list of data
  init = init.alpha,      # initial value
  chains = 3,             # number of Markov chains
  # warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = parallel::detectCores(), # number of cores (could use one per chain)
  refresh = 1000             # no progress shown
))

posterior <- extract(fit1)

# plot(fit1, plotfun = "trace", pars = c("lp__"), nrow = 3)
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("Scenario A2") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/appendix_",n,"_traceplot_scA.pdf"), width=22, height = 3)

# tau.samples <- summary(fit1, par=c("tau"), probs = c(0.05,0.5, 0.95))$summary
theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary

newgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
newgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
newalpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary

gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q2 <- gamma.samples[,1]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q2 <- theta.samples[,5]
theta.q3 <- theta.samples[,6]

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
  geom_errorbar(aes(ymin = l, ymax = u), alpha = 0.4, width = 4, linewidth = 1.2) +
  # geom_point(aes(y=true), size =4, color ="red")+
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_point(size = 4) + ylab("") + xlab("" ) + 
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

data.smooth <- data.frame("x"=seq(0,1,length.out = n),
                          "true" = as.vector(g.new),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]")),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=true, colour = "True"), linewidth=2, linetype=2) + 
  # geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1.8) + 
  geom_line(aes(y=post.mean, colour = "Posterior Median"), linewidth=1.8) + 
  ylab("") + xlab(expression(c)) + 
  facet_grid(covariates ~ ., scales = "free_x", switch = "y", 
              labeller = label_parsed) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) + #ylim(-1.3, 1.3) + 
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position="none",
        plot.margin = margin(0,0,0,-20),
        strip.text.y = element_text(size = 38, colour = "black", angle = 0, face = "bold.italic"),
        strip.placement = "outside",
        axis.title.x = element_text(size = 45),
        axis.text = element_text(size=30))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_smooth_sc1-wi.pdf"), width=12.5, height = 15)

data.linear <- data.frame("x"=seq(0,1, length.out = n),
                          "true" = as.vector(l.new),
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
  # geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  geom_line(aes(y=post.mean, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab(expression(c)) +
  facet_wrap(covariates ~ ., scales = "free_x", nrow = 5,
             labeller = label_parsed, strip.position = "left") + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
         fill = guide_legend(order = 1)) + #ylim(-0.55, 0.5) + 
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


data.nonlinear <- data.frame("x"=seq(0,1, length.out=n),
                             "true" = as.vector(nl.new),
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

data.scenario <- data.frame("x" = seq(0,1, length.out = n),
                            # "constant" = newx,
                            "true" = (alp.new),
                            "post.mean" = (newalpha.samples[,1]),
                            "q2" = (newalpha.samples[,5]),
                            "q1" = (newalpha.samples[,4]),
                            "q3" = (newalpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,...,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2, linetype=2) + 
  # geom_line(aes(y=q2, col = "Posterior Median"), linewidth=1.5) +
  geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=1.5) +
  scale_color_manual(values=c("steelblue", "red")) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  theme_minimal(base_size = 40) + 
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 30))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_test_sc1-wi.pdf"), width=10, height = 7.78)

mcmc.alpha <- posterior$alpha
len <- dim(mcmc.alpha)[1]
r <- matrix(, nrow = n, ncol = 30)
T <- 30
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y.origin[i], u, alpha = mcmc.alpha[round(runif(1,1,len)),i]))
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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_qqplot_sc1-wi.pdf"), width=10, height = 7.78)
library(loo)
fit.log.lik <- extract_log_lik(fit1)
loo(fit.log.lik, is_method = "sis", cores = 2)
