library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(rmutil)
library(qqboxplot)

set.seed(123)

n <- 15000
psi.origin <- psi <- 10
threshold <- 0.95
p <- 5

# Function to generate Gaussian copula
C <- diag(p)
x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))

f1 <- function(x) {-0.4 * sin(2 * pi * (x-1.5)^2)*(x-1.5)^3}
f5 <- function(x) {0.1* cos(3 * pi * (x-1.4)^2)*(x-1.4)^2}

make.nl <- function(x, raw_y) {
  fit <- lm(raw_y ~ x)
  
  return(list(
    nl = residuals(fit), 
    slope = coef(fit)[["x"]],
    intercept = coef(fit)[["(Intercept)"]]
  ))
}
theta.origin <- c(1, 0, -0.3, 0, 0, 0.5)

f1.nl <- f1(x.origin[,1])
f5.nl <- f5(x.origin[,5])
f1.l <- theta.origin[2]*x.origin[,1]
f2.l <- theta.origin[3]*x.origin[,2]
f5.l <- theta.origin[6]*x.origin[,5]

eta_lin <-  f1.l + f2.l + f5.l
eta_nonlin <- f1.nl + f5.nl
eta <- theta.origin[1] + eta_lin + eta_nonlin
alp.origin <- exp(eta)

y.origin <- rt(n, df = alp.origin)

u <- quantile(y.origin, threshold)
excess.index <- which(y.origin>u)
x.origin <- x.origin[excess.index,]
y.origin <- y.origin[excess.index]
n <- length(y.origin)

f1.hidden <- make.nl(x.origin[,1], f1(x.origin[,1]))
f5.hidden <- make.nl(x.origin[,5], f5(x.origin[,5]))
theta.adjusted <- c(theta.origin[1] + f1.hidden$intercept + f5.hidden$intercept,
                    theta.origin[2] + f1.hidden$slope,
                    theta.origin[3],
                    theta.origin[4],
                    theta.origin[5],
                    theta.origin[6] + f5.hidden$slope)
newx <- seq(0,1,length.out = n)
xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))
g1.nl <- f1(newx) - (f1.hidden$intercept + f1.hidden$slope*newx)
g5.nl <- f5(newx) - (f5.hidden$intercept + f5.hidden$slope*newx)

g1.l <- theta.adjusted[2]*newx
g2.l <- theta.adjusted[3]*newx
g5.l <- theta.adjusted[6]*newx
g1 <- g1.l + g1.nl
g5 <- g5.l + g5.nl
eta.g <- rep(theta.adjusted[1], n) + g1 + g2.l + g5
alp.new <- as.vector(exp(theta.origin[1] + xholder %*% theta.origin[-1] + f1(newx) + f5(newx)))
# alp.new <- exp(eta.g)
plot(newx, alp.new)

grid.zero <- rep(0, n)
g.new <- c(g1, g2.l, grid.zero, grid.zero, g5)
l.new <- c(g1.l, g2.l, grid.zero, grid.zero, g5.l)
nl.new <- c(g1.nl, grid.zero, grid.zero, grid.zero, g5.nl)


bs.linear <- model.matrix(~ ., data = data.frame(x.origin))
psi <- psi - 2
# Initialize storage
group.map <- c()
Z.list <- list()        # Stores the final non-linear design matrices
scale_stats_list <- list() 
projection_coefs_list <- list() #
spec_decomp_list <- list() # Store eigen-decomp info for prediction
qr_list <- list()          # Store QR info for prediction
sm_spec_list <- list()     # Store smooth objects
keep_cols_list <- list()

covariates <- colnames(data.frame(x.origin))
for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  x_vec <- x.origin[, i]
  X_lin <- model.matrix(~ x_vec) 
  sm_spec <- smoothCon(s(x_vec, bs = "tp", k = psi + 2), 
                      data = data.frame(x_vec = x_vec), 
                      knots = NULL)[[1]]
  
  X_raw <- sm_spec$X
  S     <- sm_spec$S[[1]] 
  
  eig <- eigen(S, symmetric = TRUE)
  max_lambda <- max(eig$values)
  tol <- max_lambda * 1e-6  # Relative threshold (Robust)
  
  pos_idx <- which(eig$values > tol)
  
  if(length(pos_idx) == 0) stop("No wiggles found. Check data scaling.")
  
  U_pen <- eig$vectors[, pos_idx]       
  Lambda_pen <- diag(eig$values[pos_idx]) 
  
  Z_spectral <- X_raw %*% U_pen %*% solve(sqrt(Lambda_pen))
  qr_lin <- qr(X_lin)
  Q_lin <- qr.Q(qr_lin)
  R_lin <- qr.R(qr_lin)
  
  Gamma_Q <- t(Q_lin) %*% Z_spectral
  Gamma_Original <- backsolve(R_lin, Gamma_Q)
  Z_orth <- Z_spectral - X_lin %*% Gamma_Original
  keep_cols <- colSums(Z_orth^2) > 1e-9
  Z_final <- Z_orth[, keep_cols, drop = FALSE]
  train_scale <- apply(Z_final, 2, sd)
  Z_final <- scale(Z_final, center = FALSE, scale = train_scale)
  
  # Store Results
  Z.list[[i]] <- Z_final
  group.map <- c(group.map, rep(i, ncol(Z_final)))
  
  # Store Transforms
  spec_decomp_list[[i]] <- list(U_pen = U_pen, Lambda_sqrt_inv = solve(sqrt(Lambda_pen)))
  projection_coefs_list[[i]] <- Gamma_Original # Store the X-basis coefs
  keep_cols_list[[i]] <- keep_cols
  scale_stats_list[[i]] <- train_scale         # Store the scale
  sm_spec_list[[i]] <- sm_spec
}

bs.nonlinear <- do.call(cbind, Z.list)


colnames(xholder) <- covariates
grid_Z_list <- list()

for (i in seq_along(covariates)) {
  grid_df  <- data.frame(x_vec = newx)
  X_lin_grid <- model.matrix(~ x_vec, data = grid_df)
  X_raw_grid <- PredictMat(sm_spec_list[[i]], grid_df)
  
  decomp <- spec_decomp_list[[i]]
  Z_spectral_grid <- X_raw_grid %*% decomp$U_pen %*% decomp$Lambda_sqrt_inv
  Z_orth_grid <- Z_spectral_grid - X_lin_grid %*% projection_coefs_list[[i]]
  Z_final_grid <- Z_orth_grid[, keep_cols_list[[i]], drop = FALSE]
  Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = scale_stats_list[[i]])
  
  grid_Z_list[[i]] <- Z_final_grid
}

xholder.nonlinear <- do.call(cbind, grid_Z_list)
xholder.linear <- model.matrix(~ ., data = data.frame(xholder))
Z_scales <- unlist(scale_stats_list)
X_means <- colMeans(bs.linear[,-1])
X_sd   <- apply(bs.linear[,-1], 2, sd)
bs.linear[,-1] <- scale(bs.linear[,-1], center = X_means, scale = X_sd)


# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,2)]) %*% bs.nonlinear[,c((1):(psi))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,3)]) %*% bs.nonlinear[,c((psi+1):(psi*2))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,4)]) %*% bs.nonlinear[,c((psi*2+1):(psi*3))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,5)]) %*% bs.nonlinear[,c((psi*3+1):(psi*4))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,6)]) %*% bs.nonlinear[,c((psi*4+1):(psi*5))]), "\n")


model.stan <- "// Stan model for BLAST Student-t Samples
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
  vector[p] X_sd;
  vector[(psi*p)] Z_scales;
}

parameters {
  vector[(p+1)] theta; // linear predictor
  array[p] vector[psi] gamma_raw;
  array[p] real <lower=0> lambda1; // lasso penalty //array[p] real <lower=0> 
  array[p] real <lower=0> lambda2; // lambda2 group lasso penalty
  array[p] real <lower=0> tau;
}

transformed parameters {
  array[n] real <lower=0> alpha; // covariate-adjusted tail index
  
  array[p] vector[psi] gamma;
  {
    matrix[n, p] gsmooth; // linear component
    for (j in 1:p){
      for (k in 1:psi){
        int idx = (j-1)*psi + k;
        gamma[j][k] = gamma_raw[j][k] * sqrt(tau[j]) * Z_scales[idx];
      }; 
      gsmooth[,j] = bsLinear[,j] * theta[j+1] + bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    
    for (i in 1:n){
      alpha[i] = exp(theta[1] + sum(gsmooth[i,]));
    };
  }
}

model {
  // likelihood
  for (i in 1:n){
   target += student_t_lpdf(y[i] | alpha[i], 0, 1);
   target += -1*log(1-student_t_cdf(u, alpha[i], 0, 1));
  }
  // target += pareto_lpdf(y | rep_vector(u, n), alpha);
  target += normal_lpdf(theta[1] | 0, 100);
  for (j in 1:p){
    target += gamma_lpdf(lambda1[j] | 2, 1); 
    target += gamma_lpdf(lambda2[j] | 1e-3, 1e-3);
    target += double_exponential_lpdf(theta[(j+1)] | 0, 1/(lambda1[j]));
    target += gamma_lpdf(tau[j] | atau, square(lambda2[j])*0.5);
    target += std_normal_lpdf(gamma_raw[j]);
  }
}

generated quantities {
  // Used in Posterior predictive check
  array[n] real <lower=0> gridalpha; // new tail index
  matrix[n, p] gridgnl; // nonlinear component
  matrix[n, p] gridgl; // linear component
  matrix[n, p] gridgsmooth; // linear component
  vector[(p+1)] theta_origin;

  for (j in 1:p){
    theta_origin[j+1] = theta[j+1] / X_sd[j];
  }
  theta_origin[1] = theta[1] - dot_product(X_means, theta_origin[2:(p+1)]);
  for (j in 1:p){
      gridgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
      gridgl[,j] = xholderLinear[,j] * theta_origin[j+1];
      gridgsmooth[,j] = gridgl[,j] + gridgnl[,j];
  };

  for (i in 1:n){
      gridalpha[i] = exp(theta_origin[1] + sum(gridgsmooth[i,]));
  };
}
"

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                  atau = ((psi+1)/2), X_means = X_means, X_sd=X_sd,
                  bsLinear = bs.linear[,-1], bsNonlinear = bs.nonlinear,
                  xholderLinear = xholder.linear[,-1], xholderNonlinear = xholder.nonlinear, Z_scales = Z_scales)

init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p,psi)),
                        theta = rep(-0.1, (p+1)), tau = rep(0.1, p), 
                        lambda1 = rep(0.1, p), lambda2 = rep(1, p)),
                  list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p,psi)),
                        theta = rep(0.05, (p+1)), tau = rep(0.5, p),
                        lambda1 = rep(0.5, p), lambda2 = rep(0.2, p)),
                  list(gamma_raw = array(rep(0.1, (psi*p)), dim=c(p,psi)),
                        theta = rep(0.05, (p+1)), tau = rep(0.2, p),
                        lambda1 = rep(0.2, p), lambda2 = rep(0.5, p)))

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

bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("Scenario C") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/appendix_",n,"_traceplot_scC.pdf"), width=22, height = 3)

# tau.samples <- summary(fit1, par=c("tau"), probs = c(0.05,0.5, 0.95))$summary
theta.samples <- summary(fit1, par=c("theta_origin"), probs = c(0.05,0.5, 0.95))$summary
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
                       "true" = theta.adjusted,
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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_theta_scC.pdf"), width=10, height = 7.78)

df.gamma <- data.frame("seq" = seq(1, (psi*p)), 
                       "m" = as.vector(gamma.q2),
                       "l" = as.vector(gamma.q1),
                       "u" = as.vector(gamma.q3))
df.gamma$covariate <- factor(rep(seq(1, 1 + nrow(df.gamma) %/% psi), each = psi, length.out = nrow(df.gamma)))
df.gamma$labels <- factor(1:(psi*p))
ggplot(df.gamma, aes(x =labels, y = m, color = covariate)) + 
  geom_errorbar(aes(ymin = l, ymax = u), alpha = 0.4, width = 4, linewidth = 1.2) +
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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_gamma_scC.pdf"), width=10, height = 7.78)

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

data.smooth <- data.frame("x"=newx,
                          "true" = g.new,
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]")),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=true, colour = "True"), linewidth=2, linetype=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1.5) + 
  ylab("") + xlab(expression(c)) + ylim(-2.5, 2.5) +
  facet_grid(covariates ~ ., scales = "free_x", switch = "y", 
              labeller = label_parsed) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position="none",
        plot.margin = margin(0,0,0,-20),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 45),
        axis.text.x = element_text(size=30))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_smooth_scC.pdf"), width=11, height = 15)

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
  geom_line(aes(y=true, colour = "True"), linewidth=2, linetype=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1.5) + 
  # geom_line(aes(y=post.mean, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab(expression(c)) + ylim(-2.5, 2.5) + 
  facet_grid(covariates ~ ., scales = "free_x", switch = "y",
              labeller = label_parsed) +  
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
         fill = guide_legend(order = 1)) + 
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position="none",
        plot.margin = margin(0,0,0,-20),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 45),
        axis.text.x = element_text(size=30))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_linear_scC.pdf"), width=11, height = 15)


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
  geom_line(aes(y=true, colour = "True"), linewidth=2, linetype=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1.5) + 
  ylab("") + xlab(expression(c)) + ylim(-2.5, 2.5) + 
  facet_grid(covariates ~ ., scales = "free_x", switch = "y",
              labeller = label_parsed) +  
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue", "red")) + 
  guides(color = guide_legend(order = 2), 
         fill = guide_legend(order = 1)) + 
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position="none",
        plot.margin = margin(0,0,0,-20),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 45),
        axis.text.x = element_text(size=30))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_nonlinear_scC.pdf"), width=11, height = 15)

data.scenario <- data.frame("x" = newx,
                            "true" = (alp.new),
                            "post.mean" = (newalpha.samples[,1]),
                            "q2" = (newalpha.samples[,5]),
                            "q1" = (newalpha.samples[,4]),
                            "q3" = (newalpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab("") + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2, linetype=2) + ylim(0, 21.2) + 
  geom_line(aes(y=q2, col = "Posterior Median"), linewidth=1.5) +
  # geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=1.5) +
  scale_color_manual(values=c("steelblue", "red")) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  theme_minimal(base_size = 40) + 
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 30))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_scC.pdf"), width=9.5, height = 7.78)

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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_qqplot_scC.pdf"), width=10, height = 7.78)
