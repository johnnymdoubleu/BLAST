# library(npreg)
library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(qqboxplot)
# library(qgam)
library(evgam)

set.seed(1001)
n <- 10000
psi <- 10
threshold_prob <- 0.95
p <- 5
C <- diag(p)
x.random <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
time.seq <- 1:n
period <- 365 
x.season <- (time.seq %% period) / period 

# Convert continuous season to Factor for the 'by' argument
season_code_full <- cut(x.season, breaks = c(-0.1, 0.25, 0.5, 0.75, 1.1), labels = c("1","2","3","4"))



f.season.scale <- function(x) {
  return(1 + 0.6 * sin(2 * pi * x) + 0.4 * cos(2 * pi * x)) 
}

x.origin <- cbind(x.random, x.season)
covariates <- colnames(data.frame(x.origin))[1:p] 

# True Functions
f2 <- function(x) {1.5 * sin(2 * pi * x^2)*x^3}
f3 <- function(x) {-1.5 * cos(3 * pi * x^2)*x^2}

theta.origin <- c(1, 0, 0.8, -0.4, 0, 0) # Coefs for 5 random vars + intercept
alp.origin <- as.vector(exp(
  rep(theta.origin[1],n) + 
  x.origin[,1:5] %*% theta.origin[-1] + 
  f2(x.origin[,2]) + 
  f3(x.origin[,3]) 
  # f.season.shape(x.origin[,6]) is removed
))

y.noise <- rPareto(n, rep(1,n), alpha = alp.origin)
y.origin <- y.noise + f.season.scale(x.origin[,6])

# Dynamic Thresholding
evgam.df <- data.frame(
  y = y.origin,
  sin.time = sin(2 * pi * x.season),
  cos.time = cos(2 * pi * x.season)
)
evgam.cov <- y ~ 1 + cos.time + sin.time
ald.cov.fit <- evgam(evgam.cov, data = evgam.df, family = "ald", ald.args=list(tau = threshold_prob))
u.vec <- (predict(ald.cov.fit)$location)

excess.index <- which(y.origin > u.vec)
x.origin <- x.origin[excess.index,]
y.origin <- y.origin[excess.index]
u <- u.vec[excess.index]
season_code_full <- season_code_full[excess.index]
n <- length(y.origin)

make.nl <- function(x, raw_y) {
  fit <- lm(raw_y ~ x)
  
  return(list(
    nl = residuals(fit), 
    slope = coef(fit)[["x"]],
    intercept = coef(fit)[["(Intercept)"]]
  ))
}

newx <- seq(0,1,length.out = n)
xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))

season_levels <- c("1", "2", "3", "4")
true_alpha_list <- list()
true_smooth_list <- list()

f2.hidden <- make.nl(x.origin[,2], f2(x.origin[,2]))
f3.hidden <- make.nl(x.origin[,3], f3(x.origin[,3]))

# --- MODIFIED: Ground Truth Reconstruction (Constant across seasons) ---
grid_n <- n
grid_x <- seq(0, 1, length.out = grid_n)
raw_f2 <- f2(grid_x)
raw_f3 <- f3(grid_x)
theta.adjusted <- c(theta.origin[1] + f2.hidden$intercept + f3.hidden$intercept,
                    theta.origin[2], # Theta for X1 (0)
                    theta.origin[3] + f2.hidden$slope, # Theta for X2
                    theta.origin[4] + f3.hidden$slope, # Theta for X3
                    theta.origin[5], # Theta for X4 (0)
                    theta.origin[6]) # Theta for X5 (0)

g1.total <- theta.adjusted[2] * grid_x 

g2.nl <- raw_f2 - (f2.hidden$intercept + f2.hidden$slope * grid_x)
g2.l  <- theta.adjusted[3] * grid_x
g2.total <- g2.l + g2.nl

g3.nl <- raw_f3 - (f3.hidden$intercept + f3.hidden$slope * grid_x)
g3.l  <- theta.adjusted[4] * grid_x
g3.total <- g3.l + g3.nl

g4.total <- theta.adjusted[5] * grid_x
g5.total <- theta.adjusted[6] * grid_x

true_alpha <- exp(theta.adjusted[1] + g2.total + g3.total)

for(s_label in season_levels) {
  true_alpha_list[[s_label]] <- data.frame(
    x = grid_x, 
    alpha = true_alpha, 
    season = s_label
  )
  
  true_smooth_list[[s_label]] <- rbind(
    data.frame(x=grid_x, y=g1.total, variable="g1", season=s_label),
    data.frame(x=grid_x, y=g2.total, variable="g2", season=s_label),
    data.frame(x=grid_x, y=g3.total, variable="g3", season=s_label),
    data.frame(x=grid_x, y=g4.total, variable="g4", season=s_label),
    data.frame(x=grid_x, y=g5.total, variable="g5", season=s_label)
  )
}

df_true_alpha <- do.call(rbind, true_alpha_list)
df_true_smooth <- do.call(rbind, true_smooth_list)
df_true_smooth$global <- c(g1.total, g2.total, g3.total, g4.total, g5.total)

fwi.basis_list <- list()
linear_basis_list <- list()
nonlinear_basis_list <- list()
model_data <- list() 

fwi.scaled <- x.origin 
fwi.origin <- x.origin 

for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  x_vec <- fwi.scaled[, i]
  x_true <- fwi.origin[,i]
  
  model_data[[var_name]] <- list() 
  sm_list <- smoothCon(mgcv::s(x_vec, by=season_code_full, bs = "tp", k = psi + 2), 
                       data = data.frame(x_vec = x_vec, season_code_full = season_code_full), 
                       knots = NULL)
  
  for (j in 1:length(sm_list)) {
    season_label <- levels(season_code_full)[j]
    sm_spec <- sm_list[[j]] 
    mask <- as.numeric(season_code_full == season_label)
    vec_intercept <- mask
    vec_slope <- x_vec * mask
    
    lin_col_name <- paste(var_name, "S", season_label, sep="_")
    linear_basis_list[[lin_col_name]] <- matrix(vec_slope, ncol=1, dimnames=list(NULL, lin_col_name))
    
    fwi.basis_list[[lin_col_name]] <- matrix(x_true*mask, ncol=1, dimnames=list(NULL, lin_col_name))
    X_lin_local <- cbind(vec_intercept, vec_slope)
    
    qr_lin <- qr(X_lin_local)
    Q_lin <- qr.Q(qr_lin)
    R_lin <- qr.R(qr_lin)
    
    # 2. Spectral Decomposition
    X_raw <- sm_spec$X
    S     <- sm_spec$S[[1]] 
    
    eig <- eigen(S, symmetric = TRUE)
    # Adjusted tolerance slightly for robustness
    pos_idx <- which(eig$values > max(eig$values) * 1e-8)
    
    if(length(pos_idx) == 0) next 
    U_pen <- eig$vectors[, pos_idx]       
    Lambda_pen <- diag(eig$values[pos_idx]) 
    Lambda_sqrt_inv <- solve(sqrt(Lambda_pen))
    Z_spectral <- X_raw %*% U_pen %*% Lambda_sqrt_inv
    Gamma_Q <- t(Q_lin) %*% Z_spectral
    Gamma_Original <- backsolve(R_lin, Gamma_Q)
    Z_orth <- Z_spectral - X_lin_local %*% Gamma_Original
    keep_cols <- colSums(Z_orth^2) > 1e-9
    Z_final <- Z_orth[, keep_cols, drop = FALSE]
    
    train_scale <- apply(Z_final, 2, sd)
    train_scale[train_scale < 1e-12] <- 1 
    Z_final <- scale(Z_final, center = FALSE, scale = train_scale)
    Z_final <- Z_final * mask
    
    col_names <- paste(var_name, "S", season_label, 1:ncol(Z_final), sep="_")
    colnames(Z_final) <- col_names
    nonlinear_basis_list[[paste(var_name, season_label, sep="_")]] <- Z_final
    model_data[[var_name]][[season_label]] <- list(
      sm_spec = sm_spec, 
      U_pen = U_pen, 
      Lambda_sqrt_inv = Lambda_sqrt_inv, 
      projection_coefs = Gamma_Original, 
      keep_cols = keep_cols, 
      scale_stats = train_scale
    )
  }
}

# Create Final Training Matrices
true.linear <- do.call(cbind, fwi.basis_list)
fwi.linear <- bs.linear <- do.call(cbind, linear_basis_list) 
bs.nonlinear <- do.call(cbind, nonlinear_basis_list)

X_means <- colMeans(bs.linear)
X_sd   <- apply(bs.linear, 2, sd)
Z_scales <- rep(1, ncol(bs.nonlinear)) 
bs.linear_check <- cbind(rep(1,n), bs.linear)
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,2)]) %*% bs.nonlinear[,c(1:psi)]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,2,3,4,5)]) %*% bs.nonlinear[,c((1):(4*psi))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,6,7,8,9)]) %*% bs.nonlinear[,c((4*(psi+1)):(4*(psi*2)))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,10,11,12,13)]) %*% bs.nonlinear[,c((4*(psi*2+1)):(4*(psi*3)))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,14,15,16,17)]) %*% bs.nonlinear[,c((4*(psi*3+1)):(4*(psi*4)))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,18,19,20,21)]) %*% bs.nonlinear[,c((4*(psi*4+1)):(4*(psi*5)))]), "\n")

newx <- seq(0, 1, length.out = n)
xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))
grid_Z_list <- list()
grid_L_list <- list()

for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  for (season_lbl in levels(season_code_full)) {
    params <- model_data[[var_name]][[season_lbl]]
    grid_L_list[[paste(var_name, season_lbl)]] <- matrix(newx, ncol=1)
    grid_df <- data.frame(
      x_vec = newx, 
      season_code_full = factor(season_lbl, levels = levels(season_code_full))
    )
    
    X_raw_grid <- PredictMat(params$sm_spec, grid_df)
    X_lin_grid <- cbind(rep(1, n), newx)
    Z_spectral_grid <- X_raw_grid %*% params$U_pen %*% params$Lambda_sqrt_inv
    Z_orth_grid <- Z_spectral_grid - X_lin_grid %*% params$projection_coefs
    Z_final_grid <- Z_orth_grid[, params$keep_cols, drop = FALSE]
    Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)
    
    grid_Z_list[[paste(var_name, season_lbl)]] <- Z_final_grid
  }
}

xholder.linear <- do.call(cbind, grid_L_list)
xholder.nonlinear <- do.call(cbind, grid_Z_list)


scale_stats_list <- list()
for (var_name in names(model_data)) {
  for (season_label in names(model_data[[var_name]])) {
    stats <- model_data[[var_name]][[season_label]]$scale_stats
    scale_stats_list[[paste(var_name, "S", season_label, sep="_")]] <- stats
  }
}

Z_scales <- unlist(scale_stats_list)
X_sd   <- apply(bs.linear, 2, sd)
bs.linear <- scale(bs.linear, center = FALSE, scale = X_sd)

model.stan <- "// Stan model for BLAST Pareto Samples
data {
  int <lower=1> n; // Sample size
  int <lower=1> p; // regression coefficient size
  int <lower=1> psi; // splines coefficient size
  int <lower=1> n_seasons;
  matrix[n, (p*n_seasons)] bsLinear; // fwi dataset
  matrix[n, (psi*p*n_seasons)] bsNonlinear; // thin plate splines basis
  matrix[n, (p*n_seasons)] xholderLinear; // fwi dataset
  matrix[n, (psi*p*n_seasons)] xholderNonlinear; // thin plate splines basis
  vector[n] u; // large threshold value
  vector[n] y; // extreme response
  real <lower=0> atau;
  vector[(psi*p*n_seasons)] Z_scales;
  vector[p * n_seasons] X_sd;
}

parameters {
  real theta0;
  vector[p * n_seasons] theta; // linear predictor
  array[p, n_seasons] vector[psi] gamma_raw;
  array[p, n_seasons] real <lower=0> lambda1; // lasso penalty //array[p] real <lower=0> 
  array[p, n_seasons] real <lower=0> lambda2; // lambda2 group lasso penalty
  array[p, n_seasons] real <lower=0> tau;
}

transformed parameters {
  vector[n] alpha;        // Tail index
  array[p, n_seasons] vector[psi] gamma; // Scaled spline coefficients
  
  {
    vector[n] lin_predictor;
    lin_predictor = rep_vector(theta0, n); 
    for (j in 1:p){
      for (s in 1:n_seasons){
        int lin_idx = (j-1) * n_seasons + s;
        int nl_start = (lin_idx-1) * psi + 1; 
        int nl_end   = lin_idx * psi;
        
        for (k in 1:psi){
           int flat_idx = nl_start + k - 1;
           gamma[j,s][k] = gamma_raw[j,s][k] * sqrt(tau[j,s]);
        }
        
        lin_predictor += col(bsLinear, lin_idx) * theta[lin_idx] 
                       + block(bsNonlinear, 1, nl_start, n, psi) * gamma[j,s];
      }
    }
    alpha = exp(lin_predictor);
  }
}

model {
  // likelihood
  target += pareto_lpdf(y | u, alpha);
  target += normal_lpdf(theta0 | 0, 10);
  for (j in 1:p){
    for (s in 1:n_seasons){
      int idx = (j-1) * n_seasons + s;
      target += gamma_lpdf(lambda1[j,s] | 1, 1); 
      target += gamma_lpdf(lambda2[j,s] | 1e-2, 1e-2);
      target += double_exponential_lpdf(theta[idx] | 0, 1 / lambda1[j,s]);
      target += gamma_lpdf(tau[j,s] | atau, square(lambda2[j,s]) * 0.5);
      target += std_normal_lpdf(gamma_raw[j,s]);
    }
  }
}

generated quantities {
  matrix[n, n_seasons] alphaseason;
  vector[n] log_lik;
  matrix[n, p * n_seasons] gridgl;  // Linear component
  matrix[n, p * n_seasons] gridgnl; // Non-linear component
  matrix[n, p * n_seasons] gridgsmooth; // Non-linear component

  vector[p * n_seasons] theta_origin = theta ./ X_sd;

  for (j in 1:p){
    for (s in 1:n_seasons){
      int lin_idx = (j-1) * n_seasons + s;
      int nl_start = (lin_idx-1) * psi + 1;
      
      gridgl[, lin_idx] = col(xholderLinear, lin_idx) * theta_origin[lin_idx];
      gridgnl[, lin_idx] = block(xholderNonlinear, 1, nl_start, n, psi) * gamma[j,s];
    }
  }
  gridgsmooth = gridgl + gridgnl;
  
  {
    for (s in 1:n_seasons) {
      vector[n] current_season_pred = rep_vector(theta0, n);
      for (j in 1:p) {
        int lin_idx = (j-1) * n_seasons + s;
        current_season_pred += gridgsmooth[, lin_idx];
      }
      alphaseason[, s] = exp(current_season_pred);
    }
  }
  
  for (i in 1:n){
    log_lik[i] = pareto_lpdf(y[i] | u[i], alpha[i]);
  }
}
"


data.stan <- list(
  y = as.vector(y.origin), 
  u = u, 
  n = n, 
  p = p, 
  n_seasons = 4, 
  psi = psi, 
  bsLinear = bs.linear, 
  bsNonlinear = bs.nonlinear,
  xholderLinear = xholder.linear, 
  xholderNonlinear = xholder.nonlinear, 
  atau = ((psi+1)/2), 
  Z_scales = Z_scales,
  X_sd = X_sd
)

init.alpha <- list(
  list(
    theta0 = 0.1,                                  
    theta = rep(-0.1, p * 4),
    gamma_raw = array(rep(0.1, p * 4 * psi), dim = c(p, 4, psi)), 
    tau = array(rep(0.1, p * 4), dim = c(p, 4)),
    lambda1 = array(rep(0.1, p * 4), dim = c(p, 4)), 
    lambda2 = array(rep(1, p * 4), dim = c(p, 4))
  ),
  list(
    theta0 = 0.05,
    theta = rep(0.05, p * 4),
    gamma_raw = array(rep(-0.1, p * 4 * psi), dim = c(p, 4, psi)),
    tau = array(rep(0.2, p * 4), dim = c(p, 4)),
    lambda1 = array(rep(0.2, p * 4), dim = c(p, 4)), 
    lambda2 = array(rep(0.5, p * 4), dim = c(p, 4))
  ),
  list(
    theta0 = -0.1,
    theta = rep(0.1, p * 4),
    gamma_raw = array(rep(0.05, p * 4 * psi), dim = c(p, 4, psi)),
    tau = array(rep(0.1, p * 4), dim = c(p, 4)),
    lambda1 = array(rep(0.1, p * 4), dim = c(p, 4)), 
    lambda2 = array(rep(0.1, p * 4), dim = c(p, 4))
  )
)

fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value 
    chains = 3,             # number of Markov chains
    iter = 3000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 1500           # no progress shown
)

posterior <- extract(fit1)
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("Scenario A") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/appendix_",n,"_traceplot_scA.pdf"), width=22, height = 3)

theta.samples <- summary(fit1, par=c("theta0", "theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary

newgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
newgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
season.samples <- summary(fit1, par=c("alphaseason"), probs = c(0.05,0.5, 0.95))$summary

gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q2 <- gamma.samples[,1]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q2 <- theta.samples[,5]
theta.q3 <- theta.samples[,6]

# df.theta <- data.frame("seq" = seq(1, (p+1)),
#                        "true" = theta.adjusted,
#                        "m" = theta.q2,
#                        "l" = theta.q1,
#                        "u" = theta.q3)
# df.theta$covariate <- factor(0:p)
# df.theta$labels <- factor(0:p)
# ggplot(df.theta, aes(x = covariate, y=m, color = covariate)) + ylab("") + xlab('') +
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_point(size = 5) + geom_point(aes(y = true), color="red", size = 4) +
#   geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, linewidth =1.2) + 
#   scale_x_discrete(labels = c(expression(bold(theta[0])),
#                               expression(bold(theta[1])),
#                               expression(bold(theta[2])),
#                               expression(bold(theta[3])),
#                               expression(bold(theta[4])),
#                               expression(bold(theta[5])))) + 
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         legend.text.align = 0,
#         legend.title = element_blank(),
#         legend.text = element_text(size=25),
#         legend.margin=margin(0,0,0,-10),
#         legend.box.margin=margin(-10,0,-10,0),
#         plot.margin = margin(0,0,0,-20),
#         axis.text.x = element_text(hjust=0.35),
#         axis.text = element_text(size = 28))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_theta_sc1-wi.pdf"), width=10, height = 7.78)

# df.gamma <- data.frame("seq" = seq(1, (psi*p)), 
#                       #  "true" = as.vector(gamma.origin),
#                        "m" = as.vector(gamma.q2),
#                        "l" = as.vector(gamma.q1),
#                        "u" = as.vector(gamma.q3))
# df.gamma$covariate <- factor(rep(seq(1, 1 + nrow(df.gamma) %/% psi), each = psi, length.out = nrow(df.gamma)))
# df.gamma$labels <- factor(1:(psi*p))
# ggplot(df.gamma, aes(x =labels, y = m, color = covariate)) + 
#   geom_errorbar(aes(ymin = l, ymax = u), alpha = 0.4, width = 4, linewidth = 1.2) +
#   # geom_point(aes(y=true), size =4, color ="red")+
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_point(size = 4) + ylab("") + xlab("" ) + 
#   scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+7), 
#                    label = c(expression(bold(gamma[1])), 
#                              expression(bold(gamma[2])), 
#                              expression(bold(gamma[3])), 
#                              expression(bold(gamma[4])), 
#                              expression(bold(gamma[5]))),
#                    expand=c(0,3)) +
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#         legend.title = element_blank(),
#         legend.text = element_text(size=25),
#         legend.margin=margin(0,0,0,-10),
#         legend.box.margin=margin(-10,0,-10,0),
#         plot.margin = margin(0,0,0,-20),
#         axis.text.x = element_text(hjust=0.5),
#         axis.text = element_text(size = 28))
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

plot.alpha <- df_true_alpha[order(df_true_alpha$season, df_true_alpha$x), ]

n_grid <- length(newx)
n_vars <- 5
n_seasons <- 4
total_curves <- n_vars * n_seasons # 20 columns in the gridgsmooth matrix
var_labels <- rep(c("g1", "g2", "g3", "g4", "g5"), each = n_grid * n_seasons)
season_labels <- rep(rep(c("1", "2", "3", "4"), each = n_grid), times = n_vars)

data.smooth <- data.frame(
  "x" = rep(newx, total_curves),
  "post.mean" = as.vector(matrix(newgsmooth.samples[,1], nrow = n_grid, byrow=TRUE)), 
  "q1" = as.vector(matrix(newgsmooth.samples[,4], nrow = n_grid, byrow=TRUE)),
  "q2" = as.vector(matrix(newgsmooth.samples[,5], nrow = n_grid, byrow=TRUE)),
  "q3" = as.vector(matrix(newgsmooth.samples[,6], nrow = n_grid, byrow=TRUE)),
  "variable" = factor(var_labels, levels = c("g1", "g2", "g3", "g4", "g5")),
  "season" = factor(season_labels)
)

plot_data <- merge(data.smooth, 
                   df_true_smooth[, c("x", "season", "variable", "y")], 
                   by = c("x", "season", "variable"),
                   sort = FALSE)
names(plot_data)[names(plot_data) == 'y'] <- 'true'

season_levels <- levels(plot_data$season)

for(s in season_levels) {
  season_data <- subset(plot_data, season == s)
  plt <- ggplot(season_data, aes(x=x)) + 
    geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 1) + 
    geom_ribbon(aes(ymin = q1, ymax = q3), fill = "steelblue", alpha = 0.2) +
    geom_line(aes(y=true, colour = "True"), linewidth=1, linetype="dashed") + 
    geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
    facet_grid(variable ~ ., scales = "free_y", switch = "y", 
               labeller = label_parsed) + 
    scale_color_manual(name = "", values = c("True" = "red", "Posterior Median" = "steelblue")) +
    
    ylab("") + xlab("c") +    
    theme_minimal(base_size = 20) +
    theme(legend.position = "none",
          strip.text.y = element_text(angle = 0, face = "bold.italic"),
          strip.placement = "outside",
          plot.title = element_text(hjust = 0.5))
  print(plt)
}



alpha.mean <- as.vector(matrix(season.samples[,1], nrow = n, byrow=TRUE))
alpha.q1 <- as.vector(matrix(season.samples[,4], nrow = n, byrow=TRUE))
alpha.q2 <- as.vector(matrix(season.samples[,5], nrow = n, byrow=TRUE))
alpha.q3 <- as.vector(matrix(season.samples[,6], nrow = n, byrow=TRUE))

data.alpha <- data.frame("x" = newx,
                          "true" = df_true_alpha$alpha,
                          "post.mean" = alpha.mean,
                          "q1" = alpha.q1,
                          "q2" = alpha.q2,
                          "q3" = alpha.q3)
seasons <- c("Winter", "Spring", "Summer", "Autumn")

grid.plts <- list()
for(i in 1:4){
  grid.plt <- ggplot(data = data.frame(data.alpha[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=true), color = "red", linewidth=1, linetype=2) + 
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  # geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab(expression(alpha(c,ldots,c))) + xlab(seasons[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  theme_minimal(base_size = 20) +
                  theme(legend.position = "none",
                        plot.margin = margin(5, 5, 5, 5),
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        axis.text = element_text(size = 18),
                        axis.title.x = element_text(size = 22))
  grid.plts[[i]] <- grid.plt
}

# marrangeGrob(grobs = grid.plts, nrow = 2, ncol = 2)
final_plot <- patchwork::wrap_plots(grid.plts, nrow = 2, ncol = 2)
print(final_plot)

g_post <- posterior$gridgsmooth 
n_iter <- dim(g_post)[1]
n_grid <- dim(g_post)[2]
# g_global_post <- array(NA, dim = c(n_iter, n_grid, p))
# for(v in 1:p) {
#   col_indices <- ((v - 1) * n_seasons + 1):(v * n_seasons)
#   g_global_post[, , v] <- apply(g_post[, , col_indices], c(1, 2), mean)
# }

g_4d <- array(g_post, dim = c(n_iter, n_grid, n_seasons, p))
g_permuted <- aperm(g_4d, c(3, 1, 2, 4))
g_global_post <- colMeans(g_permuted, dims = 1)
global.mean <- apply(g_global_post, c(2, 3), mean)
global.q1   <- apply(g_global_post, c(2, 3), quantile, prob=0.05)
global.q2   <- apply(g_global_post, c(2, 3), median)
global.q3   <- apply(g_global_post, c(2, 3), quantile, prob=0.95)




data.smooth <- data.frame("x"=newx,
                          "true" =df_true_smooth$global,
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
  # geom_line(aes(y=post.mean, colour = "Posterior Median"), linewidth=1.8) + 
  ylab("") + xlab(expression(c)) +
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
        strip.text.y = element_text(size = 38, colour = "black", angle = 0, face = "bold.italic"),
        strip.placement = "outside",
        axis.title.x = element_text(size = 45),
        axis.text = element_text(size=30))

# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_smooth_scA.pdf"), width=12.5, height = 15)

# data.linear <- data.frame("x"=seq(0,1, length.out = n),
#                           "true" = as.vector(l.new),
#                           "post.mean" = as.vector(g.linear.mean),
#                           "q1" = as.vector(g.linear.q1),
#                           "q2" = as.vector(g.linear.q2),
#                           "q3" = as.vector(g.linear.q3),
#                           "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]")),
#                           "fakelab" = rep(1, (p*n)),
#                           "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

# ggplot(data.linear, aes(x=x, group=interaction(covariates, replicate))) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#   geom_line(aes(y=true, colour = "True"), linewidth=2, linetype=2) + 
#   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1.5) + 
#   # geom_line(aes(y=post.mean, colour = "Posterior Median"), linewidth=1) + 
#   ylab("") + xlab(expression(c)) +
#   facet_grid(covariates ~ ., scales = "free_x", switch = "y",
#               labeller = label_parsed) +  
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values=c("steelblue", "red")) + 
#   guides(color = guide_legend(order = 2), 
#          fill = guide_legend(order = 1)) + 
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#           plot.margin = margin(0,0,0,-20),
#           strip.text = element_blank(),
#           axis.title.x = element_text(size = 45),  
#           axis.text = element_text(size = 30))

# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_linear_scA.pdf"), width=12.5, height = 15)


# data.nonlinear <- data.frame("x"=seq(0,1, length.out=n),
#                              "true" = as.vector(nl.new),
#                              "post.mean" = as.vector(g.nonlinear.mean),
#                              "q1" = as.vector(g.nonlinear.q1),
#                              "q2" = as.vector(g.nonlinear.q2),
#                              "q3" = as.vector(g.nonlinear.q3),
#                              "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]")),
#                              "fakelab" = rep(1, (p*n)),
#                              "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

# ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#   geom_line(aes(y=true, colour = "True"), linewidth=2, linetype=2) + 
#   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1.5) + 
#   ylab("") + xlab(expression(c)) +  
#   facet_grid(covariates ~ ., scales = "free_x", switch = "y",
#               labeller = label_parsed) +  
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values=c("steelblue", "red")) + 
#   guides(color = guide_legend(order = 2), 
#          fill = guide_legend(order = 1)) + 
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#           plot.margin = margin(0,0,0,-20),
#           strip.text = element_blank(),
#           axis.title.x = element_text(size = 45),  
#           axis.text = element_text(size = 30))

# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_nonlinear_scA.pdf"), width=12.5, height = 15)

# data.scenario <- data.frame("x" = newx,
#                             "true" = (alp.new),
#                             "post.mean" = (newalpha.samples[,1]),
#                             "q2" = (newalpha.samples[,5]),
#                             "q1" = (newalpha.samples[,4]),
#                             "q3" = (newalpha.samples[,6]))

# ggplot(data.scenario, aes(x=x)) + 
#   ylab(expression(alpha(c,...,c))) + xlab(expression(c)) + labs(col = "") +
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#   geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2, linetype=2) + #ylim(0, 20) + 
#   geom_line(aes(y=q2, col = "Posterior Median"), linewidth=1.5) +
#   # geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=1.5) +
#   scale_color_manual(values=c("steelblue", "red")) + 
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   theme_minimal(base_size = 40) + 
#   theme(legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(size = 30))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_scA.pdf"), width=10, height = 7.78)

mcmc.alpha <- posterior$alpha
len <- dim(mcmc.alpha)[1]
r <- matrix(, nrow = n, ncol = 30)
T <- 30
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y.origin[i], u[i], alpha = mcmc.alpha[round(runif(1,1,len)),i]))
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
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_qqplot_scA.pdf"), width=10, height = 7.78)
