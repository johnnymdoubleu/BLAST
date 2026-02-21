# library(npreg)
library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(qqboxplot)
# library(qgam)
library(evgam)

# set.seed(1001)
n <- 10000
psi <- 10
threshold <- 0.95
p <- 5
# C <- matrix(0.2, nrow = p, ncol = p)
# diag(C) <- 1
C <- diag(p)

time.seq <- 1:n
period <- 365 
x.season <- time.seq / period 

# Convert continuous season to Factor for the 'by' argument
season_code_full <- cut((time.seq%% period / period), breaks = c(-0.1, 0.25, 0.5, 0.75, 1.1), labels = c(1,2,3,4))
seasons <- c("Winter", "Spring", "Summer", "Autumn")
f.season.scale <- function(x) {
  return(1.5  - 0.8 * sin(2 * pi * x) - .6 * cos(2 * pi * x)) 
}

x.random <- pnorm(matrix(rnorm(n*p, sd=0.5), ncol = p) %*% chol(C) + f.season.scale(x.season))
x.origin <- cbind(x.random, season_code_full)
# plot(x.origin[,1])
covariates <- colnames(data.frame(x.origin))[1:p] 

x.linear <- matrix(0, nrow = n, ncol = p * 4)
col_names <- c()

for (j in 1:p) {
  for (s in 1:4) {
    idx <- (j - 1) * 4 + s
    # Only keep the covariate value for the active season
    mask <- as.numeric(season_code_full == s)
    x.linear[, idx] <- x.origin[, j] * mask
    col_names <- c(col_names, paste0("V", j, "_S", s))
  }
}
colnames(x.linear) <- col_names

theta0.origin <- 0.2 # Global Intercept

theta.linear.origin <- matrix(c(
  # S1,   S2,   S3,   S4
   0.0,  0.0,  0.0,  0.0,  # X1: Pure noise, completely 0 everywhere
   0.2,  0.2,  0.2,  0.2,  # X2: Linear effect gets stronger throughout the year
   0.3,  0.3,  -0.3,  -0.3,  # X3: Constant negative effect across all seasons
   0.0,  0.0,  -0.5,  -0.5,  # X4: Only active in Summer & Autumn in linear base
   0.0,  0.0,  0.0,  0.0   # X5: Pure noise, completely 0 everywhere
), nrow = p, ncol = 4, byrow = TRUE)

f2 <- function(x) { .2 * sin(2 * pi * x^2) * x }
f3 <- function(x) { -.2 * cos(3 * pi * x^2) * x}
alp.linear.predictor <- rep(theta0.origin, n)

for (i in 1:n) {
  s <- as.numeric(season_code_full[i])
  x_vec <- x.origin[,c(1:p)]  
  alp.linear.predictor[i] <- alp.linear.predictor[i] + sum(x.linear[i, ] * theta.linear.origin[, s])
  alp.linear.predictor[i] <- alp.linear.predictor[i] + f2(x_vec[2])
  alp.linear.predictor[i] <- alp.linear.predictor[i] + f3(x_vec[3])
}

alp.origin <- exp(alp.linear.predictor)

# theta.origin <- c(0.4, 0, 0.1, 0.8, 0, 0) # Coefs for 5 random vars + intercept
# alp.origin <- as.vector(exp(
#   rep(theta.origin[1],n) + 
#   x.origin[,1:p] %*% theta.origin[-1] + 
#   f2(x.origin[,2]) + 
#   f3(x.origin[,3])
# ))

y.noise <- rPareto(n, rep(1,n), alpha = alp.origin)
y.origin <- y.noise * f.season.scale(x.season)
plot(y.origin)


evgam.df <- data.frame(
  y = y.origin,
  sin.time = sin(2 * pi * time.seq / 365),
  cos.time = cos(2 * pi * time.seq / 365)
)
evgam.cov <- y ~ 1 + cos.time + sin.time
ald.cov.fit <- evgam(evgam.cov, data = evgam.df, family = "ald", ald.args=list(tau = threshold))
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
x.origin <- data.frame(x.origin)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
x.scaled <- sapply(x.origin[,c(1:p)], FUN = range01)
# x.scaled <- x.origin
x.linear <- matrix(0, nrow = n, ncol = p * 4)
col_names <- c()
for (j in 1:p) {
  for (s in 1:4) {
    idx <- (j - 1) * 4 + s
    mask <- as.numeric(season_code_full == s)
    x.linear[, idx] <- x.scaled[, j] * mask
    col_names <- c(col_names, paste0("V", j, "_S", s))
  }
}
colnames(x.linear) <- col_names

x.minmax <- sapply(x.origin[,c(1:p)], function(x) max(x)-min(x))
x.min <- sapply(x.origin[,c(1:p)], function(x) min(x))

linear_basis_list <- list()
nonlinear_basis_list <- list()
model_data <- list() 
Z_scales <- c()
psi <- psi - 2
for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  x_vec <- x.scaled[, i]
  
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
    X_lin_local <- cbind(vec_intercept, vec_slope)
    
    qr_lin <- qr(X_lin_local)
    Q_lin <- qr.Q(qr_lin)
    R_lin <- qr.R(qr_lin)
    X_raw <- sm_spec$X
    S     <- sm_spec$S[[1]] 
    
    eig <- eigen(S, symmetric = TRUE)
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
    Z_scales <- c(Z_scales, train_scale)
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

bs.linear <- do.call(cbind, linear_basis_list) 
bs.nonlinear <- do.call(cbind, nonlinear_basis_list)
# bs.linear_check <- cbind(rep(1,n), bs.linear)
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,2,3,4,5)]) %*% bs.nonlinear[,c((1):(4*psi))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,6,7,8,9)]) %*% bs.nonlinear[,c((4*(psi+1)):(4*(psi*2)))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,10,11,12,13)]) %*% bs.nonlinear[,c((4*(psi*2+1)):(4*(psi*3)))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,14,15,16,17)]) %*% bs.nonlinear[,c((4*(psi*3+1)):(4*(psi*4)))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,18,19,20,21)]) %*% bs.nonlinear[,c((4*(psi*4+1)):(4*(psi*5)))]), "\n")

grid_linear_list <- list()
grid_nonlinear_list <- list()
raw_grid_list <- list()

for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  global.min <- x.min[i]
  global.range <- x.minmax[i]
  for (season_label in names(model_data[[var_name]])) {
    params <- model_data[[var_name]][[season_label]]
    lin_col_name <- paste(var_name, "S", season_label, sep="_")
    season_idx <- which(season_code_full == season_label)
    raw_season_vals <- x.origin[season_idx, i]
    season_min <- min(raw_season_vals, na.rm = TRUE)
    season_max <- max(raw_season_vals, na.rm = TRUE)
    raw_grid <- seq(season_min, season_max, length.out = n)
    raw_grid_list[[lin_col_name]] <- raw_grid
    grid.vals <- (raw_grid - global.min) / global.range
    grid_linear_list[[lin_col_name]] <- matrix(grid.vals, ncol=1, dimnames=list(NULL, lin_col_name))
    
    pred_df <- data.frame(
      x_vec = grid.vals,
      season_code_full = factor(season_label, levels = levels(season_code_full))
    )
    
    X_grid.vals <- PredictMat(params$sm_spec, pred_df)
    Z_spectral_grid <- X_grid.vals %*% params$U_pen %*% params$Lambda_sqrt_inv
    intercept_grid <- rep(1, length(grid.vals))
    slope_grid     <- grid.vals
    X_lin_grid_local <- cbind(intercept_grid, slope_grid)
    Z_orth_grid <- Z_spectral_grid - X_lin_grid_local %*% params$projection_coefs
    Z_final_grid <- Z_orth_grid[, params$keep_cols, drop = FALSE]
    Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)
    col_names <- paste(var_name, "S", season_label, 1:ncol(Z_final_grid), sep="_")
    colnames(Z_final_grid) <- col_names
    grid_nonlinear_list[[paste(var_name, season_label, sep="_")]] <- Z_final_grid
  }
}

xholder.linear <- do.call(cbind, grid_linear_list)
xholder.nonlinear <- do.call(cbind, grid_nonlinear_list)

theta.adjusted.mat <- theta.linear.origin
hidden_intercept_sum <- 0
hidden_funcs <- list()
for (s in 1:4) {
  season_idx <- which(season_code_full == s)
  x_season <- x.origin[season_idx, ]
  hidden_funcs[[s]] <- list()
  hidden_f2 <- make.nl(x_season[, 2], f2(x_season[, 2]))
  theta.adjusted.mat[2, s] <- theta.adjusted.mat[2, s] + hidden_f2$slope
  hidden_intercept_sum <- hidden_intercept_sum + (hidden_f2$intercept * length(season_idx))
  hidden_funcs[[s]]$f2 <- hidden_f2
  hidden_f3 <- make.nl(x_season[, 3], f3(x_season[, 3]))
  theta.adjusted.mat[3, s] <- theta.adjusted.mat[3, s] + hidden_f3$slope
  hidden_intercept_sum <- hidden_intercept_sum + (hidden_f3$intercept * length(season_idx))
  hidden_funcs[[s]]$f3 <- hidden_f3
}
theta0.adjusted <- theta0.origin + (hidden_intercept_sum / n)

true_gridgsmooth <- matrix(0, nrow = n, ncol = p)
true_gsmoothseason <- matrix(0, nrow = n, ncol = p * 4)

for (j in 1:p) {
  for (s in 1:4) {
    lin_idx <- (j - 1) * 4 + s
    col_name <- paste0("V", j, "_S_", s)
    x_grid <- raw_grid_list[[col_name]]
    g_l <- theta.adjusted.mat[j, s] * x_grid
    g_nl <- rep(0, n)
    if (j == 2) {
      g_nl <- f2(x_grid) - (hidden_funcs[[s]]$f2$intercept + hidden_funcs[[s]]$f2$slope * x_grid)
    } 
    if (j == 3) {
      g_nl <- f3(x_grid) - (hidden_funcs[[s]]$f3$intercept + hidden_funcs[[s]]$f3$slope * x_grid)
    }
    true_gsmoothseason[, lin_idx] <- g_l + g_nl
    true_gridgsmooth[, j] <- true_gridgsmooth[, j] + true_gsmoothseason[, lin_idx]
  }
}

true_gridalpha <- exp(theta0.adjusted + rowSums(true_gridgsmooth))
plot(true_gridalpha)
df_true_gridgsmooth <- do.call(rbind, lapply(1:p, function(j) {
  data.frame(
    x_index = seq(0, 1, length.out = n), 
    y = true_gridgsmooth[, j],
    variable = paste0("Covariate_", j)
  )
}))

true_alphaseason <- matrix(0, nrow = n, ncol = 4)

for (s in 1:4) {
  current_season_pred <- rep(theta0.adjusted, n)
  for (j in 1:p) {
    lin_idx <- (j - 1) * 4 + s
    current_season_pred <- current_season_pred + true_gsmoothseason[, lin_idx]
  }
  true_alphaseason[, s] <- exp(current_season_pred)
}


df_true_alphaseason <- do.call(rbind, lapply(1:4, function(s) {
  data.frame(
    x_index = seq(0, 1, length.out = n), 
    alpha = true_alphaseason[, s],
    season = seasons[s]
  )
}))

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
    // vector[n] lin_predictor;
    alpha = rep_vector(theta0, n); 
    for (j in 1:p){
      for (s in 1:n_seasons){
        int lin_idx = (j-1) * n_seasons + s;
        int nl_start = (lin_idx-1) * psi + 1; 
        
        for (k in 1:psi){
          int flat_idx = nl_start + k - 1;
          gamma[j,s][k] = gamma_raw[j,s][k] * sqrt(tau[j,s]) * Z_scales[flat_idx];
        }
        alpha += col(bsLinear, lin_idx) * theta[lin_idx];
        alpha += block(bsNonlinear, 1, nl_start, n, psi) * gamma[j,s];
      }
    }
    alpha = exp(alpha);
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
  // vector[n] gridalpha;
  vector[n] log_lik;
  // matrix[n, p] gridgsmooth;               // Combined seasonal effects per covariate
  matrix[n, p * n_seasons] gsmoothseason; // Individual season-specific smooths

  vector[p * n_seasons] theta_origin = theta ./ X_sd;
  // gridgsmooth = rep_matrix(0.0, n, p);
  // vector[n] grid_vec = rep_vector(theta0, n);
  for (j in 1:p) {
    for (s in 1:n_seasons) {
      int lin_idx = (j-1) * n_seasons + s;
      int nl_start = (lin_idx-1) * psi + 1;
      gsmoothseason[, lin_idx] = col(xholderLinear, lin_idx) * theta_origin[lin_idx] 
                               + block(xholderNonlinear, 1, nl_start, n, psi) * gamma[j,s];
      // gridgsmooth[, j] += gsmoothseason[, lin_idx];
    }
    // grid_vec += gridgsmooth[, j];
  }
  // gridalpha = exp(grid_vec);

  for (s in 1:n_seasons) {
    vector[n] current_season_pred = rep_vector(theta0, n);
    for (j in 1:p) {
      int lin_idx = (j-1) * n_seasons + s;
      current_season_pred += gsmoothseason[, lin_idx];
    }
    alphaseason[, s] = exp(current_season_pred);
  }

  for (i in 1:n) {
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
    iter = 2000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 1000           # no progress shown
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
# newgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
# newgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
season.samples <- summary(fit1, par=c("alphaseason"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary

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

# g.linear.mean <- as.vector(matrix(newgl.samples[,1], nrow = n, byrow=TRUE))
# g.linear.q1 <- as.vector(matrix(newgl.samples[,4], nrow = n, byrow=TRUE))
# g.linear.q2 <- as.vector(matrix(newgl.samples[,5], nrow = n, byrow=TRUE))
# g.linear.q3 <- as.vector(matrix(newgl.samples[,6], nrow = n, byrow=TRUE))
# g.nonlinear.mean <- as.vector(matrix(newgnl.samples[,1], nrow = n, byrow=TRUE))
# g.nonlinear.q1 <- as.vector(matrix(newgnl.samples[,4], nrow = n, byrow=TRUE))
# g.nonlinear.q2 <- as.vector(matrix(newgnl.samples[,5], nrow = n, byrow=TRUE))
# g.nonlinear.q3 <- as.vector(matrix(newgnl.samples[,6], nrow = n, byrow=TRUE))
g.smooth.mean <- as.vector(matrix(newgsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(newgsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(newgsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(newgsmooth.samples[,6], nrow = n, byrow=TRUE))


est_alphaseason_median  <- apply(posterior$alphaseason, c(2, 3), median)
est_alphaseason_lower <- apply(posterior$alphaseason, c(2, 3), quantile, probs = 0.05)
est_alphaseason_upper <- apply(posterior$alphaseason, c(2, 3), quantile, probs = 0.95)

# Format into a plotting dataframe
df_est_alphaseason <- do.call(rbind, lapply(1:4, function(s) {
  data.frame(
    x_index = seq(0, 1, length.out = n),
    alpha_est = est_alphaseason_median[, s],
    lower_ci = est_alphaseason_lower[, s],
    upper_ci = est_alphaseason_upper[, s],
    season = seasons[s]
  )
}))
plot_data_alphaseason <- merge(df_est_alphaseason, df_true_alphaseason, by = c("season", "x_index"))

grid.plts <- list()
for(i in 1:4){
  grid.plt <- ggplot(plot_data_alphaseason[((((i-1)*n)+1):(i*n)),], aes(x = x_index)) +
          geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "steelblue", alpha = 0.3) +
          geom_line(aes(y = alpha_est, color = "Estimated (Posterior Mean)"), linewidth = 1) +
          geom_line(aes(y = alpha, color = "True (DGP)"), linewidth = 1, linetype = "dashed") +
          scale_color_manual(values = c("Estimated (Posterior Mean)" = "steelblue", "True (DGP)" = "red")) +
          labs(
            x = seasons[i],
            y = expression(alpha(c,ldots,c))
          ) +
          theme_minimal(base_size = 20) +
          theme(legend.position = "none",
                plot.margin = margin(5, 5, 5, 5),
                plot.title = element_text(hjust = 0.5, face = "bold"),
                axis.text = element_text(size = 18),
                axis.title.x = element_text(size = 22))  
  
  grid.plts[[i]] <- grid.plt
}

final_plot <- patchwork::wrap_plots(grid.plts, nrow = 2, ncol = 2)
print(final_plot)


# 2. Extract summary from Stan (already contains all seasons concatenated)
# Stan orders these by Season 1, then Season 2, etc., based on your alphaseason matrix
alpha_sum <- as.data.frame(season.samples)
alpha.mean <- as.vector(matrix(season.samples[,1], nrow = n, byrow=TRUE))
alpha.q1 <- as.vector(matrix(season.samples[,4], nrow = n, byrow=TRUE))
alpha.q2 <- as.vector(matrix(season.samples[,5], nrow = n, byrow=TRUE))
alpha.q3 <- as.vector(matrix(season.samples[,6], nrow = n, byrow=TRUE))
grid.plts <- list()
current_start <- 1

for(i in 1:length(season_names)){
  n_seas <- season_counts[i]
  current_end <- current_start + n_seas - 1
  
  plot_df <- data.frame(
    x = newx, # Mapping to a 0-1 grid for the plot
    true = true_gridalpha[((((i-1)*n)+1):(i*n))],
    q2 = alpha.q2[((((i-1)*n)+1):(i*n))],
    q1 = alpha.q1[((((i-1)*n)+1):(i*n))],
    q3 = alpha.q3[((((i-1)*n)+1):(i*n))]
  )
  
  grid.plt <- ggplot(plot_df, aes(x=x)) + 
    geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
    geom_line(aes(y=true), color = "red", linewidth=1, linetype=2) + 
    geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
    ylab(expression(alpha(c,ldots,c))) + 
    xlab(labels_pretty[i]) +
    scale_fill_manual(values=c("steelblue"), name = "") + 
    scale_color_manual(values=c("steelblue")) +
    theme_minimal(base_size = 20) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 18),
          axis.title.x = element_text(size = 22))
  grid.plts <- list()
  grid.plts[[i]] <- grid.plt
  current_start <- current_end + 1
}

final_plot <- patchwork::wrap_plots(grid.plts, nrow = 2, ncol = 2)
print(final_plot)




data.alpha <- data.frame("x" = newx,
                          "true" = df_true_alphaseason$alpha,
                          "post.mean" = alpha.mean,
                          "q1" = alpha.q1,
                          "q2" = alpha.q2,
                          "q3" = alpha.q3)


grid.plts <- list()
for(i in 1:4){
  grid.plt <- ggplot(data = data.frame(data.alpha[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  # geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=true), color = "red", linewidth=1, linetype=2) + 
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
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





gsmooth_samples <- posterior$gridgsmooth 
summed_effects_list <- list()

for (j in 1:p) {
  col_indices <- ((j-1) * 4 + 1) : (j * 4)
  if(length(col_indices) > 1) {
    var_j_summed <- apply(gsmooth_samples[, , col_indices], c(1, 2), sum)
  } else {
    var_j_summed <- gsmooth_samples[, , col_indices]
  }
  
  # Calculate summary statistics for the summed effect
  summed_effects_list[[j]] <- data.frame(
    x = newx,
    covariate = paste0("g", j),
    q2 = apply(var_j_summed, 2, median),
    q1 = apply(var_j_summed, 2, quantile, probs = 0.05),
    q3 = apply(var_j_summed, 2, quantile, probs = 0.95)
  )
}

# 3. Concatenate them in order
df.smooth <- do.call(rbind, summed_effects_list)


data.smooth <- data.frame("x"=newx,
                          "true" =df_true_smooth$global[1:(n*p)],
                          "q2"        = df.smooth$q2,
                          "q1"        = df.smooth$q1,
                          "q3"        = df.smooth$q3,
                          "covariates" = gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]")),
                          "replicate" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=true, colour = "True"), linewidth=2, linetype=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1.5) + 
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

est_gridalpha_median  <- apply(posterior$gridalpha, 2, median)
est_gridalpha_lower <- apply(posterior$gridalpha, 2, quantile, probs = 0.05)
est_gridalpha_upper <- apply(posterior$gridalpha, 2, quantile, probs = 0.95)

data.scenario <- data.frame("x" = newx,
                            "true" = true_gridalpha,
                            "q2" = est_gridalpha_median,
                            "q1" = est_gridalpha_lower,
                            "q3" = est_gridalpha_upper)
                            # "post.median" = (alpha.samples[,5]),
                            # "q1" = (alpha.samples[,4]),
                            # "q3" = (alpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,...,c))) + xlab(expression(c)) +
  # geom_ribbon(aes(ymin = q1, ymax = q3), fill="steelblue", alpha = 0.2) +
  geom_line(aes(y=true), color = "red", linewidth=1, linetype =2) +
  geom_line(aes(y=q2), color = "steelblue", linewidth=1) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))


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
