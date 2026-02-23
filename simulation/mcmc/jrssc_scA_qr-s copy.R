# library(npreg)
library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(qqboxplot)
# library(qgam)
library(evgam)

n <- 5000
psi <- 10
threshold <- 0.95
p <- 5
C <- diag(p)
psi <- psi - 2

time.seq <- 1:n
period <- 365 
x.season <- time.seq / period 

# Convert continuous season to Factor for the 'by' argument
season_code_full <- cut((time.seq%% period / period), breaks = c(-0.1, 0.25, 0.5, 0.75, 1.1), labels = c(1,2,3,4))
seasons <- c("Winter", "Spring", "Summer", "Autumn")

x_raw <- matrix(rnorm(n*p), ncol = p) %*% chol(C)
for(j in 1:p) {
  x_raw[,j] <- x_raw[,j] + 0.8 - 0.8 * sin(2 * pi * x.season) + (j*0.2) 
}
x.random <- 1-pnorm(x_raw)
x.origin <- cbind(x.random, season_code_full)
plot(x.origin[,3])
covariates <- colnames(data.frame(x.origin))[1:p] 
x.linear <- c()
for (i in seq_along(covariates)) {
  x_vec <- x.origin[, i]
  for (j in 1:4) {
    season_label <- levels(season_code_full)[j]
    mask <- as.numeric(season_code_full == season_label)
    vec_intercept <- mask
    vec_slope <- x_vec * mask
    x.linear <- cbind(x.linear, matrix(vec_slope, ncol = 1))
  }
}

# x.linear <- do.call(cbind, linear_basis_list)

set.seed(7071)
true.theta0 <- 0.7          # Increased baseline (was 0.7)
theta.linear.origin <- matrix(c(
  # S1,    S2,    S3,    S4
   0.0,   0.0,   0.0,   0.0,   # X1: Noise
   0.8,   0.6,   0.4,   0.3,   # X2: Stronger slopes (was 0.2,0.15,0.05,0.1)
   0.4,   0.4,  -0.2,  -0.2,   # X3: Stronger effects (was 0.1,0.1,-0.05,-0.05)
   0.0,   0.0,  -0.6,  -0.6,   # X4: Stronger Summer drop (was -0.15)
   0.0,   0.0,   0.0,   0.0    # X5: Noise
), nrow = p, ncol = 4, byrow = TRUE)

true.theta <- c(t(theta.linear.origin))

# =====================================================
# SCALED NONLINEAR FUNCTIONS (reduced amplitude)
# =====================================================
nl_func_X2 <- function(x) { 0.8 * sin(2 * pi * x) + 0.6 * (x - 0.5)^2 }  # Amplified (was 0.3,0.2)
nl_func_X3 <- function(x) {-0.5 * cos(3 * pi * x) - 0.3 * x^3 }          # Amplified (was 0.15,0.1)
true.theta <- c(t(theta.linear.origin))  # Keep linear parts
eta <- rep(true.theta0, n) + as.vector(x.linear %*% true.theta)
eta <- eta + nl_func_X2(x.origin[,2])
eta <- eta + nl_func_X3(x.origin[,3])
# mask_X2 <- gamma.mask[((2-1)*4*psi + 1) : (2*4*psi)] == 1  # Only X2 seasons 2,3
# mask_X3 <- gamma.mask[((3-1)*4*psi + 1) : (3*4*psi)] == 1  # Only X3 seasons 2,3

# # Apply nonlinear to relevant seasons (S2,S3 for X2/X3)
# for(season in 2:3) {
#   eta[mask_season] <- eta[mask_season] + nl_func_X2(x.origin[mask_season, 2])
#   eta[mask_season] <- eta[mask_season] + nl_func_X3(x.origin[mask_season, 3])
# }

alp.origin <- exp(eta)
y.noise <- rPareto(n, rep(1, n), alpha = alp.origin)
f.season.scale <- function(x) {
  return(1.5 - 0.8 * sin(2 * pi * x) - .6 * cos(2 * pi * x)) 
}
y.origin <- y.noise * f.season.scale(x.season)
plot(y.origin)

evgam.df <- data.frame(
  y = (y.origin),
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

x.origin <- data.frame(x.origin)

linear_basis_list <- list()
nonlinear_basis_list <- list()
model_data <- list() 

for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  x_vec <- x.origin[, i]  # or full data for generation
  
  model_data[[var_name]] <- list() 
  sm_list <- smoothCon(mgcv::s(x_vec, by=season_code_full, bs = "tp", k = psi + 2), 
                       data = data.frame(x_vec = x_vec, season_code_full = season_code_full))
  
  for (j in 1:length(sm_list)) {
    season_label <- levels(season_code_full)[j]
    sm_spec <- sm_list[[j]] 
    mask <- as.numeric(season_code_full == season_label)
    
    # FIXED QR: Subset to ACTIVE season rows only
    mask_idx <- which(mask == 1)
    x_vec_season <- x_vec[mask_idx]
    X_lin_season <- cbind(rep(1, length(mask_idx)), x_vec_season)
    
    qr_lin <- qr(X_lin_season)
    
    # Full matrices for spectral basis
    X_raw <- sm_spec$X
    S <- sm_spec$S[[1]] 
    eig <- eigen(S, symmetric = TRUE)
    pos_idx <- 1:psi  # Or adaptive for generation
    
    U_pen <- eig$vectors[, pos_idx]      
    Lambda_sqrt_inv <- diag(1/sqrt(eig$values[pos_idx]))
    Z_spectral <- X_raw %*% U_pen %*% Lambda_sqrt_inv
    
    # FIXED: Stable projection using SEASON-ONLY QR
    Gamma_Original <- backsolve(qr.R(qr_lin), t(qr.Q(qr_lin)) %*% Z_spectral[mask_idx, , drop = FALSE])
    
    # Project back to FULL data
    X_lin_full <- cbind(mask, x_vec * mask)
    Z_orth <- Z_spectral - X_lin_full %*% Gamma_Original
    
    keep_cols <- colSums(Z_orth^2) > 1e-9
    Z_final <- Z_orth[, keep_cols, drop = FALSE]

    # Centering/scaling on season data
    train_means <- apply(Z_final[mask == 1, , drop = FALSE], 2, mean)
    Z_final <- sweep(Z_final, 2, train_means, "-")
    
    train_scale <- apply(Z_final[mask == 1, , drop = FALSE], 2, sd)
    train_scale[train_scale < 1e-12] <- 1 
    Z_final <- scale(Z_final, center = FALSE, scale = train_scale) * mask
    
    # Store for grid prediction
    model_data[[var_name]][[season_label]] <- list(
      sm_spec = sm_spec, U_pen = U_pen, Lambda_sqrt_inv = Lambda_sqrt_inv,
      projection_coefs = Gamma_Original, keep_cols = keep_cols,
      center_stats = train_means, scale_stats = train_scale
    )
    
    # STANDARDIZED names
    lin_name <- paste0(var_name, "_S", season_label)
    nl_name <- paste(var_name, season_label, sep = "_")
    linear_basis_list[[lin_name]] <- matrix(x_vec * mask, ncol = 1)
    nonlinear_basis_list[[nl_name]] <- Z_final
  }
}

bs.linear <- do.call(cbind, linear_basis_list)
bs.nonlinear <- do.call(cbind, nonlinear_basis_list)

# grid_linear_list <- list()
# grid_nonlinear_list <- list()
# grid_n <- 200

# for (i in seq_along(covariates)) {
#   var_name <- covariates[i]
#   for (season_label in names(model_data[[var_name]])) {
#     params <- model_data[[var_name]][[season_label]]
#     lin_col_name <- paste(var_name, "S", season_label, sep="_")
    
#     season.idx <- which(season_code_full == season_label)
#     raw.vals <- x.origin[season.idx, i]
#     raw.grid <- seq(min(raw.vals), max(raw.vals), length.out = grid_n)
#     # raw.grid <- seq(0,1,length.out = grid_n)
#     grid_linear_list[[lin_col_name]] <- matrix(raw.grid, ncol=1, dimnames=list(NULL, lin_col_name))
    
#     pred_df <- data.frame(
#       x_vec = raw.grid,
#       season_code_full = factor(season_label, levels = levels(season_code_full))
#     )
    
#     X_raw.grid <- PredictMat(params$sm_spec, pred_df)
#     Z_spectral_grid <- X_raw.grid %*% params$U_pen %*% params$Lambda_sqrt_inv
#     X_lin_grid_local <- cbind(rep(1, grid_n), raw.grid)
#     Z_orth_grid <- Z_spectral_grid - X_lin_grid_local %*% params$projection_coefs
    
#     Z_final_grid <- Z_orth_grid[, params$keep_cols, drop = FALSE]
#     Z_final_grid <- sweep(Z_final_grid, 2, params$center_stats, "-")
#     Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)    
#     # Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)
    
#     col_names <- paste(var_name, "S", season_label, 1:ncol(Z_final_grid), sep="_")
#     colnames(Z_final_grid) <- col_names
#     grid_nonlinear_list[[paste(var_name, season_label, sep="_")]] <- Z_final_grid
#   }
# }

# xholder.linear <- do.call(cbind, grid_linear_list)
# xholder.nonlinear <- do.call(cbind, grid_nonlinear_list)
grid_linear_list <- list()
grid_nonlinear_list <- list()
grid_n <- 200
newx <- seq(0, 1, length.out = grid_n)
for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  for (season_label in names(model_data[[var_name]])) {
    params <- model_data[[var_name]][[season_label]]
    lin_col_name <- paste(var_name, "S", season_label, sep = "_")
    
    season.idx <- which(season_code_full == season_label)
    raw.vals <- x.origin[season.idx, i]
    raw.grid <- seq(min(raw.vals), max(raw.vals), length.out = grid_n)
    grid_linear_list[[lin_col_name]] <- matrix(raw.grid, ncol = 1, dimnames = list(NULL, lin_col_name))
    
    # FIXED: Pure season prediction (NO mask)
    pred_df <- data.frame(
      x_vec = raw.grid,
      season_code_full = factor(season_label, levels = levels(season_code_full))
    )
    
    X_raw.grid <- PredictMat(params$sm_spec, pred_df)
    Z_spectral_grid <- X_raw.grid %*% params$U_pen %*% params$Lambda_sqrt_inv
    
    X_lin_grid <- cbind(rep(1, grid_n), raw.grid)  # Pure season linear space
    Z_orth_grid <- Z_spectral_grid - X_lin_grid %*% params$projection_coefs
    
    Z_final_grid <- Z_orth_grid[, params$keep_cols, drop = FALSE]
    Z_final_grid <- sweep(Z_final_grid, 2, params$center_stats, "-")
    Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)
    # NO * mask - grid represents 100% season data!
    
    colnames(Z_final_grid) <- paste(var_name, season_label, 1:ncol(Z_final_grid), sep = "_")
    grid_nonlinear_list[[paste(var_name, season_label, sep = "_")]] <- Z_final_grid
  }
}

xholder.linear <- do.call(cbind, grid_linear_list)
xholder.nonlinear <- do.call(cbind, grid_nonlinear_list)



bs.linear_check <- cbind(rep(1,n), bs.linear)
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,2,3,4,5)]) %*% bs.nonlinear[,c((1):(4*psi))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,6,7,8,9)]) %*% bs.nonlinear[,c((4*(psi+1)):(4*(psi*2)))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,10,11,12,13)]) %*% bs.nonlinear[,c((4*(psi*2+1)):(4*(psi*3)))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,14,15,16,17)]) %*% bs.nonlinear[,c((4*(psi*3+1)):(4*(psi*4)))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear_check[,c(1,18,19,20,21)]) %*% bs.nonlinear[,c((4*(psi*4+1)):(4*(psi*5)))]), "\n")

true_alphaseason <- matrix(0, nrow = grid_n, ncol = 4)
true_gsmoothseason <- matrix(0, nrow = grid_n, ncol = p * 4)

for (s in 1:4) {
  current_season_pred <- rep(true.theta0, grid_n)
  for (j in 1:p) {
    lin_col_name <- paste(covariates[j], "S", s, sep="_")
    lin_idx <- (j-1)*4 + s
    current_season_pred <- current_season_pred + xholder.linear[, lin_col_name] * true.theta[lin_idx]
  }
  x_vals_X2 <- xholder.linear[, paste("V2_S", s, sep="_")]  # X2's OWN grid for this season
  current_season_pred <- current_season_pred + nl_func_X2(x_vals_X2)
  x_vals_X3 <- xholder.linear[, paste("V3_S", s, sep="_")]  # X3's OWN grid for this season
  current_season_pred <- current_season_pred + nl_func_X3(x_vals_X3)  
  true_alphaseason[, s] <- exp(current_season_pred)
}



true_gridgsmooth <- matrix(0, nrow = grid_n, ncol = p)

for (j in 1:p) {
  x_global <- newx  # Universal [0,1] grid
  
  # Linear average
  season_cols <- ((j-1)*4 + 1):(j*4)
  linear_avg <- rowMeans(xholder.linear[, season_cols, drop=FALSE] * true.theta[season_cols])
  nl_avg <- rep(0, grid_n)
  if (j == 2) nl_avg <- nl_func_X2(x_global)
  if (j == 3) nl_avg <- nl_func_X3(x_global)
  
  true_gridgsmooth[, j] <- linear_avg + nl_avg
}

df_true_gridgsmooth <- do.call(rbind, lapply(1:p, function(j) {
  data.frame(
    x_index = xholder.linear[,(j-1)*4+1], 
    y = true_gridgsmooth[, j],
    variable = paste0("Covariate_", j)
  )
}))

# true_alphaseason <- matrix(0, nrow = grid_n, ncol = 4)

# for (s in 1:4) {
#   current_season_pred <- rep(true.theta0, grid_n)
#   for (j in 1:p) {
#     lin_idx <- (j - 1) * 4 + s
#     current_season_pred <- current_season_pred + true_gsmoothseason[, lin_idx]
#   }
#   true_alphaseason[, s] <- exp(current_season_pred)
# }

true_gridalpha <- rowMeans(true_alphaseason)
plot(true_gridalpha)

X_sd   <- apply(bs.linear, 2, sd)
bs.linear <- scale(bs.linear, center = FALSE, scale = X_sd)
true.theta.scaled <- true.theta * X_sd


model.stan <- "
// Stan model for Seasonal Pareto Samples
data {
  int <lower=1> n; 
  int <lower=1> grid_n;
  int <lower=1> p; 
  int <lower=1> psi; 
  int <lower=1> n_seasons;
  matrix[n, (p * n_seasons)] bsLinear; 
  matrix[n, (psi * p * n_seasons)] bsNonlinear; 
  matrix[grid_n, (p * n_seasons)] xholderLinear; 
  matrix[grid_n, (psi * p * n_seasons)] xholderNonlinear; 
  vector[n] u; 
  vector[n] y; 
  real <lower=0> atau;
  // vector[(psi*p*n_seasons)] Z_scales;
  vector[p * n_seasons] X_sd;
}

parameters {
  real theta0; // Seasonal intercepts
  vector[p * n_seasons] theta; 
  array[p, n_seasons] vector[psi] gamma_raw;
  array[p, n_seasons] real <lower=0> lambda1; 
  array[p, n_seasons] real <lower=0> lambda2; 
  array[p, n_seasons] real <lower=0> tau;
}

transformed parameters {
  vector[n] alpha;         
  array[p, n_seasons] vector[psi] gamma; 
  
  // 1. Initialize with seasonal intercepts based on the data index
  alpha = rep_vector(theta0, n);

  // 2. Add linear and nonlinear effects
  for (j in 1:p) {
    for (s in 1:n_seasons) {
      int lin_idx = (j-1) * n_seasons + s;
      int nl_start = (lin_idx-1) * psi + 1; 
      
      gamma[j,s] = gamma_raw[j,s] * sqrt(tau[j,s]);
      alpha += col(bsLinear, lin_idx) * theta[lin_idx];
      alpha += block(bsNonlinear, 1, nl_start, n, psi) * gamma[j,s];
    }
  }
  alpha = exp(alpha);
}

model {
  // Likelihood
  target += pareto_lpdf(y | u, alpha);

  // Priors
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
  matrix[grid_n, n_seasons] alphaseason;
  vector[grid_n] gridalpha; 
  matrix[grid_n, p] gridgsmooth; // Marginal smooths per covariate              
  matrix[grid_n, p * n_seasons] gsmoothseason; 
  vector[p*n_seasons] theta_origin = theta ./ X_sd;

  // Compute smooth effects per variable per season
  for (j in 1:p) {
    gridgsmooth[, j] = rep_vector(0, grid_n);
    for (s in 1:n_seasons) {
      int lin_idx = (j-1) * n_seasons + s;
      int nl_start = (lin_idx-1) * psi + 1;
      
      gsmoothseason[, lin_idx] = col(xholderLinear, lin_idx) * theta[lin_idx] 
                                + block(xholderNonlinear, 1, nl_start, grid_n, psi) * gamma[j,s];
      
      // Accumulate for marginal plot (average effect of X_j)
      gridgsmooth[, j] += gsmoothseason[, lin_idx] / n_seasons;
    }
  }

  // Compute Seasonal Alphas
  for (s in 1:n_seasons) {
    vector[grid_n] current_season_pred = rep_vector(theta0, grid_n);
    for (j in 1:p) {
      int lin_idx = (j-1) * n_seasons + s;
      current_season_pred += gsmoothseason[, lin_idx];
    }
    alphaseason[, s] = exp(current_season_pred);
  }

  // Marginal Alpha (Average alpha across the year)
  gridalpha = rep_vector(0.0, grid_n);
  for (s in 1:n_seasons) {
    gridalpha += alphaseason[, s] / n_seasons;
  }
}
"

data.stan <- list(
  y = as.vector(y.origin), 
  u = u, 
  n = n, 
  p = p, 
  # season_code = as.numeric(season_code_full),
  n_seasons = 4, 
  psi = psi, 
  bsLinear = bs.linear, 
  bsNonlinear = bs.nonlinear,
  xholderLinear = xholder.linear, 
  xholderNonlinear = xholder.nonlinear, 
  atau = ((psi+1)/2), 
  # Z_scales = Z_scales,
  X_sd = X_sd
)

init.alpha <- list(
  list(
    theta0 = 0.1,
    theta = rep(-0.1, p * 4),
    # theta0_seasonal = rep(0.1, 4),
    gamma_raw = array(rep(0.1, p * 4 * psi), dim = c(p, 4, psi)), 
    tau = array(rep(0.1, p * 4), dim = c(p, 4)),
    lambda1 = array(rep(0.1, p * 4), dim = c(p, 4)), 
    lambda2 = array(rep(1, p * 4), dim = c(p, 4))
  ),
  list(
    theta0 = 0.05,
    theta = rep(0.05, p * 4),
    # theta0_seasonal = rep(0.2, 4),
    gamma_raw = array(rep(-0.1, p * 4 * psi), dim = c(p, 4, psi)),
    tau = array(rep(0.2, p * 4), dim = c(p, 4)),
    lambda1 = array(rep(0.2, p * 4), dim = c(p, 4)), 
    lambda2 = array(rep(0.5, p * 4), dim = c(p, 4))
  ),
  list(
    theta0 = 0.3,
    theta = rep(0.1, p * 4),
    # theta0_seasonal = rep(0.3, 4),
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
# newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
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
# g.smooth.mean <- as.vector(matrix(newgsmooth.samples[,1], nrow = n, byrow=TRUE))
# g.smooth.q1 <- as.vector(matrix(newgsmooth.samples[,4], nrow = n, byrow=TRUE))
# g.smooth.q2 <- as.vector(matrix(newgsmooth.samples[,5], nrow = n, byrow=TRUE))
# g.smooth.q3 <- as.vector(matrix(newgsmooth.samples[,6], nrow = n, byrow=TRUE))


est_alphaseason_median  <- apply(posterior$alphaseason, c(2, 3), median)
est_alphaseason_lower <- apply(posterior$alphaseason, c(2, 3), quantile, probs = 0.05)
est_alphaseason_upper <- apply(posterior$alphaseason, c(2, 3), quantile, probs = 0.95)

# Format into a plotting dataframe
df_est_alphaseason <- do.call(rbind, lapply(1:4, function(s) {
  data.frame(
    x_index = seq(0, 1, length.out = grid_n),
    alpha_est = est_alphaseason_median[, s],
    lower_ci = est_alphaseason_lower[, s],
    upper_ci = est_alphaseason_upper[, s],
    alpha = true_alphaseason[,s],
    season = seasons[s]
  )
}))
# plot_data_alphaseason <- merge(df_est_alphaseason, df_true_alphaseason, by = c("season", "x_index"))

grid.plts <- list()
for(i in 1:4){
  grid.plt <- ggplot(df_est_alphaseason[((((i-1)*grid_n)+1):(i*grid_n)),], aes(x = x_index)) +
          geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "steelblue", alpha = 0.3) +
          geom_line(aes(y = alpha_est), color = "steelblue", linewidth = 1) +
          geom_line(aes(y = alpha), color = "red", linewidth = 1, linetype = "dashed") +
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

gsmooth_samples <- posterior$gridgsmooth 
summed_effects_list <- list()

for (j in 1:p) {
  col_indices <- ((j-1) * 4 + 1) : (j * 4)
  if(length(col_indices) > 1) {
    var_j_summed <- apply(gsmooth_samples[, , col_indices], c(1, 2), sum)
  } else {
    var_j_summed <- gsmooth_samples[, , col_indices]
  }

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
                          "true" =df_true_smooth$global[1:(grid_n*p)],
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
