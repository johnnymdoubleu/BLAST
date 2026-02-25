library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(qqboxplot)
library(dplyr)
# library(qgam)
library(evgam)
library(dplyr)

n <- 10000
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

# x.origin <- pnorm(matrix(rnorm(n * p), ncol = p) %*% C) # matrix(runif(n*p),ncol=p) 
x.origin <- matrix(0, nrow = n, ncol = p)

# for (j in 1:p) {
#   phase_shift <- j * (2 * pi / p) 
#   seasonal_trend <- 0.5 + 0.1 * sin(2 * pi * time.seq / period + phase_shift) #+ 0.25 * cos(2 * pi * time.seq / period + phase_shift)
#   uniform_noise <- runif(n, min = -0.4, max = 0.4)
#   x.origin[, j] <- seasonal_trend + uniform_noise
# }

# for (j in 1:p) {
#   phase_shift <- j * (2 * pi / p)
#   sin_val <- sin(2 * pi * time.seq / period + phase_shift)
#   amp_vector <- ifelse(sin_val > 0, 0.64, 0.24)
#   seasonal_trend <- 0.2 + (amp_vector * sin_val)
#   uniform_noise <- runif(n, min = -0.15, max = 0.15)
#   x.origin[, j] <- pmax(pmin(seasonal_trend + uniform_noise, 1), 0)  # clip to [0,1]
# }

for (j in 1:p) {
  phase_shift <- j * (2 * pi / p) 
  base_wave <- cos(phase_shift) * cos(2 * pi * time.seq / period) + 
               sin(phase_shift) * sin(2 * pi * time.seq / period)
  asymmetric_wave <- ((base_wave + 1) / 2)^1.5
  seasonal_trend <- 0.5 + 0.1 * asymmetric_wave
  uniform_noise <- runif(n, min = -0.5, max = 0.4)
  x.origin[, j] <- seasonal_trend + uniform_noise
}

plot(x.origin[,3])
covariates <- colnames(data.frame(x.origin))[1:p]

f2 <- function(x) {1.5 * sin(2 * pi * x^2)*x^3}
f3 <- function(x) {-1.5 * cos(3 * pi * x^2)*x^2}
theta.origin <- c(1, 0, 0.8, -0.8, 0, 0) 

alp.origin <- exp(rep(theta.origin[1],n) + x.origin%*%theta.origin[-1] + f2(x.origin[,2]) + f3(x.origin[,3]))
y.noise <- rPareto(n, rep(1, n), alpha = alp.origin)
f.season.scale <- function(t){
  return(2.5 - .8 * sin(2 * pi * t / 365) - .6 * cos(2 * pi * t/365)) 
}
y.origin <- y.noise * f.season.scale(time.seq)
# plot(y.origin)

evgam.df <- data.frame(
  y = log(y.origin),
  sin.time = sin(2 * pi * time.seq / 365),
  cos.time = cos(2 * pi * time.seq / 365)
)
evgam.cov <- y ~ 1 + cos.time + sin.time
ald.cov.fit <- evgam(evgam.cov, data = evgam.df, family = "ald", ald.args=list(tau = threshold))
u.vec <- exp(predict(ald.cov.fit)$location)

excess.index <- which(y.origin > u.vec)
x.origin <- x.origin[excess.index,]
y.origin <- y.origin[excess.index]
u <- u.vec[excess.index]
season_code_full <- season_code_full[excess.index]
n <- length(y.origin)
plot(x.origin[,1])
plot(y.origin)
x.origin <- data.frame(x.origin)
# ── Training basis construction ─────────────────────────────────────────────

season_indicator_list <- list()   # 4 season intercepts (no global)
linear_basis_list     <- list()   # season-varying slopes (unscaled)
nonlinear_basis_list  <- list()
model_data            <- list()

for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  x_vec    <- x.origin[, i]

  model_data[[var_name]] <- list()
  sm_list <- smoothCon(
    mgcv::s(x_vec, by = season_code_full, bs = "tp", k = psi + 2),
    data = data.frame(x_vec = x_vec, season_code_full = season_code_full)
  )

  for (j in seq_along(sm_list)) {
    season_label <- levels(season_code_full)[j]
    sm_spec      <- sm_list[[j]]
    mask         <- as.numeric(season_code_full == season_label)
    mask_idx     <- which(mask == 1)
    x_vec_season <- x_vec[mask_idx]

    # Season indicator built once (same across all covariates)
    if (i == 1) {
      ind_key <- paste0("intercept_S", season_label)
      season_indicator_list[[ind_key]] <- matrix(mask, ncol = 1)
    }

    # Within-season linear space for orthogonalization
    X_lin_season <- cbind(rep(1, length(mask_idx)), x_vec_season)
    qr_lin       <- qr(X_lin_season)

    # Spectral reparameterization
    X_raw  <- sm_spec$X
    S      <- sm_spec$S[[1]]
    eig    <- eigen(S, symmetric = TRUE)
    pos_idx         <- 1:psi
    U_pen           <- eig$vectors[, pos_idx]
    Lambda_sqrt_inv <- diag(1 / sqrt(pmax(eig$values[pos_idx], 1e-10)))  # guard near-zero
    Z_spectral      <- X_raw %*% U_pen %*% Lambda_sqrt_inv

    # Project out {mask, x_j * mask}
    Gamma_Original <- backsolve(
      qr.R(qr_lin),
      t(qr.Q(qr_lin)) %*% Z_spectral[mask_idx, , drop = FALSE]
    )
    X_lin_full <- cbind(mask, x_vec * mask)
    Z_orth     <- Z_spectral - X_lin_full %*% Gamma_Original

    # Filter using within-season rows only (avoids spurious near-zero noise elsewhere)
    keep_cols <- colSums(Z_orth[mask_idx, , drop = FALSE]^2) > 1e-9
    Z_final   <- Z_orth[, keep_cols, drop = FALSE]

    train_scale <- apply(Z_final[mask_idx, , drop = FALSE], 2, sd)
    train_scale[train_scale < 1e-12] <- 1
    Z_final <- scale(Z_final, center = FALSE, scale = train_scale) * mask

    model_data[[var_name]][[season_label]] <- list(
      sm_spec          = sm_spec,
      U_pen            = U_pen,
      Lambda_sqrt_inv  = Lambda_sqrt_inv,
      projection_coefs = Gamma_Original,
      keep_cols        = keep_cols,
      scale_stats      = train_scale
    )

    col_key <- paste0(var_name, "_S", season_label)   # consistent key: "X1_S1"
    linear_basis_list[[col_key]]    <- matrix(x_vec * mask, ncol = 1)
    nonlinear_basis_list[[col_key]] <- Z_final
  }
}

# ── Assemble training matrices ───────────────────────────────────────────────

bs.season    <- do.call(cbind, season_indicator_list)  # n × 4  (no global intercept)
bs.linear    <- do.call(cbind, linear_basis_list)      # n × 4p (unscaled)
bs.nonlinear <- do.call(cbind, nonlinear_basis_list)   # n × Σ|K^{j,s}|

# Scale bs.linear for Stan sampling; store stats for Stan back-transformation
linear_scale_stats <- apply(bs.linear, 2, function(col) {
  vals <- col[col != 0]
  if (length(vals) > 1) c(mean = mean(vals), sd = sd(vals)) else c(mean = 0, sd = 1)
})
X_means <- linear_scale_stats["mean", ]   # pass both to Stan data block
X_sds   <- linear_scale_stats["sd",   ]

bs.linear.scaled <- do.call(cbind, lapply(seq_len(ncol(bs.linear)), function(k) {
  col <- bs.linear[, k]
  nz  <- col != 0
  col[nz] <- (col[nz] - X_means[k]) / X_sds[k]
  col
}))
colnames(bs.linear.scaled) <- colnames(bs.linear)

# Full linear design (4 season intercepts + 4p scaled slopes)
# No global intercept — avoids dummy variable trap
bs.linear.full <- cbind(bs.season, bs.linear.scaled)  # n × 4(p+1)

for (j_cov in seq_along(covariates)) {
  var_name <- covariates[j_cov]
  nl_colnames <- grep(paste0("^", var_name, "_S"), colnames(bs.nonlinear), value = TRUE)
  lin_colnames <- grep(paste0("^", var_name, "_S"), colnames(bs.linear.scaled), value = TRUE)
  nl_cols <- which(colnames(bs.nonlinear) %in% nl_colnames)
  lin_cols <- which(colnames(bs.linear.scaled) %in% lin_colnames)
  lin_block <- cbind(bs.season, bs.linear.scaled[, lin_cols, drop = FALSE])
  nl_block  <- bs.nonlinear[, nl_cols, drop = FALSE]
  
  cat(sprintf("Ortho check [%s]: %.10e\n", var_name, sum(abs(t(lin_block) %*% nl_block))))
}


# ── Grid basis construction ──────────────────────────────────────────────────

grid_linear_list    <- list()
grid_nonlinear_list <- list()
grid_n <- 100

for (i in seq_along(covariates)) {
  var_name <- covariates[i]

  for (season_label in names(model_data[[var_name]])) {
    params  <- model_data[[var_name]][[season_label]]
    col_key <- paste0(var_name, "_S", season_label)  # matches training key

    season_idx <- which(season_code_full == season_label)
    raw_vals   <- x.origin[season_idx, i]
    raw_grid   <- seq(min(raw_vals), max(raw_vals), length.out = grid_n)

    # Raw (unscaled) linear grid — Stan applies (x - mean) / sd internally
    grid_linear_list[[col_key]] <- matrix(raw_grid, ncol = 1,
                                          dimnames = list(NULL, col_key))

    # Nonlinear grid: identical pipeline to training, no * mask needed
    pred_df <- data.frame(
      x_vec            = raw_grid,
      season_code_full = factor(rep(season_label, grid_n),
                                levels = levels(season_code_full))
    )
    X_raw_grid      <- PredictMat(params$sm_spec, pred_df)
    Z_spectral_grid <- X_raw_grid %*% params$U_pen %*% params$Lambda_sqrt_inv

    X_lin_grid  <- cbind(rep(1, grid_n), raw_grid)
    Z_orth_grid <- Z_spectral_grid - X_lin_grid %*% params$projection_coefs

    Z_final_grid <- Z_orth_grid[, params$keep_cols, drop = FALSE]
    Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)

    colnames(Z_final_grid) <- paste0(col_key, "_", seq_len(ncol(Z_final_grid)))
    grid_nonlinear_list[[col_key]] <- Z_final_grid
  }
}

# ── Assemble grid matrices ───────────────────────────────────────────────────

xholder.linear    <- do.call(cbind, grid_linear_list)    # grid_n × 4p  (raw, unscaled)
xholder.nonlinear <- do.call(cbind, grid_nonlinear_list) # grid_n × Σ|K^{j,s}|


scale_stats_list <- list()
for (var_name in names(model_data)) {
  for (season_label in names(model_data[[var_name]])) {
    stats <- model_data[[var_name]][[season_label]]$scale_stats
    scale_stats_list[[paste(var_name, "S", season_label, sep="_")]] <- stats
  }
}

Z_scales <- unlist(scale_stats_list)

stats <- apply(bs.linear, 2, function(x) {
  vals <- x[x != 0]
  if(length(vals) > 1) {
    return(c(mean = mean(vals), sd = sd(vals)))
  } else {
    return(c(mean = 0, sd = 1))
  }
})

X_means <- stats["mean", ]
X_sd   <- stats["sd", ]
bs.linear_scaled <- as.data.frame(lapply(1:ncol(bs.linear), function(j) {
  col <- bs.linear[, j]
  mask <- col != 0
  col[mask] <- (col[mask] - X_means[j]) / X_sd[j]
  return(col)
})) %>% as.matrix()
bs.linear <- bs.linear_scaled

model.stan <- "
// Stan model for Seasonally Varying Pareto Tail Index
data {
  int<lower=1> n; 
  int<lower=1> grid_n;
  int<lower=1> p; 
  int<lower=1> psi; 
  int<lower=1> n_seasons;
  matrix[n, (p+1) * n_seasons] bsLinear;     // intercepts + slopes (scaled)
  matrix[n, p * n_seasons * psi] bsNonlinear;                           // dynamic cols: Σ|K^{j,s}|
  matrix[grid_n, p * n_seasons] xholderLinear;       // slopes only (raw)
  matrix[grid_n, p * n_seasons * psi] xholderNonlinear;                 // dynamic cols
  vector[n] u; 
  vector<lower=u>[n] y; 
  real<lower=0> atau;
  vector[p * n_seasons] X_mean;                      // slope centering stats
  vector[p * n_seasons] X_sd;                        // slope scaling stats
  vector[p * n_seasons * psi] Z_scales;                                 // dynamic length
}

parameters {
  vector[n_seasons] beta0;                           // season intercepts
  vector[p * n_seasons] theta_raw;                   // scaled slopes
  array[p, n_seasons] real <lower=0> tau;        // smoothness per smooth
  array[p, n_seasons] vector[psi] gamma_raw;  // dynamic
  array[p] real<lower=0> lambda1;                    // Lasso shrinkage (linear)
  array[p] real<lower=0> lambda2;                    // Ridge shrinkage (nonlinear)
}

transformed parameters {
  vector[n] alpha;
  array[p, n_seasons] vector[psi] gamma;

  // Linear predictor: season intercepts + scaled slopes + nonlinear
  alpha = bsLinear[, 1:n_seasons] * beta0 + bsLinear[, (n_seasons+1):((p+1)*n_seasons)] * theta_raw;

  // Nonlinear effects (shrinkage-scaled)
  for (j in 1:p) {
    for (s in 1:n_seasons) {
      int nl_offset = ((j-1) * n_seasons + s - 1) * psi + 1;      
      for (k in 1:psi) {
        gamma[j,s][k] = gamma_raw[j, s][k] * sqrt(tau[j, s]) * Z_scales[nl_offset + k - 1];
      }
      alpha += block(bsNonlinear, 1, nl_offset, n, psi) * gamma[j,s];
    }
  }
  alpha = exp(alpha);
}

model {
  // Likelihood
  target += pareto_lpdf(y | u, alpha);

  // Priors
  target += normal_lpdf(beta0 | 0, 10);
  for (j in 1:p) {
    target += gamma_lpdf(lambda1[j] | 1e-1, 1e-1); 
    target += gamma_lpdf(lambda2[j] | 1e-2, 1e-2);      
    for (s in 1:n_seasons) {
      int idx = (j-1) * n_seasons + s;
      target += double_exponential_lpdf(theta_raw[idx] | 0, 1 / lambda1[j]);
      target += gamma_lpdf(tau[j, s] | atau, square(lambda2[j]) * 0.5);
      target += std_normal_lpdf(gamma_raw[j, s]);
    }
  }
}

generated quantities {
  matrix[grid_n, n_seasons] alphaseason;
  vector[grid_n] gridalpha; 
  matrix[grid_n, p] gridgsmooth;
  matrix[grid_n, p * n_seasons] gsmoothseason; 
  vector[p * n_seasons] theta = theta_raw .* X_sd;

  // Nonlinear smooths
  for (j in 1:p) {
    gridgsmooth[, j] = rep_vector(0, grid_n);
    for (s in 1:n_seasons) {
      int col_idx = (j-1) * n_seasons + s;
      int nl_offset = (col_idx - 1) * psi + 1;
      
      gsmoothseason[, col_idx] = col(xholderLinear, col_idx) * theta[col_idx] 
                                + block(xholderNonlinear, 1, nl_offset, grid_n, psi) * gamma[j,s];
      gridgsmooth[, j] += gsmoothseason[,col_idx] / n_seasons;
    }
  }

  
  // Per-season alpha (intercept + season-specific smooth)
  {
    for (s in 1:n_seasons) {
      vector[grid_n] current_season_pred = rep_vector(beta0[s], grid_n);
      for (j in 1:p) {
        int lin_idx = (j-1) * n_seasons + s;
        current_season_pred += gsmoothseason[, lin_idx];
      }
      alphaseason[, s] = exp(current_season_pred);
    }
  }
  
  // Global alpha (mean intercept + average smooth)
  {
    vector[grid_n] global_eta = rep_vector(mean(beta0), grid_n);
    for (j in 1:p) {
      global_eta += gridgsmooth[, j];
    }
    gridalpha = exp(global_eta);
  }
}

"
data.stan <- list(
  n              = nrow(bs.linear.full),
  grid_n         = nrow(xholder.linear),
  p              = length(covariates),
  psi            = psi,
  n_seasons      = 4,
  bsLinear       = bs.linear.full,                    # n × 4(p+1): intercepts + slopes
  bsNonlinear    = bs.nonlinear,                      # n × Σ|K^{j,s}|
  xholderLinear  = xholder.linear,                    # grid_n × 4p: slopes only (raw)
  xholderNonlinear = xholder.nonlinear,              # grid_n × Σ|K^{j,s}|
  u              = u,
  y              = y.origin,
  atau           = 1.0,                               # hyperparameter for tau prior
  X_mean         = X_means,                           # 4p: linear slope centering stats
  X_sd           = X_sds,                             # 4p: linear slope scaling stats
  # Compute Z_scales dynamically from model_data
  Z_scales       = Z_scales                          # assuming already computed as before
)

init.stan <- function(chain_id = 1) {
  # Chain-specific initialization for 3 chains
  # set.seed(1234 + chain_id)  # reproducible but different per chain
  
  list(
    beta0     = rnorm(4, 0, 0.5),                           # season intercepts
    theta_raw = rnorm(p * 4, 0, 0.2),                       # scaled slopes
    gamma_raw = array(
      rnorm(p * 4 * psi, 0, 0.1), 
      dim = c(p, 4, psi)
    ),
    lambda1   = rep(runif(p, 0.5, 1.5), 1),          # linear shrinkage
    lambda2   = rep(runif(p, 0.05, 0.2), 1),         # nonlinear shrinkage
    tau       = array(
      runif(p * 4, 0.5, 2.0),
      dim = c(p, 4)
    )
  )
}

# Usage for 3 chains:
# fit <- stan(file = "model.stan", data = data.stan, 
#             chains = 3, chain_id = 1:3,
#             init = lapply(1:3, init.stan))



fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    # named list of data
    chain_id = 1:3,
    init = lapply(1:3, init.stan),
    # init = init.stan,      # initial value 
    chains = 3,             # number of Markov chains
    iter = 2000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 1000           # no progress shown
)

posterior <- extract(fit1)
# bayesplot::color_scheme_set("mix-blue-red")
# bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("Scenario A") +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/appendix_",n,"_traceplot_scA.pdf"), width=22, height = 3)
beta0.samples <- summary(fit1, par=c("beta0"), probs = c(0.05,0.5, 0.95))$summary
theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
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

# true_eta_season <- matrix(0, nrow = grid_n, ncol = 4)
# true_alphaseason <- matrix(0, nrow = grid_n, ncol = 4)
# theta_0 <- theta.origin[1]
# theta_2 <- theta.origin[3]
# theta_3 <- theta.origin[4]

# for (s in 1:4) {
#   x2_seq <- xholder.linear[, paste0("X2_S_", s)]
#   x3_seq <- xholder.linear[, paste0("X3_S_", s)]

#   true_eta_season[, s] <- theta_0 + 
#                           (theta_2 * x2_seq + f2(x2_seq)) + 
#                           (theta_3 * x3_seq + f3(x3_seq))
  
#   true_alphaseason[, s] <- exp(true_eta_season[, s])
# }

# true_global_eta <- rowMeans(true_eta_season)
# true_gridalpha <- exp(true_global_eta)

# ── True seasonal alphas (matching new basis structure) ─────────────────────

theta_0 <- theta.origin[1]      # global intercept 0.7
theta_2 <- theta.origin[3]      # X2 coefficient 0.8  
theta_3 <- theta.origin[4]      # X3 coefficient -0.8
theta_other <- theta.origin[c(2,5,6)]  # 0 for X1, X4, X5

# DGP nonlinear functions (defined earlier in your script)
# f2(x) = 0.7 * sin(2π x^2) * x^3
# f3(x) = -0.7 * cos(3π x^2) * x^2

true_eta_season <- matrix(0, nrow = grid_n, ncol = 4)
true_alphaseason <- matrix(0, nrow = grid_n, ncol = 4)

for (s in 1:4) {
  # Extract season-specific covariate ranges from xholder.linear
  # New naming convention: "X1_S1", "X2_S1", etc. (no extra "_")
  x1_seq <- xholder.linear[, paste0("X1_S", s)]
  x2_seq <- xholder.linear[, paste0("X2_S", s)]
  x3_seq <- xholder.linear[, paste0("X3_S", s)]
  x4_seq <- xholder.linear[, paste0("X4_S", s)]
  x5_seq <- xholder.linear[, paste0("X5_S", s)]
  
  # True DGP: GLOBAL effects (no season variation in alpha)
  # theta_0 + sum_j theta_j * x_j + f2(x2) + f3(x3)
  true_eta_season[, s] <- theta_0/4 + 
    theta_other[1] * x1_seq +
    (theta_2 * x2_seq + f2(x2_seq)) + 
    (theta_3 * x3_seq + f3(x3_seq)) +
    theta_other[2] * x4_seq +
    theta_other[3] * x5_seq
  
  true_alphaseason[, s] <- exp(true_eta_season[, s])
}

# Global truth: average across seasons (identical to any column since global)
true_global_eta <- rowMeans(true_eta_season)  # or just true_eta_season[,1]
true_gridalpha <- exp(true_global_eta)

# ── Diagnostic: true season intercepts should all be theta_0 ─────────────────
true_season_intercepts <- colMeans(true_eta_season)



est_gridalpha_median  <- apply(posterior$gridalpha, 2, median)
est_gridalpha_lower <- apply(posterior$gridalpha, 2, quantile, probs = 0.05)
est_gridalpha_upper <- apply(posterior$gridalpha, 2, quantile, probs = 0.95)

data.scenario <- data.frame("x" = seq(0,1,length.out = grid_n),
                            "true" = true_gridalpha,
                            "q2" = est_gridalpha_median,
                            "q1" = est_gridalpha_lower,
                            "q3" = est_gridalpha_upper)
                            # "post.median" = (alpha.samples[,5]),
                            # "q1" = (alpha.samples[,4]),
                            # "q3" = (alpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,...,c))) + xlab(expression(c)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), fill="steelblue", alpha = 0.2) +
  geom_line(aes(y=true), color = "red", linewidth=1, linetype =2) +
  geom_line(aes(y=q2), color = "steelblue", linewidth=1) +
  theme_minimal(base_size = 30) + #ylim(0, 5)+
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))


est_alphaseason_median  <- apply(posterior$alphaseason, c(2, 3), median)
est_alphaseason_lower <- apply(posterior$alphaseason, c(2, 3), quantile, probs = 0.05)
est_alphaseason_upper <- apply(posterior$alphaseason, c(2, 3), quantile, probs = 0.95)

# Format into a plotting dataframe
df_est_alphaseason <- do.call(rbind, lapply(1:4, function(s) {
  data.frame(
    x_index = xholder.linear[,(j-1)*4+1],#seq(0, 1, length.out = grid_n),
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
          ) + #ylim(0,10) +
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




plot_data_list <- list()

for (s in 1:4) {
  x2_seq <- xholder.linear[, paste0("X2_S", s)]
  x3_seq <- xholder.linear[, paste0("X3_S", s)]
  
  alpha_true <- outer(x2_seq, x3_seq, function(x2, x3) {
    eta <- 1 + (0.8 * x2 + f2(x2)) + (-0.8 * x3 + f3(x3))
    exp(eta)
  })
  
  idx_X2 <- (2 - 1) * 4 + s
  idx_X3 <- (3 - 1) * 4 + s
  
  g2_median <- apply(posterior$gsmoothseason[, , idx_X2], 2, median)
  g3_median <- apply(posterior$gsmoothseason[, , idx_X3], 2, median)
  theta0_post <- median(posterior$beta0[,s]) # Or post$theta0_origin if centered
  
  alpha_post <- outer(1:grid_n, 1:grid_n, function(i, j) {
    exp(theta0_post + g2_median[i] + g3_median[j])
  })
  
  grid_df <- expand.grid(X2_idx = 1:grid_n, X3_idx = 1:grid_n)
  grid_df$X2 <- x2_seq[grid_df$X2_idx]
  grid_df$X3 <- x3_seq[grid_df$X3_idx]
  
  # Flatten the matrices into the dataframe
  grid_df$True_Alpha <- as.vector(alpha_true)
  grid_df$Post_Alpha <- as.vector(alpha_post)
  grid_df$Season <- seasons[s]
  
  plot_data_list[[s]] <- grid_df
}

# 2. Combine all seasons into one data frame
full_df <- do.call(rbind, plot_data_list)

long_df <- full_df %>%
  tidyr::pivot_longer(cols = c(True_Alpha, Post_Alpha), 
               names_to = "Model", 
               values_to = "Alpha") %>%
  mutate(Model = case_when(
    Model == "True_Alpha" ~ "TRUE",
    Model == "Post_Alpha" ~ "BLAST"
  ))

ggplot(long_df, aes(x = X3, y = X2, z = Alpha)) +
  geom_raster(aes(fill = Alpha)) +
  geom_tile(aes(fill = Alpha)) +
  geom_contour(color = "white", alpha = 0.4, bins = 25) +
  facet_grid(Season ~ Model) +
  # scale_fill_gradient(low="red", high="steelblue") +
  scale_fill_viridis_c(option = "D") + 
  theme_minimal() +
  labs(
    x = "X2",
    y = "X3",
    fill = expression(alpha(c,ldots,c))
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

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



T <- 500
len    <- nrow(posterior$alpha)
posterior_idx <- sample(len, T, replace = TRUE)
alpha_sub <- posterior$alpha[posterior_idx, ] 
y.matrix <- matrix(y.origin, nrow = T, ncol = n, byrow = TRUE)
u.matrix <- matrix(u, nrow = T, ncol = n, byrow = TRUE)
r_vec <- qnorm(pPareto(y.matrix, u.matrix, alpha = alpha_sub))
r_mat <- matrix(r_vec, nrow = T, ncol = n)
quantile_prob <- ppoints(n)
grid          <- qnorm(quantile_prob)
traj <- t(apply(r_mat, 1, sort))
bands <- matrixStats::colQuantiles(traj, probs = c(0.05, 0.5, 0.95))

qqplot_df <- data.frame(
  grid    = grid, 
  l_band  = bands[, 1], 
  trajhat = bands[, 2], 
  u_band  = bands[, 3],
  season = factor(season_code_full, label=seasons)
)

ggplot(qqplot_df, aes(x = grid)) + 
  geom_ribbon(aes(ymin = l_band, ymax = u_band), 
              fill = "steelblue", alpha = 0.3) + 
  geom_line(aes(y = trajhat), 
            color = "steelblue", linetype = "dashed", linewidth = 1.2) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 1) + 
  facet_wrap(~season) + 
  # coord_fixed(xlim = c(0, 4), ylim = c(0, 4)) +
  coord_fixed(xlim = c(-3, 3), ylim = c(-3, 3)) +
  labs(x = "", y = "Sample Quantiles") +
  theme_minimal(base_size = 30) +
  theme(axis.text = element_text(size = 20))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_qqplot.pdf"), width=10, height = 7.78)


posterior_idx <- sample(len, T, replace = TRUE)
y.matrix <- matrix(y.origin, nrow = T, ncol = n, byrow = TRUE)
u.matrix <- matrix(u, nrow = T, ncol = n, byrow = TRUE)
p_matrix <- 1 - (u.matrix/ y.matrix)^(posterior$alpha[posterior_idx,])
p_mean <- colMeans(p_matrix)
rp_integrated <- qnorm(p_mean)
p_q1 <- apply(p_matrix, 2, quantile, probs = 0.25)
p_q3 <- apply(p_matrix, 2, quantile, probs = 0.75)
rp_q1 <- qnorm(p_q1)
rp_q3 <- qnorm(p_q3)

rp <- data.frame(rp=as.numeric(rp_integrated), rp_q1 = rp_q1, rp_q3 = rp_q3, group = factor("residuals"), season = factor(season_code_full, label=seasons))
library(qqboxplot)
ggplot(data = rp) + 
  geom_qqboxplot(aes(y=rp, group = group), notch=FALSE, varwidth=FALSE, reference_dist="norm", width = 0.15, qq.colour = "steelblue")+
  labs(x = "", y = "Residuals") + 
  ylim(-4,4) + xlim(-.2,.2)+
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30))

ggplot(data = rp, aes(x = season, y = rp, group = season)) + 
  geom_qqboxplot(
    notch = FALSE, 
    varwidth = TRUE, 
    reference_dist = "norm", 
    width = 0.5, 
    qq.colour = "steelblue",
    fatten = 1
  ) +
  coord_cartesian(ylim = c(-4, 4)) + 
  labs(
    x = "QQ boxplot", 
    y = "Residuals"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text = element_text(size = 20),
    legend.position = "none"
  )



# Set up a 2x2 plotting grid to see all 4 seasons at once
# par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))

# # True parameter values based on your theta.origin = c(1, 0, 0.8, -0.8, 0, 0)
# theta_0 <- 1
# theta_2 <- 0.8
# theta_3 <- -0.8
# library(plotly)
# plot_list <- list()
# for (s in 1:4) {
#   # 1. Extract the specific 200-length grid for X2 and X3 for this specific season
#   x2_seq <- xholder.linear[, paste0("X2_S_", s)]
#   x3_seq <- xholder.linear[, paste0("X3_S_", s)]
  
#   # 2. Create the 200x200 surface matrix using outer()
#   # This calculates true alpha for every single combination of X2 and X3
#   alpha_matrix <- outer(x2_seq, x3_seq, function(x2, x3) {
#     # eta = intercept + linear X2 + nonlinear X2 + linear X3 + nonlinear X3
#     eta <- theta_0 + (theta_2 * x2 + f2(x2)) + (theta_3 * x3 + f3(x3))
#     exp(eta)
#   })
#   fig <- plot_ly(x = ~x2_seq, y = ~x3_seq)
  
#   # 1. Add the True Surface (Blue scale)
#   # add_surface(z = ~alpha_matrix, 
#   #             colorscale = "Blues", 
#   #             opacity = 0.9, # Slightly transparent
#   #             name = "True Alpha",
#   #             colorbar = list(title = "True", x = 1.0)) %>%
#   if (s == 1) {
#     fig <- fig %>%
#       add_surface(z = ~alpha_matrix, colorscale = "steelblue", opacity = 0.9, 
#                   name = "True Alpha", 
#                   colorbar = list(title = "True", x = 1.0)) #%>%
#       # add_surface(z = ~alpha_post_matrix, colorscale = "Reds", opacity = 0.7, 
#                   # name = "Posterior Median", 
#                   # colorbar = list(title = "Posterior", x = 1.15))
#   } else {
#     fig <- fig %>%
#       add_surface(z = ~alpha_matrix, colorscale = "steeblue", opacity = 0.9, 
#                   name = "True Alpha", showscale = FALSE)# %>%
#       # add_surface(z = ~alpha_post_matrix, colorscale = "Reds", opacity = 0.7, 
#                   # name = "Posterior Median", showscale = FALSE)
#   }
#   # 2. Add the Posterior Median Surface (Red/Orange scale)
#   # add_surface(z = ~alpha_post_matrix, 
#   #             colorscale = "Reds", 
#   #             opacity = 0.7, # More transparent so you can see the true surface beneath
#   #             name = "Posterior Median",
#   #             colorbar = list(title = "Posterior", x = 1.15)) %>%
  
#   # 3. Clean up the layout and axis labels
#   fig <- fig %>% layout(
#     title = seasons[s],
#     scene = list(
#       xaxis = list(title = "X2"),
#       yaxis = list(title = "X3"),
#       zaxis = list(title = "alpha(c,...,c)")
#     ),
#     annotations = list(
#       list(
#         text = seasons[s],
#         x = 0.5, y = 1.05, # Position title slightly above the plot
#         xref = "paper", yref = "paper",
#         showarrow = FALSE,
#         font = list(size = 16)
#       )
#     )
#   )
#   plot_list[[s]] <- fig
# }

# plot_list[[2]]
# fig
# final_dashboard <- subplot(plot_list, nrows = 2, margin = 0.05)
# final_dashboard

# Display the interactive plot


# par(mfrow = c(1, 1))

# true_alphaseason <- matrix(0, nrow = grid_n, ncol = 4)
# true_gsmoothseason <- matrix(0, nrow = grid_n, ncol = p * 4)

# for (s in 1:4) {
#   current_season_pred <- rep(true.theta0, grid_n)
#   for (j in 1:p) {
#     lin_col_name <- paste(covariates[j], "S", s, sep="_")
#     lin_idx <- (j-1)*4 + s
#     current_season_pred <- current_season_pred + xholder.linear[, lin_col_name] * true.theta[lin_idx]
#   }
#   x_vals_X2 <- xholder.linear[, paste("V2_S", s, sep="_")]  # X2's OWN grid for this season
#   current_season_pred <- current_season_pred + nl_func_X2(x_vals_X2)
#   x_vals_X3 <- xholder.linear[, paste("V3_S", s, sep="_")]  # X3's OWN grid for this season
#   current_season_pred <- current_season_pred + nl_func_X3(x_vals_X3)  
#   true_alphaseason[, s] <- exp(current_season_pred)
# }


# true_gridgsmooth <- matrix(0, nrow = grid_n, ncol = p)

# for (j in 1:p) {
#   smooth_sum <- rep(0, grid_n)
  
#   for (s in 1:4) {
#     lin_col <- paste(covariates[j], "S", s, sep = "_")
#     x_s <- xholder.linear[, lin_col]  # Season-specific grid
    
#     # Linear: evaluated at season-specific grid
#     smooth_sum <- smooth_sum + x_s * true.theta[(j-1)*4 + s]
    
#     # Nonlinear: evaluated at SAME season-specific grid
#     if (j == 2) smooth_sum <- smooth_sum + nl_func_X2(x_s)
#     if (j == 3) smooth_sum <- smooth_sum + nl_func_X3(x_s)
#   }
  
#   # Average across seasons (matching Stan's / n_seasons)
#   true_gridgsmooth[, j] <- smooth_sum / 4
# }

# df_true_gridgsmooth <- do.call(rbind, lapply(1:p, function(j) {
#   data.frame(
#     x_index = xholder.linear[, (j-1)*4 + 1], 
#     y = true_gridgsmooth[, j],
#     variable = paste0("Covariate_", j)
#   )
# }))

# true_global_eta <- rep(true.theta0, grid_n)
# for (j in 1:p) {
#   true_global_eta <- true_global_eta + true_gridgsmooth[, j]
# }

# true_gridalpha <- exp(true_global_eta)
# plot(true_gridalpha)
