library(VGAM)
library(mgcv)
library(Pareto)
suppressMessages(library(tidyverse))
library(readxl)
library(gridExtra)
library(corrplot)
library(rstan)
library(loo)
library(qqboxplot)
library(ggdensity)
library(ggforce)
library(ggdist)
library(evgam)
options(mc.cores = parallel::detectCores())

# Structure of the FWI System
#DSR : Daily Severity Rating
#FWI : Fire Weather Index
#BUI : Buildup Index
#ISI : Initial Spread Index
#FFMC : Fine FUel Moisture Code
#DMC : Duff Moisture Code
#DC : Drought Code

setwd("C:/Users/Johnny Lee/Documents/GitHub")
df <- read_excel("./BLAST/application/AADiarioAnual.xlsx", col_types = c("date", rep("numeric",40)))
df.long <- gather(df, condition, measurement, "1980":"2019", factor_key=TRUE)
missing.values <- which(!is.na(df.long$measurement))

Y <- df.long$measurement[!is.na(df.long$measurement)]
# Reduced psi to tighten credible bands slightly via less extreme flexibility
psi.origin <- psi <- 8
threshold <- 0.95

multiplesheets <- function(fname) {
    setwd("C:/Users/Johnny Lee/Documents/GitHub")
    sheets <- excel_sheets(fname)
    tibble <- lapply(sheets, function(x) read_excel(fname, sheet = x, col_types = c("date", rep("numeric", 41))))
    data_frame <- lapply(tibble, as.data.frame)
    names(data_frame) <- sheets
    return(data_frame)
}
setwd("C:/Users/Johnny Lee/Documents/GitHub")
path <- "./BLAST/application/DadosDiariosPT_FWI.xlsx"

cov <- multiplesheets(path)
fwi.scaled <- fwi.index <- data.frame(DSR = double(length(Y)),
                                      FWI = double(length(Y)),
                                      BUI = double(length(Y)),
                                      ISI = double(length(Y)),
                                      FFMC = double(length(Y)),
                                      DMC = double(length(Y)),
                                      DC = double(length(Y)),
                                      stringsAsFactors = FALSE)
for(i in 1:length(cov)){
    cov.long <- gather(cov[[i]][,1:41], condition, measurement, "1980":"2019", factor_key=TRUE)
    fwi.index[,i] <- cov.long$measurement[missing.values]
    fwi.scaled[,i] <- cov.long$measurement[missing.values]
}

fwi.scaled$time <- fwi.index$time <- seq(1,length(Y), length.out=length(Y))
fwi.scaled$sea <- fwi.index$sea <- fwi.index$time %% 365 / 365
fwi.scaled$cos.time <- fwi.index$cos.time <- cos(2*pi*seq(1,length(Y), length.out=length(Y))/365)
# FIX: Corrected variable assignment here
fwi.scaled$sin.time <- fwi.index$sin.time <- sin(2*pi*seq(1,length(Y), length.out=length(Y))/365) 

fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                          levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)

# NOTE: The detrending GAM loop has been completely removed to preserve absolute physical covariates.

above.0 <- which(Y>0)
Y <- Y[above.0]
fwi.scaled <- fwi.scaled[above.0,]
fwi.index <- fwi.index[above.0,]

# Calculate varying threshold u.vec using harmonics
qr.df <- data.frame(y = (Y), fwi.scaled)
evgam.cov <- y ~ 1 + cos.time + sin.time 
ald.cov.fit <- evgam(evgam.cov, data = qr.df, family = "ald", ald.args=list(tau = threshold))
u.vec <- (predict(ald.cov.fit)$location)

excess <- which(Y>u.vec)
y <- Y[excess]
u <- u.vec[excess]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Extract the covariates we want (5 FWI variables + 2 Harmonics)
fwi.cols <- c("BUI", "ISI", "FFMC", "DMC", "DC")
harm.cols <- c("cos.time", "sin.time")
all_linear_vars <- fwi.scaled[excess, c(fwi.cols, harm.cols)]

# Min/Max for unscaling later
X_minmax <- sapply(all_linear_vars, function(x) max(x)-min(x))
X_min <- sapply(all_linear_vars, function(x) min(x))

# Scale ALL linear terms to 0-1 for fair Lasso penalization
all_scaled <- as.data.frame(sapply(all_linear_vars, FUN = range01))

# Separate the FWI components for the nonlinear QR projection
fwi.scaled.only <- all_scaled[, 1:length(fwi.cols)]

n <- dim(all_scaled)[[1]]
p <- length(fwi.cols)      # 5 Nonlinear covariates
p_lin <- ncol(all_scaled)  # 7 Linear covariates (5 FWI + 2 Harmonics)

fwi.origin <- data.frame(all_linear_vars, BA=y)
max.fwi <- fwi.origin[which.max(y),]
fwi.grid <- data.frame(lapply(all_linear_vars, function(x) seq(min(x), max(x), length.out = n)))

bs.linear <- model.matrix(~ ., data = all_scaled)[,-1]

psi <- psi - 2
group.map <- c()
Z.list <- list()        
scale_stats_list <- list() 
projection_coefs_list <- list() 
spec_decomp_list <- list() 
qr_list <- list()          
sm_spec_list <- list()     
keep_cols_list <- list()

covariates <- fwi.cols

# Apply QR Spline loop ONLY to the 5 FWI Variables
for (i in seq_along(covariates)) {
    var_name <- covariates[i]
    x_vec <- fwi.scaled.only[, i]
    X_lin <- model.matrix(~ x_vec) 
    sm_spec <- smoothCon(mgcv::s(x_vec, bs = "tp", k = psi + 2), 
                         data = data.frame(x_vec = x_vec), 
                         knots = NULL)[[1]]
    
    X_raw <- sm_spec$X
    S     <- sm_spec$S[[1]] 
    
    eig <- eigen(S, symmetric = TRUE)
    max_lambda <- max(eig$values)
    tol <- max_lambda * 1e-8  
    
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
    train_scale[train_scale < 1e-12] <- 1 
    Z_final <- scale(Z_final, center = FALSE, scale = train_scale)

    Z.list[[i]] <- Z_final
    group.map <- c(group.map, rep(i, ncol(Z_final)))
    
    spec_decomp_list[[i]] <- list(U_pen = U_pen, Lambda_sqrt_inv = solve(sqrt(Lambda_pen)))
    projection_coefs_list[[i]] <- Gamma_Original 
    keep_cols_list[[i]] <- keep_cols
    scale_stats_list[[i]] <- train_scale         
    sm_spec_list[[i]] <- sm_spec
}

bs.nonlinear <- do.call(cbind, Z.list)

grid.n <- n
newx <- seq(0, 1, length.out = grid.n)
xholder_fwi <- do.call(cbind, lapply(1:p, function(j) {newx}))
colnames(xholder_fwi) <- fwi.cols

# 1. Define your meaningful day of the year (e.g., August 15 is day 227)
peak_day <- 227 

raw_cos_peak <- cos(2 * pi * peak_day / 365)
raw_sin_peak <- sin(2 * pi * peak_day / 365)
# scaled_cos_peak <- (raw_cos_peak + 1) / 2
# scaled_sin_peak <- (raw_sin_peak + 1) / 2
xholder_harm <- matrix(
  c(rep(raw_cos_peak, grid.n), rep(raw_sin_peak, grid.n)), 
  ncol = 2
)
colnames(xholder_harm) <- harm.cols
xholder_all <- cbind(xholder_fwi, xholder_harm)

grid_Z_list <- list()

for (i in seq_along(covariates)) {
    grid_df  <- data.frame(x_vec = xholder_fwi[,i])
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
xholder.linear <- model.matrix(~ ., data = data.frame(xholder_all))[,-1]
Z_scales <- unlist(scale_stats_list)

grid_Z_list <- list()

for (i in seq_along(covariates)) {
    scaled_x_grid <- (fwi.grid[,i] - X_min[i]) / (X_minmax[i])
    grid_df  <- data.frame(x_vec = scaled_x_grid)
    X_lin_grid <- model.matrix(~ x_vec, data = grid_df)
    X_raw_grid <- PredictMat(sm_spec_list[[i]], grid_df)
    
    decomp <- spec_decomp_list[[i]]
    Z_spectral_grid <- X_raw_grid %*% decomp$U_pen %*% decomp$Lambda_sqrt_inv
    Z_orth_grid <- Z_spectral_grid - X_lin_grid %*% projection_coefs_list[[i]]
    Z_final_grid <- Z_orth_grid[, keep_cols_list[[i]], drop = FALSE]
    Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = scale_stats_list[[i]])
    grid_Z_list[[i]] <- Z_final_grid
}

xgrid.nonlinear <- do.call(cbind, grid_Z_list)
xgrid.linear <- model.matrix(~ ., data = data.frame(fwi.grid))[,-1]
xgrid.linear <- sweep(xgrid.linear, 2, X_min, "-")

X_means <- colMeans(bs.linear)
X_sd   <- apply(bs.linear, 2, sd)
bs.linear <- scale(bs.linear, center = X_means, scale = X_sd)


model.stan <- "// Stan model for BLAST Pareto Samples
data {
  int <lower=1> n; 
  int <lower=1> grid_n; 
  int <lower=1> p;       // Non-linear components (5 FWI)
  int <lower=1> p_lin;   // Total linear components (5 FWI + 2 Harmonics)
  int <lower=1> psi; 
  vector<lower=0>[n] u; 
  matrix[n, p_lin] bsLinear; 
  matrix[n, (psi*p)] bsNonlinear; 
  matrix[grid_n, p_lin] xholderLinear; 
  matrix[grid_n, (psi*p)] xholderNonlinear;     
  matrix[n, p_lin] gridL; 
  matrix[n, (psi*p)] gridNL;       
  vector<lower=0>[n] y; 
  real <lower=0> atau;
  vector[p_lin] X_minmax;
  vector[p_lin] X_min;
  vector[p_lin] X_means;
  vector[p_lin] X_sd;
}

parameters {
  vector[(p_lin+1)] theta; // linear predictor (intercept + p_lin)
  array[p] vector[psi] gamma_raw;
  array[p_lin] real <lower=0> lambda1; // lasso penalty for ALL linear terms
  array[p] real <lower=0> lambda2; // lambda2 group lasso penalty (FWI splines only)
  array[p] real <lower=0> tau;
}

transformed parameters {
  vector[n] alpha; 
  array[p] vector[psi] gamma;
  {
    vector[n] eta = rep_vector(theta[1], n);
    
    // Add linear terms
    for (k in 1:p_lin) {
       eta += col(bsLinear, k) * theta[k+1];
    }
    
    // Add non-linear terms
    for (j in 1:p){
      gamma[j] = gamma_raw[j] * sqrt(tau[j]);
      eta += block(bsNonlinear, 1, ((j - 1) * psi + 1), n, psi) * gamma[j];
    };
    alpha = exp(eta);
  }
}

model {
  // likelihood
  target += pareto_lpdf(y | u, alpha);
  
  // TIGHTENED PRIOR to fix extreme uncertainty near 0
  target += normal_lpdf(theta[1] | 0, 2); 
  
  // Stronger regularization to squash non-informative covariates
  target += gamma_lpdf(lambda1 | 2, 0.1); 
  target += gamma_lpdf(lambda2 | 2, 0.1);    
  
  for (k in 1:p_lin){
    target += double_exponential_lpdf(theta[(k+1)] | 0, 1/(lambda1[k]));
  }
  
  for (j in 1:p){
    target += gamma_lpdf(tau[j] | atau, square(lambda2[j])*0.5);
    target += std_normal_lpdf(gamma_raw[j]);
  }
}

generated quantities {
  vector[grid_n] gridalpha; 
  matrix[grid_n, p] gridgsmooth; // Keep as 'p' for plotting just the FWI variables
  matrix[n, p] fwismooth;
  vector[n] log_lik;
  vector[n] time;
  vector[p_lin+1] theta_origin;
  vector[p_lin] theta_fwi;

  theta_origin[2:(p_lin+1)] = theta[2:(p_lin+1)] ./ X_sd;
  theta_origin[1] = theta[1] - dot_product(X_means, theta_origin[2:(p_lin+1)]);
  theta_fwi = theta_origin[2:(p_lin+1)] ./ X_minmax;

  {
    vector[grid_n] pred = rep_vector(theta_origin[1], grid_n);
    // Add Harmonics base effect to overall prediction
    for (h in (p+1):p_lin) {
       pred += col(xholderLinear, h) * theta_origin[h+1]; 
    }
    
    for (j in 1:p){
      int nl_start = (j - 1) * psi + 1;
      fwismooth[,j] = block(gridNL,1, nl_start, n, psi) * gamma[j] + col(gridL,j) * theta_fwi[j];
      gridgsmooth[,j] = col(xholderLinear, j) * theta_origin[j+1]  + block(xholderNonlinear,1, nl_start, grid_n, psi) * gamma[j];
      pred += gridgsmooth[,j];
    }
    time = col(bsLinear, p+1) * theta[p_lin] + col(bsLinear, p_lin) * theta[p_lin+1];
    gridalpha = exp(pred);
  }

  for (i in 1:n){
    log_lik[i] = pareto_lpdf(y[i] | u[i], alpha[i]);
  };
}
"

data.stan <- list(y = as.vector(y), u = u, p = p, p_lin = p_lin, n= n, psi = psi,
                  atau = ((psi+1)/2), xholderNonlinear = xholder.nonlinear, 
                  bsLinear = bs.linear, bsNonlinear = bs.nonlinear, X_min=X_min,
                  xholderLinear = xholder.linear, X_minmax = X_minmax, 
                  X_means = X_means, X_sd = X_sd, grid_n = grid.n,
                  gridL = xgrid.linear, gridNL = xgrid.nonlinear)

# Make sure init.alpha accounts for p_lin vs p where appropriate
init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p, psi)), 
                        theta = rep(-0.1, (p_lin+1)), tau = rep(0.1, p), 
                        lambda1 = rep(0.1, p_lin), lambda2 = rep(1, p)),
                   list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p, psi)),
                        theta = rep(0.05, (p_lin+1)), tau = rep(0.2, p),
                        lambda1 = rep(2, p_lin), lambda2 = rep(5, p)),
                   list(gamma_raw= array(rep(0.1, (psi*p)), dim=c(p, psi)),
                        theta = rep(0.1, (p_lin+1)), tau = rep(0.1, p), 
                        lambda1 = rep(0.1, p_lin), lambda2 = rep(0.1, p)))

fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    
    init = init.alpha,      
    chains = 3,             
    iter = 4000,            
    cores = parallel::detectCores(), 
    refresh = 2000          
)

# saveRDS(fit1, file=paste0("./BLAST/application/",Sys.Date(),"_stanfit.rds"))
# readRDS(file=paste0("./BLAST/application/2024-11-27","_stanfit.rds"))
posterior <- rstan::extract(fit1)
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))

theta.samples <- summary(fit1, par=c("theta_origin"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
fwismooth.samples <- summary(fit1, par=c("fwismooth"), probs = c(0.05, 0.5, 0.95))$summary
# gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
# gridgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
origin.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary
time.samples <- summary(fit1, par=c("time"), probs = c(0.05,0.5, 0.95))$summary
# yrep <- summary(fit1, par=c("yrep"), probs = c(0.05,0.5, 0.95))$summary
# f.samples <- summary(fit1, par=c("f"), probs = c(0.05,0.5, 0.95))$summary
# loglik.samples <- summary(fit1, par=c("log_lik"), probs = c(0.05,0.5, 0.95))$summary

# MCMCvis::MCMCplot(fit1, params = 'theta')
# MCMCvis::MCMCplot(fit1, params = "gamma")


# post.samples <- as.matrix(fit1)
# theta_samples <- post.samples[, grepl("^theta\\[", colnames(post.samples))]
# gamma_samples <- post.samples[, grepl("^gamma\\[", colnames(post.samples))]
# num_draws <- nrow(post.samples)
# p <- length(X_sd)
# psi <- ncol(gamma_samples) / p
# theta_fwi_samples <- matrix(NA, nrow = num_draws, ncol = p)
# fwismooth_samples <- array(NA, dim = c(num_draws, n, p))

# for (s in 1:num_draws) {
#   curr_theta <- theta_samples[s, ]
#   theta_orig_slopes <- curr_theta[2:(p+1)] / X_sd
#   theta_fwi_s <- theta_orig_slopes / fwi.minmax
#   for (j in 1:p) {
#     start_idx <- ((j - 1) * psi) + 1
#     end_idx   <- j * psi
#     curr_gamma_j <- as.matrix(gamma_samples[s, start_idx:end_idx])
    
#     nl_part <- xgrid.nonlinear[, start_idx:end_idx] %*% curr_gamma_j
#     l_part  <- (xgrid.linear[, j+1]-fwi.min[j]) * theta_fwi_s[j]
    
#     fwismooth_samples[s, , j] <- as.vector(nl_part + l_part)
#   }
# }



g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = grid.n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = grid.n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = grid.n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = grid.n, byrow=TRUE))
# g.linear.mean <- as.vector(matrix(gridgl.samples[,1], nrow = grid.n, byrow=TRUE))
# g.linear.q1 <- as.vector(matrix(gridgl.samples[,4], nrow = grid.n, byrow=TRUE))
# g.linear.q2 <- as.vector(matrix(gridgl.samples[,5], nrow = grid.n, byrow=TRUE))
# g.linear.q3 <- as.vector(matrix(gridgl.samples[,6], nrow = grid.n, byrow=TRUE))
# g.nonlinear.mean <- as.vector(matrix(gridgnl.samples[,1], nrow = grid.n, byrow=TRUE))
# g.nonlinear.q1 <- as.vector(matrix(gridgnl.samples[,4], nrow = grid.n, byrow=TRUE))
# g.nonlinear.q2 <- as.vector(matrix(gridgnl.samples[,5], nrow = grid.n, byrow=TRUE))
# g.nonlinear.q3 <- as.vector(matrix(gridgnl.samples[,6], nrow = grid.n, byrow=TRUE))
fwi.smooth.mean <- as.vector(matrix(fwismooth.samples[,1], nrow = n, byrow=TRUE))
fwi.smooth.q1 <- as.vector(matrix(fwismooth.samples[,4], nrow = n, byrow=TRUE))
fwi.smooth.q2 <- as.vector(matrix(fwismooth.samples[,5], nrow = n, byrow=TRUE))
fwi.smooth.q3 <- as.vector(matrix(fwismooth.samples[,6], nrow = n, byrow=TRUE))

# fwi.smooth.mean <- as.vector(apply(fwismooth_samples, c(2, 3), mean))
# fwi.smooth.q2 <- as.vector(apply(fwismooth_samples, c(2, 3), median))
# fwi.smooth.q1 <- as.vector(apply(fwismooth_samples, c(2, 3), quantile, probs = 0.05))
# fwi.smooth.q3 <- as.vector(apply(fwismooth_samples, c(2, 3), quantile, probs = 0.95))

g.min.samples <- min(gsmooth.samples[,4])
g.max.samples <- max(gsmooth.samples[,6])
# fwi.smooth <- as.data.frame(apply(fwismooth_samples, c(2, 3), quantile, probs = 0.05))
fwi.smooth <- as.data.frame(matrix(fwismooth.samples[,4], nrow = n, byrow=TRUE))
fwi.min.samples <- sapply(fwi.smooth, min)
# data.smooth <- data.frame("x"= as.vector(xholder),
#                           "post.mean" = as.vector(g.smooth.mean),
#                           "q1" = as.vector(g.smooth.q1),
#                           "q2" = as.vector(g.smooth.q2),
#                           "q3" = as.vector(g.smooth.q3),
#                           "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
#                           "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

# ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#   ylab("") + xlab("") +
#   facet_grid(covariates ~ ., scales = "free",
#               labeller = label_parsed) + 
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values=c("steelblue")) + 
#   guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) + 
#           # ylim(-3.5, 3.5) +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#           plot.margin = margin(0,0,0,-20),
#           # strip.text = element_blank(),
#           axis.text = element_text(size = 20))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_smooth.pdf"), width=12.5, height = 15)
# xholder <- as.data.frame(xholder)
# colnames(xholder) <- colnames(fwi.scaled)[1:p]
simul.data <- data.frame(BA = y-u, fwi.scaled.only, 
                          sin.time = fwi.scaled$sin.time[excess], 
                          cos.time = fwi.scaled$cos.time[excess])#fwi.origin[c(1:p)])

gam.scale <- list(
  # 1. Formula for the log-Scale parameter
  BA ~ sin.time + cos.time + 
       s(BUI, bs = "ts", k = 10) + 
       s(ISI, bs = "ts", k = 10) + 
       s(FFMC, bs = "ts", k = 10) +
       s(DMC, bs = "ts", k = 10) + 
       s(DC, bs = "ts", k = 10),
       
  # 2. Formula for the Shape parameter (xi)
     ~ sin.time + cos.time + 
       s(BUI, bs = "ts", k = 10) + 
       s(ISI, bs = "ts", k = 10) + 
       s(FFMC, bs = "ts", k = 10) +
       s(DMC, bs = "ts", k = 10) +
       s(DC, bs = "ts", k = 10)
)
evgam.fit.scale <- evgam::evgam(gam.scale, data = simul.data, family = "gpd")
pred.scale <- predict(evgam.fit.scale, newdata = as.data.frame(xholder_all), type="response", se.fit = TRUE)
xi.pred.scale <-pred.scale$fitted$shape
# xi.se.scale <- pred.scale$se.fit$shape
# xi.low.scale <- xi.pred.scale - (1.96 * xi.se.scale)
# xi.high.scale <- xi.pred.scale + (1.96 * xi.se.scale)
# alpha.pred.scale <- 1/xi.pred.scale

# xholder.basis.scale <- predict(evgam.fit.scale, newdata = xholder, type= "lpmatrix")$shape
# psi <- 15
# xi.coef.scale <- tail(evgam.fit.scale$coefficients, (psi-1)*p)
# gamma.xi.scale <- matrix(xi.coef.scale, ncol = p)
# alpha.nonlinear.scale <- xi.nonlinear.scale <- matrix(, nrow = n, ncol = p)
# bs.nonlinear.scale <- xholder.basis.scale[,c(2:((psi-1)*p+1))]
# for(j in 1:p){
#   xi.nonlinear.scale[,j] <- bs.nonlinear.scale[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.scale[,j]
#   alpha.nonlinear.scale[,j] <- 1/(xi.nonlinear.scale[,j])
# }

gam.1 <- list(BA ~ 1,
                ~ cos.time + sin.time + 
                  s(BUI, bs = "ts", k = 10) + 
                  s(ISI, bs = "ts", k = 10) + 
                  s(FFMC, bs = "ts", k = 10) +
                  s(DMC, bs = "ts", k = 10) +
                  s(DC, bs = "ts", k = 10))
evgam.fit.1 <- evgam::evgam(gam.1, data = simul.data, family = "gpd")
pred.1 <- predict(evgam.fit.1, newdata = as.data.frame(xholder_all), type="response", se.fit = TRUE)
xi.pred.1 <-pred.1$fitted$shape
# xi.se.1 <- pred.1$se.fit$shape
# xi.low.1 <- xi.pred.1 - (1.96 * xi.se.1)
# xi.high.1 <- xi.pred.1 + (1.96 * xi.se.1)
# alpha.pred.1 <- 1/xi.pred.1

# xholder.basis.1 <- predict(evgam.fit.1, newdata = xholder, type= "lpmatrix")$shape
# xi.coef.1 <- tail(evgam.fit.1$coefficients, (psi-1)*p)
# gamma.xi.1 <- matrix(xi.coef.1, ncol = p)
# alpha.nonlinear.1 <- xi.nonlinear.1 <- matrix(, nrow = n, ncol = p)
# bs.nonlinear.1 <- xholder.basis.1[,c(2:((psi-1)*p+1))]
# for(j in 1:p){
#   xi.nonlinear.1[,j] <- bs.nonlinear.1[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.1[,j]
#   alpha.nonlinear.1[,j] <- 1/(xi.nonlinear.1[,j])
# }
# simul.data <- data.frame(BA = y, fwi.scaled[,c(1:p)])
# vgam.fit.scale <- vgam(BA ~ sm.ps(BUI, ps.int = 28) + 
#                            sm.ps(ISI, ps.int = 28) + 
#                            sm.ps(FFMC, ps.int = 28) + 
#                            sm.ps(DMC, ps.int = 28) + 
#                            sm.ps(DC, ps.int = 28),
#                         data = simul.data,
#                         family = gpd(threshold= u,
#                                       lshape="loglink",
#                                       zero = NULL),
#                         trace = TRUE,
#                         control = vgam.control(maxit = 200))
  # fitted.linear <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "link")
  # fitted.terms <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "terms")
  # vgam.xi.scale <- exp(fitted.linear[,2])
  # vgam.sigma.scale <- exp(fitted.linear[,1])

# vgam.fit.1 <- vgam(BA ~ sm.ps(BUI, ps.int = 28) + 
#                        sm.ps(ISI, ps.int = 28) + 
#                        sm.ps(FFMC, ps.int = 28) + 
#                        sm.ps(DMC, ps.int = 28) + 
#                        sm.ps(DC, ps.int = 28),
#                       data = simul.data,
#                       family = gpd(threshold= u,
#                                     lshape="loglink",
#                                     zero = 1),
#                       trace = TRUE,
#                       control = vgam.control(maxit = 200))
# fitted.linear <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "link")
# fitted.terms <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "terms")
# vgam.xi.1 <- exp(fitted.linear[,2])
# vgam.sigma.1 <- exp(fitted.linear[,1])


data.scenario <- data.frame("x" = newx,
                            "post.mean" = (alpha.samples[,1]),
                            "post.median" = (alpha.samples[,5]),
                            "q1" = (alpha.samples[,4]),
                            "q3" = (alpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(bold(c),bold(t)))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue")) + ylim(0, 10) +
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_alpha.pdf"), width=10, height = 7.78)

xi.scenario <- data.frame("x" = newx,
                            "post.mean" = (1/alpha.samples[,1]),
                            "post.median" = (1/alpha.samples[,5]),
                            "q1" = (1/alpha.samples[,4]),
                            "q3" = (1/alpha.samples[,6]),
                            "evgam.1" = xi.pred.1,
                            # "evgam.1.q1" = xi.low.1,
                            # "evgam.1.q3" = xi.high.1,
                            # "vgam.1" = vgam.xi.1,
                            # "vgam.scale" = vgam.xi.scale,
                            # "evgam.scale.q1" = xi.low.scale,
                            # "evgam.scale.q3" = xi.high.scale,
                            "evgam.scale" = xi.pred.scale)

ggplot(xi.scenario, aes(x=x)) + 
  ylab(expression(xi(c,...,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  # geom_ribbon(aes(ymin = evgam.scale.q1, ymax = evgam.scale.q3), fill= "orange", alpha = 0.2) +
  geom_line(aes(y=evgam.scale), colour = "orange", linewidth=1, linetype=3) +
  # geom_ribbon(aes(ymin = evgam.1.q1, ymax = evgam.1.q3), fill= "purple", alpha = 0.2) +
  geom_line(aes(y=evgam.1), colour = "purple", linewidth=1, linetype=3) +
  # geom_line(aes(y=vgam.scale), colour = "orange", linewidth=1, linetype=4) +
  # geom_line(aes(y=vgam.1), colour = "purple", linewidth=1, linetype=4) +  
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))



T <- 100
len <- dim(posterior$alpha)[1]
posterior_idx <- sample(1:len, T, replace = TRUE)
r <- sapply(1:T, function(t) {
  alpha_t <- posterior$alpha[posterior_idx[t], ]  # Extract vector of length n
  qnorm(pPareto(y, u, alpha = alpha_t))           # Vectorized across all y
})

quantile.prob <- ppoints(n)
grid <- qnorm(quantile.prob)
traj <- t(apply(r, 2, quantile, probs = quantile.prob, type = 2))
# r <- matrix(, nrow = n, ncol = T)
# for(i in 1:n){
#   for(t in 1:T){
#     r[i, t] <- qnorm(pPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i]))
#   }
# }

# traj <- matrix(NA, nrow = T, ncol = lgrid)
# for (t in 1:T){
#   traj[t, ] <- quantile(r[, t], quantile.prob, type = 2)
# }

l.band <- apply(traj, 2, quantile, prob = 0.025)
trajhat <- apply(traj, 2, quantile, prob = 0.5)
u.band <- apply(traj, 2, quantile, prob = 0.975)
qqplot.df <- data.frame(grid = grid, 
                        l.band = l.band,  
                        trajhat = trajhat,
                        u.band = u.band)
ggplot(data = ) + 
  geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
              fill = "steelblue",
              alpha = 0.4, linetype = "dashed") + 
  geom_line(aes(x = grid, y = trajhat), colour = "steelblue", linetype = "dashed", linewidth = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.2) + 
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
  theme_minimal(base_size = 30) +
  theme(axis.text = element_text(size = 20)) + 
  coord_fixed(xlim = c(-3, 3),
              ylim = c(-3, 3))

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_qqplot.pdf"), width=10, height = 7.78)
# save(loglik.samples, data.smooth, data.scenario, qqplot.df, file="./BLAST/application/blast_1.Rdata")

rp <-c()
for(i in 1:n){
  # rp[i] <- rPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i])
  rp[i] <- qnorm(pPareto(y[i], u[i], alpha = posterior$alpha[round(runif(1,1,len)),i]))
}
rp <- data.frame(rp, group = rep("residuals", n))

ggplot(data = rp) + 
  # geom_qqboxplot(aes(factor(group, levels=c("residuals")), y=rp), notch=FALSE, varwidth=TRUE, reference_dist="norm")+ 
  geom_qqboxplot(aes(y=rp), notch=FALSE, varwidth=FALSE, reference_dist="norm", width = 0.15, qq.colour = "steelblue")+
  labs(x = "", y = "Residuals") + ylim(-4,4) + xlim(-.2,.2)+
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_qqboxplot.pdf"), width = 10, height = 7.78)
             
cat("Finished Running")

# relative_eff(exp(fit.log.lik))
#https://discourse.mc-staqan.org/t/four-questions-about-information-criteria-cross-validation-and-hmc-in-relation-to-a-manuscript-review/13841/3


data.smooth <- data.frame("x" = as.vector(as.matrix(xholder.linear[,1:p])),
                          "true" = as.vector(as.matrix(fwi.scaled.only)),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "covariates" = gl(p, grid.n, (p*grid.n), labels = names(fwi.scaled)))


grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*grid.n)+1):(i*grid.n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  # ylim(g.min.samples, g.max.samples) +
                  ylim(-3, 3) +
                  theme_minimal(base_size = 20) +
                  theme(legend.position = "none",
                        plot.margin = margin(5, 5, 5, 5),
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        axis.text = element_text(size = 18),
                        axis.title.x = element_text(size = 22))
                  # theme_minimal(base_size = 30) +
                  # theme(legend.position = "none",
                  #         plot.margin = margin(0,0,0,-20),
                  #         axis.text = element_text(size = 35),
                  #         axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled.only[which.max(y),i], y=g.min.samples, color = "red", size = 7)
}
marrangeGrob(grobs = grid.plts, nrow = 1, ncol = 5, top = NULL)
# grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_DC.pdf"), grid.plts[[7]], width=10, height = 7.78)


# data.linear <- data.frame("x" = as.vector(as.matrix(xholder)),
#                           "true" = as.vector(as.matrix(fwi.scaled)),
#                           "post.mean" = as.vector(g.linear.mean),
#                           "q1" = as.vector(g.linear.q1),
#                           "q2" = as.vector(g.linear.q2),
#                           "q3" = as.vector(g.linear.q3),
#                           "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


# grid.plts <- list()
# for(i in 1:p){
#   grid.plt <- ggplot(data = data.frame(data.linear[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
#                   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#                   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#                   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#                   geom_rug(aes(x=true, y=q2), sides = "b") +
#                   ylab("") + xlab(names(fwi.scaled)[i]) +
#                   scale_fill_manual(values=c("steelblue"), name = "") + 
#                   scale_color_manual(values=c("steelblue")) +
#                   ylim(g.min.samples, g.max.samples) +
#                   theme_minimal(base_size = 30) +
#                   theme(legend.position = "none",
#                           plot.margin = margin(0,0,0,-20),
#                           axis.text = element_text(size = 35),
#                           axis.title.x = element_text(size = 45))
#   grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=g.min.samples, color = "red", size = 7)
# }

# grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)


# data.nonlinear <- data.frame("x" = as.vector(as.matrix(xholder)),
#                           "true" = as.vector(as.matrix(fwi.scaled)),
#                           "post.mean" = as.vector(g.nonlinear.mean),
#                           "q1" = as.vector(g.nonlinear.q1),
#                           "q2" = as.vector(g.nonlinear.q2),
#                           "q3" = as.vector(g.nonlinear.q3),
#                           "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


# grid.plts <- list()
# for(i in 1:p){
#   grid.plt <- ggplot(data = data.frame(data.nonlinear[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
#                   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#                   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#                   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#                   geom_rug(aes(x=true, y=q2), sides = "b") +
#                   ylab("") + xlab(names(fwi.scaled)[i]) +
#                   scale_fill_manual(values=c("steelblue"), name = "") + 
#                   scale_color_manual(values=c("steelblue")) +
#                   ylim(g.min.samples, g.max.samples) +
#                   theme_minimal(base_size = 30) +
#                   theme(legend.position = "none",
#                           plot.margin = margin(0,0,0,-20),
#                           axis.text = element_text(size = 35),
#                           axis.title.x = element_text(size = 45))
#   grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=g.min.samples, color = "red", size = 7)
# }

# grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

data.smooth <- data.frame("x" = as.vector(as.matrix(fwi.grid[,1:p])),
                          "true" = as.vector(as.matrix(fwi.origin[,c(1:p)])),
                          "post.mean" = as.vector(fwi.smooth.mean),
                          "q1" = as.vector(fwi.smooth.q1),
                          "q2" = as.vector(fwi.smooth.q2),
                          "q3" = as.vector(fwi.smooth.q3),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) + ylim(-3.3,3.3) +
                  # ylim(min(fwi.min.samples), 7) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.grid[which.max(y),i], y=-3.3, color = "red", size = 7)
}
grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)
# time.df <- data.frame(q1 = as.vector(time.samples[,4]), 
#                       q2 = as.vector(time.samples[,5]), 
#                       q3 = as.vector(time.samples[,6]), 
#                       x=fwi.scaled[excess,"time"])
# grid.plts[[p+1]] <- ggplot(data = time.df, aes(x=x)) + 
#                       geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#                       geom_ribbon(aes(ymin = q1, ymax = q3), fill= "steelblue", alpha = 0.2) +
#                       geom_line(aes(y=q2), color = "steelblue", linewidth=1) + 
#                       # geom_rug(aes(x=true, y=q2), sides = "b") +
#                       ylab("") + xlab("Year") + ylim(-1,1) +
#                       theme_minimal(base_size = 30) +
#                       theme(legend.position = "none",
#                               plot.margin = margin(0,0,0,-20),
#                               axis.text = element_text(size = 35),
#                               axis.title.x = element_text(size = 45))


summary(fit1, par=c("theta_fwi"), probs = c(0.05,0.5, 0.95))$summary
# saveRDS(data.smooth, file="./BLAST/application/figures/comparison/full_stanfit.rds")

#Predictive Distribution check
# y.container <- as.data.frame(matrix(, nrow = n, ncol = 0))  
# random.alpha.idx <- floor(runif(100, 1, ncol(t(posterior$f))))
# for(i in random.alpha.idx){
#   y.container <- cbind(y.container, log(t(posterior$f)[,i]))
# }
# colnames(y.container) <- paste("col", 1:100, sep="")
# y.container$x <- seq(1,n)
# y.container$logy <- log(y)
# plt <- ggplot(data = y.container, aes(x = logy)) + ylab("Density") + xlab("log(Burned Area)") + labs(col = "")

# for(i in names(y.container)){
#   plt <- plt + geom_density(aes(x=.data[[i]]), color = "slategray1", alpha = 0.1, linewidht = 0.7)
# }

# print(plt + geom_density(aes(x=logy), color = "steelblue", linewidth = 2) +
#         theme_minimal(base_size = 30) + ylim(0, 1.25) + xlim(7.5,30) +
#         theme(legend.position = "none",
#                 axis.text = element_text(size = 35)))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_BLAST_predictive_distribution.pdf"), width=10, height = 7.78)


# extreme.container <- as.data.frame(matrix(, nrow = n, ncol = 3000))
# for(i in 1:3000){
#   extreme.container[,i] <- density(log(posterior$f[i,]), n=n)$y
# }
# extreme.container <- cbind(extreme.container, t(apply(extreme.container[,1:3000], 1, quantile, c(0.05, .5, .95))))
# colnames(extreme.container)[(dim(extreme.container)[2]-2):(dim(extreme.container)[2])] <- c("q1","q2","q3")
# colnames(extreme.container)[1:3000] <- as.character(1:3000)
# extreme.container$mean <- rowMeans(extreme.container[,1:3000])
# extreme.container$y <- seq(0, 30, length.out = n)
# extreme.container <- as.data.frame(extreme.container)


# plt <- ggplot(data = extreme.container, aes(x = y)) + xlab("log(Burned Area)") + ylab("Density")+
#         geom_line(aes(y=mean), colour = "steelblue", linewidth = 1.5) +
#         geom_ribbon(aes(ymin = q1, ymax = q3), fill = "steelblue", alpha = 0.2) + 
#         theme_minimal(base_size = 30) + 
#         theme(legend.position = "none",
#               axis.title = element_text(size = 30))
# d <- ggplot_build(plt)$data[[1]]
# print(plt + 
#         geom_segment(x=12.44009, xend=12.44009, 
#               y=0, yend=approx(x = d$x, y = d$y, xout = 12.4409)$y,
#               colour="red", linewidth=1.2, linetype = "dotted"))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_post_generative.pdf"), width = 10, height = 7.78)

# random.alpha.idx <- floor(runif(1000, 1, ncol(t(posterior$alpha))))
# ev.y1 <- ev.y2 <- as.data.frame(matrix(, nrow = 1, ncol = 0))
# ev.alpha.single <- c()  
# for(i in random.alpha.idx){
#   ev.y1 <- rbind(ev.y1, as.numeric(posterior$yrep[i]))
# }
# ev.y1 <- as.data.frame(log(ev.y1))
# ev.y1$logy <- max(log(y))
# colnames(ev.y1) <- c("yrep", "logy")
# ev.y1$group <- rep("15th Oct 2017",1000)
# ggplot(data=ev.y, aes(x=yrep, y = group)) +
#   ylab("") + 
#   xlab("log(Burnt Area)") + labs(col = "") +  
#   stat_slab(scale = 0.6, colour = "steelblue", fill=NA, slab_linewidth = 1.5, trim = FALSE, expand = TRUE, density = "unbounded", subguide="outside", justification = -0.01) +
#   # stat_spike(aes(linetype = after_stat(at)), at = c("median"), scale=0.7)+
#   stat_dotsinterval(subguide = 'integer', side = "bottom", scale = 0.6, slab_linewidth = NA, position = "dodge") +
#   # geom_point(position = position_jitter(seed = 1, height = 0.05), alpha = 0.1) +  
#   # geom_boxplot(width = 0.2, notch = TRUE, alpha = 0.25, outlier.color = NA) +
#   geom_vline(xintercept = log(max(y)), linetype="dashed", color = "red",) +
#   # geom_label(aes(log(max(y)), 1), label = "Target Length", show.legend = FALSE)+
#   geom_vline(xintercept = log(y[133]), linetype="dashed", color = "black",) +
#   # geom_label(aes(log(y[133]), 1), label = "Target Length", show.legend = FALSE)+
#   theme_minimal(base_size = 30) +  
#   theme(legend.position = "none",
#         plot.margin = margin(0,0,0,25),
#         axis.text.y = element_text(angle = 90, size = 15, vjust = 15, hjust = 0.5),
#         axis.title = element_text(size = 30)) +
#         annotate(x=(log(max(y))+2), y= 0.1, label = "15th Oct 2017", geom="label") +
#         annotate(x=(log(y[133])-2), y= 0.1, label = "18th Jun 2017", geom="label")
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_BLAST_two_generative.pdf"), width = 10, height = 7.78)

# plt <- ggplot(data = ev.y1, aes(x = yrep)) + ylab("Density") + xlab("log(Burned Area)") + labs(col = "") +
#   geom_density(color = "steelblue", linewidth = 1.2) + 
#   geom_rug(alpha = 0.1) + 
#   xlim(5.5, 40) +
#   theme_minimal(base_size = 30) +  
#   theme(legend.position = "none",
#         axis.title = element_text(size = 30))

# d <- ggplot_build(plt)$data[[1]]
# print(plt + geom_area(data = subset(d, x>12.44009), aes(x=x,y=y), fill = "slategray1", alpha = 0.5) +
#         geom_segment(x=12.44009, xend=12.44009, 
#               y=0, yend=approx(x = d$x, y = d$y, xout = 12.4409)$y,
#               colour="red", linewidth=1.2, linetype = "dotted"))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_BLAST_generative.pdf"), width = 10, height = 7.78)

# library(ismev)
# gpd.fit(y-u, u)

# fit.log.lik <- extract_log_lik(fit1)
# constraint.elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 2)
# save(constraint.elpd.loo, constraint.waic, file = (paste0("./BLAST/application/BLAST_constraint_",Sys.Date(),"_",psi,"_",floor(threshold*100),"quantile_IC.Rdata")))


# x_data <- sort(Y[Y > 0], decreasing = TRUE)
# n <- length(x_data)
# k_range_pick <- 15:500
# pick_res <- data.frame(k = k_range_pick, xi = NA, lower = NA, upper = NA)

# for(i in seq_along(k_range_pick)) {
#   k <- k_range_pick[i]
  
#   # Quantiles: X_{n-k+1}, X_{n-2k+1}, X_{n-4k+1}
#   q1 <- x_data[k]
#   q2 <- x_data[2*k]
#   q3 <- x_data[4*k]
  
#   # 1. Point Estimate
#   # Formula: (1/ln2) * ln( (q1-q2) / (q2-q3) )
#   xi_hat <- (1/log(2)) * log((q1 - q2) / (q2 - q3))
  
#   # 2. Asymptotic Variance Formula for Pickands
#   # Var = [ xi^2 * (2^(2xi+1) + 1) ] / [ (2 * (2^xi - 1) * ln2)^2 ]
#   # Handle the xi=0 case to prevent division by zero (limit is approx 3.24)
#   if(abs(xi_hat) < 1e-6) {
#     asy_var <- 3.24 # Approx limit
#   } else {
#     num <- (xi_hat^2) * ((2^(2*xi_hat + 1)) + 1)
#     den <- (2 * (2^xi_hat - 1) * log(2))^2
#     asy_var <- num / den
#   }
  
#   # Standard Error = sqrt(Var / k)
#   se <- sqrt(asy_var / k)
  
#   pick_res$xi[i] <- xi_hat
#   pick_res$lower[i] <- xi_hat - 1.96 * se
#   pick_res$upper[i] <- xi_hat + 1.96 * se
# }

# # --- 3. MOMENT (DedH) ESTIMATOR ---
# # Constraint: Uses all k up to n-1
# k_range_mom <- 15:500
# mom_res <- data.frame(k = k_range_mom, xi = NA, lower = NA, upper = NA)

# for(i in seq_along(k_range_mom)) {
#   k <- k_range_mom[i]
#   y_k <- x_data[1:k]
  
#   # 1. Point Estimate
#   log_excess <- log(y_k) - log(x_data[k+1])
#   M1 <- mean(log_excess)
#   M2 <- mean(log_excess^2)
  
#   xi_hat <- M1 + 1 - 0.5 * (1 - (M1^2 / M2))^(-1)
  
#   # 2. Asymptotic Variance (for xi > 0)
#   # Var = 1 + xi^2
#   asy_var <- 1 + xi_hat^2
#   se <- sqrt(asy_var / k)
  
#   mom_res$xi[i] <- xi_hat
#   mom_res$lower[i] <- xi_hat - 1.96 * se
#   mom_res$upper[i] <- xi_hat + 1.96 * se
# }


# p1 <- ggplot(pick_res, aes(x=k, y=xi)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), fill="darkgreen", alpha=0.2) +
#   geom_line(color="darkgreen", linewidth=1) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_cartesian(ylim=c(-0.15, 1.5)) + # Zoom to relevant range
#   labs(title="Pickands Estimator", 
#        y="Extreme Value Index") +
#   theme_minimal(base_size = 30) + ylim(-10,10) +
#   theme(legend.position = "none",
#           axis.text = element_text(size = 35),
#           axis.title.x = element_text(size = 45))

# p2 <- ggplot(mom_res, aes(x=k, y=xi)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), fill="purple", alpha=0.2) +
#   geom_line(color="purple", linewidth=1) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_cartesian(ylim=c(-0.15, 1.5)) + 
#   labs(title="Moment Estimator", y="") +
#   theme_minimal(base_size = 30) + ylim(0, 2) +
#   theme(legend.position = "none",
#           axis.text = element_text(size = 35),
#           axis.text.y = element_blank(),
#           axis.title.x = element_text(size = 45))

# grid.plt <- grid.arrange(p1, p2, nrow=1)
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_heavytail.pdf"), grid.plt, width=22, height = 7.78)
fit.log.lik <- extract_log_lik(fit1)
elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 4)
elpd.loo
# save(elpd.loo, file = (paste0("./BLAST/application/BLAST_full_",Sys.Date(),"_",psi+2,"_",floor(threshold*100),"quantile_IC.Rdata")))
