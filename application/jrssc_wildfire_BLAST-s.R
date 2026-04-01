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
library(forecast)
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

#NAs on Feb 29 most years, and Feb 14, 1999
#considering the case of leap year, the missing values are the 29th of Feb
#Thus, each year consist of 366 data with either 1 or 0 missing value.
Y <- df.long$measurement[!is.na(df.long$measurement)]
psi.origin <- psi <- 10
threshold <- 0.95

multiplesheets <- function(fname) {
    setwd("C:/Users/Johnny Lee/Documents/GitHub")
    # getting info about all excel sheets
    sheets <- excel_sheets(fname)
    tibble <- lapply(sheets, function(x) read_excel(fname, sheet = x, col_types = c("date", rep("numeric", 41))))
    data_frame <- lapply(tibble, as.data.frame)
    # assigning names to data frames
    names(data_frame) <- sheets
    return(data_frame)
}
setwd("C:/Users/Johnny Lee/Documents/GitHub")
path <- "./BLAST/application/DadosDiariosPT_FWI.xlsx"
# importing fire weather index
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
    fwi.scaled[,i] <- fwi.index[,i] <- cov.long$measurement[missing.values]
}

# era5 <- read_excel("./BLAST/application/ERA_5.xlsx")
# era5 <- era5[era5$year>1979,]
# era5 <- era5[!(era5$year == 1999 & era5$month == 2 & era5$day == 14), ]
# fwi.index$ERA5 <- fwi.scaled$ERA5 <- as.numeric(era5$ERA_5)
fwi.scaled$time <- fwi.index$time <- seq(1,length(Y), length.out=length(Y))
fwi.scaled$date <- as.Date(fwi.scaled$time - 1, origin = "1980-01-01")
# head(fwi.scaled[, c("time", "date")])
fwi.scaled$sea <- fwi.index$sea <- fwi.index$time %% 365.25 / 365.25
fwi.scaled$cos.time <- fwi.index$cos.time <- cos(2*pi*seq(1,length(Y), length.out=length(Y))/365.25)
fwi.scaled$sin.time <- fwi.index$sin.time <- sin(2*pi*seq(1,length(Y), length.out=length(Y))/365.25)
fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)

# time_arima <- seq(1, length(Y))
# xreg_trend_seasonal <- cbind(
#   trend = time_arima,
#   cos_season = cos(2 * pi * time_arima / 365.25),
#   sin_season = sin(2 * pi * time_arima / 365.25)
# )

# fit.list <- list()
# for (j in 1:7) {
#   y_ts <- ts(fwi.scaled[, j], frequency = 365.25)
#   fit.list[[j]] <- forecast::auto.arima(
#     y_ts,
#     xreg = xreg_trend_seasonal,
#     seasonal = FALSE,
#     stepwise = TRUE,
#     approximation = FALSE
#   )
#   fwi.index[, j] <- fwi.scaled[, j] <- as.numeric(residuals(fit.list[[j]]))
# }

xreg.season <- cbind(
  trend = c(1:length(Y)),
  cos_season = cos(2 * pi * c(1:length(Y)) / 365.25),
  sin_season = sin(2 * pi * c(1:length(Y)) / 365.25)
)

fit.list <- list()
x.detrended <- matrix(nrow = length(Y), ncol = 7)
for (j in 1:7) {
  y_ts <- ts(fwi.scaled[, j], frequency = 365.25) 
  fit.list[[j]] <- auto.arima(y_ts, seasonal = FALSE, xreg = xreg.season, stepwise = TRUE, approximation = FALSE, max.d = 0)
  x.detrended[, j] <- as.numeric(residuals(fit.list[[j]]))
}
fwi.index[,1:7] <- fwi.scaled[, 1:7] <- x.detrended
# acf(fwi.index$BUI)
# acf(fwi.index$ISI)
# acf(fwi.index$FFMC)
# acf(fwi.index$DMC)
# acf(fwi.index$DC)

above.0 <- which(Y > 0)
Y_pos <- Y[above.0]
fwi_pos <- fwi.scaled[above.0, ]
# Y[!above.0] <- 1e-5
# Y_pos <- Y
# fwi_pos <- fwi.scaled
# pca_result <- prcomp(fwi_pos[,3:7], center = TRUE, scale. = TRUE)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# fwi_pos[,1:7] <- as.data.frame(sapply(fwi_pos[,1:7], FUN = range01))
# qr.df <- data.frame(y = log(Y_pos), pca_result$x, cos.time = fwi_pos$cos.time, sin.time = fwi_pos$sin.time) #fwi_pos)
qr.df <- data.frame(y = log(Y_pos), (fwi_pos[,1:7]), cos.time = fwi_pos$cos.time, sin.time = fwi_pos$sin.time, time = fwi_pos$sea)
# evgam.cov <- y ~ 1 + cos.time + sin.time + s(PC1, k=15) + s(PC2, k=15) + s(PC3, k=15) + s(PC4, k=15) + s(PC5, k=15)
# evgam.cov <- y ~ s(time, bs="cc", k=5) + s(BUI, bs="ts", k = 5) + s(ISI, bs="ts", k = 5) + s(FFMC, bs="ts", k = 5) + s(DMC, bs="ts", k = 5) + s(DC, bs="ts", k = 5)
s.cov <- c(3:7)
evgam.cov <- as.formula(paste0("y ~ cos.time + sin.time + ", paste0("s(", colnames(fwi_pos[,s.cov]), ", k = ", psi+2, ", bs='" ,"ts')", collapse = " + ")))

qr.fit <- evgam(evgam.cov, data = qr.df, family = "ald", ald.args=list(tau = threshold))

# qr.lin <- evgam(y ~ cos.time + sin.time + BUI + ISI + FFMC + DMC + DC, data = qr.df, family = "ald", ald.args=list(tau = threshold))
# qr.cov <- evgam(y ~ s(BUI, bs = "ts", k = 30) + s(ISI, bs = "ts", k = 30) + s(FFMC, bs = "ts", k = 30) + s(DMC, bs = "ts", k = 30) + s(DC, bs = "ts", k = 30), data = qr.df, family = "ald", ald.args=list(tau = threshold))
# qr.time <- evgam(y ~ cos.time + sin.time, data = qr.df, family = "ald", ald.args=list(tau = threshold))
# qr.null <- evgam(y ~ 1, data = qr.df, family = "ald", ald.args=list(tau = threshold))
u.vec <- exp(predict(qr.fit)$location)
# save(u.c, qr.fit, file = paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_qr-c.Rdata"))
# load("./BLAST/application/figures/2026-03-25_pareto_qr-ct.Rdata")
# u.vec <- u.ct
# qr.fit <- quantreg::rq(y ~ 1 + cos.time + sin.time + BUI + ISI + FFMC + DMC + DC, data = qr.df, tau = threshold)
# u.vec <- exp(predict(qr.fit))  # threshold on raw scale for Y_pos
# AIC(qr.fit, qr.cov, qr.time, qr.null)
# BIC(qr.fit, qr.cov, qr.time, qr.null)

plot(fwi.scaled[above.0,"date"], log(Y_pos))
lines(fwi.scaled[above.0,"date"], log(u.vec), type = "l", col = "red")


excess.idx <- which(Y_pos > u.vec)
y <- Y_pos[excess.idx]
u <- u.vec[excess.idx]
# fwi.scaled <- fwi_pos
fwi.01 <- fwi_pos[excess.idx, s.cov] # BUI to DC

plot(fwi_pos[excess.idx, "date"], log(y), xlab= "Year")
lines(fwi_pos[excess.idx, "date"], log(u), type = "l", col = "red", xlab= "Year")

fwi.cols <- colnames(fwi_pos[,s.cov])
fwi.max <- fwi.01[which.max(y),]
X_minmax <- sapply(fwi.01, function(x) max(x)-min(x))
X_min <- sapply(fwi.01, function(x) min(x))
fwi.01 <- as.data.frame(sapply(fwi.01, FUN = range01))
grid.n <- n <- nrow(fwi.01)
p <- ncol(fwi.01)
psi <- psi.origin - 2 # Adjusted basis dimension

fwi.grid <- data.frame(lapply(fwi_pos[,s.cov], function(x) seq(min(x), max(x), length.out =grid.n)))



Z.list <- list()
projection_coefs_list <- list()
spec_decomp_list <- list()
sm_spec_list <- list()
scale_stats_list <- list()
keep_cols_list <- list()

for (i in 1:p) {
  x_vec <- as.numeric(fwi.01[, i])
  X_lin <- model.matrix(~ x_vec)
  sm_spec <- smoothCon(s(x_vec, bs = "tp", k = psi + 2), data = data.frame(x_vec))[[1]]
  
  S <- sm_spec$S[[1]]
  eig <- eigen(S, symmetric = TRUE)
  pos_idx <- which(eig$values > max(eig$values) * 1e-8)
  
  U_pen <- eig$vectors[, pos_idx]
  Lambda_inv <- solve(sqrt(diag(eig$values[pos_idx])))
  
  Z_spectral <- sm_spec$X %*% U_pen %*% Lambda_inv
  qr_lin <- qr(X_lin)
  Gamma_Original <- backsolve(qr.R(qr_lin), t(qr.Q(qr_lin)) %*% Z_spectral)
  Z_orth <- Z_spectral - X_lin %*% Gamma_Original
  keep_cols <- colSums(Z_orth^2) > 1e-9
  Z_final <- Z_orth[, keep_cols, drop = FALSE]

  z_scale <- apply(Z_final, 2, sd)
  Z_final <- scale(Z_final, center = FALSE, scale = z_scale)
  
  Z.list[[i]] <- Z_final
  spec_decomp_list[[i]] <- list(U = U_pen, L_inv = Lambda_inv)
  projection_coefs_list[[i]] <- Gamma_Original
  scale_stats_list[[i]] <- z_scale
  keep_cols_list[[i]] <- keep_cols
  sm_spec_list[[i]] <- sm_spec
}

bs.linear <- as.matrix(fwi.01)
bs.nonlinear <- do.call(cbind, Z.list)


newx <- seq(0, 1, length.out = grid.n) # Standard Deviation units
xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))
colnames(xholder) <- fwi.cols
xholder.linear <- matrix(rep(newx, p), ncol = p)

grid_Z_list <- list()

for (i in 1:p) {
  grid_df <- data.frame(x_vec = newx)
  X_raw_grid <- PredictMat(sm_spec_list[[i]], grid_df)
  Z_spectral_grid <- X_raw_grid %*% spec_decomp_list[[i]]$U %*% spec_decomp_list[[i]]$L_inv
  X_lin_grid <- model.matrix(~ x_vec, data = grid_df)
  Z_orth_grid <- Z_spectral_grid - X_lin_grid %*% projection_coefs_list[[i]]
  Z_final_grid <- Z_orth_grid[, keep_cols_list[[i]], drop = FALSE]
  grid_Z_list[[i]] <- scale(Z_final_grid, center = FALSE, scale = scale_stats_list[[i]])
}
xholder.nonlinear <- do.call(cbind, grid_Z_list)

grid_Z_list <- list()

for (i in 1:p) {
  scaled_x_grid <- (fwi.grid[,i] - X_min[i]) / (X_minmax[i])
  grid_df  <- data.frame(x_vec = scaled_x_grid)
  X_raw_grid <- PredictMat(sm_spec_list[[i]], grid_df)
  Z_spectral_grid <- X_raw_grid %*% spec_decomp_list[[i]]$U %*% spec_decomp_list[[i]]$L_inv
  X_lin_grid <- model.matrix(~ x_vec, data = grid_df)
  Z_orth_grid <- Z_spectral_grid - X_lin_grid %*% projection_coefs_list[[i]]
  Z_final_grid <- Z_orth_grid[, keep_cols_list[[i]], drop = FALSE]
  grid_Z_list[[i]] <- scale(Z_final_grid, center = FALSE, scale = scale_stats_list[[i]])
}
xgrid.linear <- fwi.grid
xgrid.nonlinear <- do.call(cbind, grid_Z_list)

X_means <- colMeans(bs.linear)
X_sd   <- apply(bs.linear, 2, sd)
bs.linear <- scale(bs.linear, center = X_means, scale = X_sd)

model.stan <- "
data {
  int n; int grid_n; int p; int psi;
  vector[n] u; 
  vector[n] y;
  matrix[n, p] bsLinear; matrix[n, p*psi] bsNonlinear;
  matrix[n, p] gridL; matrix[n, p*psi] gridNL;
  matrix[grid_n, p] xholderLinear; matrix[grid_n, p*psi] xholderNonlinear;
  real atau; vector[p] X_means; vector[p] X_sd; vector[p] X_minmax;
}
parameters {
  vector[p+1] theta;
  array[p] vector[psi] gamma_raw;
  vector<lower=0>[p] lambda1; 
  vector<lower=0>[p] lambda2;
  vector<lower=0>[p] tau;
}

transformed parameters {
  vector[n] alpha;
  array[p] vector[psi] gamma;
  {
    vector[n] eta = rep_vector(theta[1], n);
    for (j in 1:p) {
      gamma[j] = gamma_raw[j] * sqrt(tau[j]);
      eta += col(bsLinear, j) * theta[j+1] + block(bsNonlinear, 1, (j-1)*psi + 1, n, psi) * gamma[j];
    }
    alpha = exp(eta);
  }
}

model {
  target += pareto_lpdf(y | u, alpha);
  target += normal_lpdf(theta[1] | 0, 10);
  target += gamma_lpdf(lambda1 | 1e-2, 1e-2);
  // target += exponential_lpdf(lambda1 | 0.1); 
  target += gamma_lpdf(lambda2 | 1e-2, 1e-2);
  for (j in 1:p) {
    target += double_exponential_lpdf(theta[j+1] | 0, 1/lambda1[j]);
    target += gamma_lpdf(tau[j] | atau, 0.5 * square(lambda2[j]));
    target += std_normal_lpdf(gamma_raw[j]);
  }
}

generated quantities {
  vector[n] log_lik;
  vector[grid_n] gridalpha;
  matrix[grid_n, p] gridgsmooth;
  vector[p+1] theta_origin;
  vector[p] theta_fwi;

  theta_origin[2:(p+1)] = theta[2:(p+1)] ./ X_sd;
  theta_origin[1] = theta[1] - dot_product(X_means, theta_origin[2:(p+1)]);
  theta_fwi[1:p] = theta_origin[2:(p+1)] ./ X_minmax;


  for (i in 1:n) log_lik[i] = pareto_lpdf(y[i] | u[i], alpha[i]);
  {
    vector[grid_n] pred = rep_vector(theta_origin[1], grid_n);
    for (j in 1:p) {
      gridgsmooth[,j] = col(xholderLinear, j) * theta_origin[j+1] + block(xholderNonlinear, 1, (j-1)*psi + 1, grid_n, psi) * gamma[j]; // gridgsmooth[,j] = col(gridL, j) * theta_fwi[j] + block(gridNL, 1, (j-1)*psi + 1, grid_n, psi) * gamma[j]
      pred += col(xholderLinear, j) * theta_origin[j+1] + block(xholderNonlinear, 1, (j-1)*psi + 1, grid_n, psi) * gamma[j];
    }
    gridalpha = exp(pred);
  }
}
"

data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi,
                  atau = ((psi+1)/2), xholderNonlinear = xholder.nonlinear, 
                  bsLinear = bs.linear, bsNonlinear = bs.nonlinear, #X_min=fwi.min,
                  xholderLinear = xholder.linear, X_minmax = X_minmax, 
                  X_means = X_means, X_sd = X_sd, grid_n = grid.n,
                  gridL = xgrid.linear, gridNL = xgrid.nonlinear)

init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p, psi)), 
                        theta = rep(-0.1, (p+1)), tau = rep(0.1, p), 
                        lambda1 = rep(0.1,p), lambda2 = rep(1, p)),
                   list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p, psi)),
                        theta = rep(0.05, (p+1)), tau = rep(0.2, p),
                        lambda1 = rep(2,p), lambda2 = rep(5, p)),
                   list(gamma_raw= array(rep(0.1, (psi*p)), dim=c(p, psi)),
                        theta = rep(0.1, (p+1)), tau = rep(0.1, p), 
                        lambda1 = rep(0.1,p), lambda2 = rep(0.1, p)))

blast_model <- stan_model(model_code = model.stan)
map_fit <- optimizing(
  blast_model,
  data = data.stan, 
  init = init.alpha[[1]],
  as_vector = FALSE, # Returns a structured list
  algorithm = "BFGS",
  iter = 10000        # Give it plenty of iterations to converge
)

init.alpha <- list(list(gamma_raw= array(map_fit$par$gamma_raw, dim=c(p, psi)), 
                        theta = map_fit$par$theta, tau = map_fit$par$tau, 
                        lambda1 = map_fit$par$lambda1, lambda2 = map_fit$par$lambda2),
                    list(gamma_raw= array(map_fit$par$gamma_raw, dim=c(p, psi)), 
                          theta = map_fit$par$theta, tau = map_fit$par$tau, 
                          lambda1 = map_fit$par$lambda1, lambda2 = map_fit$par$lambda2),
                    list(gamma_raw= array(map_fit$par$gamma_raw, dim=c(p, psi)), 
                    theta = map_fit$par$theta, tau = map_fit$par$tau, 
                    lambda1 = map_fit$par$lambda1, lambda2 = map_fit$par$lambda2))

fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value 
    chains = 3,             # number of Markov chains
    iter = 4000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 2000           # no progress shown
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
# fwismooth.samples <- summary(fit1, par=c("fwismooth"), probs = c(0.05, 0.5, 0.95))$summary
# gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
# gridgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
# origin.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary
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
# fwi.smooth.mean <- as.vector(matrix(fwismooth.samples[,1], nrow = n, byrow=TRUE))
# fwi.smooth.q1 <- as.vector(matrix(fwismooth.samples[,4], nrow = n, byrow=TRUE))
# fwi.smooth.q2 <- as.vector(matrix(fwismooth.samples[,5], nrow = n, byrow=TRUE))
# fwi.smooth.q3 <- as.vector(matrix(fwismooth.samples[,6], nrow = n, byrow=TRUE))
# fwi.smooth.mean <- as.vector(apply(fwismooth_samples, c(2, 3), mean))
# fwi.smooth.q2 <- as.vector(apply(fwismooth_samples, c(2, 3), median))
# fwi.smooth.q1 <- as.vector(apply(fwismooth_samples, c(2, 3), quantile, probs = 0.05))
# fwi.smooth.q3 <- as.vector(apply(fwismooth_samples, c(2, 3), quantile, probs = 0.95))

g.min.samples <- min(gsmooth.samples[,4])
g.max.samples <- max(gsmooth.samples[,6])
# fwi.smooth <- as.data.frame(apply(fwismooth_samples, c(2, 3), quantile, probs = 0.05))
# fwi.smooth <- as.data.frame(matrix(fwismooth.samples[,4], nrow = n, byrow=TRUE))
# fwi.min.samples <- sapply(fwi.smooth, min)
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
xholder <- as.data.frame(xholder)
# colnames(xholder) <- colnames(fwi.scaled)[1:p]
simul.data <- data.frame(BA = y-u, fwi.01)#fwi.origin[c(1:p)])
# gam.scale <- list(BA ~ s(DSR, bs = "tp", k = 30) + 
#                       s(FWI, bs = "tp", k = 30) + 
#                       s(BUI, bs = "tp", k = 30) + 
#                       s(ISI, bs = "tp", k = 30) + 
#                       s(FFMC, bs = "tp", k = 30) +
#                       s(DMC, bs = "tp", k = 30) + 
#                       s(DC, bs = "tp", k = 30),
#                     ~ s(DSR, bs = "tp", k = 30) + 
#                       s(FWI, bs = "tp", k = 30) + 
#                       s(BUI, bs = "tp", k = 30) + 
#                       s(ISI, bs = "tp", k = 30) + 
#                       s(FFMC, bs = "tp", k = 30) +
#                       s(DMC, bs = "tp", k = 30) +
#                       s(DC, bs = "tp", k = 30))
# gam.scale <- list(BA ~ s(BUI, bs = "ts", k = 30) + 
#                       s(ISI, bs = "ts", k = 30) + 
#                       s(FFMC, bs = "ts", k = 30) +
#                       s(DMC, bs = "ts", k = 30) + 
#                       s(DC, bs = "ts", k = 30),
#                     ~ s(BUI, bs = "ts", k = 30) + 
#                       s(ISI, bs = "ts", k = 30) + 
#                       s(FFMC, bs = "ts", k = 30) +
#                       s(DMC, bs = "ts", k = 30) +
#                       s(DC, bs = "ts", k = 30))
# evgam.fit.scale <- evgam::evgam(gam.scale, data = simul.data, family = "gpd")
# pred.scale <- predict(evgam.fit.scale, newdata = xholder, type="response", se.fit = TRUE)
# xi.pred.scale <-pred.scale$fitted$shape
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
                ~ s(BUI, bs = "ts", k = 10) + 
                  s(ISI, bs = "ts", k = 10) + 
                  s(FFMC, bs = "ts", k = 10) +
                  s(DMC, bs = "ts", k = 10) +
                  s(DC, bs = "ts", k = 10))
evgam.fit.1 <- evgam::evgam(gam.1, data = simul.data, family = "gpd")
pred.1 <- predict(evgam.fit.1, newdata = xholder, type="response", se.fit = TRUE)
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


data.scenario <- data.frame("x" = seq(0,1,length.out=grid.n),
                            "post.mean" = (alpha.samples[,1]),
                            "post.median" = (alpha.samples[,5]),
                            "q1" = (alpha.samples[,4]),
                            "q3" = (alpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,ldots,c))) + xlab(expression(c)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), fill="steelblue", alpha = 0.2) +
  geom_line(aes(y=post.median), colour = "steelblue", linewidth=1) +
  theme_minimal(base_size = 30) + scale_y_log10() +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_alpha.pdf"), width=10, height = 7.78)

xi.scenario <- data.frame("x" = seq(0,1,length.out=grid.n),
                            "post.mean" = (1/alpha.samples[,1]),
                            "post.median" = (1/alpha.samples[,5]),
                            "q1" = (1/alpha.samples[,4]),
                            "q3" = (1/alpha.samples[,6]))
                            # "evgam.1" = xi.pred.1,
                            # "evgam.1.q1" = xi.low.1,
                            # "evgam.1.q3" = xi.high.1,
                            # "vgam.1" = vgam.xi.1,
                            # "vgam.scale" = vgam.xi.scale,
                            # "evgam.scale.q1" = xi.low.scale,
                            # "evgam.scale.q3" = xi.high.scale,
                            # "evgam.scale" = xi.pred.scale)

ggplot(xi.scenario, aes(x=x)) + 
  ylab(expression(xi(c,ldots,c))) + xlab(expression(c)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), fill="steelblue", alpha = 0.2) +
  geom_line(aes(y=post.median), color = "steelblue", linewidth=1) +
  # geom_ribbon(aes(ymin = evgam.scale.q1, ymax = evgam.scale.q3), fill= "orange", alpha = 0.2) +
  # geom_line(aes(y=evgam.scale), colour = "orange", linewidth=1, linetype=3) +
  # geom_ribbon(aes(ymin = evgam.1.q1, ymax = evgam.1.q3), fill= "purple", alpha = 0.2) +
  # geom_line(aes(y=evgam.1), colour = "purple", linewidth=1, linetype=3) +
  # geom_line(aes(y=vgam.scale), colour = "orange", linewidth=1, linetype=4) +
  # geom_line(aes(y=vgam.1), colour = "purple", linewidth=1, linetype=4) +  
  theme_minimal(base_size = 30) + scale_y_log10() +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_xi.pdf"), width=10, height = 7.78)


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

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_qqplot.pdf"), width=10, height = 7.78)
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


data.smooth <- data.frame("x" = as.vector(as.matrix(fwi.grid)),
                          "true" = as.vector(as.matrix(fwi_pos[excess.idx,s.cov])),
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
                  ylab("") + xlab(names(fwi.01)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  # ylim(g.min.samples, g.max.samples) +
                  # ylim(-10, 10) +
                  theme_minimal(base_size = 20) +
                  theme(legend.position = "none",
                        plot.margin = margin(5, 5, 5, 5),
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        axis.text.y = element_text(size = 18),
                        axis.text.x = element_text(size = 12),
                        axis.title.x = element_text(size = 22))
                  # theme_minimal(base_size = 30) +
                  # theme(legend.position = "none",
                  #         plot.margin = margin(0,0,0,-20),
                  #         axis.text = element_text(size = 35),
                  #         axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.max[,i], y=g.min.samples, color = "red", size = 7)
}
marrangeGrob(grobs = grid.plts, nrow = 1, ncol = p, top = NULL)
# grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_cov.pdf"), marrangeGrob(grobs = grid.plts, nrow = 1, ncol = 5, top = NULL), width=10, height = 7.78)


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
