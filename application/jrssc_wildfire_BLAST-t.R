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
library(mboost)
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
psi.origin <- psi <- 30
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
    fwi.index[,i] <- cov.long$measurement[missing.values]
    fwi.scaled[,i] <- cov.long$measurement[missing.values]
}

# era5 <- read_excel("./BLAST/application/ERA_5.xlsx")
# era5 <- era5[era5$year>1979,]
# era5 <- era5[!(era5$year == 1999 & era5$month == 2 & era5$day == 14), ]
# fwi.index$ERA5 <- fwi.scaled$ERA5 <- as.numeric(era5$ERA_5)
fwi.scaled$time <- fwi.index$time <- c(1:length(Y))
fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)



above.0 <- which(Y>0)
Y <- Y[above.0]
fwi.scaled <- fwi.scaled[above.0,]
fwi.index <- fwi.index[above.0,]

fwi.season <- fwi.index %>%
  mutate(
    year = as.numeric(as.character(year)), 
    Month_Num = match(month, month.abb), 
    
    season = case_when(
      Month_Num %in% c(12, 1, 2) ~ "Winter",
      Month_Num %in% c(3, 4, 5) ~ "Spring",
      Month_Num %in% c(6, 7, 8) ~ "Summer",
      Month_Num %in% c(9, 10, 11) ~ "Autumn"
    ),

    SeasonYear = ifelse(Month_Num == 12, year + 1, year)
  )
fwi.season <- fwi.season %>% 
  mutate(Label = paste(SeasonYear, season)) %>%
  mutate(Label = factor(Label, levels = unique(Label)))

fwi.season$code <- factor(fwi.season$season, levels=c("Winter", "Spring", "Summer", "Autumn"), labels = c(1,2,3,4))
fwi.season$winter <- ifelse(fwi.season$season == "Winter", 1, 0)
fwi.season$spring <- ifelse(fwi.season$season == "Spring", 1, 0)
fwi.season$summer <- ifelse(fwi.season$season == "Summer", 1, 0)
fwi.season$autumn <- ifelse(fwi.season$season == "Autumn", 1, 0)
fwi.season$BA <- Y
fwi.season$log.BA <- log(fwi.season$BA)
fwi.season$cos.time <- cos(2*pi*fwi.season$time / 365)
fwi.season$sin.time <- sin(2*pi*fwi.season$time / 365)

fwi.scaled <- fwi.df <- fwi.season[,c(3,4,5,6,7,16)]
that.season <- which(fwi.scaled$code==2)
fwi.scaled <- fwi.scaled[that.season,]
Y <- Y[that.season]

# load(paste0("./BLAST/application/qr-",threshold*1000,"-t.Rdata"))

u <- quantile(Y, threshold)
# preds <- preds.evgam
excess <- which(Y>u)
y <- Y[excess]
u <- rep(u, length(y))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

p <- dim(fwi.scaled)[[2]]-1
fwi.scaled <- as.data.frame(sapply(fwi.scaled[excess,c(1:p)], FUN = range01))
# fwi.df <- fwi.df[excess,]
# fwi.scaled <- fwi.scaled[which(fwi.df$code==1),]
n <- dim(fwi.scaled)[[1]]

# fwi.origin <- data.frame(, BA=y)
# max.fwi <- fwi.origin[which.max(y),]
# fwi.grid <- data.frame(lapply(fwi.origin[,c(1:p)], function(x) seq(min(x), max(x), length.out = nrow(fwi.scaled))))
# fwi.minmax <- sapply(fwi.origin[,c(1:p)], function(x) max(x)-min(x))
# fwi.min <- sapply(fwi.origin[,c(1:p)], function(x) min(x))

bs.linear <- model.matrix(~ ., data = data.frame(fwi.scaled))
psi <- psi - 2
group.map <- c()
Z.list <- list()        # Stores the final non-linear design matrices
scale_stats_list <- list() 
projection_coefs_list <- list() #
spec_decomp_list <- list() # Store eigen-decomp info for prediction
qr_list <- list()          # Store QR info for prediction
sm_spec_list <- list()     # Store smooth objects
keep_cols_list <- list()

covariates <- colnames(fwi.scaled)
for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  x_vec <- fwi.scaled[, i]
  X_lin <- model.matrix(~ x_vec) 
  sm_spec <- smoothCon(mgcv::s(x_vec, bs = "tp", k = psi + 2), 
                      data = data.frame(x_vec = x_vec), 
                      knots = NULL)[[1]]
  
  X_raw <- sm_spec$X
  S     <- sm_spec$S[[1]] 
  
  eig <- eigen(S, symmetric = TRUE)
  max_lambda <- max(eig$values)
  tol <- max_lambda * 1e-8  # Relative threshold (Robust)
  
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
  
  # Store Transforms
  spec_decomp_list[[i]] <- list(U_pen = U_pen, Lambda_sqrt_inv = solve(sqrt(Lambda_pen)))
  projection_coefs_list[[i]] <- Gamma_Original # Store the X-basis coefs
  keep_cols_list[[i]] <- keep_cols
  scale_stats_list[[i]] <- train_scale         # Store the scale
  sm_spec_list[[i]] <- sm_spec
}

bs.nonlinear <- do.call(cbind, Z.list)

# xholder <- do.call(cbind, lapply(1:p, function(j) {seq(min(fwi.scaled[,j]),max(fwi.scaled[,j]), length.out = n)}))
newx <- seq(0, 1, length.out = n)
xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))

colnames(xholder) <- covariates
grid_Z_list <- list()

for (i in seq_along(covariates)) {
  grid_df  <- data.frame(x_vec = xholder[,i])
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

# grid_Z_list <- list()

# for (i in seq_along(covariates)) {
#   scaled_x_grid <- (fwi.grid[,i] - fwi.min[i]) / (fwi.minmax[i])
#   grid_df  <- data.frame(x_vec = scaled_x_grid)
#   X_lin_grid <- model.matrix(~ x_vec, data = grid_df)
#   X_raw_grid <- PredictMat(sm_spec_list[[i]], grid_df)
  
#   decomp <- spec_decomp_list[[i]]
#   Z_spectral_grid <- X_raw_grid %*% decomp$U_pen %*% decomp$Lambda_sqrt_inv
#   Z_orth_grid <- Z_spectral_grid - X_lin_grid %*% projection_coefs_list[[i]]
#   Z_final_grid <- Z_orth_grid[, keep_cols_list[[i]], drop = FALSE]
#   Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = scale_stats_list[[i]])
#   grid_Z_list[[i]] <- Z_final_grid
# }

# xgrid.nonlinear <- do.call(cbind, grid_Z_list)
# xgrid.linear <- model.matrix(~ ., data = data.frame(fwi.grid))
# xgrid.linear[,-1] <- sweep(xgrid.linear[, -1], 2, fwi.min, "-")
X_means <- colMeans(bs.linear[,-1])
X_sd   <- apply(bs.linear[,-1], 2, sd)
bs.linear[,-1] <- scale(bs.linear[,-1], center = X_means, scale = X_sd)

model.stan <- "// Stan model for BLAST Pareto Samples
data {
  int <lower=1> n; // Sample size
  int <lower=1> p; // regression coefficient size
  int <lower=1> psi; // splines coefficient size
  array[n] real <lower=0> u; // large threshold value
  matrix[n, p] bsLinear; // fwi dataset
  matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
  matrix[n, p] xholderLinear; // fwi dataset
  matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
  array[n] real <lower=0> y; // extreme response
  real <lower=0> atau;
  vector[(psi*p)] Z_scales;
  // vector[p] X_minmax;
  // vector[p] X_min;
  vector[p] X_means;
  vector[p] X_sd;
}

parameters {
  vector[(p+1)] theta; // linear predictor
  array[p] vector[psi] gamma_raw;
  array[p] real <lower=0> lambda1; // lasso penalty //array[p] real <lower=0> 
  array[p] real <lower=0> lambda2; // lambda2 group lasso penalty
  array[p] real <lower=0> tau;
}

transformed parameters {
  vector[n] alpha; // covariate-adjusted tail index
  
  array[p] vector[psi] gamma;
  {
    vector[n] eta = rep_vector(theta[1],n);
    for (j in 1:p){
      for (k in 1:psi){
        int idx = (j-1)*psi + k;
        gamma[j][k] = gamma_raw[j][k] * sqrt(tau[j]) * Z_scales[idx];
      }; 
      eta += bsLinear[,j] * theta[j+1] + bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    alpha = exp(eta);
  }
}

model {
  // likelihood
  target += pareto_lpdf(y | u, alpha);
  target += normal_lpdf(theta[1] | 0, 10);
  for (j in 1:p){
    target += gamma_lpdf(lambda1[j] | 1,1); 
    target += gamma_lpdf(lambda2[j] | 1e-2, 1e-2);  
    target += double_exponential_lpdf(theta[(j+1)] | 0, 1/(lambda1[j]));
    target += gamma_lpdf(tau[j] | atau, square(lambda2[j])*0.5);
    target += std_normal_lpdf(gamma_raw[j]);
  }
}

generated quantities {
  // Used in Posterior predictive check
  vector[n] gridalpha; // new tail index
  matrix[n, p] gridgnl; // nonlinear component
  matrix[n, p] gridgl; // linear component
  matrix[n, p] gridgsmooth; // linear component
  vector[n] log_lik;
  vector[p+1] theta_origin;
  // vector[p] theta_fwi;

  for (j in 1:p){
    theta_origin[j+1] = theta[j+1] / X_sd[j];
    // theta_fwi[j] = theta_origin[j+1] / X_minmax[j];
  }
  theta_origin[1] = theta[1] - dot_product(X_means, theta_origin[2:(p+1)]);
  for (j in 1:p){
    gridgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    gridgl[,j] = xholderLinear[,j] * theta_origin[j+1];
  };
  gridgsmooth = gridgl + gridgnl;
  {
  vector[n] pred = rep_vector(theta_origin[1], n);
    for(j in 1:p){
      pred += gridgsmooth[,j];      
    }
    gridalpha = exp(pred);
  }


  for (i in 1:n){
    log_lik[i] = pareto_lpdf(y[i] | u[i], alpha[i]);
  };
}
"


data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, Z_scales= Z_scales,
                  atau = ((psi+1)/2), xholderNonlinear = xholder.nonlinear, 
                  bsLinear = bs.linear[,-1], bsNonlinear = bs.nonlinear, #X_min=fwi.min,
                  xholderLinear = xholder.linear[,-1], #X_minmax = fwi.minmax, 
                  X_means = X_means, X_sd = X_sd)#,
                  #gridL = xgrid.linear[,-1], gridNL = xgrid.nonlinear)

init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p, psi)), 
                        theta = rep(-0.1, (p+1)), tau = rep(0.1, p), 
                        lambda1 = rep(0.1,p), lambda2 = rep(1, p)),
                   list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p, psi)),
                        theta = rep(0.05, (p+1)), tau = rep(0.2, p),
                        lambda1 = rep(0.2,p), lambda2 = rep(0.5, p)),
                   list(gamma_raw= array(rep(0.1, (psi*p)), dim=c(p, psi)),
                        theta = rep(0.1, (p+1)), tau = rep(0.1, p), 
                        lambda1 = rep(0.1,p), lambda2 = rep(0.1, p)))

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

# saveRDS(fit1, file=paste0("./BLAST/application/",Sys.Date(),"_stanfit.rds"))
# readRDS(file=paste0("./BLAST/application/2024-11-27","_stanfit.rds"))
posterior <- rstan::extract(fit1)
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))

theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
# fwismooth.samples <- summary(fit1, par=c("fwismooth"), probs = c(0.05, 0.5, 0.95))$summary
gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
gridgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
origin.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
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



g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))
g.linear.mean <- as.vector(matrix(gridgl.samples[,1], nrow = n, byrow=TRUE))
g.linear.q1 <- as.vector(matrix(gridgl.samples[,4], nrow = n, byrow=TRUE))
g.linear.q2 <- as.vector(matrix(gridgl.samples[,5], nrow = n, byrow=TRUE))
g.linear.q3 <- as.vector(matrix(gridgl.samples[,6], nrow = n, byrow=TRUE))
g.nonlinear.mean <- as.vector(matrix(gridgnl.samples[,1], nrow = n, byrow=TRUE))
g.nonlinear.q1 <- as.vector(matrix(gridgnl.samples[,4], nrow = n, byrow=TRUE))
g.nonlinear.q2 <- as.vector(matrix(gridgnl.samples[,5], nrow = n, byrow=TRUE))
g.nonlinear.q3 <- as.vector(matrix(gridgnl.samples[,6], nrow = n, byrow=TRUE))
# fwi.smooth.mean <- as.vector(matrix(fwismooth.samples[,1], nrow = n, byrow=TRUE))
# fwi.smooth.q1 <- as.vector(matrix(fwismooth.samples[,4], nrow = n, byrow=TRUE))
# fwi.smooth.q2 <- as.vector(matrix(fwismooth.samples[,5], nrow = n, byrow=TRUE))
# fwi.smooth.q3 <- as.vector(matrix(fwismooth.samples[,6], nrow = n, byrow=TRUE))
g.min.samples <- min(gsmooth.samples[,4])
g.max.samples <- max(gsmooth.samples[,6])
# fwi.smooth <- as.data.frame(matrix(fwismooth.samples[,4], nrow = n, byrow=TRUE))
# fwi.min.samples <- sapply(fwi.smooth, min)


xholder <- as.data.frame(xholder)
colnames(xholder) <- colnames(fwi.scaled)[1:p]
simul.data <- data.frame(BA = y-u, fwi.scaled[,c(1:p)])
gam.scale <- list(BA ~ s(BUI, bs = "ts", k = 30) + 
                      s(ISI, bs = "ts", k = 30) + 
                      s(FFMC, bs = "ts", k = 30) +
                      s(DMC, bs = "ts", k = 30) + 
                      s(DC, bs = "ts", k = 30),
                    ~ s(BUI, bs = "ts", k = 30) + 
                      s(ISI, bs = "ts", k = 30) + 
                      s(FFMC, bs = "ts", k = 30) +
                      s(DMC, bs = "ts", k = 30) +
                      s(DC, bs = "ts", k = 30))
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
                ~ s(BUI, bs = "ts", k = 30) + 
                  s(ISI, bs = "ts", k = 30) + 
                  s(FFMC, bs = "ts", k = 30) +
                  s(DMC, bs = "ts", k = 30) +
                  s(DC, bs = "ts", k = 30))
# evgam.fit.1 <- evgam::evgam(gam.1, data = simul.data, family = "gpd")
# pred.1 <- predict(evgam.fit.1, newdata = xholder, type="response", se.fit = TRUE)
# xi.pred.1 <-pred.1$fitted$shape
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
  ylab(expression(alpha(c,...,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_alpha.pdf"), width=10, height = 7.78)

# xi.scenario <- data.frame("x" = newx,
#                             "post.mean" = (1/alpha.samples[,1]),
#                             "post.median" = (1/alpha.samples[,5]),
#                             "q1" = (1/alpha.samples[,4]),
#                             "q3" = (1/alpha.samples[,6]),
#                             "evgam.1" = xi.pred.1,
#                             # "evgam.1.q1" = xi.low.1,
#                             # "evgam.1.q3" = xi.high.1,
#                             # "vgam.1" = vgam.xi.1,
#                             # "vgam.scale" = vgam.xi.scale,
#                             # "evgam.scale.q1" = xi.low.scale,
#                             # "evgam.scale.q3" = xi.high.scale,
#                             "evgam.scale" = xi.pred.scale)

# ggplot(xi.scenario, aes(x=x)) + 
#   ylab(expression(xi(c,...,c))) + xlab(expression(c)) + labs(col = "") +
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
#   geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
#   # geom_ribbon(aes(ymin = evgam.scale.q1, ymax = evgam.scale.q3), fill= "orange", alpha = 0.2) +
#   geom_line(aes(y=evgam.scale), colour = "orange", linewidth=1, linetype=3) +
#   # geom_ribbon(aes(ymin = evgam.1.q1, ymax = evgam.1.q3), fill= "purple", alpha = 0.2) +
#   geom_line(aes(y=evgam.1), colour = "purple", linewidth=1, linetype=3) +
#   # geom_line(aes(y=vgam.scale), colour = "orange", linewidth=1, linetype=4) +
#   # geom_line(aes(y=vgam.1), colour = "purple", linewidth=1, linetype=4) +  
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values = c("steelblue")) + 
#   guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(size = 20))



T <- 100
len    <- nrow(posterior$alpha)
posterior_idx <- sample(len, T, replace = TRUE)
alpha_sub <- posterior$alpha[posterior_idx, ] 
y.matrix <- matrix(y, nrow = T, ncol = n, byrow = TRUE)
u.matrix <- matrix(u, nrow = T, ncol = n, byrow = TRUE)
r_vec <- qnorm(pPareto(y.matrix, u.matrix, alpha = alpha_sub))
r_mat <- matrix(r_vec, nrow = T, ncol = n)
quantile_prob <- ppoints(n)
grid          <- qnorm(quantile_prob)
traj <- t(apply(r_mat, 1, sort))
bands <- matrixStats::colQuantiles(traj, probs = c(0.05, 0.5, 0.95))

qqplot_df <- data.frame(
  grid    = grid, 
  l.band  = bands[, 1], 
  trajhat = bands[, 2], 
  u.band  = bands[, 3]
)
ggplot(data = qqplot_df) + 
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

posterior_idx <- sample(len, T, replace = TRUE)
y.matrix <- matrix(y, nrow = T, ncol = n, byrow = TRUE)
u.matrix <- matrix(u, nrow = T, ncol = n, byrow = TRUE)
p_matrix <- 1 - (u.matrix/ y.matrix)^(posterior$alpha[posterior_idx,])
p_mean <- colMeans(p_matrix)
rp_integrated <- qnorm(p_mean)

rp <- data.frame(rp=as.numeric(rp_integrated), group = factor("residuals"))

ggplot(data = rp) +  
  geom_qqboxplot(aes(y=rp), notch=FALSE, varwidth=FALSE, reference_dist="norm", width = 0.15, qq.colour = "steelblue")+
  labs(x = "", y = "Residuals") + ylim(-4,4) + xlim(-.2,.2)+
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_qqboxplot.pdf"), width = 10, height = 7.78)
             
cat("Finished Running")

data.smooth <- data.frame("x" = as.vector(as.matrix(xholder)),
                          "true" = as.vector(as.matrix(fwi.scaled)),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
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
                  scale_color_manual(values=c("steelblue")) +
                  ylim(-10, 10) +
                  # theme_minimal(base_size = 30) +
                  # theme(legend.position = "none",
                  #         plot.margin = margin(0,0,0,-20),
                  #         axis.text = element_text(size = 35),
                  #         axis.title.x = element_text(size = 45))
                  theme_minimal(base_size = 20) +
                  theme(legend.position = "none",
                        plot.margin = margin(5, 5, 5, 5),
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        axis.text = element_text(size = 18),
                        axis.title.x = element_text(size = 22))                  
  grid.plts[[i]] <- grid.plt #+ annotate("point", x= fwi.scaled[which.max(y),i], y=g.min.samples, color = "red", size = 7)
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

fit.log.lik <- extract_log_lik(fit1)
elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 4)
elpd.loo
# save(elpd.loo, file = (paste0("./BLAST/application/BLAST_full_",Sys.Date(),"_",psi+2,"_",floor(threshold*100),"quantile_IC.Rdata")))
