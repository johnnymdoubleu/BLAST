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
threshold <- 0.975

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

fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)



# fwi.index <- fwi.index[which(Y>1),]
# load("./BLAST/application/quant-t.Rdata")
load("./BLAST/application/quant-t_10.Rdata")
# load("./BLAST/application/qgam_5_10.Rdata")
preds <- predict(quant.fit)
# u <- rep(quantile(Y, threshold),ceiling(nrow(fwi.index)*(1-threshold)))
# excess <- which(Y>u)
excess <- which(Y>preds)
u <- preds[excess]
# excess <- which(fwi.dd$excess==TRUE)
# u <- fwi.dd$origin_Model_Smooth_975[excess]
y <- Y[excess]

# fwi.scaled <- data.frame(fwi.scaled[which(Y>u),])
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# range01 <- function(x){(x)/max(x)}
fwi.scaled <- as.data.frame(sapply(fwi.scaled[excess,c(-1,-2)], FUN = range01))

n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]


fwi.origin <- data.frame(fwi.index[excess,c(-1,-2)], BA=y)
max.fwi <- fwi.origin[which.max(y),]
fwi.grid <- data.frame(lapply(fwi.origin[,c(1:p)], function(x) seq(min(x), max(x), length.out = nrow(fwi.scaled))))
fwi.minmax <- sapply(fwi.origin[,c(1:p)], function(x) max(x)-min(x))
fwi.min <- sapply(fwi.origin[,c(1:p)], function(x) min(x))
# ggplot(fwi.origin, aes(x=DSR, y=FFMC)) + 
#   geom_point(aes(colour = BA), size= 2.5) + 
#   scale_colour_stepsn(colours = c("slategray1", "red"), labels=function(x) format(x, big.mark = ",", scientific = TRUE), breaks=c(0.1e5, 0.5e5, 1e5, 2e5)) +
#   geom_density2d(colour="steelblue", linewidth = 1.3) + 
#   geom_mark_circle(aes(x = max.fwi$DSR, y = max.fwi$FFMC, label = "15th Oct 2017"), con.type = "straight",
#                    radius = unit(2.5, "mm"), color = "steelblue", size = 1, 
#                    con.colour = "steelblue", con.cap = unit(0, "mm"),
#                    label.colour = "steelblue", label.buffer = unit(5, "mm"),
#                    label.fill = "transparent")  +
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 30),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         strip.text = element_blank(),
#         axis.title = element_text(size = 30))
# ggsave("./BLAST/application/figures/extremeviz.pdf", width = 10, height = 7.78)

# ggplot(fwi.origin, aes(x=as.numeric(year), y=log(BA), color = BA)) + 
#   ylab("Hectares (log)") + xlab("Time (years)") + 
#   geom_point(size= 2.5, alpha = 0.5) + 
#   scale_colour_stepsn(colours = c("slategray1", "red"), labels=function(y) format(y, big.mark = ",", scientific = TRUE), 
#   breaks = quantile(fwi.origin$BA, probs = seq(0,1,length.out = 20))) + 
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 30),
#         legend.position = "none",
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         strip.text = element_blank(),
#         axis.title = element_text(size = 30))

# ggsave("./BLAST/application/figures/hectareslog.pdf", width = 10, height = 7.78)

# ggplot(fwi.origin, aes(x=as.numeric(year))) + 
#   ylab("Density") + xlab("Time (years)") + 
#   geom_histogram(aes(y = after_stat(density)), fill = "steelblue", color = "gray", alpha = .2) +
#   geom_rug() +
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 30),
#         legend.position = "none",
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         strip.text = element_blank(),
#         axis.title = element_text(size = 30))
# M <- cor(fwi.origin[,c(1:7)])
# corrplot(M, order = 'AOE', type = 'upper', tl.pos = 'tp')
# corrplot(M, add = TRUE, type = 'lower', method = 'number', order = 'AOE',
        #  col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n')
# ggsave("./BLAST/application/figures/intensityfn.pdf", width = 10, height = 7.78)
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
  sm_spec <- smoothCon(s(x_vec, bs = "tp", k = psi + 2), 
                      data = data.frame(x_vec = x_vec), 
                      knots = NULL)[[1]]
  
  X_raw <- sm_spec$X
  S     <- sm_spec$S[[1]] 
  
  eig <- eigen(S, symmetric = TRUE)
  max_lambda <- max(eig$values)
  tol <- max_lambda * 1e-7  # Relative threshold (Robust)
  
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
  # if(ncol(Z_final) < psi){
  #   n.pad <- psi - ncol(Z_final)
  #   zero.pad <- matrix(0, nrow = nrow(Z_final), ncol = n.pad)
  #   Z_final <- cbind(Z_final, zero.pad)    
  # }  
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

grid_Z_list <- list()

for (i in seq_along(covariates)) {
  scaled_x_grid <- (fwi.grid[,i] - fwi.min[i]) / (fwi.minmax[i])
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
# fwi.grid <- as.data.frame(sapply(fwi.grid, FUN = range01))
xgrid.linear <- model.matrix(~ ., data = data.frame(fwi.grid))
xgrid.linear[,-1] <- sweep(xgrid.linear[, -1], 2, fwi.min, "-")
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
  matrix[n, p] gridL; // fwi dataset
  matrix[n, (psi*p)] gridNL; // thin plate splines basis      
  array[n] real <lower=u> y; // extreme response
  real <lower=0> atau;
  vector[(psi*p)] Z_scales;
  vector[p] X_minmax;
  vector[p] X_min;
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
  array[n] real <lower=0> gridalpha; // new tail index
  matrix[n, p] gridgnl; // nonlinear component
  matrix[n, p] gridgl; // linear component
  matrix[n, p] gridgsmooth; // linear component 
  matrix[n, p] fwismooth;
  vector[n] log_lik;
  vector[p+1] theta_origin;
  vector[p] theta_fwi;

  for (j in 1:p){
    theta_origin[j+1] = theta[j+1] / X_sd[j];
    theta_fwi[j] = theta_origin[j+1] / X_minmax[j];
  }
  theta_origin[1] = theta[1] - dot_product(X_means, theta_origin[2:(p+1)]);
  for (j in 1:p){
    gridgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    gridgl[,j] = xholderLinear[,j] * theta_origin[j+1];
    gridgsmooth[,j] = gridgl[,j] + gridgnl[,j];
    fwismooth[,j] = gridNL[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j] + gridL[,j] * theta_fwi[j];
  };

  for (i in 1:n){
    gridalpha[i] = exp(theta_origin[1] + sum(gridgsmooth[i,]));
    log_lik[i] = pareto_lpdf(y[i] | u[i], alpha[i]);
  };
}
"

  # for (j in 1:p){
  #   theta_origin[j+1] = theta[j+1] / (X_max[j]-X_min[j]);
  # }
  # theta_origin[1] = theta[1] - dot_product(X_means, theta_origin[2:(p+1)]);


# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,2)]) %*% bs.nonlinear[,c((1):(psi))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,3)]) %*% bs.nonlinear[,c((psi+1):(psi*2))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,4)]) %*% bs.nonlinear[,c((psi*2+1):(psi*3))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,5)]) %*% bs.nonlinear[,c((psi*3+1):(psi*4))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,6)]) %*% bs.nonlinear[,c((psi*4+1):(psi*5))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,7)]) %*% bs.nonlinear[,c((psi*5+1):(psi*6))]), "\n")
# cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,8)]) %*% bs.nonlinear[,c((psi*6+1):(psi*7))]), "\n")

data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, Z_scales= Z_scales,
                  atau = ((psi+1)/2), xholderNonlinear = xholder.nonlinear, 
                  bsLinear = bs.linear[,-1], bsNonlinear = bs.nonlinear, X_min=fwi.min,
                  xholderLinear = xholder.linear[,-1], X_minmax = fwi.minmax, 
                  X_means = X_means, X_sd = X_sd,
                  gridL = xgrid.linear[,-1], gridNL = xgrid.nonlinear)

init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p, psi)), #rho = 1, 
                        theta = rep(-0.1, (p+1)), tau = rep(0.1, p), #sigma_ridge = rep(0.1, p),
                        lambda1 = rep(0.1,p), lambda2 = rep(1, p)),
                   list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p, psi)),# rho = 1,
                        theta = rep(0.05, (p+1)), tau = rep(0.2, p), #sigma_ridge = rep(0.2, p), 
                        lambda1 = rep(2,p), lambda2 = rep(5, p)),
                   list(gamma_raw= array(rep(0.1, (psi*p)), dim=c(p, psi)), #rho = 1,
                        theta = rep(0.1, (p+1)), tau = rep(0.1, p), #sigma_ridge = rep(0.1, p),
                        lambda1 = rep(0.1,p), lambda2 = rep(0.1, p)))

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
# bayesplot::color_scheme_set("mix-blue-red")
# bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("") +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(size = 18))

theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
fwismooth.samples <- summary(fit1, par=c("fwismooth"), probs = c(0.05, 0.5, 0.95))$summary
gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
gridgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
origin.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary
# yrep <- summary(fit1, par=c("yrep"), probs = c(0.05,0.5, 0.95))$summary
# f.samples <- summary(fit1, par=c("f"), probs = c(0.05,0.5, 0.95))$summary
# loglik.samples <- summary(fit1, par=c("log_lik"), probs = c(0.05,0.5, 0.95))$summary

# MCMCvis::MCMCplot(fit1, params = 'theta')
# MCMCvis::MCMCplot(fit1, params = "gamma")


post.samples <- as.matrix(fit1)
theta_samples <- post.samples[, grepl("^theta\\[", colnames(post.samples))]
gamma_samples <- post.samples[, grepl("^gamma\\[", colnames(post.samples))]
num_draws <- nrow(post.samples)
p <- length(X_sd)
psi <- ncol(gamma_samples) / p
theta_fwi_samples <- matrix(NA, nrow = num_draws, ncol = p)
fwismooth_samples <- array(NA, dim = c(num_draws, n, p))

for (s in 1:num_draws) {
  curr_theta <- theta_samples[s, ]
  theta_orig_slopes <- curr_theta[2:(p+1)] / X_sd
  theta_fwi_s <- theta_orig_slopes / fwi.minmax
  for (j in 1:p) {
    start_idx <- ((j - 1) * psi) + 1
    end_idx   <- j * psi
    curr_gamma_j <- as.matrix(gamma_samples[s, start_idx:end_idx])
    
    nl_part <- xgrid.nonlinear[, start_idx:end_idx] %*% curr_gamma_j
    l_part  <- (xgrid.linear[, j+1]-fwi.min[j]) * theta_fwi_s[j]
    
    fwismooth_samples[s, , j] <- as.vector(nl_part + l_part)
  }
}



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
xholder <- as.data.frame(xholder)
colnames(xholder) <- colnames(fwi.scaled)[1:p]
simul.data <- data.frame(BA = y-u, bs.linear[,-1])#fwi.origin[c(1:p)])
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
gam.scale <- list(BA ~ s(BUI, bs = "tp", k = 30) + 
                      s(ISI, bs = "tp", k = 30) + 
                      s(FFMC, bs = "tp", k = 30) +
                      s(DMC, bs = "tp", k = 30) + 
                      s(DC, bs = "tp", k = 30),
                    ~ s(BUI, bs = "tp", k = 30) + 
                      s(ISI, bs = "tp", k = 30) + 
                      s(FFMC, bs = "tp", k = 30) +
                      s(DMC, bs = "tp", k = 30) +
                      s(DC, bs = "tp", k = 30))
evgam.fit.scale <- evgam::evgam(gam.scale, data = simul.data, family = "gpd")
pred.scale <- predict(evgam.fit.scale, newdata = xholder, type="response", se.fit = TRUE)
xi.pred.scale <-pred.scale$fitted$shape
xi.se.scale <- pred.scale$se.fit$shape
xi.low.scale <- xi.pred.scale - (1.96 * xi.se.scale)
xi.high.scale <- xi.pred.scale + (1.96 * xi.se.scale)
alpha.pred.scale <- 1/xi.pred.scale

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

# gam.1 <- list(BA ~ 1,
#                 ~ s(DSR, bs = "tp", k = 30) + 
#                     s(FWI, bs = "tp", k = 30) + 
#                     s(BUI, bs = "tp", k = 30) + 
#                     s(ISI, bs = "tp", k = 30) + 
#                     s(FFMC, bs = "tp", k = 30) +
#                     s(DMC, bs = "tp", k = 30) + 
#                     s(DC, bs = "tp", k = 30)) 

gam.1 <- list(BA ~ 1,
                ~ s(BUI, bs = "tp", k = 30) + 
                  s(ISI, bs = "tp", k = 30) + 
                  s(FFMC, bs = "tp", k = 30) +
                  s(DMC, bs = "tp", k = 30) +
                  s(DC, bs = "tp", k = 30))
evgam.fit.1 <- evgam::evgam(gam.1, data = simul.data, family = "gpd")
pred.1 <- predict(evgam.fit.1, newdata = xholder, type="response", se.fit = TRUE)
xi.pred.1 <-pred.1$fitted$shape
xi.se.1 <- pred.1$se.fit$shape
xi.low.1 <- xi.pred.1 - (1.96 * xi.se.1)
xi.high.1 <- xi.pred.1 + (1.96 * xi.se.1)
alpha.pred.1 <- 1/xi.pred.1

# xholder.basis.1 <- predict(evgam.fit.1, newdata = xholder, type= "lpmatrix")$shape
# xi.coef.1 <- tail(evgam.fit.1$coefficients, (psi-1)*p)
# gamma.xi.1 <- matrix(xi.coef.1, ncol = p)
# alpha.nonlinear.1 <- xi.nonlinear.1 <- matrix(, nrow = n, ncol = p)
# bs.nonlinear.1 <- xholder.basis.1[,c(2:((psi-1)*p+1))]
# for(j in 1:p){
#   xi.nonlinear.1[,j] <- bs.nonlinear.1[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.1[,j]
#   alpha.nonlinear.1[,j] <- 1/(xi.nonlinear.1[,j])
# }

data.scenario <- data.frame("x" = newx,
                            "post.mean" = (alpha.samples[,1]),
                            "post.median" = (alpha.samples[,5]),
                            "q1" = (alpha.samples[,4]),
                            "q3" = (alpha.samples[,6]),
                            "post.mean.org" = (origin.samples[,1]),
                            "post.median.org" = (origin.samples[,5]),
                            "q1.org" = (origin.samples[,4]),
                            "q3.org" = (origin.samples[,6]))

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

xi.scenario <- data.frame("x" = newx,
                            "post.mean" = (1/alpha.samples[,1]),
                            "post.median" = (1/alpha.samples[,5]),
                            "q1" = (1/alpha.samples[,4]),
                            "q3" = (1/alpha.samples[,6]),
                            "evgam.1" = xi.pred.1,
                            "evgam.1.q1" = xi.low.1,
                            "evgam.1.q3" = xi.high.1,
                            "evgam.scale.q1" = xi.low.scale,
                            "evgam.scale.q3" = xi.high.scale,
                            "evgam.scale" = xi.pred.scale)

ggplot(xi.scenario, aes(x=x)) + 
  ylab(expression(xi(c,...,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  geom_ribbon(aes(ymin = evgam.scale.q1, ymax = evgam.scale.q3), fill= "orange", alpha = 0.2) +
  geom_line(aes(y=evgam.scale), colour = "orange", linewidth=1, linetype=2) +
  geom_ribbon(aes(ymin = evgam.1.q1, ymax = evgam.1.q3), fill= "purple", alpha = 0.2) +
  geom_line(aes(y=evgam.1), colour = "purple", linewidth=1, linetype=3) +  
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
  rp[i] <- qnorm(pPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i]))
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
                  ylim(g.min.samples, g.max.samples) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=g.min.samples, color = "red", size = 7)
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_DC.pdf"), grid.plts[[7]], width=10, height = 7.78)


data.linear <- data.frame("x" = as.vector(as.matrix(xholder)),
                          "true" = as.vector(as.matrix(fwi.scaled)),
                          "post.mean" = as.vector(g.linear.mean),
                          "q1" = as.vector(g.linear.q1),
                          "q2" = as.vector(g.linear.q2),
                          "q3" = as.vector(g.linear.q3),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(data.linear[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  ylim(g.min.samples, g.max.samples) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=g.min.samples, color = "red", size = 7)
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)


data.nonlinear <- data.frame("x" = as.vector(as.matrix(xholder)),
                          "true" = as.vector(as.matrix(fwi.scaled)),
                          "post.mean" = as.vector(g.nonlinear.mean),
                          "q1" = as.vector(g.nonlinear.q1),
                          "q2" = as.vector(g.nonlinear.q2),
                          "q3" = as.vector(g.nonlinear.q3),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(data.nonlinear[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  ylim(g.min.samples, g.max.samples) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=g.min.samples, color = "red", size = 7)
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

data.smooth <- data.frame("x" = as.vector(as.matrix(fwi.grid)),
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
                  scale_color_manual(values=c("steelblue")) +
                  # ylim(fwi.min.samples[i], 7) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.grid[which.max(y),i], y=fwi.min.samples[i], color = "red", size = 7)
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

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
loo(fit.log.lik, is_method = "sis", cores = 4)

