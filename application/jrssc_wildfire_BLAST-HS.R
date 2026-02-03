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

model.stan <- "// Stan model for BHST Pareto Samples
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
  array[n] real <lower=1> y; // extreme response
  vector[(psi*p)] Z_scales;
  vector[p] X_minmax;
  vector[p] X_means;
  vector[p] X_sd;
}

parameters {
  vector[(p+1)] theta; // linear predictor
  array[p] vector[psi] gamma_raw;
  array[p] real <lower=0> lambda1; // lasso penalty //array[p] real <lower=0> 
  array[p] real <lower=0> lambda2; // lambda2 group lasso penalty
  real <lower=0> t1;
  real <lower=0> b1;
  real <lower=0> t2;
  real <lower=0> b2;
  array[p] real <lower=0> lt1;
  array[p] real <lower=0> lt2;
}

transformed parameters {
  array[n] real <lower=0> alpha; // covariate-adjusted tail index
  
  array[p] vector[psi] gamma;
  {
    matrix[n, p] gsmooth; // linear component
    for (j in 1:p){
      for (k in 1:psi){
        int idx = (j-1)*psi + k;
        gamma[j][k] = gamma_raw[j][k] * sqrt(lambda2[j]) * sqrt(t2) * Z_scales[idx];
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
  target += inv_gamma_lpdf(t1 | 0.5, 1/b1); 
  target += inv_gamma_lpdf(b1 | 0.5, 1); 
  target += inv_gamma_lpdf(t2 | 0.5, 1/b2);
  target += inv_gamma_lpdf(b2 | 0.5, 1);
  for (j in 1:p){ 
    target += inv_gamma_lpdf(lambda1[j] | 0.5, 1/lt1[j]);
    target += cauchy_lpdf(lt1[j] | 0, 1);
    target += normal_lpdf(theta[(j+1)] | 0, t1 * lambda1[j]);
    target += inv_gamma_lpdf(lambda2[j] | 0.5, 1/lt2[j]);
    target += cauchy_lpdf(lt2[j] | 0, 1);
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


data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, Z_scales= Z_scales, 
                  xholderNonlinear = xholder.nonlinear, 
                  bsLinear = bs.linear[,-1], bsNonlinear = bs.nonlinear, 
                  xholderLinear = xholder.linear[,-1], X_minmax = fwi.minmax, 
                  X_means = X_means, X_sd = X_sd,
                  gridL = xgrid.linear[,-1], gridNL = xgrid.nonlinear)

init.alpha <- list(list(gamma_raw = array(rep(1, (psi*p)), dim=c(p,psi)),
                        theta = rep(0, (p+1)), 
                        lambda1 = rep(0.1, p), lambda2 = rep(1, p), 
                        t1 = 0.01, t2 = 0.1, b1 = 0.1, b2 = 0.1,
                        lt1 = rep(0.1, p), lt2 = rep(0.1, p)),
                    list(gamma_raw = array(rep(2, (psi*p)), dim=c(p,psi)),
                        theta = rep(0, (p+1)), 
                        lambda1 = rep(0.5, p), lambda2 = rep(5, p), 
                        t1 = 1, t2 = 1, b1 = 0.1, b2 = 0.01,
                        lt1 = rep(0.1, p), lt2 = rep(0.1, p)),
                    list(gamma_raw = array(rep(-0.5, (psi*p)), dim=c(p,psi))),
                        theta = rep(0.2, (p+1)),
                        lambda1 = rep(1, p), lambda2 = rep(2, p), 
                        t1 = 1, t2 = 0.01, b1 = 0.01, b2 = 0.01,
                        lt1 = rep(1, p), lt2 = rep(0.01, p))

fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value 
    chains = 3,             # number of Markov chains
    iter = 4000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 1000           # no progress shown
)

posterior <- rstan::extract(fit1)
# bayesplot::color_scheme_set("mix-blue-red")
# bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("") +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(size = 18))

theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2", "t1", "t2", "lt1", "lt2", "b1", "b2"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
fwismooth.samples <- summary(fit1, par=c("fwismooth"), probs = c(0.05, 0.5, 0.95))$summary
gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
gridgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary
# yrep <- summary(fit1, par=c("yrep"), probs = c(0.05,0.5, 0.95))$summary
# f.samples <- summary(fit1, par=c("f"), probs = c(0.05,0.5, 0.95))$summary
# loglik.samples <- summary(fit1, par=c("log_lik"), probs = c(0.05,0.5, 0.95))$summary
MCMCvis::MCMCplot(fit1, params = 'theta')
MCMCvis::MCMCplot(fit1, params = "gamma")

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

g.min.samples <- min(gsmooth.samples[,4])
g.max.samples <- max(gsmooth.samples[,6])
fwi.smooth <- as.data.frame(matrix(fwismooth.samples[,4], nrow = n, byrow=TRUE))
fwi.min.samples <- sapply(fwi.smooth, min)



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


data.linear <- data.frame("x" = as.vector(xholder),
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


data.nonlinear <- data.frame("x" = as.vector(xholder),
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

fit.log.lik <- extract_log_lik(fit1)

elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 4)
save(elpd.loo, file = (paste0("./BLAST/application/BHST_full_",Sys.Date(),"_",psi+2,"_",floor(threshold*100),"quantile_IC.Rdata")))
