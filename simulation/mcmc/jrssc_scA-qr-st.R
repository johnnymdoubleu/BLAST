library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(qqboxplot)
library(evgam)

n <- 10000
psi <- 10
threshold <- 0.95
p <- 5
p_lin <- 7 # 5 covariates + 2 harmonics
psi <- psi - 2

time.seq <- 1:n
period <- 365 
x.season <- time.seq %% period / period 

seasons <- c("Winter", "Spring", "Summer", "Autumn")
x.origin <- matrix(0, nrow = n, ncol = p)

for (j in 1:p) {
  phase_shift <- j * (2 * pi / p) 
  seasonal_trend <- 0.5 + 0.1 * sin(2 * pi * time.seq / period + phase_shift) 
  uniform_noise <- runif(n, min = -0.4, max = 0.4)
  x.origin[, j] <- seasonal_trend + uniform_noise
  # x.origin[,j] <- runif(n, 0,1)
}

colnames(x.origin) <- paste0("X", 1:p)

# --- DGP MODIFICATION: ADD HARMONICS DIRECTLY TO TRUE ALPHA ---
sin.time.full <- sin(2 * pi * x.season)
cos.time.full <- cos(2 * pi * x.season)

f2 <- function(x) {-.7 * sin(2 * pi * x^2)*(x-0.5)}
f5 <- function(x) {-.7 * cos(3 * pi * x^2)*x}

theta.origin <- c(0.7, 0, 0.8, 0, 0, -0.8) 
theta.harm <- c(-0.5, 0.3)                 # True effects for sin and cos

alp.origin <- exp(rep(theta.origin[1],n) + 
                  x.origin %*% theta.origin[-1] + 
                  f2(x.origin[,2]) + f5(x.origin[,5]) + 
                  (theta.harm[1] * sin.time.full) + 
                  (theta.harm[2] * cos.time.full))

# Generate true Pareto responses
y.origin <- rPareto(n, rep(1, n), alpha = alp.origin)

# Dynamic Thresholding using ALD
evgam.df <- data.frame(
  y = (y.origin),
  sin.time = sin.time.full,
  cos.time = cos.time.full,
  x.origin
)
evgam.cov <- y ~ 1 + cos.time + sin.time 
ald.cov.fit <- evgam(evgam.cov, data = evgam.df, family = "ald", ald.args=list(tau = threshold))
u.vec <- (predict(ald.cov.fit)$location)
# u.vec <- 20^(1/alp.origin)
# Subset excesses
excess.index <- which(y.origin > u.vec)
x.origin <- x.origin[excess.index,]
sin.time <- sin.time.full[excess.index]
cos.time <- cos.time.full[excess.index]
y.origin <- y.origin[excess.index]
u <- u.vec[excess.index]
n_excess <- length(y.origin)

# Combine all 7 linear variables
all_linear_vars <- data.frame(x.origin, sin.time = sin.time, cos.time = cos.time)

# Scale linear covariates to [0,1]
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
X_minmax <- sapply(all_linear_vars, function(x) max(x)-min(x))
X_min <- sapply(all_linear_vars, function(x) min(x))
all_scaled <- as.data.frame(sapply(all_linear_vars, FUN = range01))
x_scaled_only <- all_scaled[, 1:p]

group.map <- c()
Z.list <- list()        
scale_stats_list <- list() 
projection_coefs_list <- list() 
spec_decomp_list <- list() 
sm_spec_list <- list()     
keep_cols_list <- list()

covariates <- colnames(x.origin)

# QR Spline decomposition applied ONLY to the p=5 covariates
for (i in seq_along(covariates)) {
  x_vec <- x_scaled_only[, i]
  X_lin <- model.matrix(~ x_vec) 
  sm_spec <- smoothCon(s(x_vec, bs = "tp", k = psi + 2), 
                       data = data.frame(x_vec = x_vec), 
                       knots = NULL)[[1]]
  
  X_raw <- sm_spec$X
  S     <- sm_spec$S[[1]] 
  
  eig <- eigen(S, symmetric = TRUE)
  tol <- max(eig$values) * 1e-6 
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
  Z.list[[i]] <- Z_final
  
  spec_decomp_list[[i]] <- list(U_pen = U_pen, Lambda_sqrt_inv = solve(sqrt(Lambda_pen)))
  projection_coefs_list[[i]] <- Gamma_Original 
  keep_cols_list[[i]] <- keep_cols
  scale_stats_list[[i]] <- train_scale         
  sm_spec_list[[i]] <- sm_spec
}

bs.nonlinear <- do.call(cbind, Z.list)
bs.linear <- model.matrix(~ ., data = all_scaled)[,-1]

grid.n <- 200
newx <- seq(0, 1, length.out = grid.n)
xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))
colnames(xholder) <- covariates

# Set Peak Season for Partial Plots (e.g., day 200)
peak_day <- 200
raw_cos_peak <- cos(2 * pi * peak_day / period)
raw_sin_peak <- sin(2 * pi * peak_day / period)
# scaled_cos_peak <- (raw_cos_peak - X_min["cos.time"]) / X_minmax["cos.time"]
# scaled_sin_peak <- (raw_sin_peak - X_min["sin.time"]) / X_minmax["sin.time"]

xholder_harm <- matrix(c(rep(raw_sin_peak, grid.n), rep(raw_cos_peak, grid.n)), ncol=2)
colnames(xholder_harm) <- c("sin.time", "cos.time")
xholder_all <- cbind(xholder, xholder_harm)

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

Z_scales <- unlist(scale_stats_list)
xholder.nonlinear <- do.call(cbind, grid_Z_list)
xholder.linear <- model.matrix(~ ., data = data.frame(xholder_all))[,-1]

# Convert scaled grid back to original scale to evaluate True Alpha Function
# xholder_orig_scale <- sweep(sweep(xholder_fwi, 2, X_minmax[1:p], "*"), 2, X_min[1:p], "+")
alp.new <- as.vector(exp(theta.origin[1] + xholder %*% theta.origin[-1] + 
                          f2(xholder[,2]) + f5(xholder[,5]) +
                          (theta.harm[1] * raw_sin_peak) + (theta.harm[2] * raw_cos_peak)))
g.new <- c(xholder[,1] * theta.origin[2], 
           xholder[,2] * theta.origin[3] + f2(xholder[,2]),
           xholder[,3] * theta.origin[4],
           xholder[,4] * theta.origin[5],
           xholder[,5] * theta.origin[6] + f5(xholder[,5]))


X_means <- colMeans(bs.linear)
X_sd   <- apply(bs.linear, 2, sd)
bs.linear <- scale(bs.linear, center = X_means, scale = X_sd)

# --- STAN MODEL (Matches the applied version) ---
model.stan <- "// Stan model for BLAST Pareto Samples
data {
  int <lower=1> n; 
  int <lower=1> grid_n;
  int <lower=1> p; 
  int <lower=1> p_lin;
  int <lower=1> psi; 
  vector[n] u; 
  matrix[n, p_lin] bsLinear; 
  matrix[n, (psi*p)] bsNonlinear; 
  matrix[grid_n, p_lin] xholderLinear; 
  matrix[grid_n, (psi*p)] xholderNonlinear;    
  vector<lower=0>[n] y; 
  real <lower=0> atau;
  vector[p_lin] X_means;
  vector[p_lin] X_sd;
}

parameters {
  vector[(p_lin+1)] theta; 
  array[p] vector[psi] gamma_raw;
  array[p_lin] real <lower=0> lambda1; 
  array[p] real <lower=0> lambda2; 
  array[p] real <lower=0> tau;
}

transformed parameters {
  vector[n] alpha; 
  array[p] vector[psi] gamma;
  {
    vector[n] eta = rep_vector(theta[1], n);
    for (k in 1:p_lin) {
       eta += col(bsLinear, k) * theta[k+1];
    }
    for (j in 1:p){
      gamma[j] = gamma_raw[j] * sqrt(tau[j]); 
      int nl_start = (j - 1) * psi + 1;
      eta += block(bsNonlinear,1, nl_start, n, psi) * gamma[j];
    };
    alpha = exp(eta);
  }
}

model {
  target += pareto_lpdf(y | u, alpha);
  target += normal_lpdf(theta[1] | 0, 10); // Tightened prior
  target += gamma_lpdf(lambda1 | 1e-1, 1e-1); 
  target += gamma_lpdf(lambda2 | 1e-2, 1e-2);
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
  matrix[grid_n, p] gridgnl; 
  matrix[grid_n, p] gridgl; 
  matrix[grid_n, p] gridgsmooth; 

  vector[p_lin] theta_origin = theta[2:(p_lin+1)] ./ X_sd;
  real theta0 = theta[1] - dot_product(X_means, theta_origin);

  {
    vector[grid_n] grideta = rep_vector(theta0, grid_n);
    // Add harmonic baseline to the alpha prediction
    for (h in (p+1):p_lin){
       grideta += col(xholderLinear, h) * theta_origin[h];
    }
    
    for (j in 1:p){
        gridgnl[,j] = block(xholderNonlinear,1, ((j - 1) * psi + 1), grid_n, psi) * gamma[j];
        gridgl[,j] = col(xholderLinear, j) * theta_origin[j];
        gridgsmooth[,j] = gridgl[,j] + gridgnl[,j];
        grideta += gridgsmooth[,j];
    };
    gridalpha = exp(grideta);
  }
}
"

data.stan <- list(y = as.vector(y.origin), u = u, p = p, p_lin = p_lin, n= n_excess, psi = psi, grid_n = grid.n,
                  atau = ((psi+1)/2), X_means = X_means, X_sd=X_sd, 
                  bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                  xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear)

init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p,psi)),
                        theta = rep(-0.1, (p_lin+1)), tau = rep(0.1, p),
                        lambda1 = rep(0.1, p_lin), lambda2 = rep(1, p)),
                   list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p,psi)),
                        theta = rep(0.05, (p_lin+1)), tau = rep(2, p),
                        lambda1 = rep(2, p_lin), lambda2 = rep(5, p)),
                   list(gamma_raw = array(rep(-0.75, (psi*p)), dim=c(p,psi)),
                        theta = rep(0.01, (p_lin+1)), tau = rep(1.5, p),
                        lambda1 = rep(0.5, p_lin), lambda2= rep(5, p)))

system.time(fit1 <- stan(
  model_code = model.stan, 
  data = data.stan,    
  init = init.alpha,      
  chains = 3,             
  iter = 2000,            
  cores = parallel::detectCores(), 
  refresh = 1000             
))

posterior <- extract(fit1)

posterior <- extract(fit1)
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("Scenario A") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/appendix_",n,"_traceplot_scA.pdf"), width=22, height = 3)
# theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
theta.samples <- summary(fit1, par=c("theta0","theta_origin"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
# newgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
# newgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
# season.samples <- summary(fit1, par=c("alphaseason"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary

# gamma.post.mean <- gamma.samples[,1]
# gamma.q1 <- gamma.samples[,4]
# gamma.q2 <- gamma.samples[,1]
# gamma.q3 <- gamma.samples[,6]
# theta.post.mean <- theta.samples[,1]
# theta.q1 <- theta.samples[,4]
# theta.q2 <- theta.samples[,5]
# theta.q3 <- theta.samples[,6]

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
g.smooth.mean <- as.vector(matrix(newgsmooth.samples[,1], nrow = grid.n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(newgsmooth.samples[,4], nrow = grid.n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(newgsmooth.samples[,5], nrow = grid.n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(newgsmooth.samples[,6], nrow = grid.n, byrow=TRUE))

newalpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary
data.scenario <- data.frame("x" = newx,
                            "true" = (alp.new),
                            "q2" = (newalpha.samples[,5]),
                            "q1" = (newalpha.samples[,4]),
                            "q3" = (newalpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,...,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y = true, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2, linetype=2) + #ylim(0, 10) +
  geom_line(aes(y=q2, col = "Posterior Median"), linewidth=1.5) +
  scale_color_manual(values=c("steelblue", "red")) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  theme_minimal(base_size = 40) + 
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 30))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_scA.pdf"), width=10, height = 7.78)

# est_gridalpha_median  <- apply(posterior$gridalpha, 2, median)
# est_gridalpha_lower <- apply(posterior$gridalpha, 2, quantile, probs = 0.05)
# est_gridalpha_upper <- apply(posterior$gridalpha, 2, quantile, probs = 0.95)

# data.scenario <- data.frame("x" = newx,
#                             "true" = true_gridalpha,
#                             "q2" = est_gridalpha_median,
#                             "q1" = est_gridalpha_lower,
#                             "q3" = est_gridalpha_upper)
#                             # "post.median" = (alpha.samples[,5]),
#                             # "q1" = (alpha.samples[,4]),
#                             # "q3" = (alpha.samples[,6]))

# ggplot(data.scenario, aes(x=x)) + 
#   ylab(expression(alpha(c,...,c))) + xlab(expression(c)) +
#   geom_ribbon(aes(ymin = q1, ymax = q3), fill="steelblue", alpha = 0.2) +
#   geom_line(aes(y=true), color = "red", linewidth=1, linetype =2) +
#   geom_line(aes(y=q2), color = "steelblue", linewidth=1) +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(size = 20))


data.smooth <- data.frame("x"=newx,
                          "true" = g.new,
                          "q2"        = g.smooth.q2,
                          "q1"        = g.smooth.q1,
                          "q3"        = g.smooth.q3,
                          "covariates" = gl(p, grid.n, (p*grid.n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]")),
                          "replicate" = gl(p, grid.n, (p*grid.n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

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

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_smooth_scA.pdf"), width=12.5, height = 15)

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



mcmc.alpha <- posterior$alpha
len <- dim(mcmc.alpha)[1]
r <- matrix(, nrow = n_excess, ncol = 30)
T <- 30
for(i in 1:n_excess){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y.origin[i], u[i], alpha = mcmc.alpha[round(runif(1,1,len)),i]))
  }
}

lgrid <- n_excess
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
