library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(qqboxplot)
library(evgam)
library(forecast)

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
season_code_full <- cut((time.seq %% period / period), breaks = c(-0.1, 0.25, 0.5, 0.75, 1.1), labels = c(1,2,3,4))
seasons <- c("Winter", "Spring", "Summer", "Autumn")

x.origin <- matrix(0, nrow = n, ncol = p)

# 1. Generate Raw Covariates with Seasonality
for (j in 1:p) {
  phase_shift <- j * (2 * pi / p) 
  # seasonal_trend <- 0.5 + 0.1 * sin(2 * pi * time.seq / period + phase_shift) 
  seasonal_trend <- 0.5 + 0.1 * sin(2 * pi * time.seq / period) 
  uniform_noise <- runif(n, min = -0.49, max = 0.49)
  x.origin[, j] <- seasonal_trend + uniform_noise
  # norm_noise <- rnorm(n)
  # x.origin[, j] <- seasonal_trend + norm_noise  
  # ar.noise <- as.numeric(arima.sim(n = n, model = list(ar = 0.6), sd = 0.1))
  # x.origin[, j] <- seasonal_trend + ar.noise
}

fit.list <- list()
for (j in 1:p) {
  y_ts <- ts(x.origin[, j], frequency = period) 
  decomp <- stl(y_ts, s.window = "periodic", robust = TRUE)
  seasonal <- as.numeric(decomp$time.series[, "seasonal"])
  x.origin[, j] <- x.origin[, j] - seasonal
  # fit.list[[j]] <- forecast::auto.arima(x.origin[,j] - seasonal, seasonal = FALSE, stepwise = TRUE)
  # x.origin[, j] <- as.numeric(residuals(fit.list[[j]]))
}

# acf(x.origin[,5])
# plot(x.origin[,5])

covariates <- colnames(data.frame(x.origin))[1:p]

f2 <- function(x) {-.7 * sin(2 * pi * x^2)*(x-0.5)}
f5 <- function(x) {-.7 * cos(3 * pi * x^2)*x}
theta.origin <- c(0.7, 0, 0.8, 0, 0, -0.8) 

alp.origin <- exp(rep(theta.origin[1],n) + x.origin %*% theta.origin[-1] + f2(x.origin[,2]) + f5(x.origin[,5]))
y.noise <- rPareto(n, rep(1, n), alpha = alp.origin)

f.season.scale <- function(t){
  return(2.5 - .8 * sin(2 * pi * t / period) - .6 * cos(2 * pi * t/period)) 
}

y.origin <- y.noise * f.season.scale(time.seq)
plot(y.origin)

evgam.df <- data.frame(
  y = (y.origin),
  sin.time = sin(2 * pi * time.seq / period),
  cos.time = cos(2 * pi * time.seq / period),
  x.origin
)

evgam.cov <- y ~ 1 + cos.time + sin.time 
ald.cov.fit <- evgam(evgam.cov, data = evgam.df, family = "ald", ald.args=list(tau = threshold))
u.vec <- (predict(ald.cov.fit)$location)
# u.vec <- exp(log(f.season.scale(time.seq)) + log(20) / alp.origin)
excess.index <- which(y.origin > u.vec)
x.origin <- as.data.frame(x.origin[excess.index,])
y.origin <- y.origin[excess.index]
u <- u.vec[excess.index]
season_code_full <- season_code_full[excess.index]  
n <- length(y.origin)

colnames(x.origin) <- paste0("X", 1:p)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# X_minmax <- sapply(x.origin, function(x) max(x)-min(x))
# X_min <- sapply(x.origin, function(x) min(x))
# x.origin <- as.data.frame(sapply(x.origin, FUN = range01))

group.map <- c()
Z.list <- list()        
scale_stats_list <- list() 
projection_coefs_list <- list() 
spec_decomp_list <- list() 
qr_list <- list()          
sm_spec_list <- list()     
keep_cols_list <- list()

covariates <- colnames(data.frame(x.origin))
for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  x_vec <- x.origin[, i]
  X_lin <- model.matrix(~ x_vec) 
  
  sm_spec <- smoothCon(s(x_vec, bs = "tp", k = psi + 2), 
                       data = data.frame(x_vec = x_vec), 
                       knots = NULL)[[1]]
  
  X_raw <- sm_spec$X
  S     <- sm_spec$S[[1]] 
  
  eig <- eigen(S, symmetric = TRUE)
  max_lambda <- max(eig$values)
  tol <- max_lambda * 1e-6 
  
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
  group.map <- c(group.map, rep(i, ncol(Z_final)))
  spec_decomp_list[[i]] <- list(U_pen = U_pen, Lambda_sqrt_inv = solve(sqrt(Lambda_pen)))
  projection_coefs_list[[i]] <- Gamma_Original 
  keep_cols_list[[i]] <- keep_cols
  scale_stats_list[[i]] <- train_scale         
  sm_spec_list[[i]] <- sm_spec
}

bs.nonlinear <- do.call(cbind, Z.list)
bs.linear <- model.matrix(~ ., data = data.frame(x.origin))[,-1]

grid.n <- 200
newx <- seq(0, 1, length.out = grid.n)
xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))
# xholder <- do.call(cbind, lapply(1:p, function(j) {seq(min(x.origin[,j]), max(x.origin[,j]), length.out = grid.n)}))
grid_Z_list <- list()

for (i in seq_along(covariates)) {
  x_vec <- xholder[,i]
  grid_df  <- data.frame(x_vec = x_vec)
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
xholder.linear <- model.matrix(~ ., data = data.frame(xholder))[,-1]
alp.new <- as.vector(exp(theta.origin[1] + xholder %*% theta.origin[-1] + f2(xholder[,2]) + f5(xholder[,5])))
g.new <- c(newx*theta.origin[2], 
           newx*theta.origin[3] + f2(newx),
           newx*theta.origin[4],
           newx*theta.origin[5],
           newx*theta.origin[6] + f5(newx))

X_means <- colMeans(bs.linear)
X_sd   <- apply(bs.linear, 2, sd)
bs.linear <- scale(bs.linear, center = X_means, scale = X_sd)

# 6. Stan Model
model.stan <- "// Stan model for BLAST Pareto Samples
data {
    int <lower=1> n; 
    int <lower=1> grid_n;
    int <lower=1> p; 
    int <lower=1> psi; 
    vector[n] u; 
    matrix[n, p] bsLinear; 
    matrix[n, (psi*p)] bsNonlinear; 
    matrix[grid_n, p] xholderLinear; 
    matrix[grid_n, (psi*p)] xholderNonlinear;     
    vector<lower=u>[n] y; 
    real <lower=0> atau;
    vector[p] X_means;
    vector[p] X_sd;
    vector[(psi*p)] Z_scales;
}

parameters {
    vector[(p+1)] theta; 
    array[p] vector[psi] gamma_raw;
    array[p] real <lower=0> lambda1;  
    array[p] real <lower=0> lambda2; 
    array[p] real <lower=0> tau;
}

transformed parameters {
    vector[n] alpha; 
    
    array[p] vector[psi] gamma;
    {
      vector[n] eta = rep_vector(theta[1], n);
      for (j in 1:p){
        gamma[j] = gamma_raw[j] * sqrt(tau[j]);
        int nl_start = (j - 1) * psi + 1;
        eta += col(bsLinear, j) * theta[j+1] + block(bsNonlinear,1, nl_start, n, psi) * gamma[j];
      };
      alpha = exp(eta);
    }
}

model {
  target += pareto_lpdf(y | u, alpha);
  target += normal_lpdf(theta[1] | 0, 10);
  target += gamma_lpdf(lambda1 | 1e-1, 1e-1); 
  target += gamma_lpdf(lambda2 | 1e-2, 1e-2);  
  for (j in 1:p){
    target += double_exponential_lpdf(theta[(j+1)] | 0, 1/(lambda1[j]));
    target += gamma_lpdf(tau[j] | atau, square(lambda2[j])*0.5);
    target += std_normal_lpdf(gamma_raw[j]);
  }
}

generated quantities {
  vector[grid_n] gridalpha; 
  matrix[grid_n, p] gridgnl; 
  matrix[grid_n, p] gridgl; 
  matrix[grid_n, p] gridgsmooth; 

  vector[p] theta_origin = theta[2:(p+1)] ./ X_sd;
  real theta0 = theta[1] - dot_product(X_means, theta_origin);

  {
    vector[grid_n] grideta = rep_vector(theta0, grid_n);
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

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, grid_n = grid.n,
                  atau = ((psi+1)/2), X_means = X_means, X_sd=X_sd, Z_scales=Z_scales,
                  bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                  xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear)

init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p,psi)),
                        theta = rep(-0.1, (p+1)), tau = rep(0.1, p),
                        lambda1 = rep(0.1,p), lambda2 = rep(1, p)),
                   list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p,psi)),
                        theta = rep(0.05, (p+1)), tau = rep(2, p),
                        lambda1 = rep(2,p), lambda2 = rep(5, p)),
                   list(gamma_raw = array(rep(-0.75, (psi*p)), dim=c(p,psi)),
                        theta = rep(0.01, (p+1)), tau = rep(1.5, p),
                        lambda1 = rep(0.5,p), lambda2= rep(5, p)))

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
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("Scenario A") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/appendix_",n,"_traceplot_scA.pdf"), width=22, height = 3)

# theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
theta.samples <- summary(fit1, par=c("theta0", "theta_origin"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
# newgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
# newgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
# season.samples <- summary(fit1, par=c("alphaseason"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary

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
data.scenario <- data.frame("x" = seq(0,1,length.out=grid.n),
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


data.smooth <- data.frame("x"=seq(0,1,length.out=grid.n),
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
