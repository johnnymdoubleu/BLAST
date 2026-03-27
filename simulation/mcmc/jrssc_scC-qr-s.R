library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(qqboxplot)
library(evgam)
library(crch)
library(forecast)

set.seed(100)
n <- n.origin <- 10000
grid.n <- 200
psi.origin <- psi <- 10
threshold <- 0.95
p <- 5

C <- diag(p)
 
f1 <- function(x) { 0.4 * (1.2 - x)^3 } 
f5 <- function(x) { 0.4 * (1.1 - x)^3 }  

time.seq <- 1:n
period <- 365 
x.season <- time.seq / period 

season_code_full <- cut((time.seq%% period / period), breaks = c(-0.1, 0.25, 0.5, 0.75, 1.1), labels = c(1,2,3,4))
seasons <- c("Winter", "Spring", "Summer", "Autumn")

make.nl <- function(x, raw_y) {
  fit <- lm(raw_y ~ x)
  
  return(list(
    nl = residuals(fit), 
    slope = coef(fit)[["x"]],
    intercept = coef(fit)[["(Intercept)"]]
  ))
}

theta.origin <- c(0.4, 0, -0.6, 0, 0, 0.4)
psi <- psi -2

x.origin <- x.origin.full <- pnorm(matrix(rnorm(n.origin*p), nrow = n.origin, ncol = p))
alp.origin <- exp(theta.origin[1] + x.origin.full%*%theta.origin[-1] + 
                  f1(x.origin.full[,1]) + f5(x.origin.full[,5]))
y.origin <- rtt(n.origin, df = alp.origin, left=0)

# for (j in 1:p) {
#   phase_shift <- j * (2 * pi / p) 
#   seasonal_trend <- 0.5 * sin(2 * pi * time.seq / period - phase_shift)
#   x.origin[,j] <- seasonal_trend + x.origin.full[,j]
# }

# xreg.season <- cbind(
#   trend = time.seq,
#   cos_season = cos(2 * pi * time.seq / 365),
#   sin_season = sin(2 * pi * time.seq / 365)
# )

# fit.list <- list()
# x.detrended <- matrix(nrow = n.origin, ncol = p)
# for (j in 1:p) {
#   y_ts <- ts(x.origin[, j], frequency = period) 
#   fit.list[[j]] <- fit <- auto.arima(y_ts, seasonal = FALSE, xreg = xreg.season, stepwise = TRUE, approximation = FALSE)
#   x.detrended[, j] <- as.numeric(residuals(fit.list[[j]]))
# }
# x.origin <- x.detrended

f.season.scale <- function(t){
  return(2.5 - .8 * sin(2 * pi * t / 365) - .6 * cos(2 * pi * t / 365)) 
}
y.origin <- y.origin * f.season.scale(time.seq)

evgam.df <- data.frame(
  y = (y.origin),
  sin.time = sin(2 * pi * time.seq / 365),
  cos.time = cos(2 * pi * time.seq / 365),
  x.origin
)

evgam.cov <- y ~ cos.time + sin.time + s(X1) + s(X2) + s(X3) + s(X4) + s(X5)
ald.cov.fit <- evgam(evgam.cov, data = evgam.df, family = "ald", ald.args=list(tau = threshold))
u.vec <- (predict(ald.cov.fit, type = "response")$location)

excess.index <- which(y.origin > u.vec)
x.origin <- x.origin[excess.index,]
y.origin <- y.origin[excess.index]
u <- u.vec[excess.index]
season_code_full <- season_code_full[excess.index]  
n <- length(y.origin)

newx <- seq(max(apply(x.origin, 2, min)), min(apply(x.origin, 2, max)), length.out = grid.n)
xholder <- do.call(cbind, lapply(1:p, function(i) {seq(min(x.origin[,i]), max(x.origin[,i]), length.out = grid.n)}))
x.grid <- do.call(cbind, lapply(1:p, function(i) {seq(0, 1, length.out = grid.n)}))  
alp.new <- as.vector(exp(theta.origin[1] + x.grid %*% theta.origin[-1] + f1(x.grid[,1]) + f5(x.grid[,5])))

f1.hidden <- make.nl(x.origin[,1], f1(x.origin[,1]))
f5.hidden <- make.nl(x.origin[,5], f5(x.origin[,5]))
theta.adjusted <- c(theta.origin[1] + f1.hidden$intercept + f5.hidden$intercept,
                    theta.origin[2] + f1.hidden$slope,
                    theta.origin[3],
                    theta.origin[4],
                    theta.origin[5],
                    theta.origin[6] + f5.hidden$slope)

g1.nl <- f1(x.grid[,1]) - (f1.hidden$intercept + f1.hidden$slope*x.grid[,1])
g5.nl <- f5(x.grid[,5]) - (f5.hidden$intercept + f5.hidden$slope*x.grid[,5])

g1.l <- theta.adjusted[2]*x.grid[,1]
g2 <- g2.l <- theta.adjusted[3]*x.grid[,2]
g5.l <- theta.adjusted[6]*x.grid[,5]
g1 <- g1.l + g1.nl
g5 <- g5.l + g5.nl
# g1 <- g1 - mean(g1)
# g2 <- g2 - mean(g2)
# g5 <- g5 - mean(g5)
alp.new <- exp(theta.adjusted[1] + g1 + g2 + g5)
grid.zero <- rep(0, grid.n)
g.new <- c(g1, g2.l, grid.zero, grid.zero, g5)
l.new <- c(g1.l, g2.l, grid.zero, grid.zero, g5.l)
nl.new <- c(g1.nl, grid.zero, grid.zero, grid.zero, g5.nl)
colnames(xholder) <- colnames(x.origin) <- covariates <- paste0("X", 1:p)


group.map <- c()
Z.list <- list()        # Stores the final non-linear design matrices
scale_stats_list <- list() 
projection_coefs_list <- list() #
spec_decomp_list <- list() # Store eigen-decomp info for prediction
qr_list <- list()          # Store QR info for prediction
sm_spec_list <- list()     # Store smooth objects
keep_cols_list <- list()

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
  projection_coefs_list[[i]] <- Gamma_Original # Store the X-basis coefs
  keep_cols_list[[i]] <- keep_cols
  scale_stats_list[[i]] <- train_scale         # Store the scale
  sm_spec_list[[i]] <- sm_spec
}

bs.nonlinear <- do.call(cbind, Z.list)
bs.linear <- model.matrix(~ ., data = data.frame(x.origin))[,-1]
Z_scales <- unlist(scale_stats_list)

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

xholder.linear <- model.matrix(~ ., data = data.frame(xholder))[,-1]
xholder.nonlinear <- do.call(cbind, grid_Z_list)

X_means <- colMeans(bs.linear)
X_sd   <- apply(bs.linear, 2, sd)
bs.linear <- scale(bs.linear, center = X_means, scale = X_sd)


model.stan <- "// Stan model for BLAST t Samples
data {
    int <lower=1> n; // Sample size
    int <lower=1> grid_n;
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    vector[n] u; // large threshold value
    matrix[n, p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[grid_n, p] xholderLinear; // fwi dataset
    matrix[grid_n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector<lower=u>[n] y; // extreme response
    real <lower=0> atau;
    vector[p] X_means;
    vector[p] X_sd;
}

parameters {
    vector[(p+1)] theta; // linear predictor
    array[p] vector[psi] gamma_raw;
    array[p] real <lower=0> lambda1; // 
    array[p] real <lower=0> lambda2; //
    array[p] real <lower=0> tau;
}

transformed parameters {
    vector[n] alpha; // covariate-adjusted tail index
    
    array[p] vector[psi] gamma;
    {
      vector[n] eta = rep_vector(theta[1], n);
      for (j in 1:p){
        gamma[j] = gamma_raw[j] * sqrt(tau[j]); 
        eta += col(bsLinear, j) * theta[j+1] + block(bsNonlinear,1, ((j - 1) * psi + 1), n, psi) * gamma[j];
      };
      
      alpha = exp(eta);
    }
}

model {
  //likelihood
  target += pareto_lpdf(y | u, alpha);
  target += normal_lpdf(theta[1] | 0, 10);
  target += gamma_lpdf(lambda1 | 1e-2, 1e-2); 
  target += gamma_lpdf(lambda2 | 1e-2, 1e-2);
  for (j in 1:p){
    target += double_exponential_lpdf(theta[(j+1)] | 0, 1/(lambda1[j]));
    target += gamma_lpdf(tau[j] | atau, square(lambda2[j])*0.5);
    target += std_normal_lpdf(gamma_raw[j]);
  }
}

generated quantities {
  // Used in Posterior predictive check
  vector[grid_n] gridalpha; // new tail index
  matrix[grid_n, p] gridgnl; // nonlinear component
  matrix[grid_n, p] gridgl; // linear component
  matrix[grid_n, p] gridgsmooth; // linear component

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
                  atau = ((psi+1)/2), X_means = X_means, X_sd=X_sd,
                  bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                  xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear)

init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p,psi)),
                        theta = rep(-0.1, (p+1)), tau = rep(0.1, p), 
                        lambda1 = rep(0.1, p), lambda2 = rep(1, p)),
                  list(gamma_raw = array(rep(0.15, (psi*p)), dim=c(p,psi)),
                        theta = rep(-0.05, (p+1)), tau = rep(0.3, p),
                        lambda1 = rep(0.2, p), lambda2 = rep(0.5, p)),
                  list(gamma_raw = array(rep(0.1, (psi*p)), dim=c(p,psi)),
                        theta = rep(0.05, (p+1)), tau = rep(0.2, p),
                        lambda1 = rep(0.2, p), lambda2 = rep(0.5, p)))

system.time(fit1 <- stan(
  model_code = model.stan,  # Stan program
  data = data.stan,    # named list of data
  init = init.alpha,      # initial value
  chains = 3,             # number of Markov chains
  iter = 2000,            # total number of iterations per chain
  cores = parallel::detectCores(), # number of cores (could use one per chain)
  refresh = 1000             # no progress shown
))

posterior <- extract(fit1)
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("Scenario B") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/appendix_",n,"_traceplot_scA.pdf"), width=22, height = 3)

theta.samples <- summary(fit1, par=c("theta0", "theta_origin"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
# newgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
# newgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
# season.samples <- summary(fit1, par=c("alphaseason"), probs = c(0.05,0.5, 0.95))$summary
# alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary

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
  xlab(expression(c)) + ylab("") +
  geom_ribbon(aes(ymin = q1, ymax = q3), fill="steelblue", alpha = 0.2) +
  geom_line(aes(y = true), color = "red", linewidth = 2, linetype=2) +
  geom_line(aes(y=q2), color = "steelblue", linewidth=1.5) +
  theme_minimal(base_size = 40) + ylim(0, 15) +
  theme(legend.position = "none",
          strip.text = element_blank(),
          axis.text.y = element_blank(),
          axis.text = element_text(size = 30))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_alpha_scC_",n.origin,"_ct.pdf"), width=9.5, height = 7.78)

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
  theme_minimal(base_size = 30) + ylim(-1.5, 1.5) + 
  theme(legend.position = "none",
          plot.margin = margin(0,0,0,-20),
          axis.text.y = element_blank(),
          strip.text = element_blank(),
          axis.title.x = element_text(size = 45),                
          axis.text = element_text(size = 30))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_smooth_scC_",n.origin,"_ct.pdf"), width=11, height = 15)

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

# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_linear_scC.pdf"), width=12.5, height = 15)


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
