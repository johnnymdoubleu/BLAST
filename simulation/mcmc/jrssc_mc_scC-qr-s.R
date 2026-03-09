library(mgcv)
library(Pareto)
suppressMessages(library(ggplot2))
library(rstan)
library(MESS)
library(evgam)
library(crch)
library(forecast)
 
# Scenario C
# array.id <- commandArgs(trailingOnly=TRUE)

total.iter <- 5

n <- n.origin <- 10000
grid.n <- 200
psi.origin <- psi <- 10
threshold <- 0.95
p <- 5

C <- diag(p)

f1 <- function(x) { 1.7 * (1.2 - x)^2 }  # Quadratic decrease: ~4.8 (at x=0) to ~0.4 (at x=1)
f5 <- function(x) { -0.4 * (1.1 - x)^2 }   # Quadratic decrease: ~3.9 (at x=0) to ~0.3 (at x=1)

time.seq <- 1:n
period <- 365 
x.season <- time.seq / period 

# Convert continuous season to Factor for the 'by' argument
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
    vector<lower=min(u)>[n] y; // extreme response
    real <lower=0> atau;
    vector[p] X_means;
    vector[p] X_sd;
    vector[(psi*p)] Z_scales;
    vector[grid_n] trueAlpha;
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
  target += student_t_lpdf(y | alpha, 0, 1) - student_t_lccdf(u | alpha, 0, 1);
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
  // Used in Posterior predictive check
  vector[grid_n] gridalpha; // new tail index
  matrix[grid_n, p] gridgnl; // nonlinear component
  matrix[grid_n, p] gridgl; // linear component
  matrix[grid_n, p] gridgsmooth; // linear component
  vector[grid_n] se;

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
  se = pow((gridalpha-trueAlpha), 2);
}
"

gridgnl.container <- as.data.frame(matrix(, nrow = (p*grid.n), ncol = total.iter))
gridgl.container <- as.data.frame(matrix(, nrow = (p*grid.n), ncol = total.iter))
gridgsmooth.container <- as.data.frame(matrix(, nrow = (p*grid.n), ncol = total.iter))
alpha.container <- as.data.frame(matrix(, nrow = grid.n, ncol = total.iter))
true.container <- as.data.frame(matrix(, nrow = grid.n, ncol = total.iter))
mise.container <- c()
# qqplot.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))

for(iter in 1:total.iter){
  n <- n.origin
  x.origin.full <- matrix(0, nrow = n.origin, ncol = p)

  for (j in 1:p) {
    # seasonal_trend <- 0.1 * sin(2 * pi * time.seq / period) 
    # uniform_noise <- runif(n, min = -0.49, max = 0.49)
    # x.origin.full[, j] <- seasonal_trend + uniform_noise
    # ar_noise <- arima.sim(model = list(ar = 0.5), n = n.origin, sd=0.5)
    # x.origin.full[,j] <-(seasonal_trend + ar_noise)
    x.origin.full[,j] <- pnorm(rnorm(n))
  }
  x.origin <- x.origin.full
  # plot(x.origin.full[,j])
  # fit.list <- list()
  # x.detrended <- matrix(nrow = n.origin, ncol = p)
  # for (j in 1:p) {
    # y_ts <- ts(x.origin.full[, j], frequency = period) 
    # decomp <- stl(y_ts, s.window = "periodic", robust = TRUE)
    # seasonal <- as.numeric(decomp$time.series[, "seasonal"])
    # x.detrended[,j] <- x.origin.full[, j] - seasonal
    # fit.list[[j]] <- forecast::auto.arima(x.origin.full[, j] - seasonal, seasonal = FALSE, stepwise = TRUE)
    # fit.list[[j]] <- fit <- auto.arima(y_ts, seasonal = TRUE, stepwise = TRUE, approximation = FALSE)
    # x.detrended[, j] <- as.numeric(residuals(fit.list[[j]]))
  # }
  # x.origin <- x.detrended
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  # X_min <- apply(x.origin, 2, min)
  # X_minmax <- sapply(x.origin, function(x) max(x)-min(x))
  # x.origin <- (sapply(as.data.frame(x.origin), FUN = range01))

  alp.origin <- exp(rep(theta.origin[1],n) + x.origin%*%theta.origin[-1] + f1(x.origin[,1]) + f5(x.origin[,5]))
  # y.noise <- rt(n, df = alp.origin)
  y.noise <- rtt(n, df = alp.origin, left=0)
  f.season.scale <- function(t){
    return(2.5 - .8 * sin(2 * pi * t / 365) - .6 * cos(2 * pi * t/365)) 
  }
  y.origin <- y.noise * f.season.scale(time.seq)

  evgam.df <- data.frame(
    y = (y.origin),
    sin.time = sin(2 * pi * time.seq / 365),
    cos.time = cos(2 * pi * time.seq / 365),
    x.origin
  )
  evgam.cov <- y ~ 1 + cos.time + sin.time #+ s(X1,k=3) + s(X2,k=3) + s(X3,k=3) + s(X4,k=3) + s(X5,k=3)
  ald.cov.fit <- evgam(evgam.cov, data = evgam.df, family = "ald", ald.args=list(tau = threshold))
  u.vec <- (predict(ald.cov.fit)$location)

  excess.index <- which(y.origin > u.vec)

  x.origin <- x.origin[excess.index,]
  y.origin <- y.origin[excess.index]
  u <- u.vec[excess.index]
  season_code_full <- season_code_full[excess.index]  
  n <- length(y.origin)
  colnames(x.origin) <- paste0("X", 1:p)
  # x.origin <- (sapply(as.data.frame(x.origin), FUN = range01))
  newx <- seq(max(apply(x.origin, 2, min)), min(apply(x.origin, 2, max)), length.out = grid.n)
  # xholder <- do.call(cbind, lapply(1:p, function(i) {seq(min(x.origin[,i]), max(x.origin[,i]), length.out = grid.n)}))  
  # newx <- seq(0.1, 0.9, length.out = grid.n)
  xholder <- do.call(cbind, lapply(1:p, function(i) {newx}))
  # newx <- seq(0, 1, length.out = grid.n)
  # xholder <- do.call(cbind, lapply(1:p, function(i) {newx}))  
  alp.new <- as.vector(exp(theta.origin[1] + xholder %*% theta.origin[-1] + 
                          f1(xholder[,1]) + f5(xholder[,5])))

  f1.hidden <- make.nl(x.origin[,1], f1(x.origin[,1]))
  f5.hidden <- make.nl(x.origin[,5], f5(x.origin[,5]))
  theta.adjusted <- c(theta.origin[1] + f1.hidden$intercept + f5.hidden$intercept,
                      theta.origin[2] + f1.hidden$slope,
                      theta.origin[3],
                      theta.origin[4],
                      theta.origin[5],
                      theta.origin[6] + f5.hidden$slope)

  g1.nl <- f1(xholder[,1]) - (f1.hidden$intercept + f1.hidden$slope*xholder[,1])
  g5.nl <- f5(xholder[,5]) - (f5.hidden$intercept + f5.hidden$slope*xholder[,5])

  g1.l <- theta.adjusted[2]*xholder[,1]
  g2.l <- theta.adjusted[3]*xholder[,2]
  g5.l <- theta.adjusted[6]*xholder[,5]
  g1 <- g1.l + g1.nl
  g5 <- g5.l + g5.nl
  eta.g <- rep(theta.adjusted[1], grid.n) + g1 + g2.l + g5
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
  Z_scales <- unlist(scale_stats_list)
  # xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))

  grid_Z_list <- list()

  for (i in seq_along(covariates)) {
    x_vec <- xholder[,i] #(seq(min(x.origin[,i]), max(x.origin[,i]), length.out = grid.n))
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

  X_means <- colMeans(x.origin)
  X_sd   <- apply(x.origin, 2, sd)
  bs.linear <- scale(bs.linear, center = X_means, scale = X_sd)

  data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, grid_n = grid.n,
                  atau = ((psi+1)/2), X_means = X_means, X_sd=X_sd, Z_scales=Z_scales,
                  bsLinear = bs.linear, bsNonlinear = bs.nonlinear, trueAlpha = alp.new,
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
    
  fit1 <- stan(
      model_code = model.stan,
      data = data.stan,    # named list of data
      init = init.alpha,      # initial value
      chains = 3,             # number of Markov chains
      iter = 2000,            # total number of iterations per chain
      cores = parallel::detectCores(), # number of cores (could use one per chain)
      refresh = 1000             # no progress shown
  )
  # posterior <- extract(fit1)
  theta.samples <- summary(fit1, par=c("theta0", "theta_origin"), probs = c(0.5))$summary
  lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.5))$summary
  gridgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.5))$summary
  gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.5))$summary
  gridgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.5))$summary
  newalpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.5))$summary
  se.samples <- summary(fit1, par=c("se"), probs = c(0.5))$summary 

  alpha.container[,iter] <- newalpha.samples[,4]
  true.container[,iter] <- alp.new

  gridgsmooth.container[,iter] <- as.vector(matrix(gridgsmooth.samples[,4], nrow = grid.n, byrow=TRUE))
  gridgl.container[,iter] <- as.vector(matrix(gridgl.samples[,4], nrow = grid.n, byrow=TRUE))
  gridgnl.container[,iter] <- as.vector(matrix(gridgnl.samples[,4], nrow = grid.n, byrow=TRUE))
  # newx <- seq(0.1, , length.out = grid.n)
  mise.container[iter] <- auc(newx, se.samples[,4], type="spline")

  # mcmc.alpha <- rstan::extract(fit1)$alpha
  # r <- matrix(, nrow = n, ncol = 30)
  # T <- 30
  # for(i in 1:n){
  #     for(t in 1:T){
  #         r[i, t] <- qnorm(pPareto(y.origin[i], u[i], alpha = mcmc.alpha[round(runif(1,1, dim(mcmc.alpha)[1])),i]))
  #     }
  # }
  # lgrid <- n
  # grid <- qnorm(ppoints(lgrid))
  # traj <- matrix(NA, nrow = T, ncol = lgrid)
  # for (t in 1:T){
  #     traj[t, ] <- quantile(r[, t], ppoints(lgrid), type = 2)
  # }
  # T <- 500
  # mcmc.alpha <- rstan::extract(fit1)$alpha
  # len    <- nrow(mcmc.alpha)
  # posterior_idx <- sample(len, T, replace = TRUE)
  # alpha_sub <- mcmc.alpha[posterior_idx, ] 
  # y.matrix <- matrix(y.origin, nrow = T, ncol = n, byrow = TRUE)
  # u.matrix <- matrix(u, nrow = T, ncol = n, byrow = TRUE)
  # r_vec <- qnorm(pPareto(y.matrix, u.matrix, alpha = alpha_sub))
  # r_mat <- matrix(r_vec, nrow = T, ncol = n)
  # quantile_prob <- ppoints(n)
  # grid          <- qnorm(quantile_prob)
  # traj <- t(apply(r_mat, 1, sort))  
  # qqplot.container[iter] <- apply(traj, 2, mean)#quantile, prob = 0.5)
}

newx <- seq(0, 1, length.out = grid.n)
xholder <- do.call(cbind, lapply(1:p, function(i) {newx}))  
alp.new <- as.vector(exp(theta.origin[1] + xholder %*% theta.origin[-1] + f1(xholder[,1]) + f5(xholder[,5])))
alpha.container$x <- seq(0, 1, length.out = grid.n)
alpha.container$true <- alp.new
# alpha.container$true <- rowMeans(true.container)
alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
alpha.container <- as.data.frame(alpha.container)

# load(paste0("./simulation/results/MC-Scenario_A/2026-02-07_",total.iter,"_MC_scC_",n.origin,".Rdata"))

plt <- ggplot(data = alpha.container, aes(x = x)) + xlab(expression(c)) + labs(col = "") + ylab(expression(alpha(c,...,c))) #+ ylab("")
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
} else{
  for(i in 50:100){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}
print(plt + 
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.8) + ylim(0,10) +
        scale_fill_manual(values=c("steelblue"), name = "") +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 40) + 
        theme(legend.position = "none",
                strip.text = element_blank(),
                axis.text = element_text(size = 30)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_alpha_scC_",n.origin,".pdf"), width=10, height = 7.78)

gridgsmooth.container$x <- seq(0, 1, length.out = grid.n)
gridgsmooth.container$true <- g.new
gridgsmooth.container$mean <- rowMeans(gridgsmooth.container[,1:total.iter])
gridgsmooth.container$covariate <- gl(p, grid.n, (p*grid.n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]"))
gridgsmooth.container <- as.data.frame(gridgsmooth.container)

plt <- ggplot(data = gridgsmooth.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(c))
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(gridgsmooth.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
} else{
  # for(i in 1:total.iter){
  for(i in 50:100){
    plt <- plt + geom_line(aes(y = .data[[names(gridgsmooth.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}
print(plt + geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + ylim(-3, 3) +
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5) + 
        facet_grid(covariate ~ ., scales = "free_x", switch = "y",
                    labeller = label_parsed) +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                plot.margin = margin(0,0,0,-20),
                strip.text.y = element_text(size = 38, colour = "black", angle = 0, face = "bold.italic"),
                strip.placement = "outside",
                axis.title.x = element_text(size = 45),                
                axis.text = element_text(size = 30)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_smooth_scA_",n.origin,".pdf"), width=12.5, height = 15)


gridgl.container$x <- newx
gridgl.container$true <- as.vector(l.new)
gridgl.container$mean <- rowMeans(gridgl.container[,1:total.iter])
gridgl.container$covariate <- gl(p, grid.n, (p*grid.n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]"))
gridgl.container <- as.data.frame(gridgl.container)

plt <- ggplot(data = gridgl.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(x))
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(gridgl.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}else{
  for(i in 1:50){
    plt <- plt + geom_line(aes(y = .data[[names(gridgl.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}
print(plt + geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.2) + 
        ylim(-3, 3)+
        facet_grid(covariate ~ ., scales = "free_x", switch = "y",
                    labeller = label_parsed) +  
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                plot.margin = margin(0,0,0,-20),
                strip.text = element_blank(),
                axis.title.x = element_text(size = 45),  
                axis.text = element_text(size = 30)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_linear_scC",n.origin,".pdf"), width=12.5, height = 7.78)                

gridgnl.container$x <- newx
gridgnl.container$true <- as.vector(nl.new)
gridgnl.container$mean <- rowMeans(gridgnl.container[,1:total.iter])
gridgnl.container$covariate <- gl(p, grid.n, (p*grid.n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]"))
gridgnl.container <- as.data.frame(gridgnl.container)

plt <- ggplot(data = gridgnl.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(x))
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(gridgnl.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}else{
  for(i in 1:50){
    plt <- plt + geom_line(aes(y = .data[[names(gridgnl.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}
print(plt + geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.2) + 
        ylim(-2.5,2.5)+
        facet_grid(covariate ~ ., scales = "free_x", switch = "y",
                    labeller = label_parsed) +  
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                plot.margin = margin(0,0,0,-20),
                strip.text = element_blank(),
                axis.title.x = element_text(size = 45),  
                axis.text = element_text(size = 30)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_nonlinear_scC",n.origin,".pdf"), width=12.5, height = 7.78)

# qqplot.container$grid <- grid
# qqplot.container$mean <- rowMeans(qqplot.container[,1:total.iter])
# plt <- ggplot(data = qqplot.container, aes(x = grid))
# if(total.iter <= 50){
#   for(i in 1:total.iter){
#     plt <- plt + geom_line(aes(y = .data[[names(qqplot.container)[i]]]), alpha = 0.05, linewidth = 0.7)
#   }
# }else{
#   for(i in 50:100){
#     plt <- plt + geom_line(aes(y = .data[[names(qqplot.container)[i]]]), alpha = 0.05, linewidth = 0.7)
#   }
# }
# print(plt + 
#         geom_line(aes(y = mean), colour = "steelblue", linewidth = 1.5, linetype = 2) + 
#         labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
#         theme_minimal(base_size = 30) +
#         theme(text = element_text(size = 20)) +
#         coord_fixed(xlim = c(-2, 2),  
#                     ylim = c(-2, 2)))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_qqplot_scC_",n.origin,".pdf"), width=10, height = 7.78)

# save(alpha.container, gridgnl.container, gridgl.container, gridgsmooth.container, mise.container, file = paste0(Sys.Date(),"_",total.iter,"_MC_scC_",n.origin,"_",array.id,".Rdata"))

mean(mise.container)
