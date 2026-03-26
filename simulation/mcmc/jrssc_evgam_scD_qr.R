library(VGAM)
library(Pareto)
suppressMessages(library(ggplot2))
library(rstan)
library(MESS)
library(evgam)
library(rmutil)
library(forecast)
library(mgcv)

# Scenario D
# array.id <- commandArgs(trailingOnly=TRUE)

total.iter <- 2
n <- n.origin <- 10000
grid.n <- 200
psi.origin <- psi <- 10
threshold <- 0.95
p <- 5

C <- diag(p)
f2 <- function(x) {-1.5 * sin(2 * pi * (x-1.1)^2)*(x-1.1)^3}
f3 <- function(x) {0.5 * cos(3 * pi * (x)^2)*(x)^2}

time.seq <- 1:n
period <- 365 
x.season <- time.seq / period 

# Convert continuous season to Factor for the 'by' argument
season_code_full <- cut((time.seq%% period / period), breaks = c(-0.1, 0.25, 0.5, 0.75, 1.1), labels = c(1,2,3,4))
seasons <- c("Winter", "Spring", "Summer", "Autumn")


theta.origin <- c(0.7, 1.2, 0, -1.2, 1.2, 0)
psi <- psi -2

model.stan <- "// Stan model for BLAST Burr Samples
data {
  int <lower=1> n; // Sample size
  int <lower=1> grid_n; //grid size
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
  vector[grid_n] trueAlpha;
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
  matrix[grid_n, p] gridgsmooth; // linear component
  vector[grid_n] se;

  vector[p] theta_origin = theta[2:(p+1)] ./ X_sd;
  real theta0 = theta[1] - dot_product(X_means, theta_origin);

  {
    vector[grid_n] grideta = rep_vector(-theta0, grid_n);
    for (j in 1:p){
      gridgsmooth[,j] = col(xholderLinear, j) * theta_origin[j] + block(xholderNonlinear,1, ((j - 1) * psi + 1), grid_n, psi) * gamma[j];
      grideta -= gridgsmooth[,j];
    };
    gridalpha = exp(grideta);
  }
  se = pow((gridalpha-trueAlpha), 2);
}
"

newgsmooth.container <- as.data.frame(matrix(, nrow = (p*grid.n), ncol = total.iter))
smooth.scale.container <- as.data.frame(matrix(, nrow = (p*grid.n), ncol = total.iter))
smooth.1.container <- as.data.frame(matrix(, nrow = (p*grid.n), ncol = total.iter))
alpha.container <- as.data.frame(matrix(, nrow = grid.n, ncol = total.iter))
true.container <- as.data.frame(matrix(, nrow = grid.n, ncol = total.iter))
vgam.scale.container <- vgam.1.container <- evgam.scale.container <- evgam.1.container <- as.data.frame(matrix(, nrow = grid.n, ncol = total.iter))
mise.vgam.scale.container <- mise.vgam.1.container <- mise.evgam.scale.container <- mise.evgam.1.container <- mise.container <- c()

for(iter in 1:total.iter){
  n <- n.origin
  x.origin <- x.origin.full <- pnorm(matrix(rnorm(n.origin*p), nrow = n.origin, ncol = p))
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  alp.origin <- exp(theta.origin[1] + x.origin.full%*%theta.origin[-1] + f2(x.origin.full[,2]) + f3(x.origin.full[,3]))

  y.origin <- NULL
  for(i in 1:n){
    y.origin[i] <- rmutil::rburr(1, m=1, s=alp.origin[i], f=1)
  }

  for (j in 1:p) {
    phase_shift <- j * (2 * pi / p) 
    seasonal_trend <- 0.5 * sin(2 * pi * time.seq / period - phase_shift)
    x.origin[,j] <- seasonal_trend + x.origin.full[,j]
  }

  xreg.season <- cbind(
    trend = time.seq,
    cos_season = cos(2 * pi * time.seq / 365),
    sin_season = sin(2 * pi * time.seq / 365)
  )

  fit.list <- list()
  x.detrended <- matrix(nrow = n.origin, ncol = p)
  for (j in 1:p) {
    y_ts <- ts(x.origin[, j], frequency = period) 
    fit.list[[j]] <- fit <- auto.arima(y_ts, seasonal = FALSE, xreg = xreg.season, stepwise = TRUE, approximation = FALSE)
    x.detrended[, j] <- as.numeric(residuals(fit.list[[j]]))
  }
  x.origin <- x.detrended 

  f.season.scale <- function(t){
    return(2.5 - .8 * sin(2 * pi * t / 365) - .6 * cos(2 * pi * t/365)) 
  }
  y.origin <- y.origin * f.season.scale(time.seq)
  evgam.df <- data.frame(
    y = log(y.origin),
    sin.time = sin(2 * pi * time.seq / 365),
    cos.time = cos(2 * pi * time.seq / 365),
    x.season = (time.seq %% period) / period,
    x.origin
  )
  evgam.cov <- y ~ 1 + cos.time + sin.time + s(X1) + s(X2) + s(X3) + s(X4) + s(X5)
  ald.cov.fit <- evgam(evgam.cov, data = evgam.df, family = "ald", ald.args=list(tau = threshold))
  u.vec <- exp(predict(ald.cov.fit)$location)

  excess.index <- which(y.origin > u.vec)
  x.origin <- data.frame(x.origin[excess.index,])
  y.origin <- y.origin[excess.index]
  u <- u.vec[excess.index]
  n <- length(y.origin)

  
  newx <- seq(max(apply(x.origin, 2, min)), min(apply(x.origin, 2, max)), length.out = grid.n)
  xholder <- do.call(cbind, lapply(1:p, function(i) {seq(min(x.origin[,i]), max(x.origin[,i]), length.out = grid.n)}))
  x.grid <- do.call(cbind, lapply(1:p, function(i) {seq(0, 1, length.out = grid.n)}))  
  alp.new <- as.vector(exp(theta.origin[1] + x.grid %*% theta.origin[-1] + f2(x.grid[,2]) + f3(x.grid[,3])))
  colnames(xholder) <- colnames(x.origin) <- covariates <- paste0("X", 1:p)

  group.map <- c()
  Z.list <- list()
  scale_stats_list <- list() 
  projection_coefs_list <- list()
  spec_decomp_list <- list() # Store eigen-decomp info for prediction
  sm_spec_list <- list()     # Store smooth objects
  keep_cols_list <- list()

  for (i in seq_along(covariates)) {
    var_name <- covariates[i]
    x_vec <- x.origin[, i]
    X_lin <- model.matrix(~ x_vec) 
    sm_spec <- smoothCon(mgcv::s(x_vec, bs = "tp", k = psi + 2), 
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
  
  data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                    X_sd = X_sd, grid_n = grid.n,
                    atau = ((psi+1)/2), X_means = X_means, trueAlpha = 1/alp.new,
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
    
  fit1 <- stan(
      model_code = model.stan,  # Stan program
      data = data.stan,    # named list of data
      init = init.alpha,      # initial value
      chains = 3,             # number of Markov chains
      iter = 2000,            # total number of iterations per chain
      cores = parallel::detectCores(), # number of cores (could use one per chain)
      refresh = 1000             # no progress shown
  )

  newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.5))$summary
  newalpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.5))$summary
  se.samples <- summary(fit1, par=c("se"), probs = c(0.5))$summary

  simul.data <- data.frame(y = y.origin, x.origin)
  vgam.fit.scale <- vgam(y ~ sm.ps(X1, ps.int = 8, outer.ok=TRUE) + 
                             sm.ps(X2, ps.int = 8, outer.ok=TRUE) + 
                             sm.ps(X3, ps.int = 8, outer.ok=TRUE) + 
                             sm.ps(X4, ps.int = 8, outer.ok=TRUE) + 
                             sm.ps(X5, ps.int = 8, outer.ok=TRUE),
                        data = simul.data,
                        family = gpd(threshold= u,
                                      lshape="loglink",
                                      zero = NULL),
                        trace = TRUE,
                        control = vgam.control(maxit = 500))
  fitted.linear <- predict(vgam.fit.scale,newdata = data.frame(xholder), type = "link")
  fitted.terms <- predict(vgam.fit.scale,newdata = data.frame(xholder), type = "terms")
  vgam.xi.scale <- exp(fitted.linear[,2])

  vgam.fit.1 <- vgam(y ~ sm.ps(X1, ps.int = 8, outer.ok=TRUE) + 
                         sm.ps(X2, ps.int = 8, outer.ok=TRUE) + 
                         sm.ps(X3, ps.int = 8, outer.ok=TRUE) + 
                         sm.ps(X4, ps.int = 8, outer.ok=TRUE) + 
                         sm.ps(X5, ps.int = 8, outer.ok=TRUE),
                        data = simul.data,
                        family = gpd(threshold= u,
                                      lshape="loglink",
                                      zero = 1),
                        trace = TRUE,
                        control = vgam.control(maxit = 500))
  fitted.linear <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "link")
  fitted.terms <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "terms")
  vgam.xi.1 <- exp(fitted.linear[,2])
  
  simul.data <- data.frame(y = y.origin-u, x.origin)
  gam.scale <- list(y ~ s(X1, bs = "ts", k = 10) + 
                        s(X2, bs = "ts", k = 10) + 
                        s(X3, bs = "ts", k = 10) + 
                        s(X4, bs = "ts", k = 10) + 
                        s(X5, bs = "ts", k = 10),
                      ~ s(X1, bs = "ts", k = 10) + 
                        s(X2, bs = "ts", k = 10) + 
                        s(X3, bs = "ts", k = 10) + 
                        s(X4, bs = "ts", k = 10) + 
                        s(X5, bs = "ts", k = 10))
  evgam.fit.scale <- evgam::evgam(gam.scale, data = simul.data, family = "gpd")
  xi.pred.scale <-predict(evgam.fit.scale, newdata = data.frame(xholder), type="response")$shape

  gam.1 <- list(y ~ 1,
                  ~ s(X1, bs = "ts", k = 10) + 
                    s(X2, bs = "ts", k = 10) + 
                    s(X3, bs = "ts", k = 10) + 
                    s(X4, bs = "ts", k = 10) + 
                    s(X5, bs = "ts", k = 10))
  evgam.fit.1 <- evgam::evgam(gam.1, data = simul.data, family = "gpd")
  xi.pred.1 <-predict(evgam.fit.1, newdata = data.frame(xholder), type="response")$shape

  vgam.1.container[,iter] <- vgam.xi.1
  vgam.scale.container[,iter] <- vgam.xi.scale
  evgam.1.container[,iter] <- xi.pred.1
  evgam.scale.container[,iter] <- xi.pred.scale
  alpha.container[,iter] <- newalpha.samples[,4]
  newgsmooth.container[,iter] <- as.vector(matrix(newgsmooth.samples[,4], nrow = grid.n, byrow=TRUE))
  true.container[,iter] <- 1/alp.new
  mise.container[iter] <- auc(newx, se.samples[,4], type="spline")
  mise.evgam.1.container[iter] <- auc(newx, (((1/alp.new)-xi.pred.1))^2, type="spline")
  mise.evgam.scale.container[iter] <- auc(newx, ((1/alp.new)-xi.pred.scale)^2,type="spline")
  mise.vgam.1.container[iter] <- auc(newx, (((1/alp.new)-vgam.xi.1))^2, type="spline")
  mise.vgam.scale.container[iter] <- auc(newx, ((1/alp.new)-vgam.xi.scale)^2,type="spline") 
}


alpha.container$x <- seq(0, 1, length.out = grid.n)
alpha.container$true <- rowMeans(true.container)
alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
alpha.container$evgam.1 <- rowMeans(evgam.1.container[,1:total.iter])
alpha.container$evgam.scale <- rowMeans(evgam.scale.container[,1:total.iter])
alpha.container$vgam.1 <- rowMeans(vgam.1.container[,1:total.iter])
alpha.container$vgam.scale <- rowMeans(vgam.scale.container[,1:total.iter])
alpha.container <- as.data.frame(alpha.container)

# save(newgsmooth.container, alpha.container, mise.container, evgam.1.container, evgam.scale.container, mise.evgam.1.container, mise.evgam.scale.container, vgam.1.container, vgam.scale.container, mise.vgam.1.container, mise.vgam.scale.container, file=paste0("evgam_mc_scD_",n.origin,"_",array.id ,".Rdata"))
# load(paste0("./simulation/results/2026-02-10_evgam_mc_scD_",(n.origin*0.05),".Rdata"))

plt <- ggplot(data = alpha.container, aes(x = x)) + xlab(expression(c)) + labs(col = "") + ylab("")
plot_limit <- min(total.iter, 50)
for(i in 1:plot_limit){
  plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.05, linewidth = 0.7)
}

print(plt +
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.8) +
        geom_line(aes(y=evgam.1), colour="purple", linewidth = 1.8, linetype = 3) +
        geom_line(aes(y=evgam.scale), colour="orange", linewidth = 1.8, linetype = 3) +
        geom_line(aes(y=vgam.1), colour="purple", linewidth = 1.8, linetype = 4) +
        geom_line(aes(y=vgam.scale), colour="orange", linewidth = 1.8, linetype = 4) +
        scale_fill_manual(values=c("steelblue"), name = "") +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 40) + ylim(-0.65, 3.55) +
        theme(legend.position = "none",
                strip.text = element_blank(),
                axis.text.y = element_blank(),
                axis.text = element_text(size = 30)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_evgam_scD_",n.origin, ".pdf"), width=9.5, height = 7.78)


# newgsmooth.container$x <- seq(0,1, length.out = n)
# newgsmooth.container$true <- as.vector(g.new)
# newgsmooth.container <- cbind(newgsmooth.container, t(apply(newgsmooth.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
# colnames(newgsmooth.container)[(dim(newgsmooth.container)[2]-2):(dim(newgsmooth.container)[2])] <- c("q1","q2","q3")
# newgsmooth.container$mean <- rowMeans(newgsmooth.container[,1:total.iter])
# newgsmooth.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]"))
# newgsmooth.container <- as.data.frame(newgsmooth.container)

# plt <- ggplot(data = newgsmooth.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(c))
# if(total.iter <= 50){
#   for(i in 1:total.iter){
#     plt <- plt + geom_line(aes(y = .data[[names(newgsmooth.container)[i]]]), alpha = 0.075, linewidth = 0.7)
#   }
# } else{
#   for(i in 1:total.iter){
#     plt <- plt + geom_line(aes(y = .data[[names(newgsmooth.container)[i]]]), alpha = 0.075, linewidth = 0.7)
#   }
# }
# print(plt + 
#         geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
#         geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5) + 
#         ylim(-1, 1) +
#         facet_grid(covariate ~ ., scales = "free_x", switch = "y",
#                     labeller = label_parsed) +        
#         scale_color_manual(values = c("steelblue", "red"))+
#         guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) +
#         theme_minimal(base_size = 30) +
#         theme(legend.position = "none",
#                 plot.margin = margin(0,0,0,-20),
#                 strip.text.y = element_text(size = 38, colour = "black", angle = 0, face = "bold.italic"),
#                 strip.placement = "outside",
#                 axis.title.x = element_text(size = 40),
#                 axis.text = element_text(size = 30)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_smooth_scD_",n.origin,".pdf"), width=11, height = 15)

cat("BLAST:   ", mean(mise.container, na.rm=TRUE), "±", sd(mise.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.container))), "\n", 
    "EVGAM:  ", mean(mise.evgam.1.container, na.rm=TRUE), "±", sd(mise.evgam.1.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.evgam.1.container))), "\n",
    "EVGAM-σ:", mean(mise.evgam.scale.container, na.rm=TRUE), "±", sd(mise.evgam.scale.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.evgam.scale.container))), "\n",
    "VGAM:   ", mean(mise.vgam.1.container, na.rm=TRUE), "±", sd(mise.vgam.1.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.vgam.1.container))), "\n",
    "VGAM-σ: ", mean(mise.vgam.scale.container, na.rm=TRUE), "±", sd(mise.vgam.scale.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.vgam.scale.container))), "\n")
