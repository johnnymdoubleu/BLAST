library(mgcv)
library(Pareto)
suppressMessages(library(tidyverse))
library(rstan)
library(MESS)

# Scenario A-2
total.iter <- 2

n <- n.origin <- 15000
psi.origin <- psi <- 10
threshold <- 0.95
p <- 5

no.theta <- 1

C <- diag(p)

raw_f2_nl <- function(x) {1.5 * sin(2 * pi * x^2)*x^3}
raw_f3_nl <- function(x) {-1.5 * cos(3 * pi * x^2)*x^2}

make.nl <- function(x, raw_y) {
  fit <- lm(raw_y ~ x)
  
  return(list(
    nl = residuals(fit), 
    slope = coef(fit)[["x"]],
    intercept = coef(fit)[["(Intercept)"]]
  ))
}

theta.origin <- c(1, -0.3, 0, 0.5, 0, 0)
psi <- psi -2

model.stan <- "// Stan model for BLAST Student-t Samples
data {
  int <lower=1> n; // Sample size
  int <lower=1> p; // regression coefficient size
  int <lower=1> psi; // splines coefficient size
  real <lower=0> u; // large threshold value
  matrix[n, p] bsLinear; // fwi dataset
  matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
  matrix[n, p] xholderLinear; // fwi dataset
  matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
  array[n] real <lower=1> y; // extreme response
  real <lower=0> atau;
  vector[p] X_means;
  vector[p] X_sd;
  array[n] real <lower=0> trueAlpha;
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
      gamma[j] = gamma_raw[j] * sqrt(tau[j]);
      gsmooth[,j] = bsLinear[,j] * theta[j+1] + bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    
    for (i in 1:n){
      alpha[i] = exp(theta[1] + sum(gsmooth[i,]));
    };
  }
}

model {
  // likelihood
  for (i in 1:n){
    target += student_t_lpdf(y[i] | alpha[i], 0, 1);
    target += -1*log(1-student_t_cdf(u, alpha[i], 0, 1));
  }
  target += normal_lpdf(theta[1] | 0, 100);

  for (j in 1:p){
    target += gamma_lpdf(lambda1[j] | 2, 1); 
    target += gamma_lpdf(lambda2[j] | 0.01, 0.01);    
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
  vector[(p+1)] theta_origin;
  vector[n] se;
  real <lower=0> mse;

  for (j in 1:p){
    theta_origin[j+1] = theta[j+1] / X_sd[j];
  }
  theta_origin[1] = theta[1] - dot_product(X_means, theta_origin[2:(p+1)]);
  for (j in 1:p){
      gridgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
      gridgl[,j] = xholderLinear[,j] * theta_origin[j+1];
      gridgsmooth[,j] = gridgl[,j] + gridgnl[,j];
  };

  for (i in 1:n){
      gridalpha[i] = exp(theta_origin[1] + sum(gridgsmooth[i,]));
      se[i] = pow((gridalpha[i]-trueAlpha[i]), 2);
  };
  mse = mean(se);
}
"

gridgnl.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
gridgl.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
gridgsmooth.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
alpha.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
mise.container <- c()
qqplot.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))

for(iter in 1:total.iter){
  n <- n.origin
  x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
  f2.nl <- raw_f2_nl(x.origin[,2])
  f3.nl <- raw_f3_nl(x.origin[,3])
  f1.l <- theta.origin[2]*x.origin[,1]
  f3.l <- theta.origin[4]*x.origin[,3]

  eta_lin <-  f1.l + f3.l
  eta_nonlin <- f2.nl + f3.nl
  eta <- theta.origin[1] + eta_lin + eta_nonlin
  alp.origin <- exp(eta)
  
  y.origin <- rt(n, df = alp.origin)
  # for(i in 1:n){
  #   y.origin[i] <- rt(1, df = alp.origin[i])
  # }
  u <- quantile(y.origin, threshold)
  excess.index <- which(y.origin>u)
  x.origin <- x.origin[excess.index,]
  y.origin <- y.origin[excess.index]
  n <- length(y.origin)

  f2.hidden <- make.nl(x.origin[,2], raw_f2_nl(x.origin[,2]))
  f3.hidden <- make.nl(x.origin[,3], raw_f3_nl(x.origin[,3]))
  theta.adjusted <- c(theta.origin[1] + f2.hidden$intercept + f3.hidden$intercept,
                      theta.origin[2],
                      theta.origin[3] + f2.hidden$slope,
                      theta.origin[4] + f3.hidden$slope,
                      theta.origin[5],
                      theta.origin[6])
  newx <- seq(0,1,length.out = n)
  g2.nl <- raw_f2_nl(newx) - (f2.hidden$intercept + f2.hidden$slope*newx)
  g3.nl <- raw_f3_nl(newx) - (f2.hidden$intercept + f2.hidden$slope*newx)
  g1.l <- theta.adjusted[2]*newx
  g2.l <- theta.adjusted[3]*newx
  g3.l <- theta.adjusted[4]*newx
  g2 <- g2.l + g2.nl
  g3 <- g3.l + g3.nl
  eta.g <- rep(theta.adjusted[1], n) +g1.l + g2 + g3
  alp.new <- exp(eta.g)
  grid.zero <- rep(0, n)
  g.new <- c(g1.l, g2, g3, grid.zero, grid.zero)
  l.new <- c(g1.l, g2.l, g3.l, grid.zero, grid.zero)
  nl.new <- c(grid.zero, g2.nl, g3.nl, grid.zero, grid.zero)

  bs.linear <- model.matrix(~ ., data = data.frame(x.origin))
  
  group.map <- c()
  Z.list <- list()        # Stores the final non-linear design matrices
  scale_stats_list <- list() 
  projection_coefs_list <- list() #
  spec_decomp_list <- list() # Store eigen-decomp info for prediction
  sm_spec_list <- list()     # Store smooth objects
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
    tol <- max_lambda * 1e-6  # Relative threshold (Robust)
    
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

  xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))
  colnames(xholder) <- covariates
  grid_Z_list <- list()

  for (i in seq_along(covariates)) {
    grid_df  <- data.frame(x_vec = newx)
    X_lin_grid <- model.matrix(~ x_vec, data = grid_df)
    X_raw_grid <- PredictMat(sm_spec_list[[i]], grid_df)
    
    decomp <- spec_decomp_list[[i]]
    Z_spectral_grid <- X_raw_grid %*% decomp$U_pen %*% decomp$Lambda_sqrt_inv
    Z_orth_grid <- Z_spectral_grid - X_lin_grid %*% projection_coefs_list[[i]]
    Z_final_grid <- Z_orth_grid[, keep_cols_list[[i]], drop = FALSE]
    Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = scale_stats_list[[i]])
    
    grid_Z_list[[i]] <- Z_final_grid
  }

  xholder.linear <- model.matrix(~ ., data = data.frame(xholder))
  xholder.nonlinear <- do.call(cbind, grid_Z_list)

  X_means <- colMeans(bs.linear[,c(2:(p+1))])
  X_sd   <- apply(bs.linear[,c(2:(p+1))], 2, sd)
  bs.linear[,c(2:(p+1))] <- scale(bs.linear[,c(2:(p+1))], center = X_means, scale = X_sd)

  
  data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, X_sd=X_sd,
                    atau = ((psi+1)/2), X_means = X_means, trueAlpha = alp.new,
                    bsLinear = bs.linear[,c(2:(p+1))], bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear[,c(2:(p+1))], xholderNonlinear = xholder.nonlinear)

  init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p,psi)),
                          theta = rep(-0.1, (p+1)), tau = rep(0.1, p), # rho = 1, 
                          lambda1 = rep(0.1, p), lambda2 = rep(1, p)),
                    list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p,psi)),
                          theta = rep(0.05, (p+1)), tau = rep(0.5, p), #rho = 1,
                          lambda1 = rep(0.5, p), lambda2 = rep(0.2, p)),
                    list(gamma_raw = array(rep(0.1, (psi*p)), dim=c(p,psi)),
                          theta = rep(0.05, (p+1)), tau = rep(0.2, p), #rho = 1,
                          lambda1 = rep(0.2, p), lambda2 = rep(0.5, p)))
    
  fit1 <- stan(
      model_code = model.stan,
      data = data.stan,    # named list of data
      init = init.alpha,      # initial value
      chains = 3,             # number of Markov chains
      # warmup = 1000,          # number of warmup iterations per chain
      iter = 2000,            # total number of iterations per chain
      cores = parallel::detectCores(), # number of cores (could use one per chain)
      refresh = 1000             # no progress shown
  )
  # posterior <- extract(fit1)
  gridgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.5))$summary
  gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.5))$summary
  gridgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.5))$summary
  newalpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.5))$summary
  se.samples <- summary(fit1, par=c("se"), probs = c(0.5))$summary 

  alpha.container[,iter] <- newalpha.samples[,1]

  gridgsmooth.container[,iter] <- as.vector(matrix(gridgsmooth.samples[,4], nrow = n, byrow=TRUE))
  gridgl.container[,iter] <- as.vector(matrix(gridgl.samples[,4], nrow = n, byrow=TRUE))
  gridgnl.container[,iter] <- as.vector(matrix(gridgnl.samples[,4], nrow = n, byrow=TRUE))
  newx <- seq(0, 1, length.out = n)
  mise.container[iter] <- auc(newx, se.samples[,1], type="spline")

  mcmc.alpha <- rstan::extract(fit1)$alpha
  r <- matrix(, nrow = n, ncol = 30)
  T <- 30
  for(i in 1:n){
      for(t in 1:T){
          r[i, t] <- qnorm(pPareto(y.origin[i], u, alpha = mcmc.alpha[round(runif(1,1, dim(mcmc.alpha)[1])),i]))
      }
  }
  lgrid <- n
  grid <- qnorm(ppoints(lgrid))
  traj <- matrix(NA, nrow = T, ncol = lgrid)
  for (t in 1:T){
      traj[t, ] <- quantile(r[, t], ppoints(lgrid), type = 2)
  }
  qqplot.container[iter] <- apply(traj, 2, mean)#quantile, prob = 0.5)
}

alpha.container$x <- newx
alpha.container$true <- alp.new
alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
alpha.container <- as.data.frame(alpha.container)

# load(paste0("./simulation/results/MC-Scenario_A/2026-01-24_",total.iter,"_MC_scC_",n.origin,".Rdata"))

plt <- ggplot(data = alpha.container, aes(x = x)) + xlab(expression(c)) + labs(col = "") + ylab("")
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
} else{
  for(i in 50:100){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}
print(plt + ylim(0, 20) +
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.8) +
        scale_fill_manual(values=c("steelblue"), name = "") +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 40) + 
        theme(legend.position = "none",
                strip.text = element_blank(),
                axis.text.y = element_blank(),
                axis.text = element_text(size = 30)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_alpha_scA_",n.origin,".pdf"), width=9.5, height = 7.78)
# n<- 750
# newx <- seq(0, 1, length.out = n)

gridgsmooth.container$x <- newx
gridgsmooth.container$true <- g.new
gridgsmooth.container$mean <- rowMeans(gridgsmooth.container[,1:total.iter])
gridgsmooth.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]"))
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
print(plt + geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + ylim(-2.3, 2.3) +
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        # geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5) + 
        facet_grid(covariate ~ ., scales = "free_x", switch = "y",
                    labeller = label_parsed) +        
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
              plot.margin = margin(0,0,0,-20),
              axis.title.x = element_text(size = 45),                
              axis.text = element_text(size = 30)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_smooth_scC_",n.origin,".pdf"), width=11, height = 15)


gridgl.container$x <- newx
gridgl.container$true <- as.vector(l.new)
gridgl.container$mean <- rowMeans(gridgl.container[,1:total.iter])
gridgl.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]", "g[6]"))
gridgl.container <- as.data.frame(gridgl.container)

plt <- ggplot(data = gridgl.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(x))
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(gridgl.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}else{
  for(i in 50:100){
    plt <- plt + geom_line(aes(y = .data[[names(gridgl.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}
print(plt + geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5) + 
        ylim(-2.3, 2.3) +
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

# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_linear_scC",n.origin,".pdf"), width=11, height = 7.78)                

gridgnl.container$x <- newx
gridgnl.container$true <- as.vector(nl.new)
gridgnl.container$mean <- rowMeans(gridgnl.container[,1:total.iter])
gridgnl.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]", "g[6]"))
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
        # geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5) + 
        ylim(-2.3, 2.3) +
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

qqplot.container$grid <- grid
qqplot.container$mean <- rowMeans(qqplot.container[,1:total.iter])
plt <- ggplot(data = qqplot.container, aes(x = grid))
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(qqplot.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}else{
  for(i in 50:100){
    plt <- plt + geom_line(aes(y = .data[[names(qqplot.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}
print(plt + 
        geom_line(aes(y = mean), colour = "steelblue", linewidth = 1.5, linetype = 2) + 
        labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
        theme_minimal(base_size = 30) +
        theme(text = element_text(size = 20)) +
        coord_fixed(xlim = c(-2, 2),  
                    ylim = c(-2, 2)))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_qqplot_scC_",n.origin,".pdf"), width=10, height = 7.78)

# save(alpha.container, gridgnl.container, gridgl.container, gridgsmooth.container, mise.container, qqplot.container, file = paste0(Sys.Date(),"_",total.iter,"_MC_scC_",n.origin,".Rdata"))

mean(mise.container)
