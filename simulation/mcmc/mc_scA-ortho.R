library(mgcv)
library(Pareto)
suppressMessages(library(tidyverse))
library(rstan)
library(MESS)

# Scenario A-2
total.iter <- 160

n <- n.origin <- 15000
psi <- 10
threshold <- 0.95
p <- 5

no.theta <- 1

C <- diag(p)
## Generate sample
f <- function(x){
  0.25*x[,2] +-cos(pi * x[,2])*x[,2] + 
  0.25*x[,3] + -0.5*sin(2*pi *x[,3])
}
theta.origin <- c(0.5, 0, 0.25, 0.25, 0, 0)
psi <- psi -2

model.stan <- "// Stan model for BLAST Pareto Samples
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
  array[n] real trueAlpha;
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
  array[n] real <lower=0> gridalpha; // new tail index
  matrix[n, p] gridgnl; // nonlinear component
  matrix[n, p] gridgl; // linear component
  matrix[n, p] gridgsmooth; // linear component  
  real theta0;
  array[p] vector[psi] gamma;
  real <lower=0> mse; // mean squared error
  array[n] real <lower=0> se; // squared error

  {
    matrix[n, p] gsmooth; // linear component

    for (j in 1:p){
      gamma[j] = gamma_raw[j] * sqrt(tau[j]);
      gsmooth[,j] = bsLinear[,j] * theta[j+1] + bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
      gridgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
      gridgl[,j] = xholderLinear[,j] * theta[j+1];
      gridgsmooth[,j] = gridgl[,j] + gridgnl[,j];
    };
    theta0 = theta[1] - dot_product(X_means, theta[2:(p+1)]);
    for (i in 1:n){
      alpha[i] = exp(theta[1] + sum(gsmooth[i,]));
      gridalpha[i] = exp(theta[1] + sum(gridgsmooth[i,]));
      se[i] = pow((gridalpha[i]-trueAlpha[i]), 2);
    };
    mse = mean(se);
  }
}

model {
  // likelihood
  for (i in 1:n){
      target += pareto_lpdf(y[i] | u, alpha[i]);
  }
  target += normal_lpdf(theta[1] | 0, 1); //
  for (j in 1:p){
      target += gamma_lpdf(lambda1[j] | 1, 1); 
      target += gamma_lpdf(lambda2[j] | 0.01, 0.01);
      target += double_exponential_lpdf(theta[(j+1)] | 0, 1/lambda1[j]);
      target += gamma_lpdf(tau[j] | atau, square(lambda2[j])*0.5);
      target += multi_normal_lpdf(gamma_raw[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)));
  }
}
"

gridgnl.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
gridgl.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
gridgsmooth.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
alpha.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
# alpha.lower.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
# alpha.upper.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
mise.container <- c()
qqplot.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))

for(iter in 1:total.iter){
  n <- n.origin
  x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
  x.origin <- data.frame(x.origin)
  alp.origin <- exp(rep(theta.origin[1], n) + f(x.origin))
  y.origin <- rPareto(n, rep(1,n), alpha = alp.origin)
  u <- quantile(y.origin, threshold)
  excess.index <- which(y.origin>u)
  x.origin <- x.origin[excess.index,]
  y.origin <- y.origin[y.origin > u]
  n <- length(y.origin)
  covariates <- colnames(x.origin) 


  bs.linear <- model.matrix(~ ., data = x.origin)
  qr.x <- base::qr(bs.linear)
  q.x <- qr.Q(qr.x)
  group.map <- c()
  Z.list <- list()
  Gamma_list <- list()      
  keep_cols_list <- list()  
  sm_spec_list <- list() 
  scale_stats_list <- list()

  for (i in seq_along(covariates)) {
    var_name <- covariates[i]
    sm_expr <- call("s", as.symbol(var_name), bs = "tp", k = (psi+2))
    sm_spec <- eval(sm_expr)
    sm <- smoothCon(sm_spec, data = x.origin, knots = NULL)[[1]]
    Z_raw <- as.matrix(sm$X)
    Gamma <- t(q.x) %*% Z_raw 
    Z_orth <- Z_raw - q.x %*% Gamma
    keep_cols <- colSums(abs(Z_orth)) > 1e-9
    Z_clean <- Z_orth[, keep_cols, drop = FALSE]
    
    col_sds <- apply(Z_clean, 2, sd)
    col_sds[col_sds < 1e-9] <- 1 # Safety
    Z_clean <- scale(Z_clean, center = FALSE, scale = col_sds)
    Z.list[[i]] <- Z_clean
    group.map <- c(group.map, rep(i, ncol(Z_clean)))
    Gamma_linear <- backsolve(qr.R(qr.x), Gamma)
    
    Gamma_list[[i]] <- Gamma_linear
    keep_cols_list[[i]] <- keep_cols
    sm_spec_list[[i]] <- sm
    scale_stats_list[[i]] <- col_sds
  }
  bs.nonlinear <- do.call(cbind, Z.list)
  grid_seq <- seq(0, 1, length.out = n)
  grid_Z_list <- list() 

  xholder <- matrix(0, nrow = n, ncol = p)
  colnames(xholder) <- covariates

  for (i in 1:p) {
    var_name <- covariates[i] 
    xholder[,i] <- grid_seq
    grid_df <- data.frame(matrix(0, nrow = n, ncol = p))
    colnames(grid_df) <- covariates # Must match x.origin
    grid_df[[var_name]] <- grid_seq 
    X_grid <- model.matrix(~ ., data = grid_df)
    Z_grid_raw <- PredictMat(sm_spec_list[[i]], grid_df)
    Gamma_curr <- Gamma_list[[i]]
    Z_grid_orth <- Z_grid_raw - X_grid %*% Gamma_curr
    keep_cols <- keep_cols_list[[i]]
    Z_grid_final <- Z_grid_orth[, keep_cols, drop = FALSE]
    Z_grid_final <- scale(Z_grid_final, center = FALSE, scale = scale_stats_list[[i]])
    
    grid_Z_list[[i]] <- Z_grid_final
  }

  xholder.linear <- model.matrix(~ ., data = data.frame(xholder))
  xholder.nonlinear <- do.call(cbind, grid_Z_list)

  alp.origin <- exp(rep(theta.origin[1], n) + f(bs.linear[,c(2:(p+1))]))
  alp.new <- exp(rep(theta.origin[1], n) + f(xholder.linear[,c(2:(p+1))]))

  X_means <- colMeans(bs.linear[,c(2:(p+1))])
  bs.linear[,c(2:(p+1))] <- scale(bs.linear[,c(2:(p+1))], center = X_means, scale = FALSE)

  
  data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), X_means = X_means, trueAlpha = alp.new,
                    bsLinear = bs.linear[,c(2:(p+1))], bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear[,c(2:(p+1))], xholderNonlinear = xholder.nonlinear)

  init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p,psi)),
                          theta = rep(-0.1, (p+1)), tau = rep(0.1, p), # rho = 1, 
                          lambda1 = rep(0.1, p), lambda2 = rep(1, p)),
                    list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p,psi)),
                          theta = rep(0.05, (p+1)), tau = rep(2, p), #rho = 1,
                          lambda1 = rep(2, p), lambda2 = rep(5, p)),
                    list(gamma_raw = array(rep(-0.75, (psi*p)), dim=c(p,psi)),
                          theta = rep(0.01, (p+1)), tau = rep(1.5, p), # rho = 1,
                          lambda1 = rep(0.5, p), lambda2= rep(5, p)))
    
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
  posterior <- extract(fit1)
  gridgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.5))$summary
  gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.5))$summary
  gridgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.5))$summary
  newalpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.5))$summary
  se.samples <- summary(fit1, par=c("se"), probs = c(0.5))$summary 

  alpha.container[,iter] <- newalpha.samples[,1]

  gridgsmooth.container[,iter] <- as.vector(matrix(gridgsmooth.samples[,1], nrow = n, byrow=TRUE))
  gridgl.container[,iter] <- as.vector(matrix(gridgl.samples[,1], nrow = n, byrow=TRUE))
  gridgnl.container[,iter] <- as.vector(matrix(gridgnl.samples[,1], nrow = n, byrow=TRUE))
  newx <- seq(0, 1, length.out = n)
  mise.container[iter] <- auc(newx, se.samples[,1], type="spline")

  mcmc.alpha <- extract(fit1)$alpha
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

alpha.container$x <- seq(0,1, length.out = n)
alpha.container$true <- alp.new
# alpha.container <- cbind(alpha.container, t(apply(alpha.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
# colnames(alpha.container)[(dim(alpha.container)[2]-2):(dim(alpha.container)[2])] <- c("q1","q2","q3")
alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
alpha.container <- as.data.frame(alpha.container)


# load(paste0("./simulation/results/MC-Scenario_A/2024-05-01_",total.iter,"_MC_scA_",n.origin,".Rdata"))

plt <- ggplot(data = alpha.container, aes(x = x)) + xlab(expression(c)) + labs(col = "") + ylab(expression(alpha(c,...,c))) #+ ylab("")
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
} else{
  # for(i in 1:total.iter){
  for(i in 50:100){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}
print(plt + 
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.8) +
        scale_fill_manual(values=c("steelblue"), name = "") +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 40) + 
        theme(legend.position = "none",
                strip.text = element_blank(),
                axis.text = element_text(size = 30)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_alpha_scA_",n.origin,".pdf"), width=10, height = 7.78)
n<- 750
newx <- seq(0, 1, length.out = n)
f2.l <- function(x) {0.25*x}
f2.nl <- function(x) {-cos(pi * x)*x}
f2 <- f2.l(newx) + f2.nl(newx)
f3.l <- function(x) {0.25*x}
f3.nl <- function(x) {-0.5*sin(2*pi *x)}
f3 <- f3.l(newx) + f3.nl(newx)
grid.zero <- rep(0, n)
g.new <- c(grid.zero, f2, f3, grid.zero, grid.zero)
l.new <- c(grid.zero, f2.l(newx), f3.l(newx), grid.zero, grid.zero)
nl.new <- c(grid.zero, f2.nl(newx), f3.nl(newx), grid.zero, grid.zero)
gridgsmooth.container$x <- newx
gridgsmooth.container$true <- g.new
# gridgsmooth.container <- cbind(gridgsmooth.container, t(apply(gridgsmooth.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
# colnames(gridgsmooth.container)[(dim(gridgsmooth.container)[2]-2):(dim(gridgsmooth.container)[2])] <- c("q1","q2","q3")
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
print(plt + 
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5) + 
        # ylim(-1, 1) +
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
gridgl.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]", "g[6]"))
gridgl.container <- as.data.frame(gridgl.container)

plt <- ggplot(data = gridgl.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(x))
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(gridgl.container)[i]]]), alpha = 0.05, linewidth = 0.7)
    # plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
  }
}else{
  for(i in 1:50){
    plt <- plt + geom_line(aes(y = .data[[names(gridgl.container)[i]]]), alpha = 0.05, linewidth = 0.7)
  }
}
print(plt + 
        geom_line(aes(y=true, col = "True"), linewidth = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5, linetype = 2) + 
        # ylim(-0.23, 0.2) +
        facet_grid(covariate ~ ., scales = "free_x", switch = "y",
                    labeller = label_parsed) +  
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                plot.margin = margin(0,0,0,-20),
                strip.text = element_blank(),
                axis.text = element_text(size = 20)))

# # ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_linear_sc1-wi.pdf"), width=10, height = 7.78)                

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
print(plt + 
        geom_line(aes(y=true, col = "True"), linewidth = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5, linetype = 2) + 
        # ylim(-0.23, 0.2) +
        facet_grid(covariate ~ ., scales = "free_x", switch = "y",
                    labeller = label_parsed) +  
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                plot.margin = margin(0,0,0,-20),
                strip.text = element_blank(),
                axis.text = element_text(size = 20)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_nonlinear_sc1-wi.pdf"), width=10, height = 7.78)

qqplot.container$grid <- grid
qqplot.container$mean <- rowMeans(qqplot.container[,1:total.iter])
plt <- ggplot(data = qqplot.container, aes(x = grid))
for(i in 1:total.iter){
# for(i in 50:100){
  plt <- plt + geom_line(aes(y = .data[[names(qqplot.container)[i]]]), alpha = 0.05, linewidth = 0.7)
}
print(plt + 
        geom_line(aes(y = mean), colour = "steelblue", linewidth = 1.5, linetype = 2) + 
        labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
        theme_minimal(base_size = 30) +
        theme(text = element_text(size = 20)) +
        coord_fixed(xlim = c(-2, 2),  
                    ylim = c(-2, 2)))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_qqplot_scA_",n.origin,".pdf"), width=10, height = 7.78)

# save(alpha.container, gridgnl.container, gridgl.container, gridgsmooth.container, mise.container, qqplot.container, file = (paste0(Sys.Date(),"_",total.iter,"_MC_scA-2_",n.origin,".Rdata")))

mean(mise.container)
