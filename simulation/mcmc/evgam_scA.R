library(mgcv)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(evgam)
library(gridExtra)

#Scenario 1
# set.seed(10)
set.seed(806)

n <- 15000
psi <- 10
threshold <- 0.95
p <- 5
no.theta <- 1
simul.no <- 50

# Function to generate Gaussian copula
C <- diag(p)
x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))

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

f2.nl <- raw_f2_nl(x.origin[,2])
f3.nl <- raw_f3_nl(x.origin[,3])
theta.origin <- c(1, 0, 0.8, -0.8, 0, 0) 
f2.l <- theta.origin[3]*x.origin[,2]
f3.l <- theta.origin[4]*x.origin[,3]
eta_lin <-  f2.l + f3.l
eta_nonlin <- f2.nl + f3.nl
psi <- psi - 2

alp.origin <- exp(rep(theta.origin[1],n) + x.origin%*%theta.origin[-1] + raw_f2_nl(x.origin[,2]) + raw_f3_nl(x.origin[,3]))
y.pre <- y.origin <- rPareto(n, rep(1,n), alpha = alp.origin)

u <- quantile(y.origin, threshold)
excess.index <- which(y.origin>u)
x.origin <- x.origin[excess.index,]
y.origin <- y.origin[excess.index]
n <- length(y.origin)
newx <- seq(0,1,length.out = n)


f2.hidden <- make.nl(x.origin[,2], raw_f2_nl(x.origin[,2]))
f3.hidden <- make.nl(x.origin[,3], raw_f3_nl(x.origin[,3]))
theta.adjusted <- c(theta.origin[1] + f2.hidden$intercept + f3.hidden$intercept,
                    theta.origin[2],
                    theta.origin[3] + f2.hidden$slope,
                    theta.origin[4] + f3.hidden$slope,
                    theta.origin[5],
                    theta.origin[6])
g2.nl <- raw_f2_nl(newx) - (f2.hidden$intercept + f2.hidden$slope*newx)
g3.nl <- raw_f3_nl(newx) - (f2.hidden$intercept + f2.hidden$slope*newx)
g2.l <- theta.adjusted[3]*newx
g3.l <- theta.adjusted[4]*newx
g2 <- g2.l + g2.nl
g3 <- g3.l + g3.nl
eta.g <- rep(theta.adjusted[1], n) + g2 + g3
alp.new <- exp(eta.g)
grid.zero <- rep(0, n)
g.new <- c(grid.zero, g2, g3, grid.zero, grid.zero)
l.new <- c(grid.zero, g2.l, g3.l, grid.zero, grid.zero)
nl.new <- c(grid.zero, g2.nl, g3.nl, grid.zero, grid.zero)


bs.linear <- model.matrix(~ ., data = data.frame(x.origin))

# Initialize storage
group.map <- c()
Z.list <- list()        # Stores the final non-linear design matrices
# bs.linear <- NULL       # Stores the strictly linear matrices
scale_stats_list <- list() 
projection_coefs_list <- list() #
spec_decomp_list <- list() # Store eigen-decomp info for prediction
qr_list <- list()          # Store QR info for prediction
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

xholder.nonlinear <- do.call(cbind, grid_Z_list)
xholder.linear <- model.matrix(~ ., data = data.frame(xholder))

X_means <- colMeans(bs.linear[,c(2:(p+1))])
X_sd   <- apply(bs.linear[,c(2:(p+1))], 2, sd)
bs.linear[,c(2:(p+1))] <- scale(bs.linear[,c(2:(p+1))], center = X_means, scale = X_sd)


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
      target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += normal_lpdf(theta[1] | 0, 100);
    for (j in 1:p){
      target += gamma_lpdf(lambda1[j] | 2, 1); 
      target += gamma_lpdf(lambda2[j] | 1e-2, 1e-2);
      target += double_exponential_lpdf(theta[(j+1)] | 0, 1/(lambda1[j]));
      target += gamma_lpdf(tau[j] | atau, square(lambda2[j])*0.5);
      target += std_normal_lpdf(gamma_raw[j]); // target += multi_normal_lpdf(gamma_raw[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)))
    }
}

generated quantities {
    // Used in Posterior predictive check    
    vector[n] log_lik;
    array[n] real <lower=0> gridalpha; // new tail index
    matrix[n, p] gridgnl; // nonlinear component
    matrix[n, p] gridgl; // linear component
    matrix[n, p] gridgsmooth; // linear component
    vector[(p+1)] theta_origin;

    for (j in 1:p){
      theta_origin[j+1] = theta[j+1] / X_sd[j];
    }
    theta_origin[1] = theta[1] - dot_product(X_means, theta_origin[2:(p+1)]);
    for(i in 1:n){
      log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
    for (j in 1:p){
        gridgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        gridgl[,j] = xholderLinear[,j] * theta_origin[j+1];
        gridgsmooth[,j] = gridgl[,j] + gridgnl[,j];
    };

    for (i in 1:n){
        gridalpha[i] = exp(theta_origin[1] + sum(gridgsmooth[i,]));
    };
}
"

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                  atau = ((psi+1)/2), X_means = X_means, X_sd=X_sd,
                  bsLinear = bs.linear[,c(2:(p+1))], bsNonlinear = bs.nonlinear,
                  xholderLinear = xholder.linear[,c(2:(p+1))], xholderNonlinear = xholder.nonlinear)

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
  model_code = model.stan,  # Stan program
  data = data.stan,    # named list of data
  init = init.alpha,      # initial value
  chains = 3,             # number of Markov chains
  # warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = parallel::detectCores(), # number of cores (could use one per chain)
  refresh = 1000             # no progress shown
))

posterior <- extract(fit1)

theta.samples <- summary(fit1, par=c("theta_origin"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c('alpha'), probs = c(0.05, 0.5, 0.95))$summary
newgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
newgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
newalpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary

g.smooth.mean <- as.vector(matrix(newgsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(newgsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(newgsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(newgsmooth.samples[,6], nrow = n, byrow=TRUE))

simul.evgam <- data.frame(y = y.origin-u, x.origin)
gam.scale <- list(y ~ s(X1, bs = "tp", k = 10) + 
                      s(X2, bs = "tp", k = 10) + 
                      s(X3, bs = "tp", k = 10) + 
                      s(X4, bs = "tp", k = 10) + 
                      s(X5, bs = "tp", k = 10),
                    ~ s(X1, bs = "tp", k = 10) + 
                      s(X2, bs = "tp", k = 10) + 
                      s(X3, bs = "tp", k = 10) + 
                      s(X4, bs = "tp", k = 10) + 
                      s(X5, bs = "tp", k = 10))
evgam.fit.scale <- evgam::evgam(gam.scale, data = simul.evgam, family = "gpd", outer = "Newton")
plot.data <- as.data.frame(evgam.fit.scale$plotdata)
xi.pred.scale <-predict(evgam.fit.scale, newdata = data.frame(xholder), type="response")$shape
alpha.pred.scale <- 1/xi.pred.scale

xholder.basis.scale <- predict(evgam.fit.scale, newdata = data.frame(xholder), type= "lpmatrix")$shape
xi.coef.scale <- tail(evgam.fit.scale$coefficients, (psi-1)*p)
gamma.xi.scale <- matrix(xi.coef.scale, ncol = p)
alpha.nonlinear.scale <- xi.nonlinear.scale <- matrix(, nrow = n, ncol = p)
bs.nonlinear.scale <- xholder.basis.scale[,c(2:((psi-1)*p+1))]
for(j in 1:p){
  xi.nonlinear.scale[,j] <- bs.nonlinear.scale[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.scale[,j]
  alpha.nonlinear.scale[,j] <- 1/(xi.nonlinear.scale[,j])
}

gam.1 <- list(y ~ 1,
                ~ s(X1, bs = "tp", k = 10) + 
                  s(X2, bs = "tp", k = 10) + 
                  s(X3, bs = "tp", k = 10) + 
                  s(X4, bs = "tp", k = 10) + 
                  s(X5, bs = "tp", k = 10))
evgam.fit.1 <- evgam::evgam(gam.1, data = simul.evgam, family = "gpd")

# save(evgam.fit.scale, evgam.fit.1, file= "./simulation/results/evgam_scA")

plot.data <- as.data.frame(evgam.fit.1$plotdata)
xi.pred.1 <-predict(evgam.fit.1, newdata = data.frame(xholder), type="response")$shape
alpha.pred.1 <- 1/xi.pred.1

xholder.basis.1 <- predict(evgam.fit.1, newdata = data.frame(xholder), type= "lpmatrix")$shape
xi.coef.1 <- tail(evgam.fit.1$coefficients, (psi-1)*p)
gamma.xi.1 <- matrix(xi.coef.1, ncol = p)
alpha.nonlinear.1 <- xi.nonlinear.1 <- matrix(, nrow = n, ncol = p)
bs.nonlinear.1 <- xholder.basis.1[,c(2:((psi-1)*p+1))]
for(j in 1:p){
  xi.nonlinear.1[,j] <- bs.nonlinear.1[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.1[,j]
  alpha.nonlinear.1[,j] <- 1/(xi.nonlinear.1[,j])
}


alpha.smooth <- data.frame("x" = as.vector(xholder),
                          "true" = as.vector(g.new),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "evgam.scale" = as.vector(alpha.nonlinear.scale),
                          "evgam.1" = as.vector(alpha.nonlinear.1),
                          "covariates" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))


grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(alpha.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=true, colour = "True"), linewidth=1, linetype=2) + 
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1) + 
                  geom_line(aes(y=evgam.1), colour = "orange", linewidth=1) + 
                  ylab("") + xlab(xlab(expression(c))) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue", "red")) +
                  ylim(-2.3, 2.3) +
                  theme_minimal(base_size = 10) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 15),
                          axis.title.x = element_text(size = 15))
  grid.plts[[i]] <- grid.plt
}

grid.arrange(grobs = grid.plts, ncol = 3, nrow = 2)

xi.smooth <- data.frame("x" = as.vector(xholder),
                          "true" = as.vector(1/exp(g.new)),
                          "post.mean" = as.vector(1/exp(g.smooth.mean)),
                          "q1" = as.vector(1/exp(g.smooth.q1)),
                          "q2" = as.vector(1/exp(g.smooth.q2)),
                          "q3" = as.vector(1/exp(g.smooth.q3)),
                          "evgam.scale" = as.vector(xi.nonlinear.scale),
                          "evgam.1" = as.vector(xi.nonlinear.1),
                          "covariates" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(xi.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=true, colour = "True"), linewidth=1, linetype = 2) + 
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1) + 
                  geom_line(aes(y=evgam.1), colour = "orange", linewidth=1) + 
                  ylab("") + xlab(expression(c)) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue", "red")) +
                  # ylim(-2.8, 2.8) + xlim(0,1) +
                  theme_minimal(base_size = 10) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 15),
                          axis.title.x = element_text(size = 15))
  grid.plts[[i]] <- grid.plt
}

grid.arrange(grobs = grid.plts, ncol = 3, nrow = 2)

alpha.scenario <- data.frame("x" = newx,
                            "constant" = newx,
                            "true" = (alp.new),
                            "post.mean" = (newalpha.samples[,1]),
                            "post.median" = (newalpha.samples[,5]),
                            "evgam.scale" = 1/xi.pred.scale,
                            "evgam.1" = 1/xi.pred.1,
                            "q1" = (newalpha.samples[,4]),
                            "q3" = (newalpha.samples[,6]))

ggplot(alpha.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,ldots,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y = true, col = "True"), linewidth = 2, linetype=2) + 
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1.5) +
  geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1.5, linetype=2) +
  geom_line(aes(y=evgam.1), colour = "orange", linewidth=1.5, linetype=2) +
  scale_color_manual(values=c("steelblue", "red")) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  theme_minimal(base_size = 30) + ylim(0,20)+
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_evgam_scA.pdf"), width=10, height = 7.78)

xi.scenario <- data.frame("x" = newx,
                            "constant" = newx,
                            "true" = 1/(alp.new),
                            "post.mean" = 1/(newalpha.samples[,1]),
                            "post.median" = 1/(newalpha.samples[,5]),
                            "evgam.scale" = xi.pred.scale,
                            "evgam.1" = xi.pred.1,
                            "q1" = 1/(newalpha.samples[,4]),
                            "q3" = 1/(newalpha.samples[,6]))

ggplot(xi.scenario, aes(x=x)) + 
  ylab(expression(xi(c,ldots,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y = true, col = "True"), linewidth = 2, linetype=2) + 
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1.5) +
  geom_line(aes(y=evgam.scale), colour = "orange", linewidth=1.5, linetype=2) +
  geom_line(aes(y=evgam.1), colour = "purple", linewidth=1.5, linetype=2) +
  scale_color_manual(values=c("steelblue", "red")) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  theme_minimal(base_size = 30) + ylim(-1.5, 3.1)+
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))



# orderred <- rev(sort(y.pre)[14325:length(y.pre)])

orderred <- sort(y.pre, decreasing = TRUE)[1:(n+1)]
n.hill <- length(orderred)
k <- 1:n.hill
loggs <- logb(orderred/u)
avesumlog <- cumsum(loggs)/k
xihat <- c(NA, (avesumlog)[2:n.hill])
ses.xi <- xihat * sqrt(k)
xx <- seq(from = n.hill, to = 2)
y.xi <- xihat[xx]
qq <- qnorm(1-(1-threshold)/2)
uu.xi <- y.xi + ses.xi[xx] * qq
ll.xi <- y.xi - ses.xi[xx] * qq

alphahat <- 1/xihat
ses.alpha <- alphahat / sqrt(k)
y.alpha <- alphahat[xx]
uu.alpha <- y.alpha + ses.alpha[xx] * qq
ll.alpha <- y.alpha - ses.alpha[xx] * qq


gamma.alpha <- (xx * y.alpha) / (xx - 1)
gamma.q1.alpha <- 1 / qgamma(0.975, shape = xx, rate = xx * y.alpha)
gamma.q3.alpha <- 1 / qgamma(0.025, shape = xx, rate = xx * y.alpha)
gamma.xi <- (xx * y.xi) / (xx - 1)
gamma.q1.xi <- 1 / qgamma(0.975, shape = xx, rate = xx * y.xi)
gamma.q3.xi <- 1 / qgamma(0.025, shape = xx, rate = xx * y.xi)

data.hill.alpha <- data.frame("k" = c(1:n),
                        "u" = sort(uu.alpha),
                        "l" = sort(ll.alpha),
                        "alpha" = sort(y.alpha),
                        "order" = xx,
                        "blast.mean" = sort(alpha.samples[,1]),
                        "blast.q1" = sort(alpha.samples[,4]),
                        "blast.q3" = sort(alpha.samples[,6]),
                        "gamma.mean" = sort(gamma.alpha),
                        "gamma.q1" = sort(gamma.q1.alpha),
                        "gamma.q3" = sort(gamma.q3.alpha),
                        "evgam.1" = sort(1/predict(evgam.fit.1, type="response")$shape),
                        "evgam.scale" = sort(1/predict(evgam.fit.scale, type="response")$shape))
ggplot(data = data.hill.alpha) + 
  geom_ribbon(aes(x = order, ymin = l, ymax = u),
              alpha = 0.2, linetype = "dashed", fill = "black") + 
  geom_line(aes(x = order, y = alpha), linewidth = 1.2, colour = "black") +
  geom_line(aes(x = order, y = blast.mean), linewidth = 1.2, colour = "steelblue") + 
  geom_ribbon(aes(x = order, ymin = blast.q1, ymax = blast.q3),
              alpha = 0.2, linetype = "dashed", fill = "steelblue") + 
  geom_line(aes(x=order, y= gamma.mean), linewidth = 1.2, colour = "brown") + 
  geom_ribbon(aes(x = order, ymin = gamma.q1, ymax = gamma.q3),
              alpha = 0.2, linetype = "dashed", fill = "brown") + 
  geom_line(aes(x=order, y= evgam.1), linewidth = 1.2, linetype=3, colour = "orange") +
  geom_line(aes(x=order, y= evgam.scale), linewidth = 1.2, linetype=3, colour = "purple") +               
  labs(x = "Order Statistics", y = expression(alpha)) +
  theme_minimal(base_size = 30) +
  theme(text = element_text(size = 30), 
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")

data.hill.xi <- data.frame("k" = c(1:n),
                        "u" = sort(uu.xi),
                        "l" = sort(ll.xi),
                        "alpha" = sort(y.xi),
                        "order" = xx,
                        "evgam.1" = sort(predict(evgam.fit.1, type="response")$shape),
                        "evgam.scale" = sort(predict(evgam.fit.scale, type="response")$shape),
                        "blast.mean" = sort(1/alpha.samples[,1]),
                        "blast.q1" = sort(1/alpha.samples[,4]),
                        "blast.q3" = sort(1/alpha.samples[,6]),
                        "gamma.mean" = sort(gamma.xi),
                        "gamma.q1" = sort(gamma.q1.xi), 
                        "gamma.q3" = sort(gamma.q3.xi))
ggplot(data = data.hill.xi) + 
  # geom_ribbon(aes(x = order, ymin = l, ymax = u),
              # alpha = 0.2, linetype = "dashed", fill = "black") + 
  # geom_line(aes(x = order, y = alpha), linewidth = 1.2, colour = "black") +
  geom_line(aes(x = order, y = blast.mean), linewidth = 1.2, colour = "steelblue") + 
  geom_ribbon(aes(x = order, ymin = blast.q1, ymax = blast.q3),
              alpha = 0.2, linetype = "dashed", fill = "steelblue") + 
  # geom_line(aes(x=order, y= gamma.mean), linewidth = 1.2, colour = "brown") + 
  # geom_ribbon(aes(x = order, ymin = gamma.q1, ymax = gamma.q3),
              # alpha = 0.2, linetype = "dashed", fill = "brown") + 
  geom_line(aes(x=order, y= evgam.1), linewidth = 1.2, linetype=3, colour = "orange") +
  geom_line(aes(x=order, y= evgam.scale), linewidth = 1.2, linetype=3, colour = "purple") +               
  labs(x = "Order Statistics", y = expression(xi)) + #ylim(0, 32) +
  theme_minimal(base_size = 30) +
  theme(text = element_text(size = 30), 
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")




hill.est <- function(data, k) {
  n <- length(data)
  if (k <= 0 | k >= n) stop("k must be between 1 and n-1")
  sorted.data <- sort(data, decreasing = TRUE)
  
  X.k <- sorted.data[k]
  X.top.k <- sorted.data[1:k]
  H.k <- (1/k) * sum(log(X.top.k) - log(X.k))
  
  return(H.k)
}


bayes.hill.est <- function(y, k) {
  # Inverse gamma = 1 / Gamma, so sample from Gamma and invert
  # 1 / rgamma(length(y), shape = k, rate = k * hill.est(y, k))
  k * hill.est(y,k) / (k-1)
}

bayes.hill.est(y.origin, 100)

