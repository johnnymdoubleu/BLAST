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

# M <- cor(fwi.origin[,c(1:p)])
# corrplot(M, order = 'AOE', type = 'upper', tl.pos = 'tp')
# corrplot(M, add = TRUE, type = 'lower', method = 'number', order = 'AOE',
#          col = 'black', diag = FALSE, tl.pos = 'n', cl.pos = 'n')
bs.linear <- model.matrix(~ ., data = data.frame(fwi.scaled))
covariates <- colnames(fwi.scaled)

# xholder <- do.call(cbind, lapply(1:p, function(j) {seq(min(fwi.scaled[,j]),max(fwi.scaled[,j]), length.out = n)}))
newx <- seq(0, 1, length.out = n)
xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))

colnames(xholder) <- covariates
xholder.linear <- model.matrix(~ ., data = data.frame(xholder))
# fwi.grid <- as.data.frame(sapply(fwi.grid, FUN = range01))
xgrid.linear <- model.matrix(~ ., data = data.frame(fwi.grid))
xgrid.linear[,-1] <- sweep(xgrid.linear[, -1], 2, fwi.min, "-")
X_means <- colMeans(bs.linear[,-1])
X_sd   <- apply(bs.linear[,-1], 2, sd)
bs.linear[,-1] <- scale(bs.linear[,-1], center = X_means, scale = X_sd)

model.stan <- "// Stan model for BLAST Linear
data {
  int <lower=1> n; // Sample size
  int <lower=1> p; // regression coefficient size
  array[n] real <lower=0> u; // large threshold value
  matrix[n, p] bsLinear; // fwi dataset
  matrix[n, p] xholderLinear; 
  matrix[n, p] gridL; // fwi dataset
  array[n] real <lower=u> y; // extreme response
  vector[p] X_means;
  vector[p] X_sd;
  vector[p] X_minmax;
}

parameters {
  vector[(p+1)] theta; // linear predictor
  array[p] real <lower=0> lambda1; // lasso penalty //array[p] real <lower=0> 
}

transformed parameters {
  array[n] real <lower=0> alpha; // covariate-adjusted tail index
  {
    matrix[n, p] gsmooth; // linear component
    for (j in 1:p){
      gsmooth[,j] = bsLinear[,j] * theta[j+1];
    };
    
    for (i in 1:n){
      alpha[i] = exp(theta[1] + sum(gsmooth[i,]));
    };
  }
}

model {
  // likelihood
  target += pareto_lpdf(y | u, alpha);
  target += std_normal_lpdf(theta);
}

generated quantities {
  // Used in Posterior predictive check
  array[n] real <lower=0> gridalpha; // new tail index
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
    gridgl[,j] = xholderLinear[,j] * theta_origin[j+1];
    gridgsmooth[,j] = gridgl[,j];
    fwismooth[,j] = gridL[,j] * theta_fwi[j];
  };

  for (i in 1:n){
    gridalpha[i] = exp(theta_origin[1] + sum(gridgsmooth[i,]));
    log_lik[i] = pareto_lpdf(y[i] | u[i], alpha[i]);
  };
}
"


data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi,
                  bsLinear = bs.linear[,-1], gridL = xgrid.linear[,-1],
                  xholderLinear = xholder.linear[,-1], 
                  X_minmax = fwi.minmax, X_means = X_means, X_sd = X_sd)

init.alpha <- list(list(theta = rep(-0.1, (p+1))),
                   list(theta = rep(0.05, (p+1))),
                   list(theta = rep(0.1, (p+1))))

fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value 
    chains = 3,             # number of Markov chains
    iter = 8000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 1000           # no progress shown
)

posterior <- rstan::extract(fit1)
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))

theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
fwismooth.samples <- summary(fit1, par=c("fwismooth"), probs = c(0.05, 0.5, 0.95))$summary
origin.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary


MCMCvis::MCMCplot(fit1, params = 'theta')

g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))
g.linear.mean <- as.vector(matrix(gridgl.samples[,1], nrow = n, byrow=TRUE))
g.linear.q1 <- as.vector(matrix(gridgl.samples[,4], nrow = n, byrow=TRUE))
g.linear.q2 <- as.vector(matrix(gridgl.samples[,5], nrow = n, byrow=TRUE))
g.linear.q3 <- as.vector(matrix(gridgl.samples[,6], nrow = n, byrow=TRUE))
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


data.smooth <- data.frame("x" = as.vector(xholder),
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

fit.log.lik <- extract_log_lik(fit1)
loo(fit.log.lik, is_method = "sis", cores = 4)
