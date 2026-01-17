library(npreg)
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
library(splines)
library(quantreg)
library(qgam)
library(mgcViz)
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
# Y[Y==0] <- 1e-5
psi <- 30

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
    fwi.scaled[,i] <- fwi.index[,i] <- cov.long$measurement[missing.values]
}

fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)
fwi.origin <- fwi.scaled

fwi.origin <- data.frame(fwi.origin, time = c(1:length(Y)), BA=Y)
# fwi.origin <- data.frame(fwi.origin[which(Y>1),], BA=Y[Y>1])
# BA.shifted <- ifelse(fwi.origin$BA == 0, 1e-5, fwi.origin$BA)
# fwi.origin$log.BA <- log(fwi.origin$BA+1)
vars <- colnames(fwi.origin)[1:7]
seq_list <- lapply(vars, function(var) seq(min(fwi.origin[[var]]), max(fwi.origin[[var]]), length.out = length(Y)))
grid.df <- as.data.frame(setNames(seq_list, vars))


# 1. Get predictions at tau = 0.5 and tau = 0.95
# fit.50 <- quantreg::rq(BA ~ bs(DSR) + bs(FWI) + bs(BUI) + bs(ISI) + bs(FFMC) + bs(DMC) + bs(DC), tau = 0.5, data = fwi.origin)

# fit.95 <- quantreg::rq(BA ~ bs(DSR) + bs(FWI) + bs(BUI) + bs(ISI) + bs(FFMC) + bs(DMC) + bs(DC), tau = 0.95, data = fwi.origin)


# pred.50 <- predict(fit.50, newdata = fwi.origin)
# pred.95 <- predict(fit.95, newdata = fwi.origin)

# par(mfrow = c(1, 2))
# plot(pred.50, fwi.origin$BA,
#      main = "Predicted 50th Percentile vs Actual BA",
#      xlab = expression(hat(F)^(-1)),
#      ylab = expression(y[i]))
#     #  xlim = c(min_val_50, max_val_50),  # Equal x and y range
#     #  ylim = c(min_val_50, max_val_50))
# abline(0, 1, col = "red", lty = 2, lwd = 2)
# plot(pred.95, fwi.origin$BA,
#      main = "Predicted 95th Percentile vs Actual BA",
#      xlab = expression(hat(F)^(-1)),
#      ylab = expression(y))
# #     #  xlim = c(min_val_95, max_val_95),  # Equal x and y range
# #     #  ylim = c(min_val_95, max_val_95))
# abline(0, 1, col = "red", lty = 2, lwd = 2)
# par(mfrow = c(1, 1))


# # qu <- exp(predict(quant.fit))-1
# summary(quant.fit)
# residuals_qr <- residuals(quant.fit)
# fitted_vals <- fitted(quant.fit)

# par(mfrow=c(1,2))

# # Residuals vs. Fitted
# plot(fitted_vals, residuals_qr, main="Residuals vs Fitted", 
# 			xlab="Fitted values", ylab="Residuals")
# abline(h=0, col="red", lty=2)
# qqnorm(residuals_qr, main = "Q-Q plot of residuals")
# qqline(residuals_qr, col="red", lty=2)
# par(mfrow=c(1,1))
# predictions <- predict(quant.fit, newdata = fwi.origin)
# pred.grid <- predict(quant.fit, newdata = grid.df)
# coverage <- mean(fwi.origin$BA <= predictions)
# cat("Target coverage level:", 0.975, "\n")
# cat("Empirical coverage:", coverage, "\n")


# For multiple quantiles
# taus <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.98, 0.99)
# fits <- quantreg::rq(BA ~ bs(DSR) + bs(FWI) + bs(BUI) + bs(ISI) + bs(FFMC) + bs(DMC) + bs(DC), tau = taus, data = fwi.origin)
# plot(summary(fits))

# colors <- rainbow(length(taus))
# par(mfrow = c(3,3))
# main_vars <- colnames(fwi.origin)[1:7]
# for (j in 1:7) {
#   plot(fwi.origin[,j], fwi.origin$BA,
#        main = paste(main_vars[j]), xlab="",ylab="")
  
#   # Add quantile lines
#   for (i in seq_along(taus)) {
#     fit_tau <- quantreg::rq(BA ~ bs(DSR) + bs(FWI) + bs(BUI) + bs(ISI) + bs(FFMC) + bs(DMC) + bs(DC), tau = c(taus)[i], data = fwi.origin)
#     pred_tau <- predict(fit_tau, newdata = grid.df)
#     lines(grid.df[,j], pred_tau, col = colors[i], lwd = 1)
#   }
# }

# par(mfrow = c(3, 3))
# # taus <- seq(0.9, 1, 0.01)[-c(1, 11)]
# taus <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.98, 0.99)
# for (i in seq_along(taus)) {
#     tau_i <- taus[i]

#     fit_tau <- quantreg::rq(BA ~ bs(DSR) + bs(FWI) + bs(BUI) + bs(ISI) + bs(FFMC) + bs(DMC) + bs(DC), tau = tau_i, data = fwi.origin)
#     pred_tau <- predict(fit_tau, newdata = fwi.origin)

#     # Bin predictions
#     nbins <- 100
#     bins <- cut(pred_tau, 
#                 breaks = quantile(pred_tau, probs = seq(0, 1, by = 1/nbins)),
#                 include.lowest = TRUE)

#     # Calculate empirical coverage per bin
#     coverage_by_bin <- tapply((fwi.origin$BA <= pred_tau), bins, mean)
#     mean_pred_bin <- tapply(pred_tau, bins, mean)

#     # Plot calibration curve
#     plot(mean_pred_bin, coverage_by_bin, 
#             main = paste("tau =", tau_i),
#             xlab = "Mean Predicted Quantile in Bin",
#             ylab = "Empirical Coverage",
#             ylim = c(0, 1),
#             xlim = c(min(pred_tau), max(pred_tau)))

#     # Perfect calibration line
#     abline(0, 1, lty = 2, col = "red", lwd = 2)
#     # Target line
#     abline(h = tau_i, lty = 3, col = "black", lwd = 1)
# }



# results <- data.frame(
#     tau = numeric(),
#     target_coverage = numeric(),
#     empirical_coverage = numeric(),
#     coverage_error = numeric(),
#     n_u = numeric()
# )

# for (i in seq_along(taus)) {
#     tau_i <- taus[i]

#     # Fit model at this quantile
#     fit_tau <- quantreg::rq(BA ~ bs(DSR) + bs(FWI) + bs(BUI) + bs(ISI) + bs(FFMC) + bs(DMC) + bs(DC), tau = tau_i, data = fwi.origin)

#     # Get predictions
#     pred_tau <- predict(fit_tau, newdata = fwi.origin)
#     hist(pred_tau, main = paste("tau =", tau_i))
#     # Calculate empirical coverage: proportion of actual values <= predicted
#     empirical_cov <- mean(fwi.origin$BA <= pred_tau, na.rm = TRUE)

#     # How many observations fall below predicted quantile
#     n_below <- sum(fwi.origin$BA <= pred_tau, na.rm = TRUE)

#     results <- rbind(results, data.frame(
#         tau = tau_i,
#         target_coverage = tau_i,
#         empirical_coverage = empirical_cov,
#         coverage_error = abs(empirical_cov - tau_i),
#         n_u = length(Y) - n_below
#     ))
# }
# results
# par(mfrow = c(1, 1))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.df <- as.data.frame(sapply(fwi.origin[, c(1:8)], FUN = range01))
fwi.df$BA <- fwi.origin$BA
# quant.fit <- qgam(BA ~ s(DSR, k = 30) + s(FWI, k = 30) + s(BUI, k = 30) + s(ISI, k = 30) + s(FFMC, k = 30) + s(DMC, k = 30) + s(DC, k = 30), qu=0.975, data = fwi.df)
# quant.fit <- qgam(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC) + s(time), qu=0.99, data = fwi.df)

# quant.u <- predict(quant.fit)

# quant.fit <- qgamV(BA ~ s(DSR, k = 30) + s(FWI, k = 30) + s(BUI, k = 30) + s(ISI, k = 30) + s(FFMC, k = 30) + s(DMC, k = 30) + s(DC, k = 30), qu=0.975, data = fwi.origin)
taus <- c(0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
library(parallel)
n.cores <- length(taus)
cl <- makeCluster(n.cores)
# quant.fit <- mqgamV(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC) + s(time), qu=taus, data = fwi.df, aQgam = list(cluster = cl, multicore = TRUE, ncores = n.cores))
quant.fit <- qgam::mqgam(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC) + s(time), qu=taus, data = fwi.df, cluster = cl, multicore = TRUE, ncores = n.cores)
stopCluster(cl)
# save(quant.fit, fwi.df, fwi.origin, taus, file="./BLAST/application/mqgam_time.Rdata")
load("./BLAST/application/mqgam_time.Rdata")

print(plot(quant.fit, allTerms = TRUE), pages = 1)
for(i in 1:length(taus)){
    print(check1D(quant.fit[[i]], fwi.df[1:8]) + l_gridQCheck1D(qu = taus[i]), pages=1)
}
# print(check1D(quant.fit[[8]], fwi.df[1:8]) + l_gridQCheck1D(qu = 0.97), pages=1)
# print(check1D(quant.fit[[9]], fwi.df[1:8]) + l_gridQCheck1D(qu = 0.98), pages=1)
# print(check1D(quant.fit[[10]], fwi.df[1:8]) + l_gridQCheck1D(qu = 0.99), pages=1)
# qdo(quant.fit, predict, newdata = )

# plot(quant.fit)
# fwi.scaled <- fwi.origin[,c(1:7)]
y <- Y[Y>quant.u]
u <- quant.u[which(Y>quant.u)]
fwi.scaled <- as.data.frame(sapply(fwi.origin[which(Y>quant.u),c(1:7)], FUN = range01))
n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]
# save(quant.fit, quant.u, fwi.scaled, y, n,p, u, file="./BLAST/application/quant975_prep.Rdata")
save(fwi.scaled, fwi.origin, fwi.df, Y, n, psi, file = "./BLAST/application/wildfire_time_prep.Rdata")
# load("./BLAST/application/quant975_prep.Rdata")

# m.gam <- mqgam(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC), data = fwi.origin, qu = taus)
# m.gam <- mqgamV(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC), data = fwi.origin, qu = taus)
# save(m.gam, file ="./BLAST/application/wildfire_mgam2.Rdata")
# load("./BLAST/application/wildfire_mgam.Rdata")

# summary(m.gam[[7]])

no.theta <- 1 #represents the no. of linear predictors for each smooth functions
newx <- seq(0, 1, length.out=n)
xholder.linear <- xholder.nonlinear <- bs.linear <- bs.nonlinear <- matrix(,nrow=n, ncol=0)
xholder <- matrix(nrow=n, ncol=p)
end.holder <- basis.holder <- matrix(, nrow = 2, ncol = 0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                      matrix(c(which.min(fwi.scaled[,i]),
                              which.max(fwi.scaled[,i])), ncol=2))
}

for(i in 1:p){
  xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
  test.knot <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  knots <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  tps <- basis.tps(fwi.scaled[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  basis.holder <- cbind(basis.holder, 
          solve(matrix(c(tps[index.holder[i,1], no.theta+1],
                  tps[index.holder[i,1], no.theta+psi],
                  tps[index.holder[i,2], no.theta+1],
                  tps[index.holder[i,2], no.theta+psi]), 
                  nrow = 2, ncol = 2)))
  end.holder <- cbind(end.holder, 
                matrix(c(tps[index.holder[i,1], no.theta+1],
                  tps[index.holder[i,1], no.theta+psi],
                  tps[index.holder[i,2], no.theta+1],
                  tps[index.holder[i,2], no.theta+psi]), 
                  nrow = 2, ncol = 2))
  # mid.holder <- cbind(mid.holder,)                  
  bs.linear <- cbind(bs.linear, tps[,1:no.theta])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
}


model.stan <- "data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    array[n] real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    array[n] real <lower=1> y; // extreme response
    real <lower=0> atau;
    matrix[2, (2*p)] basisFL;
    array[(p*2)] int indexFL;
}
parameters {
    vector[(p+1)] theta; // linear predictor
    vector[(psi-2)] gammaTemp[p]; // constraint splines coefficient from 2 to psi-1
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    array[p] real <lower=0> tau1;
    array[p] real <lower=0> tau2;
}
transformed parameters {
    array[n] real <lower=0> alpha; // covariate-adjusted tail index
    vector[psi] gamma[p]; // splines coefficient 
    real <lower=0> lambda2o;

    lambda2o=lambda2*100;
    {
      vector[2] gammaFL[p];
      matrix[n, p] gsmooth; // linear component
      
      for(j in 1:p){
          gamma[j][2:(psi-1)] = gammaTemp[j][1:(psi-2)]*100;
          gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * (bsNonlinear[indexFL[(((j-1)*2)+1):(((j-1)*2)+2)], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j]*100) * (-1);
          gamma[j][1] = gammaFL[j][1];
          gamma[j][psi] = gammaFL[j][2];
      };
      for (j in 1:p){
          gsmooth[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j] + bsLinear[,j] * theta[j+1];
      };

      for (i in 1:n){
          alpha[i] = exp(theta[1] + sum(gsmooth[i,]));
      };
    }

}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u[i], alpha[i]);
    }
    target += normal_lpdf(theta[1] | 0, 100);
    target += gamma_lpdf(lambda1 | 1, 1e-3);
    target += gamma_lpdf(lambda2o | 1, 1e-3);
    target += (2*p*log(lambda2o));
    for (j in 1:p){
        target += gamma_lpdf(tau1[j] | 1, lambda1^2*0.5);
        target += normal_lpdf(theta[(j+1)] | 0, sqrt(1/tau1[j]));
        target += gamma_lpdf(tau2[j] | atau, lambda2o^2*0.5);
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * (1/tau2[j]));
    }
}
generated quantities {
    // Used in Posterior predictive check
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgsmooth; // linear component

    for (j in 1:p){
        newgsmooth[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j] + xholderLinear[,j] * theta[j+1];
    };    

    for (i in 1:n){ 
        newalpha[i] = exp(theta[1] + sum(newgsmooth[i,]));
    };
}
"


data.stan <- list(y = as.vector(y), u = as.vector(u), p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), indexFL = as.vector(t(index.holder)),
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear, basisFL = basis.holder)

init.alpha <- list(list(gammaTemp = array(rep(0.1, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0.1, (p+1)), 
                        tau1 = rep(0.01, p),tau2 = rep(0.01, p),
                        lambda1 = 1, lambda2 = 2),
                   list(gammaTemp = array(rep(0.1, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0.2, (p+1)), 
                        tau1 = rep(0.1, p),tau2 = rep(0.1, p),
                        lambda1 = 1, lambda2 = 1),
                   list(gammaTemp = array(rep(0.1, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0.5, (p+1)), 
                        tau1 = rep(0.05, p),tau2 = rep(0.05, p),
                        lambda1 = 0.5, lambda2 = 1.5))

# stanc("C:/Users/Johnny Lee/Documents/GitHub/BLAST/application/model1.stan")
fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    chains = 3,             # number of Markov chains
    iter = 15000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 2500           # no progress shown
)

posterior <- rstan::extract(fit1)

theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
origin.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary 
alpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary
# intcpt.samples <- summary(fit1, par=c("intcpt"), probs = c(0.05, 0.5, 0.95))

g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))

alpha.smooth <- data.frame("x" = xholder[,1],
                            "post.mean" = (alpha.samples[,1]),
                            "post.median" = (alpha.samples[,5]),
                            "q1" = (alpha.samples[,4]),
                            "q3" = (alpha.samples[,6]),
                            "origin.post.mean" = (origin.samples[,1]),
                            "origin.post.median" = (origin.samples[,5]),
                            "origin.q1" = (origin.samples[,4]),
                            "origin.q3" = (origin.samples[,6]))

ggplot(alpha.smooth, aes(x=x)) + 
  ylab(expression(alpha(c,ldots,c))) + xlab(expression(c)) + labs(col = "") +
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

xi.smooth <- data.frame("x" = xholder[,1],
                            "post.mean" = 1/(alpha.samples[,1]),
                            "post.median" = 1/(alpha.samples[,5]),
                            "q1" = 1/(alpha.samples[,4]),
                            "q3" = 1/(alpha.samples[,6]))

ggplot(xi.smooth, aes(x=x)) + 
  ylab(expression(xi(c,ldots,c))) + xlab(expression(c)) + labs(col = "") +
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
                  ylim(-1, 1.1) + #max(data.smooth$q3[((((i-1)*n)+1):(i*n))])) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=-1, color = "red", size = 7)
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_DC.pdf"), grid.plts[[7]], width=10, height = 7.78)

