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

fwi.origin <- data.frame(fwi.origin, BA=Y)

vars <- colnames(fwi.origin)[1:7]
seq_list <- lapply(vars, function(var) seq(min(fwi.origin[[var]]), max(fwi.origin[[var]]), length.out = length(Y)))
grid.df <- as.data.frame(setNames(seq_list, vars))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled <- as.data.frame(sapply(fwi.scaled, FUN = range01))
fwi.df <- data.frame(fwi.scaled, BA=Y)
quant.fit <- qgam(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + 
                    s(FFMC) + s(DMC) + s(DC),
                data = fwi.df, qu = 0.975)

print(plot(quant.fit, allTerms = TRUE), pages = 1)
quant.viz <- getViz(quant.fit, nsim = 20)
check1D(quant.viz, fwi.df[,c(1:7)]) + l_gridQCheck1D(qu = 0.95)



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

# quant.fit <- qgam(BA ~ s(DSR, k = 30) + s(FWI, k = 30) + s(BUI, k = 30) + s(ISI, k = 30) + s(FFMC, k = 30) + s(DMC, k = 30) + s(DC, k = 30), qu=0.975, data = fwi.origin)
# quant.u <- predict(quant.fit, newdata = fwi.origin)
# quant.fit <- qgamV(BA ~ s(DSR, k = 30) + s(FWI, k = 30) + s(BUI, k = 30) + s(ISI, k = 30) + s(FFMC, k = 30) + s(DMC, k = 30) + s(DC, k = 30), qu=0.975, data = fwi.origin)
# quant.fit <- qgamV(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC), qu=0.99, data = fwi.df)


# qdo(quant.fit, predict, newdata = )

# save(quant.fit, quant.u, fwi.scaled, y, n,p, u, file="./BLAST/application/quant975_prep.Rdata")
# save(fwi.scaled, fwi.origin, qu, y, Y, u, n, p, psi, file = "./BLAST/application/wildfire_prep.Rdata")
# load("./BLAST/application/quant975_prep.Rdata")

# m.gam <- mqgam(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC), data = fwi.origin, qu = taus)
# m.gam <- mqgamV(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC), data = fwi.origin, qu = taus)
# save(m.gam, file ="./BLAST/application/wildfire_mgam2.Rdata")
# load("./BLAST/application/wildfire_mgam.Rdata")

