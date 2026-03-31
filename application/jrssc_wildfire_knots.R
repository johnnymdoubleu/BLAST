library(VGAM)
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
library(evgam)
library(forecast)
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
psi.origin <- psi <- 10
threshold <- 0.95

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
fwi.scaled$time <- fwi.index$time <- seq(1,length(Y), length.out=length(Y))
fwi.scaled$sea <- fwi.index$sea <- fwi.index$time %% 365.25 / 365.25
fwi.scaled$cos.time <- fwi.index$cos.time <- cos(2*pi*seq(1,length(Y), length.out=length(Y))/365.25)
fwi.scaled$sin.time <- fwi.index$sin.time <- sin(2*pi*seq(1,length(Y), length.out=length(Y))/365.25)
fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)

# time_arima <- seq(1, length(Y))
# xreg_trend_seasonal <- cbind(
#   trend = time_arima,
#   cos_season = cos(2 * pi * time_arima / 365.25),
#   sin_season = sin(2 * pi * time_arima / 365.25)
# )

# fit.list <- list()
# for (j in 1:7) {
#   y_ts <- ts(fwi.scaled[, j], frequency = 365.25)
#   fit.list[[j]] <- forecast::auto.arima(
#     y_ts,
#     xreg = xreg_trend_seasonal,
#     seasonal = FALSE,
#     stepwise = TRUE,
#     approximation = FALSE
#   )
#   fwi.index[, j] <- fwi.scaled[, j] <- as.numeric(residuals(fit.list[[j]]))
# }

xreg.season <- cbind(
  trend = c(1:length(Y)),
  cos_season = cos(2 * pi * c(1:length(Y)) / 365.25),
  sin_season = sin(2 * pi * c(1:length(Y)) / 365.25)
)

# fit.list <- list()
# x.detrended <- matrix(nrow = length(Y), ncol = 7)
# for (j in 1:7) {
#   y_ts <- ts(fwi.scaled[, j], frequency = 365.25) 
#   fit.list[[j]] <- fit <- auto.arima(y_ts, seasonal = FALSE, xreg = xreg.season, stepwise = TRUE, approximation = FALSE)
#   x.detrended[, j] <- as.numeric(residuals(fit.list[[j]]))
# }
# fwi.index[,1:7] <- fwi.scaled[, 1:7] <- x.detrended
# acf(fwi.index$BUI)
# acf(fwi.index$ISI)
# acf(fwi.index$FFMC)
# acf(fwi.index$DMC)
# acf(fwi.index$DC)

above.0 <- which(Y > 0)
Y_pos <- Y[above.0]
fwi_pos <- fwi.scaled[above.0, ]
# Y_pos <- Y
# fwi_pos <- fwi.scaled
# pca_result <- prcomp(fwi_pos[,3:7], center = TRUE, scale. = TRUE)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# fwi_pos[,1:7] <- as.data.frame(sapply(fwi_pos[,1:7], FUN = range01))
# qr.df <- data.frame(y = log(Y_pos), pca_result$x, cos.time = fwi_pos$cos.time, sin.time = fwi_pos$sin.time) #fwi_pos)
qr.df <- data.frame(y = log(Y_pos), scale(fwi_pos[,3:7]), cos.time = fwi_pos$cos.time, sin.time = fwi_pos$sin.time)
# evgam.cov <- y ~ 1 + cos.time + sin.time + s(PC1) + s(PC2) + s(PC3) + s(PC4) + s(PC5)
evgam.cov <- y ~ cos.time + sin.time + s(BUI, bs = "ts", k = 30) + s(ISI, bs = "ts", k = 30) + s(FFMC, bs = "ts", k = 30) + s(DMC, bs = "ts", k = 30) + s(DC, bs = "ts", k = 30) 

k.values <- sort(c(5, 7, 10, 15, 17, 20, 25, 27, 30, 35, 37, 40), decreasing = TRUE)
aic.results <- data.frame(k = k.values, AIC = NA)
library(foreach)
library(doParallel)

# Setup clusters (detect number of cores and use all but one)
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

# Run the loop in parallel
aic_results_list <- foreach(current.k = k.values, .packages = c("evgam"), .errorhandling = "pass") %dopar% {
  
  form_str <- paste0("y ~ 1 + cos.time + sin.time + ",
                    paste0("s(", colnames(fwi_pos[,3:7]), ", k = ", current.k, ", bs='" ,"ts')", collapse = " + "))
  
  fit <- evgam(as.formula(form_str), data = evgam.df, 
              family = "ald", ald.args = list(tau = threshold))
  
  data.frame(k = current.k, AIC = AIC(fit))
}

stopCluster(cl)

# Combine results
aic.results <- do.call(rbind, aic_results_list)
best.k <- aic.results[which.min(aic.results$AIC), ]
print(best.k)

load("./BLAST/application/2026-03-30_pareto_950_knots.Rdata")

ggplot(aic.results, aes(x=k, y = BIC)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkblue", size = 3) +
  geom_point(data = subset(aic.results, k == 10), aes(x = k, y = BIC), 
            color = "red", size = 5, shape = 1, stroke = 2) +
  labs(x = "Number of Knots", y = "AIC") +
  theme_minimal(base_size = 25) +
  theme(plot.title = element_text(face = "bold"))


# qr.fit <- evgam(evgam.cov, data = qr.df, family = "ald", ald.args=list(tau = threshold))

# qr.lin <- evgam(y ~ cos.time + sin.time + BUI + ISI + FFMC + DMC + DC, data = qr.df, family = "ald", ald.args=list(tau = threshold))
# qr.cov <- evgam(y ~ s(BUI, bs = "ts", k = 30) + s(ISI, bs = "ts", k = 30) + s(FFMC, bs = "ts", k = 30) + s(DMC, bs = "ts", k = 30) + s(DC, bs = "ts", k = 30), data = qr.df, family = "ald", ald.args=list(tau = threshold))
# qr.time <- evgam(y ~ cos.time + sin.time, data = qr.df, family = "ald", ald.args=list(tau = threshold))
# qr.null <- evgam(y ~ 1, data = qr.df, family = "ald", ald.args=list(tau = threshold))
# u.vec <- exp(predict(qr.fit)$location)
# save(u.c, qr.fit, file = paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_qr-c.Rdata"))
load("./BLAST/application/figures/2026-03-25_pareto_qr-ct.Rdata")
u.vec <- u.ct
# qr.fit <- quantreg::rq(y ~ 1 + cos.time + sin.time + BUI + ISI + FFMC + DMC + DC, data = qr.df, tau = threshold)
# u.vec <- exp(predict(qr.fit))  # threshold on raw scale for Y_pos
# AIC(qr.fit, qr.cov, qr.time, qr.null)
# BIC(qr.fit, qr.cov, qr.time, qr.null)

plot(c(1:length(Y_pos)), log(Y_pos))
lines(c(1:length(Y_pos)), log(u.vec), type = "l", col = "red")

