library(qgam)
library(mgcViz)
library(parallel)
library(evgam)

load("wildfire_prep.Rdata")
#load("wildfire_time_prep.Rdata")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled <- as.data.frame(sapply(fwi.origin[,c(1:7)], FUN = range01))
fwi.df <- data.frame(fwi.scaled[,c(-1,-2)], BA=Y)
fwi.df$time <- c(1:length(Y))
X_means <- colMeans(fwi.df[,c(1:5)])
X_sd   <- apply(fwi.df[,c(1:5)], 2, sd)
fwi.df[,c(1:5)] <- scale(fwi.df[,c(1:5)], center = X_means, scale = X_sd)
#quant.fit <- brm(bf(BA ~ s(DSR) + s(FWI) + s(BUI) + 
#                      s(ISI) + s(FFMC) + s(DMC) + 
#                      s(DC), quantile = 0.975),
#                data = fwi.df,
#                cores = 3,
#                chain = 3,
#                family = asym_laplace())
n.core <- 12
cl <- makeCluster(n.core)

#tun <- tuneLearnFast(BA ~ s(DSR, k = 30) + s(FWI, k = 30) + s(BUI, k = 30) + 
#                            s(ISI, k = 30) + s(FFMC, k = 30) + s(DMC, k = 30) + 
#                            s(DC, k = 30), data = fwi.df, qu = 0.975,
#			  cluster = cl, ncores = n.core)
			  #control = list(parallel ="multcore", ncores=n.core))
quant.fit <- qgam(BA ~ s(BUI, k=30, bs = "ts") + 
											 s(ISI, k=30, bs = "ts") + 
											 s(FFMC, k=30, bs = "ts") + 
											 s(DMC, k=30, bs = "ts") + 
											 s(DC, k=30, bs = "ts"), 
												data = fwi.df, qu = 0.96,
												cluster = cl, ncores = n.core)
 			#control = list(parallel ="multcore", ncores=n.core))
#fwi.df <- fwi.origin[which(fwi.origin$BA>1),c(1:8)]
#taus <- c(0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95)
#taus <- c(0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)
#m.gam <- mqgamV(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC) + s(time),
#		  data = fwi.df, qu = taus, aQgam = list(cluster = cl, ncores=n.cores, multicore = TRUE, paropts = opts))
stopCluster(cl)

#evgam.time <- BA ~ s(time, k = 30)
#ald.time.fit <- evgam(evgam.time, data = fwi.df, family = "ald", ald.args=list(tau = 0.975))
#evgam.cov <- BA ~ s(BUI, k = 30) + s(ISI, k = 30) + s(FFMC, k = 30) + s(DMC, k = 30) + s(DC, k = 30)
#ald.cov.fit <- evgam(evgam.cov, data = fwi.df, family = "ald", ald.args=list(tau = 0.975))

#check(quant.fit)
save(quant.fit, file="qgam_96_30_ts.Rdata")
#save(ald.time.fit, ald.cov.fit, file="quant-evgam-scaled.Rdata")
