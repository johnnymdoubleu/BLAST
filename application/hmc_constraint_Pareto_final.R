setwd("C:/Users/Johnny Lee/Documents/GitHub")
# source("./extremis/R/bltir.R")

library(npreg)
library(Pareto)
suppressMessages(library(tidyverse))
library(readxl)
library(patchwork)
library(corrplot)
library(rstan)
library(loo)
library(qqboxplot)
library(ggdensity)
library(ggforce)
library(ggdist)
library(extremis)


options(mc.cores = parallel::detectCores())

# Structure of the FWI System
#DSR : Dail Severity Rating
#FWI : Fire Weather Index
#BUI : Buildup Index
#ISI : Initial Spread Index
#FFMC : Fine FUel Moisture Code
#DMC : Duff Moisture Code
#DC : Drought Code


df <- read_excel("./BRSTIR/application/AADiarioAnual.xlsx", col_types = c("date", rep("numeric",40)))
df.long <- gather(df, condition, measurement, "1980":"2019", factor_key=TRUE)
df.long
head(df.long)
tail(df.long)
# View(df.long[is.na(df.long$measurement),])
missing.values <- which(!is.na(df.long$measurement))
df.long[which(is.na(df.long$measurement)),]
df.long[which(is.na(df.long$...1))+1,]

#NAs on Feb 29 most years, and Feb 14, 1999
#considering the case of leap year, the missing values are the 29th of Feb
#Thus, each year consist of 366 data with either 1 or 0 missing value.
Y <- df.long$measurement[!is.na(df.long$measurement)]

psi <- 30
threshold <- 0.975

# x.scale <- x.scale[which(y>quantile(y, threshold)),]
# u <- quantile(y, threshold)

multiplesheets <- function(fname) {
    setwd("C:/Users/Johnny Lee/Documents/GitHub")
    # getting info about all excel sheets
    sheets <- excel_sheets(fname)
    tibble <- lapply(sheets, function(x) read_excel(fname, sheet = x, col_types = c("date", rep("numeric", 41))))
    # print(tibble)
    data_frame <- lapply(tibble, as.data.frame)
    # assigning names to data frames
    names(data_frame) <- sheets
    return(data_frame)
}
setwd("C:/Users/Johnny Lee/Documents/GitHub")
path <- "./BRSTIR/application/DadosDiariosPT_FWI.xlsx"
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
# cov.long$ <- gather(cov$DSR[!is.na(df.long$measurement)][,1:41], )
for(i in 1:length(cov)){
    cov.long <- gather(cov[[i]][,1:41], condition, measurement, "1980":"2019", factor_key=TRUE)
    fwi.index[,i] <- cov.long$measurement[missing.values]
    fwi.scaled[,i] <- cov.long$measurement[missing.values]
}

fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)

fit <- bltir(Y, fwi.scaled, psi=20, T=5000, prior=list(a=1, b=1e-3), knot.pos="quantile", threshold=0.975)





n <- 5000; psi <- 10; threshold <- 0.95; p <- 6;

# Function to generate Gaussian copula in uniform margin
C <- diag(p)
xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow = n, ncol = 0)
x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))

basis.holder <- matrix(, nrow = 2, ncol = 0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                        matrix(c(which.min(x.origin[,i]),
                                 which.max(x.origin[,i])), ncol = 2))
}
for(i in 1:p){
  knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
  tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  bs.linear <- cbind(bs.linear, tps[,1])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-1])
  basis.holder <- cbind(basis.holder, 
                        solve(t(matrix(c(tps[index.holder[i,1], 2],
                                         tps[index.holder[i,1], 1+psi],
                                         tps[index.holder[i,2], 2],
                                         tps[index.holder[i,2], 1+psi]), 
                                       nrow = 2, ncol = 2))))
}

## Setting true parameters
gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
  for (ps in 1:psi){
    if(j %in% c(1,4,5,6)){gamma.origin[ps,j] <- 0}
    else {
      if(ps==1 || ps==psi){gamma.origin[ps,j] <- 0}
      else{gamma.origin[ps,j] <- -25}
    }
  }
}
theta.origin <- c(-0.5, 0, -0.5, -0.5, 0, 0, -0.5)

## Injecting constraints to each ends of the smooth functions
g.sub.origin <- matrix(, nrow = 2, ncol = p)
for(j in 1:p){
  g.sub.origin[,j] <- as.matrix(bs.nonlinear[index.holder[j,], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))], nrow = 2) %*% gamma.origin[(2:(psi-1)), j]
  gamma.origin[c(1,psi),j] <- -1 * basis.holder[,(((j-1)*2)+1):(((j-1)*2)+2)] %*% as.matrix(g.sub.origin[,j], nrow=2)
}

g.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){ 
  g.origin[,j] <- (bs.linear[,j] * theta.origin[j+1]) + (bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
}

y.origin <- NULL
for(i in 1:n){
  y.origin[i] <- rPareto(1, 1, alpha = exp(theta.origin[1] + sum(g.origin[i,]))) 
}

fit <- bltir(Y = y.origin, X = x.origin, knot.pos = "equal", 
              threshold = 0.95, psi = 10, T = 3000)

plot(fit, bands = TRUE, option = "alpha")
plot(fit, bands = TRUE, option = "smooth")
