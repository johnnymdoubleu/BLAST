library(VGAM)
library(mgcv)
library(evgam)
library(Pareto)
suppressMessages(library(tidyverse))
library(readxl)
library(gridExtra)
library(ggdist)
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
# u <- quantile(Y[Y>1], threshold)

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

fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)

# load("./BLAST/application/quant-evgam.Rdata")
# load("./BLAST/application/quant-t_10.Rdata")
load("./BLAST/application/qgam_95_30.Rdata")
# u <- quantile(Y, threshold)
# excess <- which(Y>u)
preds <- predict(quant.fit)
excess <- which(Y>preds)
u <- preds[excess]
# excess <- which(fwi.dd$excess==TRUE)
# u <- fwi.dd$origin_Model_Smooth_975[excess]
y <- Y[excess]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled <- as.data.frame(sapply(fwi.scaled[excess,c(-1,-2)], FUN = range01))

n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]


fwi.origin <- data.frame(fwi.index[excess,c(-1,-2)], BA=y)
max.fwi <- fwi.origin[which.max(y),]
fwi.grid <- data.frame(lapply(fwi.origin[,c(1:p)], function(x) seq(min(x), max(x), length.out = nrow(fwi.scaled))))
fwi.minmax <- sapply(fwi.origin[,c(1:p)], function(x) max(x)-min(x))
fwi.min <- sapply(fwi.origin[,c(1:p)], function(x) min(x))

simul.data <- data.frame(BA = y-u, fwi.scaled[,c(1:p)])
newx <- seq(0, 1, length.out = n)
xholder <- as.data.frame(do.call(cbind, lapply(1:p, function(j) {newx})))
colnames(xholder) <- colnames(fwi.origin[,c(1:p)])
gam.scale <- list(BA ~ s(BUI, bs = "tp", k = 30) + 
                      s(ISI, bs = "tp", k = 30) + 
                      s(FFMC, bs = "tp", k = 30) +
                      s(DMC, bs = "tp", k = 30) + 
                      s(DC, bs = "tp", k = 30),
                    ~ s(BUI, bs = "tp", k = 30) + 
                      s(ISI, bs = "tp", k = 30) + 
                      s(FFMC, bs = "tp", k = 30) +
                      s(DMC, bs = "tp", k = 30) + 
                      s(DC, bs = "tp", k = 30))
evgam.fit.scale <- evgam::evgam(gam.scale, data = simul.data, family = "gpd")
xi.pred.scale <-predict(evgam.fit.scale, newdata = xholder, type="response")$shape
sigma.pred.scale <-predict(evgam.fit.scale, newdata = xholder, type="response")$scale
alpha.pred.scale <- 1/xi.pred.scale

xholder.basis.scale <- predict(evgam.fit.scale, newdata = xholder, type= "lpmatrix")$shape
psi <- 30
xi.coef.scale <- tail(evgam.fit.scale$coefficients, (psi-1)*p)
gamma.xi.scale <- matrix(xi.coef.scale, ncol = p)
alpha.nonlinear.scale <- xi.nonlinear.scale <- matrix(, nrow = n, ncol = p)
bs.nonlinear.scale <- xholder.basis.scale[,c(2:((psi-1)*p+1))]
for(j in 1:p){
  xi.nonlinear.scale[,j] <- bs.nonlinear.scale[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.scale[,j]
  alpha.nonlinear.scale[,j] <- 1/(xi.nonlinear.scale[,j])
}

gam.1 <- list(BA ~ 1,
                ~ s(BUI, bs = "tp", k = 30) + 
                  s(ISI, bs = "tp", k = 30) + 
                  s(FFMC, bs = "tp", k = 30) +
                  s(DMC, bs = "tp", k = 30) + 
                  s(DC, bs = "tp", k = 30))
evgam.fit.1 <- evgam::evgam(gam.1, data = simul.data, family = "gpd")
xi.pred.1 <-predict(evgam.fit.1, newdata = xholder, type="response")$shape
sigma.pred.1 <-predict(evgam.fit.1, newdata = xholder, type="response")$scale
alpha.pred.1 <- 1/xi.pred.1

xholder.basis.1 <- predict(evgam.fit.1, newdata = xholder, type= "lpmatrix")$shape
xi.coef.1 <- tail(evgam.fit.1$coefficients, (psi-1)*p)
gamma.xi.1 <- matrix(xi.coef.1, ncol = p)
alpha.nonlinear.1 <- xi.nonlinear.1 <- matrix(, nrow = n, ncol = p)
bs.nonlinear.1 <- xholder.basis.1[,c(2:((psi-1)*p+1))]
for(j in 1:p){
  xi.nonlinear.1[,j] <- bs.nonlinear.1[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.1[,j]
  alpha.nonlinear.1[,j] <- 1/(xi.nonlinear.1[,j])
}
# plot(evgam.fit.1)

simul.data <- data.frame(BA = y, fwi.scaled[,c(1:p)])
vgam.fit.scale <- vgam(BA ~ sm.ps(BUI, ps.int = 28, outer.ok = TRUE) + 
                            sm.ps(ISI, ps.int = 28, outer.ok = TRUE) + 
                            sm.ps(FFMC, ps.int = 28, outer.ok = TRUE) + 
                            sm.ps(DMC, ps.int = 28, outer.ok = TRUE) + 
                            sm.ps(DC, ps.int = 28, outer.ok = TRUE),
                      data = simul.data,
                      family = gpd(threshold= u,
                                    lshape="loglink",
                                    zero = NULL),
                      trace = TRUE,
                      control = vgam.control(maxit = 200))
fitted.linear <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "link")
fitted.terms <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "terms")
vgam.xi.scale <- exp(fitted.linear[,2])
vgam.sigma.scale <- exp(fitted.linear[,1])

vgam.fit.1 <- vgam(BA ~ sm.ps(BUI, ps.int = 28, outer.ok = TRUE) + 
                        sm.ps(ISI, ps.int = 28, outer.ok = TRUE) + 
                        sm.ps(FFMC, ps.int = 28, outer.ok = TRUE) + 
                        sm.ps(DMC, ps.int = 28, outer.ok = TRUE) + 
                        sm.ps(DC, ps.int = 28, outer.ok = TRUE),
                      data = simul.data,
                      family = gpd(threshold= u,
                                    lshape="loglink",
                                    zero = 1),
                      trace = TRUE,
                      control = vgam.control(maxit = 200))
fitted.linear <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "link")
fitted.terms <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "terms")
vgam.xi.1 <- exp(fitted.linear[,2])
vgam.sigma.1 <- exp(fitted.linear[,1])

# alpha.smooth <- data.frame("x" = as.vector(as.matrix(xholder)),
#                           "evgam.scale" = as.vector(alpha.nonlinear.scale),
#                           "evgam.1" = as.vector(alpha.nonlinear.1),
#                           "covariates" = gl(p, n, (p*n), labels = colnames(xholder)))


# grid.plts <- list()
# for(i in 1:p){
#   grid.plt <- ggplot(data = data.frame(alpha.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
#                   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#                   geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1) + 
#                   geom_line(aes(y=evgam.1), colour = "orange", linewidth=1) + 
#                   ylab("") + xlab(expression(c)) +
#                   # ylim(-2.3, 2.3) +
#                   theme_minimal(base_size = 10) +
#                   theme(legend.position = "none",
#                           plot.margin = margin(0,0,0,-5),
#                           axis.text = element_text(size = 15),
#                           axis.title.x = element_text(size = 15))
#   grid.plts[[i]] <- grid.plt
# }

# grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

xi.smooth <- data.frame("x" = as.vector(as.matrix(xholder)),
                        "evgam.scale" = as.vector(xi.nonlinear.scale),
                        "evgam.1" = as.vector(xi.nonlinear.1),
                        "covariates" = gl(p, n, (p*n), labels = colnames(xholder)))

grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(xi.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1) + 
                  geom_line(aes(y=evgam.1), colour = "orange", linewidth=1) + 
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue", "red")) +
                  ylim(-0.5, 0.5) + xlim(0,1) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

alpha.scenario <- data.frame("x" = newx,
                            "evgam.scale" = alpha.pred.scale,
                            "evgam.1" = alpha.pred.1)

ggplot(alpha.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,ldots,c))) + xlab(expression(c)) + labs(col = "") +
  geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1.5, linetype=2) +
  geom_line(aes(y=evgam.1), colour = "orange", linewidth=1.5, linetype=2) +
  theme_minimal(base_size = 30) + #ylim(0,20)+
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_evgam_scA.pdf"), width=10, height = 7.78)

xi.scenario <- data.frame("x" = newx,
                          "vgam.scale" = vgam.xi.scale,
                          "vgam.1" = vgam.xi.1,
                          "evgam.scale" = xi.pred.scale,
                          "evgam.1" = xi.pred.1)

ggplot(xi.scenario, aes(x=x)) + 
  ylab(expression(xi(c,ldots,c))) + xlab(expression(c)) +
  geom_line(aes(y = evgam.scale, color = "Scale", linetype = "EVGAM"), linewidth = 1.5) +
  geom_line(aes(y = evgam.1,     color = "No Scale",     linetype = "EVGAM"), linewidth = 1.5) +
  geom_line(aes(y = vgam.scale,  color = "Scale", linetype = "VGAM"),  linewidth = 1.5) +
  geom_line(aes(y = vgam.1,      color = "No Scale",     linetype = "VGAM"),  linewidth = 1.5) +
  scale_color_manual(name = "Scenario", values = c("Scale" = "orange", "No Scale" = "purple")) +
  scale_linetype_manual(name = "Model", values = c("EVGAM" = 1, "VGAM" = 4)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.15, 0.8),
        legend.text = element_text(size = 16),
        strip.text = element_blank(),
        axis.text = element_text(size = 18))

sigma.scenario <- data.frame("x" = newx,
                            "vgam.scale" = unname(vgam.sigma.scale),
                            "vgam.1" = unname(vgam.sigma.1),
                            "evgam.scale" = sigma.pred.scale,
                            "evgam.1" = rep(sigma.pred.1,length(newx)))

ggplot(sigma.scenario, aes(x=x)) + 
  ylab(expression(sigma(c,ldots,c))) + xlab(expression(c)) +
  geom_line(aes(y = evgam.scale, color = "Scale", linetype = "EVGAM"), linewidth = 1.5) +
  geom_line(aes(y = evgam.1,     color = "No Scale",     linetype = "EVGAM"), linewidth = 1.5) +
  geom_line(aes(y = vgam.scale,  color = "Scale", linetype = "VGAM"),  linewidth = 1.5) +
  geom_line(aes(y = vgam.1,      color = "No Scale",     linetype = "VGAM"),  linewidth = 1.5) +
  scale_color_manual(name = "Scenario", values = c("Scale" = "orange", "No Scale" = "purple")) +
  scale_linetype_manual(name = "Model", values = c("EVGAM" = 1, "VGAM" = 4)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.15, 0.8),
        legend.text = element_text(size = 16),
        strip.text = element_blank(),
        axis.text = element_text(size = 18))




# save(vgam.fit.1, vgam.fit.scale, evgam.fit.1, evgam.fit.scale, file = "./BLAST/application/vgam-fit.Rdata")
load("./BLAST/application/vgam-fit.Rdata")

AIC(vgam.fit.1)
AIC(vgam.fit.scale)
anova(vgam.fit.1, vgam.fit.scale, test = "LRT", type=1)
# AIC(evgam.fit.1)
# AIC(evgam.fit.scale)

library(ggplot2)
library(gridExtra)
library(dplyr)

vars <- c("BUI", "ISI", "FFMC", "DMC", "DC")
xi_list <- list()

pred.vgam.1 <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "terms")
pred.vgam.scale <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "terms")
for (v in vars) {
  term_col  <- grep(v, colnames(pred.vgam.1), value = TRUE)[2]
  fitted.1 <- as.numeric(pred.vgam.1[, term_col])
  fitted.scale <- as.numeric(pred.vgam.scale[, term_col])
  xi_list[[v]] <- data.frame(
    x = newx,
    vgam.1 = fitted.1,#[367:(length(newx)*2)],
    vgam.scale = fitted.scale#[367:(length(newx)*2)]
  )
}

vgam.smooth <- bind_rows(xi_list)

xi.smooth <- data.frame("x" = as.vector(as.matrix(xholder)),
                        "vgam.scale" = vgam.smooth$vgam.scale,
                        "vgam.1" = vgam.smooth$vgam.1,
                        "evgam.scale" = as.vector(xi.nonlinear.scale),
                        "evgam.1" = as.vector(xi.nonlinear.1),
                        "covariates" = gl(p, n, (p*n), labels = colnames(xholder)))

grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(xi.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_line(aes(y=vgam.scale), colour = "purple", linewidth=1, linetype = 3) + 
                  geom_line(aes(y=vgam.1), colour = "orange", linewidth=1, linetype= 3) +                   
                  geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1) + 
                  geom_line(aes(y=evgam.1), colour = "orange", linewidth=1) + 
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue", "red")) +
                  # ylim(-0.5, 0.5) + xlim(0,1) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)
