library(npreg)
library(Pareto)
suppressMessages(library(tidyverse))
library(readxl)
library(gridExtra)
library(VGAM)


# Structure of the FWI System
#DSR : Daily Severity Rating
#FWI : Fire Weather Index
#BUI : Buildup Index
#ISI : Initial Spread Index
#FFMC : Fine FUel Moisture Code
#DMC : Duff Moisture Code
#DC : Drought Code


setwd("C:/Users/Johnny Lee/Documents/GitHub")
# setwd("A:/GitHub")
load("./BLAST/application/wildfire_prep.Rdata") #loading covariate-dependent thresholds
psi <- 30
u <- quantile(Y, 0.975)
y <- Y[which(Y>u)]
fwi.scaled <- fwi.origin[which(Y>u),c(1:7)]
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled <- as.data.frame(sapply(fwi.scaled, FUN = range01))
n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]
fwi.df <- data.frame(fwi.scaled, BA=y)

vgam.fit.scale <- vgam(BA ~ s(DSR) + s(FWI) + s(BUI) + s(ISI) + s(FFMC) + s(DMC) + s(DC),
                        data = fwi.df,
                        family = gpd(threshold = u,
                                      lshape="loglink",
                                      zero = NULL),
                        trace = TRUE,
                        control = vgam.control(maxit = 200))
par(mfrow = c(5, 2), mar=c(1.5,1.5,1.5,1.5))
plot(vgam.fit.scale, se = TRUE, shade = TRUE, shcol = "steelblue")
par(mfrow = c(1, 1))
fitted.linear <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "link")
fitted.terms <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "terms")
fitted.response <- predict(vgam.fit.scale,newdata = data.frame(xholder), type = "response")
vgam.xi.scale <- exp(fitted.linear[,2])

# save(m.gpd, file = "./BLAST/application/vgam_fit_all.Rdata")
# load("./BLAST/application/evgam_fit.Rdata")

no.theta <- 1 #represents the no. of linear predictors for each smooth functions
xholder.linear <- xholder.nonlinear <- matrix(,nrow=n, ncol=0)
xholder <- as.data.frame(matrix(nrow=n, ncol=p))
xholder.mat <- matrix(nrow=n, ncol=p)
for(i in 1:p){
  xholder.mat[,i] <- xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
  test.knot <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
}
names(xholder) <- names(fwi.scaled)

basis_list <- list()
for (v in names(fwi.scaled)) {
  sc <- mgcv::smoothCon(mgcv::s(x, bs="tp", k=30), data = data.frame(x = xholder[[v]]))
  basis_list[[v]] <- sc[[1]]$X
}


xholder.basis.1 <-predict(m.gpd, newdata = xholder, type = "lpmatrix")

xi.coef.1 <- tail(m.gpd$coefficients, 204)
gamma.xi.1 <- matrix(xi.coef.1[2:204], ncol = p)
f.nonlinear.1 <- matrix(, nrow = n, ncol = p)
bs.nonlinear.1 <- xholder.basis.1$shape[,c(2:204)]
for(j in 1:p){
  f.nonlinear.1[,j] <- bs.nonlinear.1[,(((j-1)*29)+1):(((j-1)*29)+29)] %*% gamma.xi.1[,j]
}

# smooth.func <- as.data.frame(m.gpd$plotdata)
# smooth.func <- as.matrix(smooth.func)
# smooth.func <- as.vector(smooth.func)


xholder.df <- data.frame(xholder)
names(xholder.df) <- names(fwi.scaled)
xi.pred.1 <-predict(m.gpd, newdata = xholder.df, type="link")$shape
alpha.pred.1 <- 1/xi.pred.1

xholder.basis.1 <- predict(m.gpd, newdata = xholder.df, type= "lpmatrix")$shape
xi.coef.1 <- tail(m.gpd$coefficients, (psi-1)*p)
gamma.xi.1 <- matrix(xi.coef.1, ncol = p)
alpha.nonlinear.1 <- xi.nonlinear.1 <- matrix(, nrow = n, ncol = p)
bs.nonlinear.1 <- xholder.basis.1[,c(2:((psi-1)*p+1))]
for(j in 1:p){
  xi.nonlinear.1[,j] <- bs.nonlinear.1[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.1[,j]
  alpha.nonlinear.1[,j] <- 1/xi.nonlinear.1[,j]
}

load("./BLAST/application/evgam_fit_all.Rdata")
xholder.basis.scale <-predict(m.gpd, newdata = xholder, type = "lpmatrix")

xi.coef.scale <- tail(m.gpd$coefficients, 204)
gamma.xi.scale <- matrix(xi.coef.1[2:204], ncol = p)
f.nonlinear.scale <- matrix(, nrow = n, ncol = p)
bs.nonlinear.scale <- xholder.basis.scale$shape[,c(2:204)]
for(j in 1:p){
  f.nonlinear.scale[,j] <- bs.nonlinear.scale[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.scale[,j]
}

xi.pred.scale <-predict(m.gpd, newdata = xholder.df, type="link")$shape
alpha.pred.scale <- 1/xi.pred.scale

xi.coef.scale <- tail(m.gpd$coefficients, (psi-1)*p)
gamma.xi.scale <- matrix(xi.coef.scale, ncol = p)
alpha.nonlinear.scale <- xi.nonlinear.scale <- matrix(, nrow = n, ncol = p)
bs.nonlinear.scale <- xholder.basis.scale$shape[,c(2:((psi-1)*p+1))]
for(j in 1:p){
  xi.nonlinear.scale[,j] <- bs.nonlinear.scale[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.scale[,j]
  alpha.nonlinear.scale[,j] <- 1/xi.nonlinear.scale[,j]
}

xi.smooth <- data.frame("x"= as.vector(xholder.mat),
                          # "fit" = smooth.func,
                          "true" = as.vector(as.matrix(fwi.scaled)),
                          "evgam.1" = as.vector(xi.nonlinear.1),
                          "evgam.scale" = as.vector(xi.nonlinear.scale),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
                          "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(xi.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_line(aes(y=evgam.1), colour = "purple", linewidth=1) + 
                  geom_line(aes(y=evgam.scale), colour = "orange", linewidth=1) + 
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  geom_rug(aes(x=true, y=evgam.1), sides = "b") + 
                  # scale_color_manual(values=c("purple")) +
                  ylim(min(xi.smooth$evgam.1[((((i-1)*n)+1):(i*n))])-0.5, max(xi.smooth$evgam.1[((((i-1)*n)+1):(i*n))]) + 0.5) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=min(xi.smooth$evgam.1[((((i-1)*n)+1):(i*n))]-0.5), color = "red", size = 7)
}
grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)



alpha.smooth <- data.frame("x"= as.vector(xholder.mat),
                          # "fit" = smooth.func,
                          "true" = as.vector(as.matrix(fwi.scaled)),
                          "evgam.1" = as.vector(alpha.nonlinear.1),
                          "evgam.scale" = as.vector(alpha.nonlinear.scale),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
                          "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(alpha.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_line(aes(y=evgam.1), colour = "purple", linewidth=1) + 
                  geom_line(aes(y=evgam.scale), colour = "orange", linewidth=1) +                   
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  geom_rug(aes(x=true, y=evgam.1), sides = "b") + 
                  # scale_color_manual(values=c("purple")) +
                  ylim(min(alpha.smooth$evgam.scale[((((i-1)*n)+1):(i*n))])-0.5, max(alpha.smooth$evgam.scale[((((i-1)*n)+1):(i*n))]) + 0.5) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=min(alpha.smooth$evgam.1[((((i-1)*n)+1):(i*n))]-0.5), color = "red", size = 7)
}
grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)


xi.scenario <- data.frame("x" = xholder[,1],
                           "evgam.1" = xi.pred.1,
                           "evgam.scale" = xi.pred.scale)#, 
                            # "post.mean" = (1/alpha.samples[,1]),
                            # "post.median" = (1/alpha.samples[,5]),
                            # "q1" = (1/alpha.samples[,4]),
                            # "q3" = (1/alpha.samples[,6]))

ggplot(xi.scenario, aes(x=x)) + 
  ylab(expression(xi(c,ldots,c))) + xlab(expression(c)) + labs(col = "") +  
  # geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  geom_line(aes(y=evgam.1), linewidth=1, color = "purple") +
  geom_line(aes(y=evgam.scale), linewidth=1, color = "orange") +
  # scale_fill_manual(values=c("steelblue"), name = "") +
  # scale_color_manual(values = c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))


alpha.scenario <- data.frame("x" = xholder[,1],
                              "evgam.1" = alpha.pred.1,
                              "evgam.scale" = alpha.pred.scale)#,
                            # "post.mean" = (alpha.samples[,1]),
                            # "post.median" = (alpha.samples[,5]),
                            # "q1" = (alpha.samples[,4]),
                            # "q3" = (alpha.samples[,6]))

ggplot(alpha.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,...,c))) + xlab(expression(c)) + labs(col = "") +
  # geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  geom_line(aes(y=evgam.1), linewidth=1, color = "purple") +
  geom_line(aes(y=evgam.scale), linewidth=1, color = "orange") +
  # scale_fill_manual(values=c("steelblue"), name = "") +
  # scale_color_manual(values = c("steelblue")) + ylim(0, 15) +
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))



theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary
yrep <- summary(fit1, par=c("yrep"), probs = c(0.05,0.5, 0.95))$summary
f.samples <- summary(fit1, par=c("f"), probs = c(0.05,0.5, 0.95))$summary

g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))

data.scenario <- data.frame("x" = xholder[,1],
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

logLik1 <- logLik(m.gpd)
# Extract number of parameters (degrees of freedom)
df1 <- attr(logLik1, "df")

# Calculate AIC for each model
AIC1 <- -2 * as.numeric(logLik1) + 2 * df1
# Calculate BIC for each model
BIC1 <- -2 * as.numeric(logLik1) + log(n) * df1

load("./BLAST/application/evgam_fit_all.Rdata")
logLik2 <- logLik(m.gpd)
df2 <- attr(logLik2, "df")
AIC2 <- -2 * as.numeric(logLik2) + 2 * df2
BIC2 <- -2 * as.numeric(logLik2) + log(n) * df2

# Likelihood ratio test for nested models (model2 nested in model1)
LR_stat <- 2 * (as.numeric(logLik2) - as.numeric(logLik1))
df_diff <- df2 - df1
p_value <- pchisq(LR_stat, df = df_diff, lower.tail = FALSE)

# Print comparison
cat("Model 1: AIC =", AIC1, "BIC =", BIC1, "\n")
cat("Model 2: AIC =", AIC2, "BIC =", BIC2, "\n")
cat("Likelihood Ratio Test Statistic:", LR_stat, "on", df_diff, "df\n")
cat("p-value:", p_value, "\n")
