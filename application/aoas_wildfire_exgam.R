library(npreg)
library(Pareto)
suppressMessages(library(tidyverse))
library(readxl)
library(gridExtra)
library(evgam)


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
# psi <- 30
# n <- dim(fwi.scaled)[[1]]
# p <- dim(fwi.scaled)[[2]]
gpd.formula <- list(
  BA ~ #1,
       s(DSR, bs = "tp", k = 30) +
       s(FWI, bs = "tp", k = 30) + 
       s(BUI, bs = "tp", k = 30) + 
       s(ISI, bs = "tp", k = 30) + 
       s(FFMC, bs = "tp", k = 30) + 
       s(DMC, bs = "tp", k = 30) + 
       s(DC, bs = "tp", k = 30),
     ~ s(DSR, bs = "tp", k = 30) +
       s(FWI, bs = "tp", k = 30) + 
       s(BUI, bs = "tp", k = 30) + 
       s(ISI, bs = "tp", k = 30) + 
       s(FFMC, bs = "tp", k = 30) + 
       s(DMC, bs = "tp", k = 30) + 
       s(DC, bs = "tp", k = 30)
)


fwi.origin <- data.frame(fwi.origin[which(Y>qu),], BA=y)
m.gpd <- evgam(gpd.formula, data = fwi.origin, family = "gpd")
save(m.gpd, file = "./BLAST/application/evgam_fit_all.Rdata")
# load("./BLAST/application/evgam_fit_all.Rdata")

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


xholder.basis <-predict(m.gpd, newdata = xholder, type = "lpmatrix")

xi.coef <- tail(m.gpd$coefficients, 204)
gamma.xi <- matrix(xi.coef[2:204], ncol = p)
f.nonlinear.new <- matrix(, nrow = n, ncol = p)
bs.nonlinear <- xholder.basis$shape[,c(2:204)]
for(j in 1:p){
  f.nonlinear.new[,j] <- bs.nonlinear[,(((j-1)*29)+1):(((j-1)*29)+29)] %*% gamma.xi[,j]
}

smooth.func <- as.data.frame(m.gpd$plotdata)
smooth.func <- as.matrix(smooth.func)
smooth.func <- as.vector(smooth.func)

equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

xi.smooth <- data.frame("x"= as.vector(xholder.mat),
                          "fit" = smooth.func,
                          "true" = as.vector(as.matrix(fwi.scaled)),
                          "smooth" = as.vector(f.nonlinear.new),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
                          "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(xi.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_line(aes(y=smooth, colour = "Posterior Median"), linewidth=1) + 
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  geom_rug(aes(x=true, y=smooth), sides = "b") + 
                  scale_color_manual(values=c("purple")) +
                  ylim(-1, max(xi.smooth$smooth)) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=-1, color = "red", size = 7)
}
grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

xholder.df <- data.frame(xholder)
names(xholder.df) <- names(fwi.scaled)
xi.pred <-predict(m.gpd, newdata = xholder.df, type="link")$shape
alpha.pred <- 1/xi.pred

xholder.basis <- predict(m.gpd, newdata = xholder.df, type= "lpmatrix")$shape
xi.coef <- tail(m.gpd$coefficients, (psi-1)*p)
gamma.xi <- matrix(xi.coef, ncol = p)
alpha.nonlinear.new <- xi.nonlinear.new <- matrix(, nrow = n, ncol = p)
bs.nonlinear <- xholder.basis[,c(2:((psi-1)*p+1))]
for(j in 1:p){
  xi.nonlinear.new[,j] <- bs.nonlinear[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi[,j]
  alpha.nonlinear.new[,j] <- 1/xi.nonlinear.new[,j]
}


xi.scenario <- data.frame("x" = xholder[,1],
                           "evgam" = xi.pred)#, 
                            # "post.mean" = (1/alpha.samples[,1]),
                            # "post.median" = (1/alpha.samples[,5]),
                            # "q1" = (1/alpha.samples[,4]),
                            # "q3" = (1/alpha.samples[,6]))

ggplot(xi.scenario, aes(x=x)) + 
  ylab(expression(xi(c,...,c))) + xlab(expression(c)) + labs(col = "") +  
  # geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  geom_line(aes(y=evgam), linewidth=1, color = "purple") +
  # scale_fill_manual(values=c("steelblue"), name = "") +
  # scale_color_manual(values = c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))


alpha.scenario <- data.frame("x" = xholder[,1],
                              "evgam" = alpha.pred)#,
                            # "post.mean" = (alpha.samples[,1]),
                            # "post.median" = (alpha.samples[,5]),
                            # "q1" = (alpha.samples[,4]),
                            # "q3" = (alpha.samples[,6]))

ggplot(alpha.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,...,c))) + xlab(expression(c)) + labs(col = "") +
  # geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  geom_line(aes(y=evgam), linewidth=1, color = "purple") +
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




# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_smooth.pdf"), width=12.5, height = 15)

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