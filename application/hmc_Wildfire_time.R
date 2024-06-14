library(npreg)
library(Pareto)
suppressMessages(library(tidyverse))
library(readxl)
# library(gridExtra)
library(patchwork)
library(corrplot)
library(rstan)
library(loo)
library(qqboxplot)
library(ggdensity)
library(ggforce)
library(ggdist)
library(splines)
# source("C:/Users/Johnny Lee/Documents/Github/extremis/R/bltir.R")

options(mc.cores = parallel::detectCores())

# Structure of the FWI System
#DSR : Dail Severity Rating
#FWI : Fire Weather Index
#BUI : Buildup Index
#ISI : Initial Spread Index
#FFMC : Fine FUel Moisture Code
#DMC : Duff Moisture Code
#DC : Drought Code


setwd("C:/Users/Johnny Lee/Documents/GitHub")
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

for(i in 1:length(cov)){
    cov.long <- gather(cov[[i]][,1:41], condition, measurement, "1980":"2019", factor_key=TRUE)
}


time <- as.Date(paste0(cov.long$condition[missing.values], substr(cov.long$...1[missing.values],5,10)), "%Y-%m-%d")


u <- quantile(Y, threshold)
y <- Y[Y>u]
t <-  as.numeric(time) / 10000
t <- t[Y>u]
x <- cbind(rep(1, length(y)), bs(t))
p <- ncol(x)
n <- nrow(x)

model.stan <- "data {
  int <lower=1> n; // Sample size
  int <lower=1> p; // regression coefficient size
  real <lower=0> u; // large threshold value
  
  vector[n] y; // extreme response
  matrix[n, p] x;
}

parameters {
  vector[p] beta;
}

transformed parameters {
  vector[n] invxi;

  invxi = 1/exp(x * beta);
}

model {
  target += pareto_lpdf(y | u, invxi);
  target += normal_lpdf(beta | 0, 1e2);
}
"

data.stan <- list(y = as.vector(y), u = u, p = p, n= n, x = x)

fit <- stan(
    model_code = model.stan,  # Stan program
    data = data.stan,    # named list of data
    init = "random",      # initial value
    # init_r = 1,
    chains = 3,             # number of Markov chains
    # warmup = 1000,          # number of warmup iterations per chain
    iter = 5000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 500           # no progress shown
)

beta.samples <- summary(fit, par=c("beta"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit, par=c("invxi"), probs = c(0.05,0.5, 0.95))$summary

posterior <- extract(fit)

data.scenario <- data.frame("x" = time[Y>u],
                            "post.mean" = (alpha.samples[,1]),
                            "post.median" = (alpha.samples[,5]),
                            "q1" = (alpha.samples[,4]),
                            "q3" = (alpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(t))) + xlab("Time") + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y = true, col = "True"), linewidth = 2) +
  # xlim(-1,1) + #ylim(0, 6.2) + 
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue")) + 
  # scale_y_log10() + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))

# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_alpha.pdf"), width=10, height = 7.78)


len <- dim(posterior$invxi)[1]
r <- matrix(, nrow = n, ncol = 100)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 100
for(i in 1:n){
  for(t in 1:T){
    # r[i, t] <- qnorm(pPareto(y[i], u, alpha = alp.x.samples[i,5]))
    r[i, t] <- qnorm(pPareto(y[i], u, alpha = posterior$invxi[round(runif(1,1,len)),i]))
    # r[i, t] <- qnorm(pPareto(y[i], u, alpha = posterior$newalpha[round(runif(1,1,len)),i]))
  }
}
lgrid <- n
grid <- qnorm(ppoints(lgrid))
# qqnorm(r[, 1])
# points(grid, quantile(r[, 1], ppoints(lgrid), type = 2), 
#     xlim = c(-3, 3), col = "red")
traj <- matrix(NA, nrow = T, ncol = lgrid)
for (t in 1:T){
  traj[t, ] <- quantile(r[, t], ppoints(lgrid), type = 2)
}
l.band <- apply(traj, 2, quantile, prob = 0.025)
trajhat <- apply(traj, 2, quantile, prob = 0.5)
u.band <- apply(traj, 2, quantile, prob = 0.975)

ggplot(data = data.frame(grid = grid, l.band = l.band, trajhat = trajhat, 
                         u.band = u.band)) + 
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
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_qqplot.pdf"), width=10, height = 7.78)
rp <-c()
for(i in 1:n){
  # rp[i] <- rPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i])
  rp[i] <- qnorm(pPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i]))
}
rp <- data.frame(rp, group = rep("residuals", n))

ggplot(data = rp) + 
  # geom_qqboxplot(aes(factor(group, levels=c("residuals")), y=rp), notch=FALSE, varwidth=TRUE, reference_dist="norm")+ 
  geom_qqboxplot(aes(y=rp), notch=FALSE, varwidth=FALSE, reference_dist="norm", width = 0.15, qq.colour="steelblue")+
  labs(x = "", y = "Residuals") + ylim(-3,3) + xlim(-.2,.2)+
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_qqboxplot.pdf"), width = 10, height = 7.78)
             
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))

cat("Finished Running")
gl.samples <- readRDS("./BRSTIR/application/2024-06-13_linear_stanfit.rds")
gl.mean <- as.vector(matrix(gl.samples[,1], nrow = n, byrow=TRUE))
# gl.q1 <- as.vector(matrix(gl.samples[,4], nrow = n, byrow=TRUE))
# gl.q2 <- as.vector(matrix(gl.samples[,5], nrow = n, byrow=TRUE))
# gl.q3 <- as.vector(matrix(gl.samples[,6], nrow = n, byrow=TRUE))


data.smooth <- data.frame("x" = as.vector(xholder),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "gl.mean" = as.vector(gl.mean),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_line(aes(y=gl.mean, colour = "linear"), linewidth=1, linetype=2) + 
                  geom_rug(aes(x=fwi.scaled[,i], y=q2), sides = "b") +
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("red","steelblue")) +
                  ylim(-7, 7) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=-7, color = "red", size = 7)
}
grid.plts[[1]] + grid.plts[[2]] +grid.plts[[3]] + grid.plts[[4]] + grid.plts[[5]] + grid.plts[[6]] + grid.plts[[7]] + plot_layout(widths = c(1,1))
# grid.arrange(grobs = grid.plts, ncol = 2, nrow = 4)
# grid.plts[[7]]
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_smooth.pdf"), width=10, height = 7.78)

#Predictive Distribution check
y.container <- as.data.frame(matrix(, nrow = n, ncol = 0))  
random.alpha.idx <- floor(runif(1000, 1, ncol(t(posterior$f))))
for(i in random.alpha.idx){
  # y.container <- cbind(y.container, log(t(posterior$f)[,i]))
  y.container <- cbind(y.container, t(posterior$f)[,i])
}

colnames(y.container) <- paste("col", 1:length(random.alpha.idx), sep="")
y.container$x <- seq(1,n)
y.container$logy <- seq(log(u), 30, length.out = n)
y.container$y <- y
y.container <- cbind(y.container, t(apply(y.container[,1:length(random.alpha.idx)], 1, quantile, c(0.05, .5, .95))))
colnames(y.container)[(dim(y.container)[2]-2):(dim(y.container)[2])] <- c("q1","q2","q3")
plt <- ggplot(data = y.container, aes(x = logy)) + ylab("Density") + xlab("log(Burned Area)") + labs(col = "")

for(i in names(y.container)[1:length(random.alpha.idx)]){
  # plt <- plt + geom_density(aes(x=.data[[i]]), color = "slategray1", alpha = 0.1, linewidht = 0.7)
  plt <- plt + geom_line(aes(x=logy, y=.data[[i]]), color = "slategray1", alpha = 0.1, linewidht = 0.7)  
}

print(plt + geom_line(aes(x=logy, y=q2), color = "steelblue", linewidth = 2) +
        geom_ribbon(aes(ymin = q1, ymax = q3), fill = "steelblue", alpha = 0.2) +
        # scale_y_log10() +
        theme_minimal(base_size = 30) + #ylim(0, 1) + #xlim(7.5,30) +
        theme(legend.position = "none",
                axis.text = element_text(size = 35)))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_BRSTIR_predictive_distribution.pdf"), width=10, height = 7.78)


# extreme.container <- as.data.frame(matrix(, nrow = n, ncol = 3000))
# for(i in 1:3000){
#   extreme.container[,i] <- density(log(posterior$f[i,]), n=n)$y
# }
# extreme.container <- cbind(extreme.container, t(apply(extreme.container[,1:3000], 1, quantile, c(0.05, .5, .95))))
# colnames(extreme.container)[(dim(extreme.container)[2]-2):(dim(extreme.container)[2])] <- c("q1","q2","q3")
# colnames(extreme.container)[1:3000] <- as.character(1:3000)
# extreme.container$mean <- rowMeans(extreme.container[,1:3000])
# extreme.container$y <- seq(0, 30, length.out = n)
# extreme.container <- as.data.frame(extreme.container)


# plt <- ggplot(data = extreme.container, aes(x = y)) + xlab("log(Burned Area)") + ylab("Density")+
#         # geom_line(aes(y=q2), colour = "steelblue", linewidth = 1.5) +
#         geom_line(aes(y=mean), colour = "steelblue", linewidth = 1.5) +
#         geom_ribbon(aes(ymin = q1, ymax = q3), fill = "steelblue", alpha = 0.2) + 
#         theme_minimal(base_size = 30) + 
#         theme(legend.position = "none",
#               axis.title = element_text(size = 30))
# d <- ggplot_build(plt)$data[[1]]
# print(plt + 
#         # geom_area(data = subset(d, x>12.44009), aes(x=x,y=y), fill = "slategray1", alpha = 0.5) +
#         geom_segment(x=12.44009, xend=12.44009, 
#               y=0, yend=approx(x = d$x, y = d$y, xout = 12.4409)$y,
              # colour="red", linewidth=1.2, linetype = "dotted"))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_post_generative.pdf"), width = 10, height = 7.78)

random.alpha.idx <- floor(runif(1000, 1, ncol(t(posterior$alpha))))
ev.y1 <- ev.y2 <- as.data.frame(matrix(, nrow = 1, ncol = 0))
# for(i ?in 1:ncol(t(posterior$theta))){
ev.alpha.single <- c()  
for(i in random.alpha.idx){
  ev.y1 <- rbind(ev.y1, as.numeric(posterior$yrep[i]))
  ev.y2 <- rbind(ev.y2, as.numeric(posterior$yrep[i]))
}
ev.y1 <- as.data.frame(log(ev.y1))
ev.y1$logy <- max(log(y))
colnames(ev.y1) <- c("yrep", "logy")
ev.y1$group <- rep("15th Oct 2017",1000)

plt <- ggplot(data = ev.y1, aes(x = yrep)) + ylab("Density") + xlab("log(Burned Area)") + labs(col = "") +
  geom_density(color = "steelblue", linewidth = 1.2) + 
  geom_rug(alpha = 0.1) + 
  # geom_point(aes(x=yrep,y=-Inf),color="steelblue", size = 3.5, alpha = 0.2) +
  xlim(5.5, 40) +
  theme_minimal(base_size = 30) +  
  theme(legend.position = "none",
        axis.title = element_text(size = 30))
  # geom_vline(xintercept = log(max(y)), linetype="dotted", color = "red",) +

d <- ggplot_build(plt)$data[[1]]
print(plt + geom_area(data = subset(d, x>12.44009), aes(x=x,y=y), fill = "slategray1", alpha = 0.5) +
        geom_segment(x=12.44009, xend=12.44009, 
              y=0, yend=approx(x = d$x, y = d$y, xout = 12.4409)$y,
              colour="red", linewidth=1.2, linetype = "dotted"))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_BRSTIR_generative.pdf"), width = 10, height = 7.78)

data.yrep <- data.frame("y" = seq(log(u), 20, length.out=n),
                          "post.mean" = f.samples[,1],
                          "q1" = f.samples[,4],
                          "q2" = f.samples[,5],
                          "q3" = f.samples[,6])

ggplot(data.yrep, aes(x=y)) + 
  ylab("Density") + xlab("log(Burned Area)") + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y = true, col = "True"), linewidth = 2) +
  xlim(7.5, 20) + ylim(0, 7e-4) + 
  geom_line(aes(y=q2, col = "Posterior Median"), linewidth=1.5) +
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue")) + 
  annotate("point", x= log(max(Y)), y=0, color = "red", size = 5) +
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))
ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_post_generative.pdf"), width = 10, height = 7.78)        


density.y <- density(ev.y1$yrep) # see ?density for parameters
# set an Avg.position value
Avg.pos <- log(max(y))
xt <- diff(density.y$x[density.y$x>Avg.pos])
library(zoo)
yt <- rollmean(density.y$y[density.y$x>Avg.pos],2)
# This gives you the area
sum(xt*yt)


fit.log.lik <- extract_log_lik(fit1)
loo(fit.log.lik, is_method = "sis", cores = 2)
fwi.loo <- loo(fit.log.lik, cores = 2)
plot(fwi.loo, label_points = TRUE)

constraint.elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 2)
constraint.waic <- waic(fit.log.lik, cores = 2)
# save(constraint.elpd.loo, constraint.waic, file = (paste0("./BRSTIR/application/BRSTIR_constraint_",Sys.Date(),"_",psi,"_",floor(threshold*100),"quantile_IC.Rdata")))




