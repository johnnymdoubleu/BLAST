library(npreg)
library(Pareto)
suppressMessages(library(tidyverse))
library(JOPS)
library(readxl)
library(gridExtra)
library(colorspace)
library(corrplot)
library(ReIns)
library(rstan)
library(loo)
library(bayesplot)
library(evir)
library(mev)
library(cmdstanr)
library(scales)
# library(ggExtra)


options(mc.cores = parallel::detectCores())

# Structure of the FWI System
#DSR : Dail Severity Rating
#FWI : Fire Weather Index
#BUI : Buildup Index
#ISI : Initial Spread Index
#FFMC : Fine FUel Moisture Code
#DMC : Duff Moisture Code
#DC : Drough Code


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
summary(Y) #total burnt area
length(Y)
psi <- 10
threshold <- 0.945
u <- quantile(Y, threshold)
y <- Y[Y>u]
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
# fwi.index$date <- 
fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)
fwi.scaled <- fwi.scaled[which(Y>u),]
# fwi.scaled <- as.data.frame(scale(fwi.scaled))

# plot((fwi.scaled[,2]), (log(y)))
# plot((fwi.scaled[,5]), (log(y)))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled <- as.data.frame(sapply(fwi.scaled, FUN = range01))
# fwi.scaled <- as.data.frame(lapply(fwi.scaled, rescale, to=c(-1,1)))

# ---------------------------------------------------------------------------
# Computing Hills Estimator plot
# orderred <- rev(sort(Y)[14462:14609])
# orderred <- rev(sort(Y)[13863:14609])
# ordered <- rev(sort(Y))
# n.hill <- length(orderred)
# k <- 1:n.hill
# loggs <- logb(orderred/u)
# avesumlog <- cumsum(loggs)/k
# # xihat <- c(NA, (avesumlog)[2:n.hill])
# xihat <- c(NA, (avesumlog-loggs)[2:n.hill])
# alphahat <- 1/xihat
# ses <- alphahat/sqrt(k)
# xx <- trunc(seq(from = n.hill, to = 15))
# y.alpha <- alphahat[xx]
# # ylabel <- alphahat
# yrange <- range(y.alpha)
# qq <- qnorm(1-(1-threshold)/2)
# uu <- y.alpha + ses[xx] * qq
# ll <- y.alpha - ses[xx] * qq
# yrange <- range(u, l)
# data.hill <- data.frame(k = c(15:n.hill),
#                         u = uu,
#                         l = ll,
#                         alpha = y.alpha,
#                         order = xx)
# ggplot(data = data.hill) + 
#   geom_ribbon(aes(x = order, ymin = l, ymax = u, fill = "confidenceband"),
#               alpha = 0.2, linetype = "dashed") + 
#   geom_line(aes(x = order, y = alpha, color="hillestimator" ), linewidth = 1.2) + 
#   xlim(15, n.hill) +
#   labs(x = "Order Statistics", y = "Tail Index") + 
#   scale_color_manual(values=c("steelblue")) + 
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   theme_minimal(base_size = 30) +
#   theme(text = element_text(size = 30), 
#         axis.text.x = element_text(angle = 0, hjust = 0.5),
#         legend.position = "none")
# ggsave("./BRSTIR/application/figures/hillestimator.pdf", width=10, height = 7.78)


# pdf(file = "./BRSTIR/application/figures/correlation.pdf")
# corrplot.mixed(cor(fwi.scaled),
#                 upper = "circle",
#                 lower = "number",
#                 addgrid.col = "black")
# dev.off()
# ggsave("./BRSTIR/application/figures/correlation.pdf", plot = replayPlot(p1), width=10, height = 7.78)
# -------------------------------------------------------------------

# ------------- Explanatory Analaysis
# first.extreme <- which(Y==max(y))
# second.extreme <- which(Y==max(y[-which.max(y)]))
# tenth.extreme <- which(Y==sort(y, decreasing = TRUE)[10])
# ggplot(fwi.index[((first.extreme):(first.extreme+12)),], aes(x=date)) +
#   geom_line(aes(y=DSR, color = "DSR"), linetype = 1) + 
#   geom_line(aes(y=FWI, color = "FWI"), linetype = 2) +
#   geom_line(aes(y=BUI, color = "BUI"), linetype = 3) +
#   geom_line(aes(y=ISI, color = "ISI"), linetype = 4) +
#   geom_line(aes(y=FFMC, color = "FFMC"), linetype = 5) + 
#   geom_line(aes(y=DMC, color = "DMC"), linetype = 6) +
#   geom_line(aes(y=DC, color = "DC"), linetype = 7)  + 
#   ylab("indices") + xlab("dates after extreme fire (sorted by burnt area)") + 
#   scale_color_manual(name = "Indices", values = c(
#     "DSR" = "darkblue", 
#     "FWI" = "red",
#     "BUI" = "green",
#     "ISI" = "yellow",
#     "FFMC" = "orange",
#     "DMC" = "purple",
#     "DC" = "skyblue")) +
#   theme(legend.position="right", 
#       legend.key.size = unit(1, 'cm'),
#       legend.text = element_text(size=20),
#       # plot.margin = margin(0,0,0,-1),
#       axis.title = element_text(size = 20))

# fwi.index[((first.extreme):(first.extreme+12)),]
# fwi.index[13682:13694,]
# fwi.index[second.extreme:(second.extreme+12),]

df.extreme <- cbind(y, fwi.scaled)
df.extreme <- as.data.frame(cbind(month = fwi.index$month[which(Y>u)], df.extreme))
# ggplot(df.extreme, aes(x=month, y=y, color=month)) + geom_point(size=6) + theme_minimal() +
#     theme(plot.title = element_text(hjust = 0.5, size = 20),
#         legend.title = element_blank(),
#         legend.text = element_text(size=20),
#         # axis.ticks.x = element_blank(),
#         axis.text.x = element_text(hjust=0.35),
#         axis.text = element_text(size = 20),
#         axis.title = element_text(size = 15))
# ggsave("./BRSTIR/application/figures/datavis.pdf", width=15)


n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]
no.theta <- 1
newx <- seq(0, 1, length.out=n)
xholder.linear <- xholder.nonlinear <- bs.linear <- bs.nonlinear <- matrix(,nrow=n, ncol=0)
xholder <- matrix(nrow=n, ncol=p)
basis.holder <- matrix(, nrow = 2, ncol =0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                      matrix(c(which.min(fwi.scaled[,i]),
                              which.max(fwi.scaled[,i])), ncol=2))
}
for(i in 1:p){
  # xholder[,i] <- seq(0, 1, length.out = n)
  # test.knot <- seq(0, 1, length.out = psi)
  # splines <- basis.tps(seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n), test.knot, m=2, rk=FALSE, intercept = TRUE)
  xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
  test.knot <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  knots <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  tps <- basis.tps(fwi.scaled[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  basis.holder <- cbind(basis.holder, 
          solve(matrix(c(splines[index.holder[i,1], no.theta+1],
                  splines[index.holder[i,1], no.theta+psi],
                  splines[index.holder[i,2], no.theta+1],
                  splines[index.holder[i,2], no.theta+psi]), 
                  nrow = 2, ncol = 2)))
  bs.linear <- cbind(bs.linear, tps[,1:no.theta])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
}




write("data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp; 
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector[n] y; // extreme response
    real <lower=0> atau;
    matrix[2, (2*p)] basisFL;
    array[p, 2] int indexFL;
}
parameters {
    vector[newp] theta; // linear predictor
    vector[(psi-2)] gamma[p]; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real sigma; //
    vector[p] tau; //real <lower=0,upper=100> rho //softplus parameter
}
transformed parameters {
    array[n] real <lower=0> alpha; // tail index
    vector[(psi-2)] gammaTemp[p]; // constraint gamma from 2 to psi-1
    vector[2] gammaFL[p];
    matrix[n, p] gnl; // nonlinear component
    matrix[2, p] subgnl;
    matrix[n, p] gfnl;
    matrix[n, p] glnl;
    matrix[n, p] gl; // linear component
    matrix[n, p] gsmooth; // linear component
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgnl; // nonlinear component
    matrix[n, p] newgfnl;
    matrix[n, p] newglnl;    
    matrix[n, p] newgl; // linear component
    matrix[n, p] newgsmooth; // linear component

    for(j in 1:p){
        gammaTemp[j] = gamma[j];
        subgnl[,j] = bsNonlinear[indexFL[1,], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j];
        gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * -1 * subgnl[,j];
    };

    for (j in 1:p){
        gnl[,j] = bsNonlinear[,(((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j]; 
        newgnl[,j] = xholderNonlinear[,(((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j];
        gl[,j] = bsLinear[,j] * theta[j+1]; 
        newgl[,j] = xholderLinear[,j] * theta[j+1];
        gfnl[,j] = bsNonlinear[,(((j-1)*psi)+1)] * gammaFL[j][1];
        glnl[,j] = bsNonlinear[,(((j-1)*psi)+psi)] * gammaFL[j][2];
        newgfnl[,j] = xholderNonlinear[,(((j-1)*psi)+1)] * gammaFL[j][1];
        newglnl[,j] = xholderNonlinear[,(((j-1)*psi)+psi)] * gammaFL[j][2];        
        gsmooth[,j] = gl[,j] + gnl[,j] + gfnl[,j] + glnl[,j]; 
        newgsmooth[,j] = newgl[,j] + newgnl[,j] + newgfnl[,j] + newglnl[,j];
    };

    for (i in 1:n){
        alpha[i] = exp(theta[1] + sum(gsmooth[i,])); 
        newalpha[i] = exp(theta[1] + sum(newgsmooth[i,]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 1, 10);
    target += gamma_lpdf(lambda2 | 0.1, 0.1);
    target += normal_lpdf(theta[1] | 0, 1);
    // target += normal_lpdf(theta[1] | log(1.2), 0.01) target += double_exponential_lpdf(theta[1] | 0, lambda1)
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    target += ((newp * log(lambda1)) + (p * (psi-2) * log(lambda2)));
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, lambda2/sqrt(2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi-2), diag_matrix(rep_vector(1, psi-2)) * sqrt(tau[j]) * sqrt(sigma));
    }
}
generated quantities {
    // Used in Posterior predictive check
    vector[n] log_lik;
    for (i in 1:n) {
        log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}
"
, "model_BRSTIR_constraint.stan")

data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi-1)/2), newp = (p+1), 
                    indexFL = index.holder, 
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear, basisFL = basis.holder)

set_cmdstan_path(path = NULL)
#> CmdStan path set to: /Users/jgabry/.cmdstan/cmdstan-2.32.2

# Create a CmdStanModel object from a Stan program,
# here using the example model that comes with CmdStan
file <- file.path(cmdstan_path(), "model.stan")

init.alpha <- list(list(gamma = array(rep(0, (psi*p)), dim=c((psi-2), p)),
                        theta = rep(0, (p+1)),
                        tau = rep(0.1, p), sigma = 0.1, 
                        lambda1 = 0.01, lambda2 = 0.01),
                  list(gamma = array(rep(0.02, (psi*p)), dim=c((psi-2), p)),
                        theta = rep(0.01, (p+1)),
                        tau = rep(0.01, p), sigma = 0.001,
                        lambda1 = 0.1, lambda2 = 0.001),
                  list(gamma = array(rep(0.01, (psi*p)), dim=c((psi-2), p)),
                        theta = rep(0.05, (p+1)),
                        tau = rep(0.01, p), sigma = 0.01,
                        lambda1 = 1, lambda2 = 0.01))

# stanc("C:/Users/Johnny Lee/Documents/GitHub/BRSTIR/application/model1.stan")
fit1 <- stan(
    file = "model_BRSTIR_constraint.stan",  # Stan program
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    # init_r = 1,
    chains = 3,             # number of Markov chains
    # warmup = 1000,          # number of warmup iterations per chain
    iter = 3000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 500           # no progress shown
)

# saveRDS(fit1, file=paste0("./BRSTIR/application/",Sys.Date(),"_stanfit.rds"))
posterior <- extract(fit1)
# str(posterior)

# print(as.mcmc(fit1), pars=c("alpha", "gamma", "intercept", "theta", "lambda1", "lambda2","lp__"), probs=c(.05,.5,.95))
# plot(fit1, plotfun = "trace", pars = c("theta"), nrow = 3)
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_theta_trace.pdf"), width=10, height = 7.78)
# plot(fit1, plotfun = "trace", pars = c("lambda1", "lambda2"), nrow = 2)
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_lambda.pdf"), width=10, height = 7.78)

# traceplot(fit1, pars = c("theta"))
# traceplot(fit1, pars = c("lambda1", "lambda2"), inc_warmup = TRUE, nrow = 2)
# fit.v2 <- as.mcmc(fit1)

# alpha.summary <- fit.v2$summary$all.chains

# alpha.summary[701:711,]

# MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
#             HPD = TRUE, xlab="theta", offset = 0.05, exact = TRUE,
#             horiz = FALSE, params = c("theta0", "theta"))
# MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
#             HPD = TRUE, xlab="gamma", offset = 0.5,
#             horiz = FALSE, params = c("gamma"))
# gg.fit <- ggs(fit.v2$samples)
# lambda.p1 <- gg.fit %>% filter(Parameter == c("lambda.1", "lambda.2")) %>% 
#   ggs_traceplot() + theme_minimal(base_size = 20) + theme(, legend.position = "none")
# lambda.p2 <- gg.fit %>% filter(Parameter == c("lambda.1", "lambda.2")) %>% 
#   ggs_density() + theme_minimal(base_size = 20)
# grid.arrange(lambda.p1, lambda.p2, ncol=2)
# # MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
# #             HPD = TRUE, xlab="lambda", offset = 0.5,
# #             horiz = FALSE, params = c("lambda.1", "lambda.2"))            
# # print(alpha.summary)
# MCMCsummary(object = fit.v2$samples, round = 3)[((dim(alpha.summary)[1]-p-2-(psi*p)):dim(alpha.summary)[1]),]

# print(MCMCtrace(object = fit.v2$samples,
#           pdf = FALSE, # no export to PDF
#           ind = TRUE, # separate density lines per chain
#           n.eff = TRUE,# add eff sample size
#           params = c("gamma")))
# print(MCMCtrace(object = fit.v2$samples,
#           pdf = FALSE, # no export to PDF
#           ind = TRUE, # separate density lines per chain
#           n.eff = TRUE,# add eff sample size
#           params = c("lambda.1", "lambda.2")))

# samples <- fit.v2$samples$chain1
# len <- dim(samples)[1]

theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
# gl.samples <- summary(fit1, par=c("newgl"), probs = c(0.05, 0.5, 0.95))$summary
# gnl.samples <- summary(fit1, par=c("newgnl"), probs = c(0.05, 0.5, 0.95))$summary
gnlfl.samples <- summary(fit1, par=c("newgfnl", "newglnl"), probs = c(0.05, 0.5, 0.95))$summary

gsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
# smooth.samples <- summary(fit1,par=c("gsmooth"), probs = c(0.05, 0.5, 0.95))$summary
alp.x.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary
# summary(fit1, par=c("sigma"), probs = c(0.05,0.5, 0.95))$summary
# summary(fit1, par=c("tau"), probs = c(0.05,0.5, 0.95))$summary


gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q2 <- gamma.samples[,5]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q2 <- theta.samples[,5]
theta.q3 <- theta.samples[,6]


df.theta <- data.frame("seq" = seq(1, (p+1)),
                        "m" = c(theta.q2),
                        "l" = c(theta.q1),
                        "u" = c(theta.q3))
df.theta$covariate <- factor(c("\u03b8",names(fwi.scaled)), levels = c("\u03b8",colnames(fwi.scaled)))
df.theta$labels <- factor(c("\u03b8",colnames(fwi.scaled)))

ggplot(df.theta, aes(x = covariate, y=m, color = covariate)) + ylab("") + xlab('') +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_point(size = 5) + 
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, linewidth =1.2) + 
  scale_x_discrete(labels = c(expression(bold(theta[0])),
                              expression(bold(theta[1])),
                              expression(bold(theta[2])),
                              expression(bold(theta[3])),
                              expression(bold(theta[4])),
                              expression(bold(theta[5])),
                              expression(bold(theta[6])),
                              expression(bold(theta[7])))) + 
  scale_color_discrete(labels = c(expression(theta[0]),colnames(fwi.scaled))) + 
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.text = element_text(size=25),
          legend.margin=margin(0,0,0,-10),
          legend.box.margin=margin(-10,0,-10,0),
          plot.margin = margin(0,0,0,-20),
          axis.text.x = element_text(hjust=0.35),
          axis.text = element_text(size = 28))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_theta.pdf"), width=10, height = 7.78)

# ggplot(data.frame(group = factor(1:(p+1)), m=theta.post.mean, l = theta.q1, u = theta.q3), 
#        aes(group)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1.4) +
#   geom_point(aes(x = group, y = m), size = 4.5) + 
#   #geom_point(aes(x = group, y = beta), shape=8, size = 4.5, col="red")+
#   geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, size = 1.2) + 
#   labs(x = "Regression coefficients", y = "") + 
#   ylim(-5,5) + 
#   scale_x_discrete(labels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC"))+
#                               #,expression(beta[8]),
#                               #expression(beta[9]))) + 
#   theme_minimal(base_size = 30) + 
#   theme(text = element_text(size = 30), 
#         axis.text.x = element_text(angle = 0, hjust = 0.5))


# df.gamma <- data.frame("seq" = seq(1, (psi*p)), 
#                   "m" = as.vector(gamma.q2),
#                   "l" = as.vector(gamma.q1),
#                   "u" = as.vector(gamma.q3))
# df.gamma$covariate <- factor(rep(names(fwi.scaled), each = psi, length.out = nrow(df.gamma)), levels = colnames(fwi.scaled))
# df.gamma$labels <- factor(1:(psi*p))
# ggplot(df.gamma, aes(x =labels, y = m, color = covariate)) + 
#   geom_errorbar(aes(ymin = l, ymax = u), alpha = 0.4, width = 4, linewidth = 1.2) +
#   geom_point(size = 4) + ylab("") + xlab("" ) + #xlim(1,(psi*p)) +
#   # geom_ribbon(aes(ymin = l, ymax = u)) +
#   # geom_point(size = 4, color = "black") + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+10), 
#                     label = c(expression(bold(gamma[1])), 
#                               expression(bold(gamma[2])), 
#                               expression(bold(gamma[3])), 
#                               expression(bold(gamma[4])), 
#                               expression(bold(gamma[5])), 
#                               expression(bold(gamma[6])), 
#                               expression(bold(gamma[7]))),
#                     expand=c(0,10)) +
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#           legend.title = element_blank(),
#           legend.text = element_text(size=25),
#           legend.margin=margin(0,0,0,-10),
#           legend.box.margin=margin(-10,0,-10,0),
#           plot.margin = margin(0,0,0,-20),
#           axis.text.x = element_text(hjust=0.5),
#           axis.text = element_text(size = 28))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_gamma.pdf"), width=10, height = 7.78)


# g.linear.mean <- as.vector(matrix(gl.samples[,1], nrow = n, byrow=TRUE))
# g.linear.q1 <- as.vector(matrix(gl.samples[,4], nrow = n, byrow=TRUE))
# g.linear.q2 <- as.vector(matrix(gl.samples[,5], nrow = n, byrow=TRUE))
# g.linear.q3 <- as.vector(matrix(gl.samples[,6], nrow = n, byrow=TRUE))
# g.nonlinear.mean <- as.vector(matrix(gnl.samples[,1], nrow = n, byrow=TRUE))
# g.nonlinear.q1 <- as.vector(matrix(gnl.samples[,4], nrow = n, byrow=TRUE))
# g.nonlinear.q2 <- as.vector(matrix(gnl.samples[,5], nrow = n, byrow=TRUE))
# g.nonlinear.q3 <- as.vector(matrix(gnl.samples[,6], nrow = n, byrow=TRUE))
g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))

equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

data.smooth <- data.frame("x"= as.vector(xholder),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          # "post.mean" = as.vector(sort(g.smooth.mean)),
                          # "q1" = as.vector(sort(g.smooth.q1)),
                          # "q2" = as.vector(sort(g.smooth.q2)),
                          # "q3" = as.vector(sort(g.smooth.q3)),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
                          # "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab("") +
  facet_grid(covariates ~ ., scales = "free",
              labeller = label_parsed) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) + 
  # scale_y_continuous(breaks=c(-3,-2,-1,1,2)) +        
          ylim(-3.5, 2.8) +
  # scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
          plot.margin = margin(0,0,0,-20),
          # strip.text = element_blank(),
          axis.text = element_text(size = 20))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_smooth.pdf"), width=12.5, height = 15)

# data.linear <- data.frame("x"= as.vector(xholder),
#                           "post.mean" = as.vector(g.linear.mean),
#                           "q1" = as.vector(g.linear.q1),
#                           "q2" = as.vector(g.linear.q2),
#                           "q3" = as.vector(g.linear.q3),
#                           "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
#                           "fakelab" = rep(1, (p*n)),
#                           "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

# ggplot(data.linear, aes(x=x, group=interaction(covariates, replicate))) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#   # geom_line(aes(y=true, colour = "True"), linewidth=2) + 
#   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#   ylab("") + xlab("") +
#   facet_grid(covariates ~ ., scales = "free",
#               labeller = label_parsed) + 
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values=c("steelblue")) + 
#   guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) + #ylim(-0.65, 0.3) +
#   scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#           plot.margin = margin(0,0,0,-20),
#           # strip.text = element_blank(),
#           axis.text = element_text(size = 20))
# # #ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_linear.pdf"), width=12.5, height = 15)


# data.nonlinear <- data.frame("x"=as.vector(xholder),
#                           "post.mean" = as.vector(g.nonlinear.mean),
#                           "q1" = as.vector(g.nonlinear.q1),
#                           "q2" = as.vector(g.nonlinear.q2),
#                           "q3" = as.vector(g.nonlinear.q3),
#                           "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
#                           "fakelab" = rep(1, (p*n)),
#                           "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

# ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#   # geom_line(aes(y=true, colour = "True"), linewidth=2) + 
#   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#   ylab("") + xlab("") +
#   facet_grid(covariates ~ ., scales = "free", #switch = "y",
#               labeller = label_parsed) + 
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values=c("steelblue")) + 
#   guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) + #ylim(-0.65, 0.3) +
#   scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#           plot.margin = margin(0,0,0,-20),
#           # strip.text = element_blank(),
#           axis.text = element_text(size = 20))
# #ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_nonlinear.pdf"), width=12.5, height = 15)

data.scenario <- data.frame("x" = seq(0, 1, length.out = n),
                            "post.mean" = (alpha.samples[,1]),
                            "post.median" = (alpha.samples[,5]),
                            # "post.median" = sort(exp(rowSums(g.q2) + rep(theta.samples[1,5], n)), decreasing = TRUE),
                            "q1" = (alpha.samples[,4]),
                            "q3" = (alpha.samples[,6]))
# data.scenario <- data.frame("x" = seq(-1, 1, length.out = n),
#                             "post.mean" = sort(alp.x.samples[,1], decreasing = TRUE),
#                             "post.median" = sort(alp.x.samples[,5], decreasing = TRUE),
#                             "q1" = sort(alp.x.samples[,4], decreasing = TRUE),
#                             "q3" = sort(alp.x.samples[,6], decreasing = TRUE))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(bold(x)))) + xlab(expression(c)) + labs(col = "") +
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
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position="none", 
        legend.key.size = unit(1, 'cm'),
        legend.text = element_text(size=20),
        plot.margin = margin(0,0,0,-1),
        strip.text = element_blank(),
        axis.title.x = element_text(size = 35))

# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_alpha.pdf"), width=10, height = 7.78)

len <- dim(posterior$alpha)[1]
r <- matrix(, nrow = n, ncol = 100)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 100
for(i in 1:n){
  for(t in 1:T){
    # r[i, t] <- qnorm(pPareto(y[i], u, alpha = alp.x.samples[i,5]))
    r[i, t] <- qnorm(pPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i]))
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
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 30)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_qqplot.pdf"), width=10, height = 7.78)
             
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))

cat("Finished Running")

# relative_eff(exp(fit.log.lik))
#https://discourse.mc-staqan.org/t/four-questions-about-information-criteria-cross-validation-and-hmc-in-relation-to-a-manuscript-review/13841/3
# y.rep <- as.matrix(fit1, pars = "y_rep")
# ppc_loo_pit_overlay(
#   y = y,
#   yrep = y.rep,
#   lw = weights(fwi.loo$psis_object)
# )

grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),], origin = fwi.scaled[,i]), aes(x=x)) + 
  # grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),], origin = fwi.index[which(Y>u),i]), aes(x=x)) +   
                  # geom_point(aes(x= origin, y=q2), alpha = 0.3) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x= origin, y=q2), sides = "b") +
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  ylim(-2.3, 1.2) +
                  # scale_y_continuous(breaks=equal_breaks(n=5, s=0.1)) + 
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt
}

grid.arrange(grobs = grid.plts, ncol = 2, nrow = 4)
# grid.plts[[7]]
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_smooth.pdf"), width=10, height = 7.78)

# Testing accuracy of estimated alpha(x)
# data.alpha <- data.frame(type=c(rep("Median", n), rep("Interval.Diff", n)),
#                           value = c(alp.x.samples[,5], abs(alp.x.samples[,6]-alp.x.samples[,4])))

# ggplot(data = data.alpha, aes(x=value, fill=type)) +
#   geom_histogram(alpha=0.4, bins = 250) +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "top",
#         plot.margin = margin(0,0,0,-20),
#         axis.text = element_text(size = 20),
#         axis.title.x = element_text(size = 15)) +
#   labs(fill="")

#Predictive Distribution check
y.container <- as.data.frame(matrix(, nrow = n, ncol = 0))
# for(i in random.alpha.idx){
random.alpha.idx <- floor(runif(100, 1, ncol(t(posterior$alpha))))
for(i in random.alpha.idx){
  # y.container <- cbind(y.container, dPareto(y, u, t(posterior$alpha)[i]))
  y.container <- cbind(y.container, log(rPareto(n, u, t(posterior$alpha)[i])))
}
colnames(y.container) <- paste("col", 1:100, sep="")
y.container$x <- seq(1,n)
y.container$logy <- log(y)
plt <- ggplot(data = y.container, aes(x = x)) + ylab("density") + xlab("log(Burnt Area)") + labs(col = "")

for(i in names(y.container)){
  # plt <- plt + geom_line(aes(y = .data[[i]]), alpha = 0.2, linewidth = 0.7)
  plt <- plt + geom_density(aes(x=.data[[i]]), color = "slategray1", alpha = 0.1, linewidht = 0.7)
}

print(plt + geom_density(aes(x=logy), color = "steelblue", linewidth = 2) +
        theme_minimal(base_size = 30) + ylim(0, 2) + xlim(6.8,25) +
        theme(legend.position = "none",
                axis.text = element_text(size = 35)))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_BRSTIR_predictive_distribution.pdf"), width=10, height = 7.78)

ev.linear <- as.vector(c(1,bs.linear[which.max(y),]))
ev.nonlinear <- (matrix(bs.nonlinear[which.max(y),], ncol = 10))
ev.y <- as.data.frame(matrix(, nrow = 1, ncol = 0))
for(i in 1:ncol(t(posterior$theta))){
  ev.alpha.single <- exp(sum(t(posterior$theta)[,i]*ev.linear) + 
                              sum(t(posterior$gamma[i,,]) %*% ev.nonlinear))
  ev.y <- cbind(ev.y, log(rPareto(1, t = u, alpha = ev.alpha.single)))
}
# colnames(ev.y) <- paste("col", 1:ncol(t(posterior$theta)), sep="")
ev.y <- as.data.frame(as.numeric(t(ev.y)))
ev.y <- as.data.frame(ev.y[!is.infinite(rowSums(ev.y)),])
ev.y$logy <- max(log(y))
colnames(ev.y) <- c("yrep")

plt <- ggplot(data = ev.y, aes(x = yrep)) + ylab("density") + xlab("log(Burnt Area)") + labs(col = "") +
  geom_density(color = "steelblue", linewidth = 1.2) + 
  geom_rug(alpha = 0.1) +
  # geom_point(aes(x=yrep,y=-Inf),color="steelblue", size = 3.5, alpha = 0.2) +
  xlim(min(ev.y$yrep), 300) +
  theme_minimal(base_size = 30) +  
  theme(legend.position = "none",
  axis.text = element_text(size = 35))
  # geom_vline(xintercept = log(max(y)), linetype="dotted", color = "red",) +

d <- ggplot_build(plt)$data[[1]]
print(plt + geom_area(data = subset(d, x>12.44009), aes(x=x,y=y), fill = "slategray1", alpha = 0.5) +
        geom_segment(x=12.44009, xend=12.44009, 
              y=0, yend=approx(x = d$x, y = d$y, xout = 12.4409)$y,
              colour="red", linewidth=1.2, linetype = "dotted"))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_BRSTIR_generative.pdf"), width=10, height = 7.78)
# scaleFUN <- function(x) sprintf("%.1f", x)

# p <- ggplot(data= ev.y, aes(x=yrep, y = logy)) +
#       geom_vline(xintercept = log(max(y)), color = "red", linewidth = 0.7, linetype = "dotted") +
#       geom_point(size = 0.8, color = "steelblue", alpha = 0.2) +
#       geom_point(aes(x=max(log(y))), size = 3.5, color = "steelblue") + 
#       ylab("True log(Burnt Area)") + xlab("") +
#       #xlab("Generated log(Burnt Area)") + 
#       theme_minimal(base_size = 30) + xlim(min(ev.y$yrep), 200) +
#       ylim(min(ev.y$yrep), 200) +
#       # scale_y_continuous(labels=scaleFUN) +  
#       theme(legend.position = "none",
#               axis.text = element_text(size = 35))
# grid.arrange(p, (plt+scale_y_reverse() + theme(axis.text.x=element_text(size=10))), nrow=2, ncol=1, widths=4, heights = c(4,1))

# library(ggExtra)
# p <- p %>% ggMarginal(margins = 'y', color="steelblue", size=10)
# print(p, newpage = TRUE)
# detach("package:ggExtra", unload = TRUE)

# library(ismev)
# gpd.fit(y, u)

fit.log.lik <- extract_log_lik(fit1)
fwi.loo <- loo(fit.log.lik, cores = 2)
plot(fwi.loo, label_points = TRUE)

brstir.elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 2)
brstir.waic <- waic(fit.log.lik, cores = 2)
save(brstir.elpd.loo, brstir.waic, file = (paste0("./BRSTIR/application/BRSTIR_",Sys.Date(),"_",floor(threshold*100),"quantile_IC.Rdata")))

loo(fit.log.lik)
