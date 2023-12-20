library(npreg)
library(Pareto)
suppressMessages(library(tidyverse))
library(JOPS)
library(readxl)
library(gridExtra)
library(colorspace)
library(corrplot)
library(ReIns)
library(evir)
library(rstan)
suppressMessages(library(coda))
library(ggmcmc)
library(loo)
library(bayesplot)
# library(R6)
# suppressMessages(library(igraph))
# library(mgcv)
library(MCMCvis)
library(cmdstanr)

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
threshold <- 0.95
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
    # cov.long[which(is.na(cov.long$measurement)),]
    # cov.long[which(is.na(cov.long$...1))+1,]
    # cov.long[which(!is.na(df.long$measurement)),]
    fwi.index[,i] <- cov.long$measurement[missing.values]
    # fwi.scaled[,i] <- cov.long$measurement[missing.values]
    fwi.scaled[,i] <- cov.long$measurement[missing.values]
}
fwi.index$date <- as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d")
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)
fwi.index$month <- factor(format(fwi.index$date,"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
# as.Date(substr(cov.long$...1[missing.values],1,10))
# fwi.index$day <- as.Date(substr(cov.long$...1[missing.values],9,10),"%d")
# with(cov.long[missing.values], paste(substr[...1, 6, 10],month,day,sep="-"))

fwi.scaled <- fwi.scaled[which(Y>u),]
fwi.scaled <- as.data.frame(scale(fwi.scaled))

# pdf(file = "./BRSTIR/application/figures/correlation.pdf")
# corrplot.mixed(cor(fwi.scaled),
#                 upper = "circle",
#                 lower = "number",
#                 addgrid.col = "black")
# dev.off()
# ggsave("./BRSTIR/application/figures/correlation.pdf", plot = replayPlot(p1), width=10, height = 7.78)
# cov$date <- as.Date(with(cov, paste(year,month,day,sep="-")),"%Y-%m-%d")
# cov$yearmon <- as.Date(with(cov, paste(year,month,sep="-")),"%Y-%m")
# special <- gather(fwi.scaled, cols, value) |> spread(cols, value) |> select(colnames(fwi.scaled))
# ggplot(gather(fwi.scaled, cols, value), aes(x = value)) + 
#        geom_histogram(binwidth = 0.1) + facet_grid(cols~.)
# ggplot(gather(fwi.index[which(Y>u),1:7], cols, value), aes(x = value)) + 
#        geom_histogram(binwidth = 2) + facet_grid(cols~.)
df.extreme <- cbind(y, fwi.scaled)
# df.extreme <- cbind(date = cov$date[which(Y>u)], df.extreme)
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

# # ggplot(as.data.frame(y),aes(x=y))+geom_histogram(aes(y=..density..), bins = 5)

# ggplot(df.extreme, aes(x=y)) +
#     geom_histogram(stat = "density", n=40, adjust=0.1, fill = "darkgrey") + 
#     xlab("Area Burnt") + 
#     # geom_histogram(aes(y=..density..), bins = 10^3) +
#     # geom_density(aes(y=..density..)) +
#     theme(plot.title = element_text(hjust = 0.5, size = 20),
#       legend.title = element_blank(),
#       legend.text = element_text(size=20),
#       # axis.ticks.x = element_blank(),
#       axis.text.x = element_text(hjust=0.35),
#       axis.text = element_text(size = 25),
#       axis.title = element_text(size = 30))

psi <- 10
n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]
no.theta <- 1
newx <- seq(0, 1, length.out=n)
xholder.linear <- xholder.nonlinear <- bs.linear <- bs.nonlinear <- matrix(,nrow=n, ncol=0)
xholder <- matrix(nrow=n, ncol=p)
for(i in 1:p){
  # xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
  # test.knot <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  # splines <- basis.tps(seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n), test.knot, m=2, rk=FALSE, intercept = TRUE)
  xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
  test.knot <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  knots <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  tps <- basis.tps(fwi.scaled[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  # tps <- mSpline(x.origin[,i], df=psi, Boundary.knots = range(x.origin[,i]), degree = 3, intercept=TRUE)
  bs.linear <- cbind(bs.linear, tps[,1:no.theta])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
}


write("// Stan model for simple linear regression
data {
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
}

parameters {
    vector[newp] theta; // linear predictor
    vector[psi] gamma[p]; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real sigma; //
    vector[p] tau;
}

transformed parameters {
    vector[n] alpha; // tail index
    vector[n] newalpha; // tail index    
    matrix[n, p] gnl; // nonlinear component
    matrix[n, p] gl; // linear component
    matrix[n, p] gsmooth; // smooth function
    matrix[n, p] newgnl; // nonlinear component
    matrix[n, p] newgl; // linear component
    matrix[n, p] newgsmooth; // smooth function    
    for (j in 1:p){
        gnl[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        gl[,j] = bsLinear[,j] * theta[j+1];
        gsmooth[,j] = gl[,j] + gnl[,j];
        newgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgl[,j] = xholderLinear[,j] * theta[j+1];
        newgsmooth[,j] = newgl[,j] + newgnl[,j];        
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
    target += inv_gamma_lpdf(sigma | 0.01, 0.01); // target += double_exponential_lpdf(theta[1] | 0, lambda1)
    target += (newp * log(lambda1) + (p * psi * log(lambda2)));
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, lambda2/sqrt(2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * sqrt(tau[j]) * sqrt(sigma));
    }
}
generated quantities {
    // Used in Posterior predictive check
    vector[n] log_lik;
    real y_rep[n] = pareto_rng(u, alpha);
    for (i in 1:n) {
        log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}
"
, "model.stan")

data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), newp = (p+1),
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear)

set_cmdstan_path(path = NULL)
#> CmdStan path set to: /Users/jgabry/.cmdstan/cmdstan-2.32.2

# Create a CmdStanModel object from a Stan program,
# here using the example model that comes with CmdStan
file <- file.path(cmdstan_path(), "model.stan")

init.alpha <- list(list(gamma = array(rep(0, (psi*p)), dim=c(psi, p)),
                        theta = rep(0, (p+1)), 
                        tau = rep(0.1, p), sigma = 0.1, 
                        lambda1 = 0.1, lambda2 = 0.1),
                  list(gamma = array(rep(0.02, (psi*p)), dim=c(psi, p)),
                        theta = rep(0.01, (p+1)), 
                        tau = rep(0.01, p), sigma = 0.001,
                        lambda1 = 0.01, lambda2 = 0.1),
                  list(gamma = array(rep(0.01, (psi*p)), dim=c(psi, p)),
                        theta = rep(0.05, (p+1)), 
                        tau = rep(0.01, p), sigma = 0.01,
                        lambda1 = 0.1, lambda2 = 0.01))

# stanc("C:/Users/Johnny Lee/Documents/GitHub/BRSTIR/application/model1.stan")
fit1 <- stan(
    file = "model.stan",  # Stan program
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    # init_r = 1,
    chains = 3,             # number of Markov chains
    warmup = 1500,          # number of warmup iterations per chain
    iter = 3000,            # total number of iterations per chain
    cores = 4,              # number of cores (could use one per chain)
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
gl.samples <- summary(fit1, par=c("newgl"), probs = c(0.05, 0.5, 0.95))$summary
gnl.samples <- summary(fit1, par=c("newgnl"), probs = c(0.05, 0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary
summary(fit1, par=c("sigma"), probs = c(0.05,0.5, 0.95))$summary
summary(fit1, par=c("tau"), probs = c(0.05,0.5, 0.95))$summary

gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q2 <- gamma.samples[,5]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q2 <- theta.samples[,5]
theta.q3 <- theta.samples[,6]


# theta.samples <- data.frame(apply(posterior$theta, 2, summary))


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
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_theta.pdf"), width=10, height = 7.78)

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


df.gamma <- data.frame("seq" = seq(1, (psi*p)), 
                  "m" = as.vector(gamma.q2),
                  "l" = as.vector(gamma.q1),
                  "u" = as.vector(gamma.q3))
df.gamma$covariate <- factor(rep(names(fwi.scaled), each = psi, length.out = nrow(df.gamma)), levels = colnames(fwi.scaled))
df.gamma$labels <- factor(1:(psi*p))
ggplot(df.gamma, aes(x =labels, y = m, color = covariate)) + 
  geom_errorbar(aes(ymin = l, ymax = u), alpha = 0.4, width = 4, linewidth = 1.2) +
  geom_point(size = 4) + ylab("") + xlab("" ) + #xlim(1,(psi*p)) +
  # geom_ribbon(aes(ymin = l, ymax = u)) +
  # geom_point(size = 4, color = "black") + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+10), 
                    label = c(expression(bold(gamma[1])), 
                              expression(bold(gamma[2])), 
                              expression(bold(gamma[3])), 
                              expression(bold(gamma[4])), 
                              expression(bold(gamma[5])), 
                              expression(bold(gamma[6])), 
                              expression(bold(gamma[7]))),
                    expand=c(0,10)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size=25),
          legend.margin=margin(0,0,0,-10),
          legend.box.margin=margin(-10,0,-10,0),
          plot.margin = margin(0,0,0,-20),
          axis.text.x = element_text(hjust=0.5),
          axis.text = element_text(size = 28))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_gamma.pdf"), width=10, height = 7.78)


g.linear.mean <- as.vector(matrix(gl.samples[,1], nrow = n, byrow=TRUE))
g.linear.q1 <- as.vector(matrix(gl.samples[,4], nrow = n, byrow=TRUE))
g.linear.q2 <- as.vector(matrix(gl.samples[,5], nrow = n, byrow=TRUE))
g.linear.q3 <- as.vector(matrix(gl.samples[,6], nrow = n, byrow=TRUE))
g.nonlinear.mean <- as.vector(matrix(gnl.samples[,1], nrow = n, byrow=TRUE))
g.nonlinear.q1 <- as.vector(matrix(gnl.samples[,4], nrow = n, byrow=TRUE))
g.nonlinear.q2 <- as.vector(matrix(gnl.samples[,5], nrow = n, byrow=TRUE))
g.nonlinear.q3 <- as.vector(matrix(gnl.samples[,6], nrow = n, byrow=TRUE))
g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))

# g.nonlinear.q2 <- g.linear.q2 <- g.q2 <- g.nonlinear.q1 <- g.linear.q1 <- g.q1 <- g.nonlinear.q3 <- g.linear.q3 <- g.q3 <- g.nonlinear.new <- g.linear.new <- g.new <- matrix(, nrow = n, ncol=p)
# for (j in 1:p){
#   g.linear.new[,j] <- xholder.linear[,j] * theta.post.mean[(j+1)]
#   g.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.post.mean, nrow=psi)[,j] 
#   g.new[1:n, j] <- g.linear.new[,j] + g.nonlinear.new[,j]
#   g.linear.q1[,j] <- xholder.linear[,j] * theta.q1[(j+1)]
#   g.nonlinear.q1[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q1, nrow=psi)[,j] 
#   g.q1[1:n, j] <- g.linear.q1[,j] + g.nonlinear.q1[,j]
#   g.linear.q2[,j] <- xholder.linear[,j] * theta.q2[(j+1)]
#   g.nonlinear.q2[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q2, nrow=psi)[,j] 
#   g.q2[1:n, j] <- g.linear.q2[,j] + g.nonlinear.q2[,j]  
#   g.linear.q3[,j] <- xholder.linear[,j] * theta.q3[(j+1)]
#   g.nonlinear.q3[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q3, nrow=psi)[,j] 
#   g.q3[1:n, j] <- g.linear.q3[,j] + g.nonlinear.q3[,j]
# }


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
                          "covariates" = gl(p, n, (p*n), labels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC")),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC")))

ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab("") +
  facet_grid(covariates ~ ., scales = "free", switch = "y",
              labeller = label_parsed) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) + #ylim(-0.65, 0.3) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
          plot.margin = margin(0,0,0,-20),
          strip.text = element_blank(),
          axis.text = element_text(size = 20))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_smooth.pdf"), width=12.5, height = 15)

data.linear <- data.frame("x"= as.vector(xholder),
                          "post.mean" = as.vector(g.linear.mean),
                          "q1" = as.vector(g.linear.q1),
                          "q2" = as.vector(g.linear.q2),
                          "q3" = as.vector(g.linear.q3),
                          "covariates" = gl(p, n, (p*n), labels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC")),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC")))

ggplot(data.linear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  # geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  # geom_line(aes(y=true, colour = "True"), linewidth=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab("") +
  facet_grid(covariates ~ ., scales = "free", switch = "y",
              labeller = label_parsed) + 
  # scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) + #ylim(-0.65, 0.3) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
          plot.margin = margin(0,0,0,-20),
          strip.text = element_blank(),
          axis.text = element_text(size = 20))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_linear.pdf"), width=12.5, height = 15)
# post.mean <- as.vector(apply(as.data.frame(matrix(alpha.summary[((n+(n*p))+1):(n+(2*n*p)),1], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q1 <- as.vector(apply(as.data.frame(matrix(alpha.summary[((n+(n*p))+1):(n+(2*n*p)),4], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q3 <- as.vector(apply(as.data.frame(matrix(alpha.summary[((n+(n*p))+1):(n+(2*n*p)),5], nrow = n, ncol = p)), 2, sort, decreasing=F))

data.nonlinear <- data.frame("x"=as.vector(xholder),
                          "post.mean" = as.vector(g.nonlinear.mean),
                          "q1" = as.vector(g.nonlinear.q1),
                          "q2" = as.vector(g.nonlinear.q2),
                          "q3" = as.vector(g.nonlinear.q3),
                          "covariates" = gl(p, n, (p*n), labels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC")),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC")))

ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  # geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  # geom_line(aes(y=true, colour = "True"), linewidth=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab("") +
  facet_grid(covariates ~ ., scales = "free", switch = "y",
              labeller = label_parsed) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) + #ylim(-0.65, 0.3) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
          plot.margin = margin(0,0,0,-20),
          strip.text = element_blank(),
          axis.text = element_text(size = 20))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_nonlinear.pdf"), width=12.5, height = 15)

data.scenario <- data.frame("x" = c(1:n),
                            "post.mean" = sort(alpha.samples[,1]),
                            "post.median" = sort(alpha.samples[,5]),
                            "q1" = sort(alpha.samples[,4]),
                            "q3" = sort(alpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(bold(x)))) + xlab(expression(x)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y = true, col = "True"), linewidth = 2) +
  ylim(0, 2500) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue")) + 
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

len <- dim(posterior$alpha)[1]
r <- matrix(, nrow = n, ncol = 30)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 30
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i]))
    # r[i, t] <- qnorm(pPareto(y[i], u, alpha = alpha.new[i]))
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
  theme(text = element_text(size = 20)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_qqplot.pdf"), width=10, height = 7.78)

# mcmc.gamma <- posterior$gamma
# gamma.container <- as.data.frame(matrix(NA, nrow = 20, ,ncol = 0))
# for(i in 5900:6000){
#     gamma.container <- cbind(gamma.container, mcmc.gamma[i,1,])    
# }
# corrplot.mixed(cor(gamma.container),
#                 upper = "circle",
#                 lower = "number",
#                 addgrid.col = "black")


# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_qqplot.pdf"), width=10, height = 7.78)              
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))

cat("Finished Running")


fwi.loo <- loo(fit1)
fwi.loo
y.psis <- fwi.loo$pointwise[,1]
print(paste("RMSE (PSIS) =",round(sqrt(mean((y-y.psis)^2)) ,2)))
print(paste("ELPD (PSIS)=",round(sum(y.psis),2)))
print(paste("ELPD (brute force)=",round(sum(y),2)))

plot(fwi.loo, label_points = TRUE)
fwi.waic <- waic(posterior$log_lik)
fwi.waic

y.rep <- as.matrix(fit1, pars = "y_rep")
ppc_loo_pit_overlay(
  y = y,
  yrep = y.rep,
  lw = weights(fwi.loo$psis_object)
)
# ppc_dens_overlay(y, y.rep[1:50, ])



# plot(xholder[,2], g.new[,c(2)])
# plot(xholder[,2], rowSums(g.new[,c(1,2,4)]))
# plot(xholder[,3], rowSums(g.new[,c(3,6)]))
