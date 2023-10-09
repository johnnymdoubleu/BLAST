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
# library(R6)
# suppressMessages(library(igraph))
# library(mgcv)
library(MCMCvis)
library(cmdstanr)
# library(ggplotify)
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
# corrplot.mixed(cor(fwi.scaled),
#                 upper = "circle",
#                 lower = "number",
#                 addgrid.col = "black")
# ggsave("./Laboratory/Application/figures/correlation.pdf", width=15)
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
# ggplot(df.extreme, aes(x=month, y=y, color=month)) + geom_point(size=6) +
#     theme(plot.title = element_text(hjust = 0.5, size = 20),
#         legend.title = element_blank(),
#         legend.text = element_text(size=20),
#         # axis.ticks.x = element_blank(),
#         axis.text.x = element_text(hjust=0.35),
#         axis.text = element_text(size = 20),
#         axis.title = element_text(size = 15))
# # ggsave("./Laboratory/Application/figures/datavis.pdf", width=15)

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

psi <- 20
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

# // priors
# lambda1 ~ gamma(1, 1.78);
# intercept ~ double_exponential(0, lambda1);
# lambda2 ~ gamma(0.1, 0.1);
# sigma ~ inv_gamma(0.01, 0.01);
# for (j in 1:p){
#     theta[j] ~ double_exponential(0, lambda1);
#     tau[j] ~ gamma(atau, square(lambda2));
#     gamma[j] ~ multi_normal(rep_vector(0, psi), diag_matrix(rep_vector(1,psi)) * tau[j] * sigma);
# }
# // likelihood
# for (i in 1:n){
#     target += pareto_lpdf(y[i] | u, alpha[i]);
# }

write("// Stan model for simple linear regression
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    vector[n] y; // extreme response
    real <lower=0> atau;
}

parameters {
    vector[p] theta; // linear predictor
    vector[psi] gamma[p];
    real <lower=-1, upper = 1> sigma; //
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real intercept; // intercept term
    vector[p] tau;
}

transformed parameters {
    vector[n] alpha; // tail index
    matrix[n, p] gsmooth; // nonlinear component
    cov_matrix[psi] covmat[p]; // covariance
    for(j in 1:p){
        covmat[j] <- diag_matrix(rep_vector(1, psi)) * tau[j] * sigma;
    }

    for (j in 1:p){
        gsmooth[,j] <- bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    }
    for (i in 1:n){
        alpha[i] <- exp(intercept + dot_product(bsLinear[i], theta) + (gsmooth[i,] *rep_vector(1, p)));
    }
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 1, 5);
    target += gamma_lpdf(lambda2 | 0.1, 0.1);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    target += double_exponential_lpdf(intercept | 0, lambda1);
    for (j in 1:p){
        target += double_exponential_lpdf(theta[j] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, (square(lambda2)/2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), covmat[j]);
    }
}
//generated quantities {} // Used in Posterior predictive check"
, "model.stan")

data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2),
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear)

set_cmdstan_path(path = NULL)
#> CmdStan path set to: /Users/jgabry/.cmdstan/cmdstan-2.32.2

# Create a CmdStanModel object from a Stan program,
# here using the example model that comes with CmdStan
file <- file.path(cmdstan_path(), "model.stan")

init.alpha <- list(list(gamma = array(rep(0,(psi*p)), dim=c(psi, p)),
                                    theta = rep(0, p), intercept = 0),
                  list(gamma = array(rep(0.02,(psi*p)), dim=c(psi, p)),
                                    theta = rep(0.02, p), intercept = 0.03))

# stanc("C:/Users/Johnny Lee/Documents/GitHub/BRSTIR/application/model1.stan")
fit1 <- stan(
    file = "model.stan",  # Stan program
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    chains = 2,             # number of Markov chains
    warmup = 1000,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = 2,              # number of cores (could use one per chain)
    refresh = 1             # no progress shown
)

saveRDS(fit1, file=paste0("./BRSTIR/application/",Sys.Date(),"_stanfit.rds"))
posterior <- extract(fit1)
str(posterior)

# print(as.mcmc(fit1), pars=c("alpha", "gamma", "intercept", "theta", "lambda1", "lambda2","lp__"), probs=c(.05,.5,.95))
plot(fit1, plotfun = "trace", pars = c("intercept", "theta"), nrow = (p+1))
plot(fit1, plotfun = "trace", pars = c("lambda1", "lambda2"), nrow = 2)

summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary

# traceplot(fit1, pars = c("theta"), nrow = (p))
# traceplot(fit1, pars = c("lambda1", "lambda2"), inc_warmup = TRUE, nrow = 2)
fit.v2 <- as.mcmc(fit1)

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


theta0.samples <- summary(fit1, par=c("intercept"), probs = c(0.05,0.5, 0.95))$summary
theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary

gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q3 <- theta.samples[,6]


# theta.samples <- data.frame(apply(posterior$theta, 2, summary))


df.theta <- data.frame("seq" = seq(1, (p+1)),
                        "m" = c(theta0.samples[1], theta.post.mean),
                        "l" = c(theta0.samples[4], theta.q1),
                        "u" = c(theta0.samples[6], theta.q3))
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
                  "m" = as.vector(gamma.post.mean),
                  "l" = as.vector(gamma.q1),
                  "u" = as.vector(gamma.q3))
df.gamma$covariate <- factor(rep(names(fwi.scaled), each = psi, length.out = nrow(df.gamma)), levels = colnames(fwi.scaled))
df.gamma$labels <- factor(1:(psi*p))
ggplot(df.gamma, aes(x =labels, y = m, color = covariate)) + 
  geom_point(size = 4) + ylab("") + xlab("" ) + #ylim(-15,15) +
  # geom_ribbon(aes(ymin = l, ymax = u)) +
  geom_errorbar(aes(ymin = l, ymax = u), width = 4, linewidth = 1.2) + 
  geom_point(size = 4, color = "black") + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+10), 
                    label = c(expression(bold(gamma[1])), 
                              expression(bold(gamma[2])), expression(bold(gamma[3])), expression(bold(gamma[4])), expression(bold(gamma[5])), expression(bold(gamma[6])), expression(bold(gamma[7]))),
                    expand=c(0,3)) +
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

g.nonlinear.q1 <- g.linear.q1 <- g.q1 <- g.nonlinear.q3 <- g.linear.q3 <- g.q3 <- g.nonlinear.new <- g.linear.new <- g.new <- matrix(, nrow = n, ncol=p)
g.smooth.q1 <- g.smooth.q3 <- g.smooth.new <- alpha.new <- NULL
for (j in 1:p){
  g.linear.new[,j] <- xholder.linear[,j] * theta.post.mean[j]
  g.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.post.mean, nrow=psi)[,j] 
  g.new[1:n, j] <- g.linear.new[,j] + g.nonlinear.new[,j]
  g.linear.q1[,j] <- xholder.linear[,j] * theta.q1[j]
  g.nonlinear.q1[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q1, nrow=psi)[,j] 
  g.q1[1:n, j] <- g.linear.q1[,j] + g.nonlinear.q1[,j]
  g.linear.q3[,j] <- xholder.linear[,j] * theta.q3[j]
  g.nonlinear.q3[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.q3, nrow=psi)[,j] 
  g.q3[1:n, j] <- g.linear.q3[,j] + g.nonlinear.q3[,j]
}

for(i in 1:n){
  g.smooth.new[i] <- theta0.samples[1] + sum(g.new[i,])
  g.smooth.q1[i] <- theta0.samples[4] + sum(g.q1[i,])
  g.smooth.q3[i] <- theta0.samples[6] + sum(g.q3[i,])
}

### Plotting linear and nonlinear components
# post.mean <- as.vector(apply(as.data.frame(matrix(alpha.summary[(n+1):(n+(n*p)),1], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q1 <- as.vector(apply(as.data.frame(matrix(alpha.summary[(n+1):(n+(n*p)),4], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q3 <- as.vector(apply(as.data.frame(matrix(alpha.summary[(n+1):(n+(n*p)),5], nrow = n, ncol = p)), 2, sort, decreasing=F))
equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

data.smooth <- data.frame("x"=c(1:n),
                          "post.mean" = as.vector(g.new),
                          "q1" = as.vector(g.q1),
                          "q3" = as.vector(g.q3),
                          "covariates" = gl(p, n, (p*n), labels = factor(names(fwi.scaled))),
                          "replicate" = gl(2, n, (p*n)))
ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Smooth Functions") +
  # geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=post.mean, colour = covariates), linewidth=2) + ylab ("") +
  facet_grid(covariates ~ ., scales = "free_y") + 
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_smooth.pdf"), width=10.5, height = 15)
data.linear <- data.frame("x"=c(1:n),
                          "post.mean" = as.vector(g.linear.new),
                          "q1" = as.vector(g.linear.q1),
                          "q3" = as.vector(g.linear.q3),
                          "covariates" = gl(p, n, (p*n), labels = factor(names(fwi.scaled))),
                          "replicate" = gl(2, n, (p*n)))
ggplot(data.linear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Linear Components") +
  # geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=post.mean, colour = covariates), linewidth=2) + ylab ("") +
  facet_grid(covariates ~ ., scales = "free_y") + 
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_linear.pdf"), width=10.5, height = 15)
# post.mean <- as.vector(apply(as.data.frame(matrix(alpha.summary[((n+(n*p))+1):(n+(2*n*p)),1], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q1 <- as.vector(apply(as.data.frame(matrix(alpha.summary[((n+(n*p))+1):(n+(2*n*p)),4], nrow = n, ncol = p)), 2, sort, decreasing=F))
# q3 <- as.vector(apply(as.data.frame(matrix(alpha.summary[((n+(n*p))+1):(n+(2*n*p)),5], nrow = n, ncol = p)), 2, sort, decreasing=F))

data.nonlinear <- data.frame("x"=c(1:n),
                          "post.mean" = as.vector(g.nonlinear.new),
                          "q1" = as.vector(g.nonlinear.q1),
                          "q3" = as.vector(g.nonlinear.q3),
                          "covariates" = gl(p, n, (p*n), labels = factor(names(fwi.scaled))),
                          "replicate" = gl(2, n, (p*n)))
ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Nonlinear Components") +
  # geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=post.mean, colour = covariates), linewidth=2) + ylab ("") +
  facet_grid(covariates ~ ., scales = "free_y") + 
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=45),
        legend.margin=margin(0,0,0,-10),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_nonlinear.pdf"), width=12.5, height = 15)

data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "post.mean" = sort(alpha.samples[,1]),
                            "post.median" = sort(alpha.samples[,5]),
                            "q1" = sort(alpha.samples[,4]),
                            "q3" = sort(alpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  # geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Nonlinear Components") +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=2) + ylab(expression(alpha(x))) + labs(col = "") + 
  ylim(0, (max(data.scenario$post.mean)+10)) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linetype=2, linewidth=2) +
  # geom_line(aes(y=chain2, col = "Chian 2"), linetype=3) +
  # facet_grid(covariates ~ .) + 
  # scale_y_continuous(breaks=c(0)) + 
  scale_color_manual(values = c("red", "#e0b430"))+
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position="top", 
        legend.key.size = unit(1, 'cm'),
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_alpha.pdf"), width=10, height = 7.78)

# f.nonlinear.old <- f.linear.old <- f.old <- matrix(, nrow = n, ncol=p)
# alpha.new <- NULL
# for (j in 1:p){
#   f.linear.old[,j] <- bs.linear[,j] * theta.post.mean[j]
#   f.nonlinear.old[,j] <- bs.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.post.mean[, j]
#   f.old[1:n, j] <- f.linear.old[,j] + f.nonlinear.old[,j]
# }
# for(i in 1:n){
#   alpha.new[i] <- exp(tail(alpha.summary, 1)[1] + sum(f.old[i]))
# }

mcmc.alpha <- fit.v2$samples$chain1[,1:n]

r <- matrix(, nrow = n, ncol = 30)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 30
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y[i], u, alpha = mcmc.alpha[round(runif(1,1,len)),i]))
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
  #geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
  #           color = "lightgrey", fill = "lightgrey",
  #          alpha = 0.4, linetype = "dashed") + 
  geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
              # color = "darkgrey", fill = "darkgrey",
              alpha = 0.4, linetype = "dashed") + 
  geom_line(aes(x = grid, y = trajhat), linetype = "dashed", linewidth = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.2) + 
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
  theme_minimal(base_size = 20) +
  theme(text = element_text(size = 20)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_qqplot.pdf"), width=10, height = 7.78)              
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))

cat("sc1_Alp Done")

alpha.summary[((dim(alpha.summary)[1]-p-2):(dim(alpha.summary)[1]-p-1)),] 

