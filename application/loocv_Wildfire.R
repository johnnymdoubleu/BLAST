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
# library(R6)
# suppressMessages(library(igraph))
# library(mgcv)
library(MCMCvis)
library(cmdstanr)

setwd("C:/Users/Johnny Lee/Documents/GitHub")
df <- read_excel("./BRSTIR/application/AADiarioAnual.xlsx", col_types = c("date", rep("numeric",40)))
df.long <- gather(df, condition, measurement, "1980":"2019", factor_key=TRUE)
df.long
head(df.long)
tail(df.long)
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
for(i in 1:length(cov)){
    cov.long <- gather(cov[[i]][,1:41], condition, measurement, "1980":"2019", factor_key=TRUE)
    fwi.index[,i] <- cov.long$measurement[missing.values]
    # fwi.scaled[,i] <- cov.long$measurement[missing.values]
    fwi.scaled[,i] <- cov.long$measurement[missing.values]
}
fwi.index$date <- as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d")
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)
fwi.index$month <- factor(format(fwi.index$date,"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))


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
  xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
  test.knot <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  knots <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  tps <- basis.tps(fwi.scaled[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  bs.linear <- cbind(bs.linear, tps[,1:no.theta])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
}

seed <- 9547
set.seed(seed)
test.set <- sort(sample.int(n, n*0.2))

stancode <- "// Stan model for simple linear regression
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp; 
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
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
    matrix[n, p] gsmooth; // nonlinear component
    for (j in 1:p){
        gsmooth[,j] <- bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    for (i in 1:n){
        alpha[i] <- exp(theta[1] + dot_product(bsLinear[i], theta[2:newp]) + (gsmooth[i,] *rep_vector(1, p)));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 1, 5);
    target += gamma_lpdf(lambda2 | 0.1, 0.1);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    target += double_exponential_lpdf(theta[1] | 0, lambda1); // target += normal_lpdf(theta[1] | 0, 0.1);
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, (square(lambda2)/2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * tau[j] * sigma);
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

data.train <- list(y = as.vector(y[-test.set]), u = u, p = p, 
                    n= (n-length(test.set)), psi = psi, 
                    atau = ((psi+1)/2), newp = (p+1),
                    bsLinear = bs.linear[-test.set,], 
                    bsNonlinear = bs.nonlinear[-test.set,])
data.test <- list(y = as.vector(y[test.set]), u = u, p = p, 
                    n= length(test.set), psi = psi, 
                    atau = ((psi+1)/2), newp = (p+1),
                    bsLinear = bs.linear[test.set, ], 
                    bsNonlinear = bs.nonlinear[test.set,])  

init.alpha <- list(list(gamma = array(rep(0,(psi*p)), dim=c(psi, p)),
                        theta = rep(0, (p+1)), 
                        tau = rep(0.01, p), sigma = 0.001, 
                        lambda1 = 0.01, lambda2 = 30),
                  list(gamma = array(rep(0.02,(psi*p)), dim=c(psi, p)),
                        theta = rep(0.1, (p+1)), 
                        tau = rep(0.01, p), sigma = 1,
                        lambda1 = 0.01, lambda2 = 30),
                  list(gamma = array(rep(-0.02, (psi*p)), dim=c(psi, p)),
                        theta = rep(-0.2, (p+1)), 
                        tau = rep(0.01, p), sigma = 1,
                        lambda1 = 0.07, lambda2 = 30))

stanmodel <- stan_model(model_code = stancode)

fit2 <- sampling(stanmodel, data = data.train, 
                  init = init.alpha,
                  iter = 4000, 
                  chains = 3, cores = 4, 
                  refresh = 0)

color_scheme_set("brightblue") # check out ?bayesplot::color_scheme_set
lambda.draws <- as.matrix(fit, pars = c("lambda1", "lambda2"))
mcmc_areas(lambda.draws, prob = 0.8) # c                  
y.rep <- as.matrix(fit, pars = "y")
ppc_dens_overlay(y, y_rep[1:50, ])

gen <- gqs(stanmodel, draws = as.matrix(fit2), data = data.test)
log.pd <- extract_log_lik(gen)

loo(log.pd)
waic(log.pd)
