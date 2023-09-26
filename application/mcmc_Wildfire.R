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
library(nimble, warn.conflicts = FALSE)
suppressMessages(library(coda))
# library(R6)
# suppressMessages(library(igraph))
# library(mgcv)
library(MCMCvis)

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
    fwi.scaled[,i] <- scale(cov.long$measurement[missing.values])
}
fwi.index$date <- as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d")
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)
fwi.index$month <- factor(format(fwi.index$date,"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
# as.Date(substr(cov.long$...1[missing.values],1,10))
# fwi.index$day <- as.Date(substr(cov.long$...1[missing.values],9,10),"%d")
# with(cov.long[missing.values], paste(substr[...1, 6, 10],month,day,sep="-"))

fwi.scaled <- fwi.scaled[which(Y>u),]
# corrplot.mixed(cor(fwi.scaled),
#                 upper = "circle",
#                 lower = "number",
#                 addgrid.col = "black")
# ggsave("./Laboratory/Application/figures/correlation.pdf", width=15)
# cov$date <- as.Date(with(cov, paste(year,month,day,sep="-")),"%Y-%m-%d")
# cov$yearmon <- as.Date(with(cov, paste(year,month,sep="-")),"%Y-%m")
df.extreme <- cbind(y, fwi.scaled)
# df.extreme <- cbind(date = cov$date[which(Y>u)], df.extreme)
df.extreme <- cbind(month = fwi.index$month[which(Y>u)], df.extreme)
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
###################### Custom Density assigned for Pareto Distribution
dpareto <- nimbleFunction(
  run = function(x = double(0), t = double(0, default=1), u = double(0), alpha = double(0), log = integer(0, default = 0)){
    log.f <- log(t)*alpha + log(alpha) - log(x/u)*(alpha)- log(x)
    returnType(double(0))
    if(log) return(log.f)
    else return(exp(log.f))
})

rpareto <- nimbleFunction(
  run = function(n = integer(0), t = double(0, default=1), u = double(0), alpha = double(0)){
  returnType(double(0))
  if(n != 1) print("rpareto only allows n = 1; using n = 1.")
  dev <- runif(1)
  return(t/(1-dev)^(1/alpha))
})

registerDistributions(list(
  dpareto = list(
    BUGSdist = "dpareto(t, u, alpha)",
    Rdist = "dpareto(t, u, alpha)",
    pqAvail = FALSE, 
    range = c(0, Inf)
    )
))

reExp = nimbleFunction(
    run = function(a = double(0)) {
        m <- 15
        if(a > m){
          ans <- exp(m) + exp(m)*(a-m)
        }
        else(
          ans <- exp(a)
        )
        returnType(double(0))
        return(ans)
    }
)

model.penalisation <- nimbleCode({
  #prior
  lambda.1 ~ dgamma(0.1, 0.1) #gamma distribution prior for lambda
  lambda.2 ~ dgamma(0.5, 0.5)
  theta.0 ~ ddexp(0, lambda.1)
  for (j in 1:p){
    theta[j] ~ ddexp(0, lambda.1)
    tau.square[j] ~ dgamma((psi+1)/2, (lambda.2^2)/2)
    sigma.square[j] ~ dinvgamma(0.01, 0.01)
  }

  for (j in 1:p){
    covm[1:psi, 1:psi, j] <- diag(psi) * (sigma.square[j] * tau.square[j])
    gamma[1:psi, j] ~ dmnorm(zero.vec[1:psi, 1], cov = covm[1:psi, 1:psi, j])
  }

  # Likelihood
  for (j in 1:p){
    g.linear[1:n, j] <- bs.linear[1:n,j] * theta[j]
    g.nonlinear[1:n, j] <- bs.nonlinear[1:n, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[1:psi, j]
    # holder.linear[1:n, j] <- xholder.linear[1:n,j] * theta[j]
    # holder.nonlinear[1:n, j] <- xholder.nonlinear[1:n, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[1:psi, j]
  }

  for (i in 1:n){
    alpha[i] <- reExp(theta.0 + sum(g.nonlinear[i, 1:p]) + sum(g.linear[i, 1:p]))
    # alpha[i] <- log(5) / log(1 + exp(theta.0 + sum(g.nonlinear[i, 1:p]) + sum(g.linear[i, 1:p])))
    # log(new.alpha[i]) <- theta.0 + sum(holder.nonlinear[i, 1:p]) + sum(holder.linear[i, 1:p])
  }
  for(i in 1:n){
    y[i] ~ dpareto(1, u, alpha[i])
    # spy[i] <- (alpha[i]*(y[i]/u)^(-1*alpha[i])*y[i]^(-1)) / C
    # ones[i] ~ dbern(spy[i])
  }
})

constant <- list(psi = psi, n = n, p = p)
init.alpha <- function() list(list(gamma = matrix(0.5, nrow = psi, ncol=p), 
                                    theta = rep(0.1, p), theta.0 = 0.1,
                                    covm = array(1, dim = c(psi,psi, p))),
                              list(gamma = matrix(0, nrow = psi, ncol=p),
                                    theta = rep(0, p), theta.0 = 0,
                                    covm = array(1, dim = c(psi,psi, p))),
                              list(gamma = matrix(-0.1, nrow = psi, ncol=p),
                                    theta = rep(0.2, p), theta.0 = 0.9,
                                    covm = array(1, dim = c(psi,psi, p))))
                            #   list(gamma = matrix(1, nrow = psi, ncol=p)))
                              # y = as.vector(y),
monitor.pred <- c("theta.0", "theta", "gamma", "alpha", "g.linear", "g.nonlinear",
                  "lambda.1", "lambda.2")
# monitor.pred <- c("covm")
data <- list(y = as.vector(y), bs.linear = bs.linear, 
              bs.nonlinear = bs.nonlinear,
              zero.vec = as.matrix(rep(0, psi)), #sigma = 0.75,
            #    new.x = xholder, new.bs.x = new.bs.x,
              u = u) #, #C = 1000,  ones = as.vector(rep(1, n)),
              #shape = 0.5, scale = 0.5)

fit.v2 <- nimbleMCMC(code = model.penalisation,
                  constants = constant,
                  data = data,
                  monitors = monitor.pred,
                  inits = init.alpha(),
                  thin = 20,
                  niter = 120000,
                  nburnin = 100000,
                  # setSeed = 300,
                  nchains = 3,
                  # WAIC = TRUE,-
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE)

alpha.summary <- fit.v2$summary$all.chains
# saveRDS(alpha.summary, file=paste0("./BRSTIR/application/",Sys.Date(),"_allChains.rds"))

# alpha.summary[701:711,]

MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
            HPD = TRUE, xlab="theta", offset = 0.05, exact = TRUE,
            horiz = FALSE, params = c("theta.0", "theta"))
MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
            HPD = TRUE, xlab="gamma", offset = 0.5,
            horiz = FALSE, params = c("gamma"))
# MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
#             HPD = TRUE, xlab="lambda", offset = 0.5,
#             horiz = FALSE, params = c("lambda.1", "lambda.2"))            
# print(alpha.summary)
# MCMCsummary(object = fit.v2$samples, round = 3)

print(MCMCtrace(object = fit.v2$samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          n.eff = TRUE,# add eff sample size
          params = c("lambda.1", "lambda.2")))

samples <- fit.v2$samples$chain1
len <- dim(samples)[1]

gamma.post.mean <- matrix(alpha.summary[((n+(2*n*p))+1):((n+(2*n*p))+(psi*p)),1], nrow = psi, ncol = p)
gamma.q1 <- matrix(alpha.summary[((n+(2*n*p))+1):((n+(2*n*p))+(psi*p)),4], nrow = psi, ncol = p)
gamma.q3 <- matrix(alpha.summary[((n+(2*n*p))+1):((n+(2*n*p))+(psi*p)),5], nrow = psi, ncol = p)
theta.post.mean <- alpha.summary[(n+(n*p*2)+(psi*p)+2+1):(n+(n*p*2)+(psi*p)+2+p),1]
theta.q1 <- alpha.summary[(n+(n*p*2)+(psi*p)+2+1):(n+(n*p*2)+(psi*p)+2+p),4]
theta.q3 <- alpha.summary[(n+(n*p*2)+(psi*p)+2+1):(n+(n*p*2)+(psi*p)+2+p),5]

df.theta <- data.frame("seq" = seq(1, (p+1)),
                        "m" = c(tail(alpha.summary, 1)[1], theta.post.mean),
                        "l" = c(tail(alpha.summary, 1)[4], theta.q1),
                        "u" = c(tail(alpha.summary, 1)[5], theta.q3))
df.theta$covariate <- factor(c("\u03b8",names(fwi.scaled)), levels = c("\u03b8","DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC"))
df.theta$labels <- factor(c("\u03b8","DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC"))

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
  scale_color_discrete(labels = c(expression(theta[0]),"DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC")) + 
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
df.gamma$covariate <- factor(rep(names(fwi.scaled), each = psi, length.out = nrow(df.gamma)), levels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC"))
df.gamma$labels <- factor(1:(psi*p))
ggplot(df.gamma, aes(x =labels, y = m, color = covariate)) + 
  geom_point(size = 4) + ylab("") + xlab("" ) + ylim(-50,75) +
  # geom_ribbon(aes(ymin = l, ymax = u)) +
  geom_errorbar(aes(ymin = l, ymax = u), width = 4, linewidth = 1.2) + 
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
  g.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.post.mean[,j] 
  g.new[1:n, j] <- g.linear.new[,j] + g.nonlinear.new[,j]
  g.linear.q1[,j] <- xholder.linear[,j] * theta.q1[j]
  g.nonlinear.q1[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.q1[,j] 
  g.q1[1:n, j] <- g.linear.q1[,j] + g.nonlinear.q1[,j]
  g.linear.q3[,j] <- xholder.linear[,j] * theta.q3[j]
  g.nonlinear.q3[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.q3[,j] 
  g.q3[1:n, j] <- g.linear.q3[,j] + g.nonlinear.q3[,j]
}

for(i in 1:n){
  g.smooth.new[i] <- tail(alpha.summary, 1)[1] + sum(g.new[i,])
  g.smooth.q1[i] <- tail(alpha.summary, 1)[4] + sum(g.q1[i,])
  g.smooth.q3[i] <- tail(alpha.summary, 1)[5] + sum(g.q3[i,])
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
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_mcmc_nonlinear.pdf"), width=12.5, height = 15)
data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "post.mean" = sort(alpha.summary[1:n,1]),
                            "q1" = sort(alpha.summary[1:n,4]),
                            "q3" = sort(alpha.summary[1:n,5]),
                            "chain1" = sort(fit.v2$summary$chain1[1:n,1]),
                            "chain2" = sort(fit.v2$summary$chain2[1:n,1]))

ggplot(data.scenario, aes(x=x)) + 
  # geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Nonlinear Components") +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) +
  geom_line(aes(y=post.mean, col = "Posterior Mean"), linewidth=2) + ylab(expression(alpha(x))) + ylim(0, 100) +
  geom_line(aes(y=chain1, col = "Chain 1"), linetype=2) +
  geom_line(aes(y=chain2, col = "Chian 2"), linetype=3) +
  # facet_grid(covariates ~ .) + 
  # scale_y_continuous(breaks=c(0)) + 
  scale_color_manual(values = c("red", "blue", "#e0b430"))+
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
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


r <- matrix(, nrow = n, ncol = 20)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 20
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y[i], u, alpha = alpha.summary[i,1]))
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
              color = "lightgrey", fill = "lightgrey",
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

