# list.of.packages <- c("nimble", "Pareto", "JOPS","tidyverse")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# library(gridExtra)
# library(colorspace)
library(npreg)
library(Pareto)
library(tidyverse)
library(JOPS)
library(readxl)
library(gridExtra)
library(colorspace)
library(corrplot)
library(JOPS)
library(nimble, warn.conflicts = FALSE)
suppressMessages(library(coda))
# library(R6)
# suppressMessages(library(igraph))
# library(mgcv)
library(MCMCvis)
suppressMessages(library(tidyverse))
# library(ggplotify)

#Scenario 1
# set.seed(10)
n <- 5000
psi <- 20
threshold <- 0.90
p <- 10
no.theta <- 1
simul.no <- 50


theta.container <- as.data.frame(matrix(, nrow = (no.theta *p), ncol= simul.no))
gamma.container <- as.data.frame(matrix(, nrow = (psi * p), ncol = simul.no))
# linear.container <- nonlinear.container <- f.container <- lapply(1:simul.no, matrix, data= NA, nrow=(n*(1-threshold)), ncol=p)
linear.container <- nonlinear.container <- f.container <- lapply(1:simul.no, data.frame)
alpha.container <- as.data.frame(matrix(, nrow=n, ncol = simul.no))

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
x.origin <- cbind(replicate(p, runif(n, 0, 1)))
for(i in 1:p){
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    # tps <- mSpline(x.origin[,i], df=psi, Boundary.knots = range(x.origin[,i]), degree = 3, intercept=TRUE)
    #   bs.x <- cbind(bs.x, tps)
    bs.linear <- cbind(bs.linear, tps[,1:no.theta])
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
}

gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
    for (ps in 1:psi){
        if(j %in% c(2,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
        else if(j==7){
            if(ps <= (psi/2)){gamma.origin[ps, j] <- 1}
            else{gamma.origin[ps, j] <- 1}
        }
        else {
            if(ps <= (psi/2)){gamma.origin[ps, j] <- 1}
            else{gamma.origin[ps, j] <- 1}
        }
    }
}

# theta.origin <- matrix(, nrow = 2, ncol = p)
# for(j in 1:p){
#     for (k in 1:2){
#         if(j %in% c(2,4,5,6,9,10)){theta.origin[k, j] <- 0}
#         else if(j==7){
#             if(k==1){theta.origin[k, j] <- 0.5}
#             else{theta.origin[k, j] <- -0.3}
#         }
#         else {
#             if(k==1){theta.origin[k,j] <- -0.2}
#             else{theta.origin[k,j] <- 0.8}
#         }
#     }
# }
theta.origin <- c(-0.1, 0.8, 0, 0.8, 0, 0, 0, -0.3, 0.8, 0, 0)

f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.linear.origin[,j] <- bs.linear[1:n, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
    f.origin[, j] <- rep(theta.origin[1], n) + f.linear.origin[,j] + f.nonlinear.origin[,j]
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp(sum(f.origin[i,]))
    y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i])
}

u <- quantile(y.origin, threshold)
x.origin <- x.origin[which(y.origin>u),]
y.origin <- y.origin[y.origin > u]
n <- length(y.origin)

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
newx <- seq(0, 1, length.out=n)
xholder <- bs.x <- matrix(, nrow = n, ncol = p)
for(i in 1:p){
    xholder[,i] <- seq(0, 1, length.out = n)
    test.knot <- seq(0, 1, length.out = psi)
    splines <- basis.tps(newx, test.knot, m=2, rk=FALSE, intercept = FALSE)
    xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
    xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
    # tps <- mSpline(x.origin[,i], df=psi, Boundary.knots = range(x.origin[,i]), degree = 3, intercept=TRUE)
    #   bs.x <- cbind(bs.x, tps)
    bs.linear <- cbind(bs.linear, tps[,1:no.theta])
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
}

f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    # f.origin[, j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j] + (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
    # f.linear.origin[,j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j]
    # f.nonlinear.origin[,j] <- (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
    f.linear.origin[,j] <- bs.linear[1:n, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    f.origin[, j] <- rep(theta.origin[1], n) + f.linear.origin[,j] + f.nonlinear.origin[,j]
}

alp.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp(sum(f.origin[i,]))
}


#### Test Set
f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
true.alpha <- alp.new <- NULL
for (j in 1:p){
    f.linear.new[,j] <- xholder.linear[, j] * theta.origin[j+1]
    f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    f.new[,j] <- rep(theta.origin[1], n) + f.linear.new[,j] + f.nonlinear.new[,j]
    # f.linear.origin[,j] <- as.matrix(xholder.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j]
    # f.nonlinear.origin[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    # f.origin[,j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
}
    # set.seed(100)
for(i in 1:n){
    # true.alpha[i] <- exp(sum(f.origin[i,]))
    alp.new[i] <- exp(sum(f.new[i,]))
}

# alp.origin <- y.origin <- NULL
# for(i in 1:n){
#   alp.origin[i] <- exp(sum(f.origin[i,]))
#   y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i])
# }


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

# nimRowSums = nimbleFunction(
#     run = function(a = double(2)) {
#         nrows <- dim(a)[1]
#         ncols <- dim(a)[2]
#         ans <- numeric(nrows)
#         for(j in 1:nrows) {
#             ans[j] <- sum(a[j, 1:ncols])
#         }
#         return(ans)
#         returnType(double(1))
#     }
# )

model.penalisation <- nimbleCode({
  #prior
  lambda.1 ~ dgamma(shape, scale) #gamma distribution prior for lambda
  lambda.2 ~ dgamma(shape, scale)
  # for(j in 1:p){
  #   lambda.1[j] ~ dgamma(shape, scale)
  #   lambda.2[j] ~ dgamma(shape, scale)
  # }
  # for (j in 1:p){
  #   theta[j] ~ ddexp(0, lambda.1)
  # }
  # lambda.0 ~ dgamma(shape, scale)
  theta.0 ~ ddexp(0, lambda.1)
  for (j in 1:p){
    theta[j] ~ ddexp(0, lambda.1)
    tau.square[j] ~ dgamma((psi+1)/2, (lambda.2^2)/2)
  }

  for (j in 1:p){
    covm[1:psi, 1:psi, j] <- diag(psi) * ((sigma^2) * sqrt(tau.square[j]))
    gamma[1:psi, j] ~ dmnorm(zero.vec[1:psi, 1], cov = covm[1:psi, 1:psi, j])
  }

  # Likelihood
  for (j in 1:p){
    g.linear[1:n, j] <- bs.linear[1:n,j] * theta[j]
    g.nonlinear[1:n, j] <- bs.nonlinear[1:n, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[1:psi, j]
    holder.linear[1:n, j] <- xholder.linear[1:n,j] * theta[j]
    holder.nonlinear[1:n, j] <- xholder.nonlinear[1:n, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[1:psi, j]
  }

#   gamma[1:psi, 1] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[1]) * I[1:psi, 1:psi]))
#   gamma[1:psi, 2] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[2]) * I[1:psi, 1:psi]))
#   gamma[1:psi, 3] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[3]) * I[1:psi, 1:psi]))
#   gamma[1:psi, 4] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[4]) * I[1:psi, 1:psi]))
#   gamma[1:psi, 5] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[5]) * I[1:psi, 1:psi]))
#   gamma[1:psi, 6] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[6]) * I[1:psi, 1:psi]))
#   gamma[1:psi, 7] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[7]) * I[1:psi, 1:psi]))
#   gamma[1:psi, 8] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[8]) * I[1:psi, 1:psi]))
#   gamma[1:psi, 9] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[9]) * I[1:psi, 1:psi]))
#   gamma[1:psi, 10] ~ dmnorm(zero.vec, cov = (sigma^2 * sqrt(tau.square[10]) * I[1:psi, 1:psi]))
#   g.linear[1:n, 1] <- bs.nonlinear[1:n, 1:psi] %*% gamma.1[1:psi]
#   g.linear[1:n, 2] <- bs.nonlinear[1:n, 1:psi] %*% gamma.2[1:psi]
#   g.linear[1:n, 3] <- bs.nonlinear[1:n, 1:psi] %*% gamma.3[1:psi]
#   g.linear[1:n, 4] <- bs.nonlinear[1:n, 1:psi] %*% gamma.4[1:psi]
#   g.linear[1:n, 5] <- bs.nonlinear[1:n, 1:psi] %*% gamma.5[1:psi]
#   g.linear[1:n, 6] <- bs.nonlinear[1:n, 1:psi] %*% gamma.6[1:psi]
#   g.linear[1:n, 7] <- bs.nonlinear[1:n, 1:psi] %*% gamma.7[1:psi]
#   g.linear[1:n, 8] <- bs.nonlinear[1:n, 1:psi] %*% gamma.8[1:psi]
#   g.linear[1:n, 9] <- bs.nonlinear[1:n, 1:psi] %*% gamma.9[1:psi]
#   g.linear[1:n, 10] <- bs.nonlinear[1:n, 1:psi] %*% gamma.10[1:psi]
  
  for (i in 1:n){
    log(alpha[i]) <- theta.0 + sum(g.nonlinear[i, 1:p]) + sum(g.linear[i, 1:p])
    log(new.alpha[i]) <- theta.0 + sum(holder.nonlinear[i, 1:p]) + sum(holder.linear[i, 1:p])
  }
  for(i in 1:n){
    y[i] ~ dpareto(1, u, alpha[i])
    # spy[i] <- (alpha[i]*(y[i]/u)^(-1*alpha[i])*y[i]^(-1)) / C
    # ones[i] ~ dbern(spy[i])
  }
})


constant <- list(psi = psi, n = n, p = p)
init.alpha <- function() list(list(gamma = matrix(0.5, nrow = psi, ncol=p), 
                                    theta = rep(0.2, p), theta.0 = 0.1,
                                    covm = array(1, dim = c(psi,psi,10))),
                              list(gamma = matrix(0, nrow = psi, ncol=p),
                                    theta = rep(0, p), theta.0 = 0,
                                    covm = array(1, dim = c(psi,psi,10))))
                            #   list(gamma = matrix(1, nrow = psi, ncol=p)))
                              # y = as.vector(y),
monitor.pred <- c("theta.0", "theta", "gamma", "alpha", "new.alpha", "lambda.1", "lambda.2")
# monitor.pred <- c("covm")
data <- list(y = as.vector(y.origin), bs.linear = bs.linear, 
              bs.nonlinear = bs.nonlinear,
              xholder.linear = xholder.linear,
              xholder.nonlinear = xholder.nonlinear,
              zero.vec = as.matrix(rep(0, psi)), sigma = 1,
            #    new.x = xholder, new.bs.x = new.bs.x,
              u = u, #C = 1000,  ones = as.vector(rep(1, n)),
              shape = 0.1, scale = 0.1)

fit.v2 <- nimbleMCMC(code = model.penalisation,
                  constants = constant,
                  data = data,
                  monitors = monitor.pred,
                  inits = init.alpha(),
                  thin = 20,
                  niter = 70000,
                  nburnin = 50000,
                  # setSeed = 300,
                  nchains = 2,
                  # WAIC = TRUE,-
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE)

# alpha.summary <- MCMCsummary(object = fit.v2$samples)
# saveRDS(fit.v2, file=paste0("./Simulation/BayesianPsplines/results/",date,"-",time,"_sc1_fit.rds"))
alpha.summary <- fit.v2$summary$all.chains
# saveRDS(alpha.summary, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_allChains.rds"))

# alpha.summary[701:702,]

MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
            HPD = TRUE, xlab="beta", offset = 0.05, exact = TRUE,
            horiz = FALSE, params = c("theta.0", "theta"))
MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
            HPD = TRUE, xlab="gamma", offset = 0.5,
            horiz = FALSE, params = c("gamma"))
# MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
#             HPD = TRUE, xlab="lambda", offset = 0.5,
#             horiz = FALSE, params = c("lambda.1", "lambda.2"))            
print(alpha.summary)
MCMCsummary(object = fit.v2$samples, round = 3)
MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2, 
            horiz = FALSE, params = "alpha", offset = 0.2, rank = TRUE)
# ggsave(filename = paste0("Simulation/BayesianPsplines/results/figures/",date,"-",time, "_sc1_gamma.pdf"),
#       plot = as.ggplot(function() MCMCplot(object = fit.v2$samples$chain1, 
#                                            object2 = fit.v2$samples$chain2,
#                                            HPD = TRUE, xlab="gamma", offset = 0.5,
#                                            horiz = FALSE, params = c("gamma"))),
#       width=7, height = 5.95)
# gelman.diag(fit.v2$samples)
# gelman.diag(fit.v2$samples, multivariate = FALSE)
# MCMCsummary(fit.v2$samples, params="alpha")


# ggsave(filename = paste0("Simulation/BayesianPsplines/results/figures/",date,"-",time, "_sc1_alpha.pdf"),
#       plot = as.ggplot(function() MCMCplot(object = fit.v2$samples$chain1, 
#                                            object2 = fit.v2$samples$chain2, 
#                                            offset = 0.5, rank = TRUE,
#                                            horiz = FALSE, params = c("alpha"))),
#       width=7, height = 5.95)

# MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2, 
#             horiz = FALSE, params = c("gamma[1, 2]", "gamma[2, 2]", "gamma[3, 2]", "gamma[4, 2]", "gamma[5, 2]", "gamma[6, 2]", "gamma[7, 2]", "gamma[8, 2]", "gamma[9, 2]", "gamma[10, 2]"), ISB = FALSE, exact = TRUE, HPD = TRUE)
# MCMCtrace(object = fit.v2$samples,
#           pdf = FALSE, # no export to PDF
#           ind = TRUE, # separate density lines per chain
#           n.eff = TRUE, ISB = FALSE, exact = TRUE,# add eff sample size
#           params = c("gamma[1, 2]", "gamma[2, 2]", "gamma[3, 2]", "gamma[4, 2]", "gamma[5, 2]", "gamma[6, 2]", "gamma[7, 2]", "gamma[8, 2]", "gamma[9, 2]", "gamma[10, 2]"))
# MCMCtrace(object = fit.v2$samples,
#           pdf = FALSE, # no export to PDF
#           ind = TRUE, # separate density lines per chain
#           n.eff = TRUE,# add eff sample size
#           params = "gamma")
# MCMCtrace(object = fit.v2$samples,
#           pdf = FALSE, # no export to PDF
#           ind = TRUE, # separate density lines per chain
#           n.eff = TRUE,# add eff sample size
#           params = "theta")
# MCMCtrace(object = fit.v2$samples,
#           pdf = FALSE, # no export to PDF
#           ind = TRUE, # separate density lines per chain
#           n.eff = TRUE,# add eff sample size
#           params = c("lambda.1", "lambda.2"))

systime <- Sys.time()
Sys.time()
systime <- chartr(":","-",systime)
date <- substr(systime, 1, 10)
time <- substr(systime, 12, 20)

beta.len <- dim(alpha.summary)[1]-n-n-(psi*(p-1))


samples <- fit.v2$samples$chain1
data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "post.mean" = sort(fit.v2$summary$chain1[1:n,1]),
                            "trueAlp" = sort(alp.new),
                            "meanAlp" = sort(fit.v2$summary$chain1[703:1202,1]),
                            "post.check" = sort(alp.origin))
# data.scenario1 <- data.frame("x"=c(1:n)) #, )

len <- dim(samples)[1]
for(i in 1:len){
  data.scenario <- cbind(data.scenario, 
                      data.frame(unname(sort(samples[i, 1:n]))))
}
for(i in 1:len){
  data.scenario <- cbind(data.scenario, 
                      data.frame(unname(sort(samples[i, 703:1202]))))
}
colnames(data.scenario) <- c("x", "constant", "post.mean",
                              "trueAlp", "meanAlp", "post.check",
                              paste("alp", 1:len, sep = ""),
                              # paste0("post.samp", 1:len))
                              paste0("post.samp.alp", 1:len))
tail(data.scenario[, c(1:10, ((dim(samples)[1]*2-5):(dim(samples)[1]*2)))], 15)
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))
#plotting all the points
plt <- ggplot(data = data.scenario, aes(x = x)) + ylab(expression(alpha(x))) + xlab(expression(x))
for(i in (dim(samples)[1] - 993):(dim(samples)[1]+6)){
  plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}

print(plt + geom_line(aes(y=post.mean, col = "Posterior Mean(Chain1)"), linewidth = 1.5) + 
        geom_line(aes(y = post.check, col=paste0("Simulated Alpha: ",n,"/",psi,"/",threshold)), linewidth = 1.5) +
        # theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        # theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        labs(col = "") + 
        scale_color_manual(values = c("#e0b430","red")) +
        theme(text = element_text(size = 27)) + 
        theme(legend.position="top", legend.key.size = unit(1, 'cm')))
# ggsave(paste0("./Simulation/BayesianPsplines/results/figures/",date,"-",time, "_sc1_Alp_samp_1.pdf"), 
#         width=14, height = 7.5)

cat("sc1_Alp Done")

plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x))) + xlab(expression(x))
for(i in ((dim(samples)[1]*2)-993):((dim(samples)[1]*2)+6)){
  # print(i)
  plt.samp <- plt.samp + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}
print(plt.samp + ylim(0, 30) +
      geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
      geom_line(aes(y = meanAlp, col = "Posterior Mean(Chain1)"), linewidth = 2.5) +
      labs(col = "") +
      scale_color_manual(values = c("#e0b430", "red"))+
      # scale_colour_manual("", 
      #               breaks = c("True Alpha", "Posterior Mean"),
      #               values = c("True Alpha"="red", "Posterior Mean"="#e0b430"),
      #               guide = "legend") +
      theme(text = element_text(size = 27)) + 
      theme(legend.position="top", legend.key.size = unit(1, 'cm')))


# ggsave(paste0("./Simulation/BayesianPsplines/results/figures/",date,"-",time, "_true_Alp_samp_1.pdf"), 
#         width=14, height = 7.5)

cat("true_Alp Done")

samples <- fit.v2$samples$chain2
data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "post.mean" = sort(fit.v2$summary$chain2[1:n,1]),
                            "trueAlp" = sort(trueAlp),
                            "meanAlp" = sort(fit.v2$summary$chain2[(n+beta.len+((p-1)*psi)+1):(n+beta.len+n+((p-1)*psi)),1]),
                            "post.check" = sort(cutoff.alp))
                            # "post.check" = sort(alp))
# data.scenario1 <- data.frame("x"=c(1:n)) #, )

len <- dim(samples)[1]
for(i in 1:len){
  data.scenario <- cbind(data.scenario, 
                      data.frame(unname(sort(samples[i, 1:n]))))
}
for(i in 1:len){
  data.scenario <- cbind(data.scenario, 
                      data.frame(unname(sort(samples[i, (n+beta.len+((p-1)*psi)+1):(n+beta.len+((p-1)*psi)+n)]))))
}
colnames(data.scenario) <- c("x", "constant", "post.mean",
                              "trueAlp", "meanAlp",
                              "post.check", paste("alp", 1:len, sep = ""),
                              # paste0("post.samp", 1:len))
                              paste0("post.samp.alp", 1:len))
tail(data.scenario[, c(1:10, ((dim(samples)[1]*2-5):(dim(samples)[1]*2)))], 15)
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))
#plotting all the points
plt <- ggplot(data = data.scenario, aes(x = x)) + ylab("alpha(x)") + xlab("")
for(i in (dim(samples)[1] - 993):(dim(samples)[1]+6)){
  plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}
print(plt + geom_line(aes(y=post.mean, col = "Posterior Mean(Chain2)"), linewidth = 1.5) + 
        geom_line(aes(y = post.check, col=paste0("Simulated Alpha: ",n,"/",psi,"/",threshold)), linewidth = 1.5) +
        # theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        # theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        labs(col = "") + 
        scale_color_manual(values = c("#e0b430","red")) +
        theme(text = element_text(size = 27)) + 
        theme(legend.position="top", legend.key.size = unit(1, 'cm')))
# ggsave(paste0("./Simulation/BayesianPsplines/results/figures/",date,"-",time, "_sc1_Alp_samp_1.pdf"), 
#         width=14, height = 7.5)

cat("sc2_Alp Done")

plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab("alpha(x)")
for(i in ((dim(samples)[1]*2)-993):((dim(samples)[1]*2)+6)){
  # print(i)
  plt.samp <- plt.samp + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}
print(plt.samp + ylim(0, 3) +
      geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
      geom_line(aes(y = meanAlp, col = "Posterior Mean(Chain2)"), linewidth = 2.5) +
      labs(col = "") +
      scale_color_manual(values = c("#e0b430", "red"))+
      # scale_colour_manual("", 
      #               breaks = c("True Alpha", "Posterior Mean"),
      #               values = c("True Alpha"="red", "Posterior Mean"="#e0b430"),
      #               guide = "legend") +
      theme(text = element_text(size = 27)) + 
      theme(legend.position="top", legend.key.size = unit(1, 'cm')))


# ggsave(paste0("./Simulation/BayesianPsplines/results/figures/",date,"-",time, "_true_Alp_samp_2.pdf"), 
#         width=14, height = 7.5)
# print(plt.samp + 
#       geom_line(aes(y = trueAlp), colour = "red", linewidth = 1.5))

cat("true_Alp Done")

samples <- fit.v2$samples$chain3
data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "post.mean" = sort(fit.v2$summary$chain3[1:n,1]),
                            "trueAlp" = sort(trueAlp),
                            "meanAlp" = sort(fit.v2$summary$chain3[(n+beta.len+((p-1)*psi)+1):(n+beta.len+n+((p-1)*psi)),1]),
                            "post.check" = sort(cutoff.alp))
                            # "post.check" = sort(alp))
# data.scenario1 <- data.frame("x"=c(1:n)) #, )

len <- dim(samples)[1]
for(i in 1:len){
  data.scenario <- cbind(data.scenario, 
                      data.frame(unname(sort(samples[i, 1:n]))))
}
for(i in 1:len){
  data.scenario <- cbind(data.scenario, 
                      data.frame(unname(sort(samples[i, (n+beta.len+((p-1)*psi)+1):(n+beta.len+((p-1)*psi)+n)]))))
}
colnames(data.scenario) <- c("x", "constant", "post.mean",
                              "trueAlp", "meanAlp",
                              "post.check", paste("alp", 1:len, sep = ""),
                              # paste0("post.samp", 1:len))
                              paste0("post.samp.alp", 1:len))
tail(data.scenario[, c(1:10, ((dim(samples)[1]*2-5):(dim(samples)[1]*2)))], 15)
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))
#plotting all the points
plt <- ggplot(data = data.scenario, aes(x = x)) + ylab("alpha(x)") + xlab("")
for(i in (dim(samples)[1] - 993):(dim(samples)[1]+6)){
  plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}
print(plt + geom_line(aes(y=post.mean, col = "Posterior Mean(Chain3)"), linewidth = 1.5) + 
        geom_line(aes(y = post.check, col=paste0("Simulated Alpha: ",n,"/",psi,"/",threshold)), linewidth = 1.5) +
        # theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        # theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        labs(col = "") + 
        scale_color_manual(values = c("#e0b430","red")) +
        theme(text = element_text(size = 27)) + 
        theme(legend.position="top", legend.key.size = unit(1, 'cm')))
# ggsave(paste0("./Simulation/BayesianPsplines/results/figures/",date,"-",time, "_sc1_Alp_samp_1.pdf"), 
#         width=14, height = 7.5)

cat("sc3_Alp Done")

plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab("alpha(x)")
for(i in ((dim(samples)[1]*2)-993):((dim(samples)[1]*2)+6)){
  # print(i)
  plt.samp <- plt.samp + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}
print(plt.samp + ylim(0, 3) +
      geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
      geom_line(aes(y = meanAlp, col = "Posterior Mean(Chain3)"), linewidth = 2.5) +
      labs(col = "") +
      scale_color_manual(values = c("#e0b430", "red"))+
      # scale_colour_manual("", 
      #               breaks = c("True Alpha", "Posterior Mean"),
      #               values = c("True Alpha"="red", "Posterior Mean"="#e0b430"),
      #               guide = "legend") +
      theme(text = element_text(size = 27)) + 
      theme(legend.position="top", legend.key.size = unit(1, 'cm')))

# ggsave(paste0("./Simulation/BayesianPsplines/results/figures/",date,"-",time, "_true_Alp_samp_2.pdf"), 
#         width=14, height = 7.5)

cat("true_Alp Done")

data.scenario <- data.frame("meanAlp" = alpha.summary[(n+beta.len+((p-1)*psi)+1):(n+beta.len+n+((p-1)*psi)),1],
                            "constant" = newx,
                            "trueAlp" = sort(trueAlp))
plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab("alpha(x)")
print(plt.samp + ylim(0, 3) +
      geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
      geom_line(aes(y = meanAlp, col = "Posterior Mean"), linewidth = 2.5) +
      labs(col = "") +
      scale_color_manual(values = c("#e0b430", "red"))+
      # scale_colour_manual("", 
      #               breaks = c("True Alpha", "Posterior Mean"),
      #               values = c("True Alpha"="red", "Posterior Mean"="#e0b430"),
      #               guide = "legend") +
      theme(text = element_text(size = 27)) + 
      theme(legend.position="top", legend.key.size = unit(1, 'cm')))

# print(plt.samp + 
#       geom_line(aes(y = trueAlp), colour = "red", linewidth = 1.5))

#--------------------------------------------------------------------------------------
#plotting bsplines basis function of covariates
# plts <- list()
# for(j in 2:p){
#   x.x <- sort(scale(x[,j]))
#   # # Make a matrix containing the B-spline basis
#   ndx <- 12
#   basis.x <- data.frame(x=x.x, y, id = as.factor(ndx))
#   deg <- 3
#   B <- bbase(x.x, min(x.x), max(x.x), nseg = ndx, bdeg = deg)
#   nb1 <- ncol(B)

#   # A basis for plotting the fit on the grid xg
#   ng <- length(x.x)
#   xg <- seq(min(x.x), max(x.x), length.out = ng)
#   Bg <- bbase(xg, min(x.x), max(x.x), nseg = ndx, bdeg = deg)

#   a <- solve(t(B) %*% B, t(B) %*% y, tol = 1e-20)
#   z <- Bg %*% a

#   # Make a matrix with B-splines scaled by coefficients
#   Bsc1 <- Bg %*% diag(c(a))
#   # Create data frames for ggplot
#   Zf1 <- data.frame(x = xg, y = z, id = as.factor(1))
#   Bf1 <- data.frame(x = rep(xg, nb1), y = as.vector(Bsc1), id = as.factor(rep(1:nb1, each = ng)))
#   # Bf1$y[abs(Bf1$y) < 0.0001] = NA
#   # Bf1 <- na.omit(Bf1)

#   # Build the graphs
#   plt <- ggplot(Bf1 , aes(x = x, y = y, group = id, colour = id)) +
#     ggtitle(paste("Basis Function of ", colnames(x)[j])) +
#     # geom_hline(yintercept = 0, linewidth = 0.3) +
#     geom_line(data = Zf1, linewidth = 1, colour = "grey") +
#     # geom_point(data = basis.x, colour = "grey60", size = 0.8, shape=1) +
#     geom_line(linewidth = 0.7) +
#     # geom_point(data = basis.x, colour = "red", size = 2, shape = 1) +
#     xlab("") + ylab ("") +
#     JOPS_theme() +
#     theme(legend.position = "none") +
#     scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

#   plts[[(j-1)]] <- plt
#   # ggsave(paste0("./Results/vert_bsplines_",j,".png"))
# }
# plts
# grid.arrange(grobs = plts, ncol = 2, nrow = 4)