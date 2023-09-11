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
# suppressMessages(library(coda))
# library(R6)
# suppressMessages(library(igraph))
# library(mgcv)
library(MCMCvis)
suppressMessages(library(tidyverse))
# library(ggplotify)

#Scenario 1
set.seed(10)
n <- 5000
psi <- 20

p <- 10
no.theta <- 2
simul.no <- 50


theta.container <- as.data.frame(matrix(, nrow = (no.theta *p), ncol= simul.no))
gamma.container <- as.data.frame(matrix(, nrow = (psi * p), ncol = simul.no))
# linear.container <- nonlinear.container <- f.container <- lapply(1:simul.no, matrix, data= NA, nrow=(n*(1-threshold)), ncol=p)
linear.container <- nonlinear.container <- f.container <- lapply(1:simul.no, data.frame)
alpha.container <- as.data.frame(matrix(, nrow=n, ncol = simul.no))
#Scenario 1
n <- 5000

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
x.origin <- cbind(replicate(p, runif(n, 0, 1)))
for(i in 1:p){
    knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
    tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = TRUE)
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

theta.origin <- matrix(, nrow = no.theta, ncol = p)
for(j in 1:p){
    for (k in 1:no.theta){
        if(j %in% c(2,4,5,6,9,10)){theta.origin[k, j] <- 0}
        else if(j==7){
            if(k==1){theta.origin[k, j] <- 0.5}
            else{theta.origin[k, j] <- -0.3}
        }
        else {
            if(k==1){theta.origin[k,j] <- -0.2}
            else{theta.origin[k,j] <- 0.8}
        }
    }
}


n <- 1000
xholder <- bs.x <- matrix(,nrow=n, ncol=0)
psi <- 20
p <- 10
gamma.origin <- matrix(, nrow= psi, ncol=p)
x.origin <- cbind(replicate(p, runif(n, 0, 1)))
# x.origin <- scale(x.origin)
newx <- seq(0, 1, length.out=n)



for(i in 1:p){
  # splines <- bbase(x.origin[, i], min(x.origin[,i]), max(x.origin[,i]), nseg = (psi-3), bdeg = 3)
#   splines <- bbase(x.scale[, i], min(x.scale[, i]), max(x.scale[, i]), nseg = (psi-3), bdeg = 3)
  test.knot <- seq(0, 1, length.out = (psi-1))
  splines <- basis.tps(newx, test.knot, m=2, rk=FALSE)
  xholder <- cbind(xholder, splines)
  knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = (psi-1))
  tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE)
  # tps <- mSpline(x.origin[,i], df=psi, Boundary.knots = range(x.origin[,i]), degree = 3, intercept=TRUE)
  bs.x <- cbind(bs.x, tps)
}

gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
    for (ps in 1:psi){
        # gamma.origin[ps, j] <- 1
        if(j %in% c(2,4,5,6,9,10)){
            gamma.origin[ps, j] <- 0
        }
        else if(j==7){
            if(ps <= (psi/2)){
                gamma.origin[ps, j] <- 1
            }
            else{
                gamma.origin[ps, j] <- 1
            }
        }
        else {
            if(ps <= (psi/2)){
                gamma.origin[ps, j] <- 1
            }
            else{
                gamma.origin[ps, j] <- 1
            }
        }
    }
}
f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.origin[, j] <- bs.x[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
  alp.origin[i] <- exp(sum(f.origin[i,]))
  y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i])
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
  # tau ~ dinvgamma(shape, scale)
  # w ~ dnorm(0, tau)
  # beta[1] ~ dnorm(0, 0.001)
  # lambda ~ dgamma(0.1, 0.1) #gamma distribution prior for lambda
  for (j in 1:p){
    lambda.1[1,j] ~ dgamma(shape, scale) 
    gamma[1, j] ~ ddexp(0, lambda.1[1,j])
    # gamma[2, j] ~ dnorm(0, 1)
    tau[j] ~ dinvgamma(shape, scale)
    # w[1,j] ~ dnorm(0, tau[j])
    # w[2,j] ~ dnorm(0, tau[j])
    # w[j] ~ dnorm(0, sqrt(tau[j]))
    for (i in 2:psi){
      w[i,j] ~ dnorm(0, tau[j])
      lambda.1[i,j] ~ dgamma(shape, scale)
      gamma[i,j] <- gamma[(i-1),j] + w[i,j]
      # gamma[i, j] <- (2*gamma[(i-1),j]) - gamma[(i-2),j] + w[i,j]
    }
    # lambda[j] ~ dgamma(0.1, 0.1) #gamma distribution prior for lambda
    # beta[j] ~ ddexp(0, lambda) #laplace distribution for beta
    # beta[j] ~ dlaplace(0, lambda)
  }
  # lambda ~ dgamma(0.1, 0.1) #gamma distribution prior for lambda
#   beta[1] ~ dnorm(0, 0.001)
  # for(j in 1:p){
  #   lambda.1[j] ~ dgamma(0.1, 0.1) #gamma distribution prior for lambda
  #   beta[j] ~ ddexp(0, lambda.1[j])
  # }
  # Likelihood

  # for (j in 1:n){
  #   g[j] <- inprod(x[j, 1:p], beta[1:p])
  #   new.g[j] <- inprod(new.x[j, 1:p], beta[1:p])
  # }

  for (j in 1:p){
    f[1:n, j] <- bs.x[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[1:psi, j] 
    new.f[1:n, j] <- new.bs.x[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[1:psi,j]
  }
  # log(alpha[1:n]) <- sum(f.sum[1:n, 1:(p-1)])
  # log(newalpha[1:n]) <- sum(new.f.sum[1:n, 1:(p-1)])
  for(i in 1:n){
    # g[i] <- inprod(x[i,1:p], beta[1:p])
    log(alpha[i]) <- sum(f[i, 1:p])
    log(newalpha[i]) <- sum(new.f[i, 1:p])
  }
  # newalpha[1:n] <- nimRowSums(new.f.sum[1:n,1:(p-1)])
c
  for(i in 1:n){
    y[i] ~ dpareto(1, u, alpha[i])
    # spy[i] <- (alpha[i]*(y[i]/u)^(-1*alpha[i])*y[i]^(-1)) / C
    # ones[i] ~ dbern(spy[i])
  }
})


constant <- list(psi = psi, n = n, p = p)
init.alpha <- function() list(list(gamma = matrix(0.5, nrow = psi, ncol=p)),
                              list(gamma = matrix(0, nrow = psi, ncol=p)),
                              list(gamma = matrix(1, nrow = psi, ncol=p)))
                              # y = as.vector(y),
monitor.pred <- c("gamma", "alpha", "newalpha")
data <- list(y = as.vector(y), bs.x = bs.x, x = x.scale,
               new.x = xholder, new.bs.x = new.bs.x,
              u = u, #C = 1000,  ones = as.vector(rep(1, n)),
              shape = 0.1, scale = 0.1)

fit.v2 <- nimbleMCMC(code = model.penalisation,
                  constants = constant,
                  data = data,
                  monitors = monitor.pred,
                  inits = init.alpha(),
                  thin = 20,
                  niter = 30000,
                  nburnin = 10000,
                  # setSeed = 300,
                  nchains = 3,
                  # WAIC = TRUE,-
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE)

# alpha.summary <- MCMCsummary(object = fit.v2$samples)
# saveRDS(fit.v2, file=paste0("./Simulation/BayesianPsplines/results/",date,"-",time,"_sc1_fit.rds"))
alpha.summary <- fit.v2$summary$all.chains
# saveRDS(alpha.summary, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_allChains.rds"))

MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
            HPD = TRUE, xlab="beta", offset = 0.05, exact = TRUE,
            horiz = FALSE, params = c("beta"))
MCMCplot(object = fit.v2$samples$chain1, object2 = fit.v2$samples$chain2,
            HPD = TRUE, xlab="gamma", offset = 0.5,
            horiz = FALSE, params = c("gamma"))
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
gelman.diag(fit.v2$samples, multivariate = FALSE)
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
#           params = "beta")

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
                            "trueAlp" = sort(trueAlp),
                            "meanAlp" = sort(fit.v2$summary$chain1[(n+beta.len+((p-1)*psi)+1):(n+beta.len+n+((p-1)*psi)),1]),
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

plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab("expression(alpha)(x)")
for(i in ((dim(samples)[1]*2)-993):((dim(samples)[1]*2)+6)){
  # print(i)
  plt.samp <- plt.samp + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}
print(plt.samp + ylim(0, 3) +
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