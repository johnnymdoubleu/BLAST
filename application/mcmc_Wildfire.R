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


# ggsave("./Laboratory/Application/figures/datahist.pdf", width=15)

    # scale_color_gradientn(colours = rainbow(5))
# ggplot(data = fwi.scaled) + 
#   geom_histogram(aes(x=DSR), bins = 30) +
#   geom_histogram(aes(x=FWI), bins = 30) +
#   geom_histogram(aes(x=BUI), bins = 30) +
#   geom_histogram(aes(x=ISI), bins = 30) +
#   geom_histogram(aes(x=FFMC), bins = 30) + 
#   geom_histogram(aes(x=DMC), bins = 30) +
#   geom_histogram(aes(x=DC), bins = 30)

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
  xholder[,i] <- seq(min(fwi.scaled), max(fwi.scaled), length.out = n)
  test.knot <- seq(min(fwi.scaled), max(fwi.scaled), length.out = psi)
  splines <- basis.tps(seq(min(fwi.scaled), max(fwi.scaled), length.out = n), test.knot, m=2, rk=FALSE, intercept = TRUE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  knots <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  tps <- basis.tps(fwi.scaled[,i], knots, m = 2, rk = FALSE, intercept = TRUE)
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
        if(a > 10){
          ans <- exp(10) + exp(10)*(a-10)
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
    sigma.square[j] ~ dinvgamma(0.1, 0.1)
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
    # alpha[i] <- log(10) / log(1 + exp(theta.0 + sum(g.nonlinear[i, 1:p]) + sum(g.linear[i, 1:p])))
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
                                    theta = rep(0.2, p), theta.0 = 0.1,
                                    covm = array(1, dim = c(psi,psi,10))),
                              list(gamma = matrix(0, nrow = psi, ncol=p),
                                    theta = rep(0, p), theta.0 = 0,
                                    covm = array(1, dim = c(psi,psi,10))))
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
                  niter = 30000,
                  nburnin = 10000,
                  # setSeed = 300,
                  nchains = 2,
                  # WAIC = TRUE,-
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE)

alpha.summary <- fit.v2$summary$all.chains
# saveRDS(alpha.summary, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_allChains.rds"))

alpha.summary[701:711,]

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
# MCMCsummary(object = fit.v2$samples, round = 3)

print(MCMCtrace(object = fit.v2$samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          n.eff = TRUE,# add eff sample size
          params = c("lambda.1", "lambda.2")))

samples <- fit.v2$samples$chain1
data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "post.mean" = sort(fit.v2$summary$chain1[1:n,1]))
# data.scenario1 <- data.frame("x"=c(1:n)) #, )

len <- dim(samples)[1]
for(i in 1:len){
  data.scenario <- cbind(data.scenario, 
                      data.frame(unname(sort(samples[i, 1:n]))))
}
# for(i in 1:len){
#   data.scenario <- cbind(data.scenario, 
#                       data.frame(unname(sort(samples[i, 703:1202]))))
# }
colnames(data.scenario) <- c("x", "constant", "post.mean",
                              paste("alp", 1:len, sep = ""))
tail(data.scenario[, c(1:10, ((dim(samples)[1]-5):(dim(samples)[1])))], 15)
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))
#plotting all the points
plt <- ggplot(data = data.scenario, aes(x = x)) + ylab(expression(alpha(x))) + xlab(expression(x))
for(i in (dim(samples)[1] - 996):(dim(samples)[1]+3)){
  plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}

print(plt + geom_line(aes(y=post.mean, col = "Posterior Mean(Chain1)"), linewidth = 1.5) + 
        # geom_line(aes(y = post.check, col=paste0("Simulated Alpha: ",n,"/",psi,"/",threshold)), linewidth = 1.5) +
        # theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        # theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        labs(col = "") + 
        scale_color_manual(values = c("#e0b430")) +
        theme(text = element_text(size = 27)) + 
        theme(legend.position="top", legend.key.size = unit(1, 'cm')))
# ggsave(paste0("./Simulation/BayesianPsplines/results/figures/",date,"-",time, "_sc1_Alp_samp_1.pdf"), 
#         width=14, height = 7.5)

cat("sc1_Alp Done")

covariates <- gl(p, n, (p*n), labels = factor(names(fwi.scaled)))
replicate <- gl(2, n, (p*n))
func.df <- data.frame(seq = seq(0,1,length.out = n),
                        # x = as.vector(apply(fwi.scaled, 2, sort, method = "quick")),
                        # new=as.vector(func.new),
                        new.linear=as.vector(fit.v2$summary$chain1[(n+1):(n+(n*p)),1]),
                        new.nonlinear=as.vector(fit.v2$summary$chain1[(n+(n*p)):((n+(n*p*2))),1]),
                        covariates=covariates, 
                        replicate=replicate,
                        x = rep(newx, p))

plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x))) + xlab(expression(x))
for(i in ((dim(samples)[1]*2)-993):((dim(samples)[1]*2)+6)){
  # print(i)
  plt.samp <- plt.samp + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
}
print(plt.samp + ylim(0, 50) +
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