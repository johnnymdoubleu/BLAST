library(npreg)
library(Pareto)
library(tidyverse)
library(JOPS)
library(readxl)
library(gridExtra)
library(colorspace)
library(corrplot)
library(ReIns)
library(evir)

# Structure of the FWI System
#DSR : Dail Severity Rating
#FWI : Fire Weather Index
#BUI : Buildup Index
#ISI : Initial Spread Index
#FFMC : Fine FUel Moisture Code
#DMC : Duff Moisture Code
#DC : Drough Code


setwd("C:/Users/Johnny Lee/Documents/GitHub")
df <- read_excel("./Laboratory/Application/AADiarioAnual.xlsx", col_types = c("date", rep("numeric",40)))
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
path <- "./Laboratory/Application/DadosDiariosPT_FWI.xlsx"
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
corrplot.mixed(cor(fwi.scaled),
                upper = "circle",
                lower = "number",
                addgrid.col = "black")
# ggsave("./Laboratory/Application/figures/correlation.pdf", width=15)
# cov$date <- as.Date(with(cov, paste(year,month,day,sep="-")),"%Y-%m-%d")
# cov$yearmon <- as.Date(with(cov, paste(year,month,sep="-")),"%Y-%m")
df.extreme <- cbind(y, fwi.scaled)
# df.extreme <- cbind(date = cov$date[which(Y>u)], df.extreme)
df.extreme <- cbind(month = fwi.index$month[which(Y>u)], df.extreme)
ggplot(df.extreme, aes(x=month, y=y, color=month)) + geom_point(size=6) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust=0.35),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 15))
ggsave("./Laboratory/Application/figures/datavis.pdf", width=15)

# ggplot(as.data.frame(y),aes(x=y))+geom_histogram(aes(y=..density..), bins = 5)

ggplot(df.extreme, aes(x=y)) +
    geom_histogram(stat = "density", n=40, adjust=0.1, fill = "darkgrey") + 
    xlab("Area Burnt") + 
    # geom_histogram(aes(y=..density..), bins = 10^3) +
    # geom_density(aes(y=..density..)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
      legend.title = element_blank(),
      legend.text = element_text(size=20),
      # axis.ticks.x = element_blank(),
      axis.text.x = element_text(hjust=0.35),
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 30))
ggsave("./Laboratory/Application/figures/datahist.pdf", width=15)

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
no.theta <- 2
newx <- seq(0,1,length.out=n)
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

# cov <- read_excel("./Laboratory/Application/WTS_ERA_5_ate2019_final(PORTUGAL).xlsx", , col_types = rep("numeric",11))
# # cov <- read_excel("WTS_ERA_5_ate2019_final(PORTUGAL).xlsx", col_types = rep("numeric",11))
# cov <- cov[-which(cov$year==1999&cov$month==2&cov$day==14),]
# dim(cov)
# # View(cov)
# head(cov)
# tail(cov)
# cov$date <- as.Date(with(cov, paste(year,month,day,sep="-")),"%Y-%m-%d")
# cov$yearmon <- as.Date(with(cov, paste(year,month,sep="-")),"%Y-%m")
# # View(summary(cov))

# x.cov <- cov[cov$year>=1980, ]
# df.whole <- cbind(y=Y, x.cov)
# #x<-cov[cov$ano>=1980, c(5:7)]
# fwi.scaled <- scale(x.cov[,c(5:11)])[which(Y>u),]
# x <- x.cov[which(Y>u), c(5:11)]
# df.extreme <- cbind(y, x)
# df.extreme <- cbind(date = cov$date[which(Y>u)], df.extreme)
# df.extreme <- cbind(month = cov$month[which(Y>u)], df.extreme)

# n <- dim(fwi.scaled)[[1]]
# p <- dim(fwi.scaled)[[2]]

# xholder <- bs.x <- matrix(,nrow=n, ncol=0)
# for(i in 1:p){
#   knots <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = (psi-1))
#   tps <- basis.tps(fwi.scaled[,i], knots, m = 2, rk = FALSE)
#   # tps <- mSpline(x.origin[,i], df=psi, Boundary.knots = range(x.origin[,i]), degree = 3, intercept=TRUE)
#   bs.x <- cbind(bs.x, tps)
# }

lambda.1 <- 0.001
lambda.2 <- 0
lambda.3 <- 0.001

log.posterior <- function(beta, y.origin){
    log.lik <- function(beta){
        exp.prime <- function(x, thres){
            if(x > thres){ans <- exp(thres) + exp(thres)*(x-thres)}
            else{ans <- exp(x)}
            return(ans)
        }
        theta <- matrix(beta[1:(no.theta*p)], nrow=no.theta)
        gamma <- matrix(beta[(no.theta*p)+1:(psi*p)], ncol=p)
        f <- matrix(, nrow=n, ncol=p)
        term <- first.term <- second.term <- NULL
        for(j in 1:p){
            # coef <- as.matrix(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])
            linear.term <- as.matrix(bs.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta[,j]
            nonlinear.term <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[,j]
            f[, j] <- linear.term + nonlinear.term
        }
        for(i in 1:n){
            first.term[i] <- sum(f[i,]) - log(y.origin[i])
            second.term[i] <- exp.prime(sum(f[i,]), thres = 10) * log(y.origin[i]/u)
            term[i] <- first.term[i] - second.term[i]
        }
        return(sum(term))
    }
    log.prior <- function(beta){
        moreau.envelope <- function(w){
            if(w < -1){ans <- -0.5 - w}
            else if (1 < w){ans <- w - 0.5}
            else {ans <- (w^2) / 2}
            return(ans)
        }
        theta <- matrix(beta[1:(no.theta*p)], nrow=no.theta)
        gamma <- matrix(beta[(no.theta*p)+1:(psi*p)], ncol=p)
        g.1 <- g.2 <- term <- third.term <- first.term <- second.term <- NULL
        for(j in 1:p){
            first.term[j] <- -1 * lambda.1 * sqrt(sum((gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])^2))
            # first.term[j] <- -1 * lambda.1 * abs(sum(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)]))
            second.term[j] <- -1 * lambda.2 * sum((gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])^2)
            third.term[j] <- -1 * lambda.3 * sum(abs(theta[,j]))
            term[j] <- first.term[j] + second.term[j] + third.term[j]
            # term[j] <- first.term[j] + second.term[j] - lambda.3 * sum(abs(beta))
        }
        return(sum(term))
    }
    return(log.lik(beta) + log.prior(beta))
}

beta.emp <- c(rep(0, no.theta*p), rep(0, p*psi))
# beta.emp <- c(as.vector(theta.origin), as.vector(gamma.origin))
beta.map <- optim(beta.emp, fn = log.posterior, #gr = grad.log.posterior, 
                  y.origin = y,
                  # method = "BFGS",
                  method = "CG",
                  # method = "SANN",
                  control = list(fnscale = -1))
# theta.map <- matrix(beta.map$par[1:(2*p)],nrow=2)
theta.map <- beta.map$par[1:(no.theta*p)]
gamma.map <- beta.map$par[(no.theta*p)+1:(psi*p)]
df.theta <- data.frame("seq" = seq(1, p),
                  theta.map = matrix(beta.map$par[1:(no.theta*p)],nrow=2)[2,])
# df.theta$covariate <- factor(rep(seq(1, 1 + nrow(df.theta) %/% no.theta), each = no.theta, length.out = nrow(df.theta)))
# df.theta$covariate <- factor(rep(names(fwi.scaled), each = no.theta, length.out = nrow(df.theta)))
df.theta$covariate <- factor(names(fwi.scaled), levels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC"))
df.theta$labels <- factor(1:p)
# ggplot(df.theta, aes(x = seq)) + 
#   geom_point(aes(y = theta.map, color = covariate), size = 1.5) + 
# #   geom_smooth(method="gam") +
#   # geom_point(data = df.theta, aes(y = theta.true, color = "true"), size = 2.5) +
#   labs(title=expression("MAP for"~theta)) + 
#   # ggtitle(expression(atop(paste("MAP vs True for ", bold(gamma))))) +
# #   annotate("text", x = seq(0, 330, length.out=10), y = -1, label = beta, colour = "red", size = 10) +   
# #   scale_color_manual(values = c("true"="black"))+
#   theme(plot.title = element_text(hjust = 0.5, size = 20))

ggplot(df.theta, aes(x = labels)) + ylab("") + 
  geom_point(aes(y = theta.map, color = covariate), size = 6) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + ylim(-0.5,0.5) + xlab('') + 
  # geom_point(aes(y = theta.true, color = "true"), size = 2.5) +
  # labs(title=expression("MAP vs True for"~theta)) + xlab("") +
  scale_x_discrete(labels = c(expression(bold(theta[1])),
                              expression(bold(theta[2])),
                              expression(bold(theta[3])),
                              expression(bold(theta[4])),
                              expression(bold(theta[5])),
                              expression(bold(theta[6])),
                              expression(bold(theta[7]))))+
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size=30),
          # axis.ticks.x = element_blank(),
          axis.text.x = element_text(hjust=0.35),
          axis.text = element_text(size = 30),
          panel.grid.minor.x = element_blank())
ggsave("./Laboratory/Application/figures/map_theta.pdf", width=10)
df <- data.frame("seq" = seq(1, (psi*p)), 
                  gamma.map)
# df$covariate <- factor(rep(seq(1, 1 + nrow(df) %/% psi), each = psi, length.out = nrow(df)))
df$covariate <- factor(rep(names(fwi.scaled), each = psi, length.out = nrow(df)), levels = c("DSR", "FWI", "BUI", "ISI", "FFMC", "DMC", "DC"))
df$labels <- factor(1:(psi*p))
ggplot(df, aes(x =labels , y = gamma.map, color = covariate)) + 
  geom_point(size = 4) + ylab("") + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("")+
  # geom_smooth(method="gam") +
  # geom_point(aes(y = gamma.true, color = "true")) + 
  # geom_line(aes(x = seq, y = true, color = "true"), linetype = 2) +
  # labs(title=expression("MAP vs True for"~gamma)) + 
  # ggtitle(expression(atop(paste("MAP vs True for ", bold(gamma))))) +
#   annotate("text", x = seq(0, 330, length.out=10), y = -1, label = beta, colour = "red", size = 10) +
  scale_x_discrete(labels = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[1])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
  expression(bold(gamma[2])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[3])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[4])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[5])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[6])),"", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", expression(bold(gamma[7])))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size=30),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(hjust=1.75),
          axis.text = element_text(size = 30),
          panel.grid.major.x = element_blank())
ggsave("./Laboratory/Application/figures/map_gamma.pdf", width=10)
# ggplot(df, aes(x = seq, y = gamma.map, color = covariate)) + 
#   geom_point() + 
#   # geom_smooth(method="gam") +
#   # geom_point(aes(x=seq, y = gamma.true, color = "true")) + 
#   # geom_line(aes(x = seq, y = true, color = "true"), linetype = 2) +
#   labs(title=expression("MAP for"~gamma)) + 
#   # ggtitle(expression(atop(paste("MAP vs True for ", bold(gamma))))) +
# #   annotate("text", x = seq(0, 330, length.out=10), y = -1, label = beta, colour = "red", size = 10) +   
#   theme(plot.title = element_text(hjust = 0.5, size = 20))

f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
newalpha <- NULL
for (j in 1:p){
  f.linear.new[,j] <- as.matrix(bs.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% matrix(theta.map, nrow=no.theta)[1:no.theta, j]
  f.nonlinear.new[,j] <- bs.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.new[1:n, j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}


# set.seed(100)
for(i in 1:n){
  temp <- sum(f.new[i,])
  newalpha[i] <- exp(temp)
}

f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
for (j in 1:p){
  f.linear.new[,j] <- as.matrix(xholder.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% matrix(theta.map, nrow=no.theta)[1:no.theta, j]
  f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.new[1:n, j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}

alpha.new <- NULL
for(i in 1:n){
  alpha.new[i] <- exp(sum(f.new[i,]))
}

plot(1/Hill(sort(Y)[13865:14609])$gamma)

plot(sort(hill(y,option="alpha", reverse = TRUE)$y))
hill(y, option = "alpha", reverse = TRUE)
hill(sort(Y)[13865:14609], option="alpha", reverse = TRUE)$y

# data.scenario <- data.frame("x" = c(1:n),
#                             "constant" = seq(0,1,length.out = n),
#                             "mapAlp" = sort(newalpha),
#                             "newAlp" = sort(alpha.new),
#                             "hill.est" = sort(hill(sort(Y)[13865:14609], option="alpha", reverse = TRUE)$y),
#                             "fakehill" = sort(1/Hill(sort(Y)[13878:14609])$gamma))

# plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x))) + xlab(expression(x[(i)]))

# print(plt.samp + 
#       # geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
#       geom_line(aes(y = hill.est, col = "Hill's Estimator:"), linewidth = 2.5) +
#       labs(col = "") +
#       geom_line(aes(y = mapAlp, col = paste0("MAP:",lambda.1,"/",lambda.3)), linewidth = 2.5) +
#       scale_color_manual(values = c("red","#e0b430"))+ 
#       theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
#       theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
#       theme(text = element_text(size = 15),
#         legend.position="bottom", legend.key.size = unit(1, 'cm'),
#         axis.text = element_text(size = 20),))
#       # theme(legend.position="bottom", legend.key.size = unit(1, 'cm')))
# ggsave("./Laboratory/Application/figures/map_alpha.pdf", width=10)
# print(plt.samp + 
#       # geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
#       geom_line(aes(y = newAlp, col = paste0("MAP Alpha:",lambda.1,"/",lambda.3)), linewidth = 2.5) +
#       labs(col = "") +
#       scale_color_manual(values = c("red"))+
#       theme(text = element_text(size = 30)) + 
#       theme(legend.position="bottom", legend.key.size = unit(1, 'cm')))


func.linear.new <- func.nonlinear.new <- func.new <- matrix(, nrow=n, ncol=0)
for (j in 1:p){
  func.new <- cbind(func.new, f.new[,j])
  func.linear.new <- cbind(func.linear.new, f.linear.new[,j])
  func.nonlinear.new <- cbind(func.nonlinear.new, f.nonlinear.new[,j])  
}

covariates <- gl(p, n, (p*n), labels = factor(names(fwi.scaled)))
replicate <- gl(2, n, (p*n))
func.df <- data.frame(seq = seq(0,1,length.out = n),
                        # x = as.vector(apply(fwi.scaled, 2, sort, method = "quick")),
                        new=as.vector(func.new),
                        new.linear=as.vector(func.linear.new),
                        new.nonlinear=as.vector(func.nonlinear.new),
                        covariates=covariates, 
                        replicate=replicate,
                        x = rep(newx, p))

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Smooth Functions") +
  geom_line(aes(y=new, colour = covariates), linewidth=2) + ylab ("") +
  facet_grid(covariates ~ .) + #ggtitle("MAP for Smooth Functions") + 
  scale_y_continuous(breaks=c(-10, 0, 10)) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
ggsave("./Laboratory/Application/figures/map_smooth.pdf", width=10.5, height = 15)
ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) +  
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Linear Component") + 
  geom_line(aes(y=new.linear, colour = covariates), linewidth=2) + ylab ("") +
  facet_grid(covariates ~ .) + #ggtitle("Linear Component of Smooth Functions") + 
  scale_y_continuous(breaks=c(0)) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.ticks = element_blank(),
        # panel.grid.minor.y = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 35))
ggsave("./Laboratory/Application/figures/map_linear.pdf", width=10, height = 15)
ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Nonlinear Component") +
  geom_line(aes(y=new.nonlinear, colour = covariates), linewidth=2) + ylab ("") +
  facet_grid(covariates ~ .) + #ggtitle("Nonlinear Component of Smooth Functions") + 
  scale_y_continuous(breaks=c(0)) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=33),
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 35))
ggsave("./Laboratory/Application/figures/map_nonlinear.pdf", width=12.5, height = 15)




# plot(sort(alp.origin))
# plot(sort(new.y), sort(y.origin))
# abline(a=0, b=1, col = "red", lty = 2)
rbind(matrix(theta.map, nrow = no.theta, ncol = p), matrix(gamma.map, nrow = psi, ncol = p))
# Randomized quantile residuals
r <- matrix(, nrow = n, ncol = 20)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 20
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y[i], u, alpha = newalpha[i]))
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
  theme(text = element_text(size = 20)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3))
ggsave("./Laboratory/Application/figures/map_qqplot.pdf", width=10)

model.df <- data.frame(y=y, fwi.scaled)
mod <- lm(log(y) ~ DSR + FWI + BUI + ISI + FFMC + DMC + DC, data= model.df)
mod.summary <- summary(mod)
mse.linear <- sqrt(mean(mod.summary$residuals^2))
beta.linear <- unname(mod$coefficients)
alpha.lr <- NULL

# for (i in 1:n){
#   alpha.lr[i] <- exp(sum(fwi.scaled[i,] * beta.linear))
# }
# y.lr <- qPareto(pPareto(y, u, alpha = alpha.lr), u, alpha = alpha.lr)
# mse.linear <- mean(sqrt((y-y.lr)^2))

# new.y <- NULL
new.y <- y.pareto <- NULL
y.pareto <- rPareto(n, u, alpha = newalpha)
# y.pareto <- qPareto(pPareto(y, 1, alpha = newalpha), 1, alpha = newalpha)
# for(i in 1:n){
# #   y.pareto[i] <- rPareto(1, 1, alpha = newalpha[i])
#     # new.y[i] <- qPareto(pPareto(y[i], 1, alpha = newalpha[i]), 1, alpha = newalpha[i])
#     new.y[i] <- qPareto(threshold, 1, alpha = newalpha[i])
# }
new.y <- qPareto(pPareto(y, u, alpha = newalpha), u, alpha = newalpha)
mse.pareto.1 <- mean(sqrt((y-y.pareto)^2))
mse.pareto.2 <- mean(sqrt((y-new.y)^2))
c(mse.pareto.1, mse.pareto.2, mse.linear)

# rPareto(, 1, alpha = newalpha)

# mean.y <- NULL
# for(i in 1:n){
#   mean.y[i] <- newalpha[i] * u / (newalpha[i]-1)
# }
# sqrt(mean((mean.y-y)^2))

r <- matrix(, nrow = n, ncol = 20)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 20
for(i in 1:n){
  for(t in 1:T){0
    r[i, t] <- qnorm(pnorm(mod.summary$residuals[i]))
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
  theme(text = element_text(size = 20)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3))
ggsave("./Laboratory/Application/figures/loglr_qqplot.pdf", width=10)

library(gap)
library(dgof)
pdf("./Laboratory/Application/figures/map_unifqqplot.pdf", width = 10, height= 7.78)
qqunif(pPareto(y, u, alpha=newalpha), logscale = FALSE, col = "black")
dev.off()
pdf("./Laboratory/Application/figures/loglr_unifqqplot.pdf", width = 10, height= 7.78)
qqunif(pnorm(mod.summary$residuals), logscale = FALSE, col = "black")
dev.off()

# qPareto(pPareto(), u, alpha=newalpha)
ks.test(pPareto(y, u, alpha=newalpha), "punif")
ks.test(pnorm(mod.summary$residuals), "punif")
ks.test(pnorm(mod$fitted.values), "punif")
