library(npreg)
library(tmvnsim)
library(reshape2)
library(splines2)
library(scales)
library(MASS)
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
library(ggmcmc)
library(MCMCvis)
# library(ggplotify)

#Scenario 1
set.seed(3)
# n <- 1000
# # beta <- c(0.2, 0.7)
# beta <- c(0.2, 0, 0.8, 0, 0, -0.1, 0, 0, 0, -0.4)
# p <- length(beta)
# x.scale <- cbind(replicate(p, runif(n, 0, 1)))
# x.scale <- scale(x.scale)


# alp.origin <- y <- NULL
# for(i in 1:n){
#   alp.origin[i] <- exp(sum(x.scale[i, ] * beta))
#   y[i] <- rPareto(1, 1, alpha = alp.origin[i])
# }
# # plot(y)


# x.origin <- x.scale
# y.origin <- y
# x.scale <- x.scale[which(y>quantile(y, threshold)),]
# u <- quantile(y, threshold)

# y <- as.matrix(y[y > quantile(y, threshold)])
# # n <- length(y)
# n <- length(y.origin)

# xholder <- new.bs.x <- bs.x <- matrix(,nrow=n, ncol=0)
# psi <- 33
# # new.bs.x <- bs.x <- array(NA, c(n, psi, p-1))
# for(i in 1:p){
#   splines <- bbase(x.origin[, i], min(x.origin[, i]), max(x.origin[, i]), nseg = (psi-3), bdeg = 3)
# #   splines <- bbase(x.scale[, i], min(x.scale[, i]), max(x.scale[, i]), nseg = (psi-3), bdeg = 3)
#   bs.x <- cbind(bs.x, splines)
# #   bs.x[, , i-1] <- splines
#   # phi <- dim(out[[1]][[1]]$X)[2]
#   # psi <- dim(splines)[2]
# }

# newx <- seq(0, 1, length.out=n)
# xholder <- cbind(xholder, rep(1,n))
# for(i in 1:p){
#   splines <- bbase(newx, min(newx), max(newx), nseg = (psi  -3), bdeg = 3)
#   new.bs.x <- cbind(new.bs.x, splines)
#   xholder <- cbind(xholder, newx)
# }

n <- 5000
xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
psi <- 20
p <- 10
threshold <- 0.90

x.origin <- cbind(replicate(p, runif(n, 0, 1)))
x.origin <- scale(x.origin)
newx <- seq(0, 1, length.out=n)
no.theta <- 2

for(i in 1:p){
  test.knot <- seq(0, 1, length.out = psi)
  splines <- basis.tps(newx, test.knot, m=2, rk=FALSE, intercept = TRUE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
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

f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.origin[, j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j] + (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
    f.linear.origin[,j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j]
    f.nonlinear.origin[,j] <- (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
  alp.origin[i] <- exp(sum(f.origin[i,]))
  y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i])
}

n <- 5000
psi <- 20
threshold <- 0.90
p <- 5
no.theta <- 1
simul.no <- 50

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)

sample_meanvector <- runif(p, 0, 1)
sample_covariance_matrix <- matrix(NA, nrow = p, ncol = p)
diag(sample_covariance_matrix) <- 1

mat_Sim <- matrix(data = NA, nrow = p, ncol = p)
U <- runif(n = p) * 0.5

for(i in 1:p)
{
  if(i %in% c(2,3))
  {
    U_Star <- pmin(U + 0.2 * runif(n = p), 0.99999)
    
  }else
  {
    U_Star <- pmin(pmax(U + sample(c(0, 1), size = p, replace = TRUE) * runif(n = p), 0.00001), 0.99999)
  }
  
  mat_Sim[, i] <- qnorm(U_Star)  
}

cor_Mat <- cor(mat_Sim)
sample_covariance_matrix <- cor_Mat * (p/2)
# diag(sample_covariance_matrix) <- 1
## create multivariate normal distribution
# x.origin <- mvrnorm(n = n, mu = rep(0,p), Sigma = sample_covariance_matrix)

C <- matrix(c(1, 0.3, 0.5, 0.3, 0.3,
              0.3, 1, 0.95, 0.4, 0.4,
              0.5, 0.95, 1, 0.5, 0.1,
              0.3, 0.4, 0.5 , 1, 0.5,
              0.3, 0.4, 0.5, 0.5, 1), nrow = p)
x.origin <- tmvnsim(n = n, k = p, lower = rep(0, p), means = rep(0, p), sigma = C)$samp


corrplot.mixed(cor(x.origin),
                upper = "circle",
                lower = "number",
                addgrid.col = "black")
newx <- seq(0, 1, length.out=n)

for(i in 1:p){
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

gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
    for (ps in 1:psi){
        if(j %in% c(1,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
        else if(j==7){
            if(ps <= (psi/2)){gamma.origin[ps, j] <- 0.01}
            else{gamma.origin[ps, j] <- 0.01}
        }
        else {
            if(ps <= (psi/2)){gamma.origin[ps, j] <- 0.01}
            else{gamma.origin[ps, j] <- 0.01}
        }
    }
}

theta.origin <- c(0, 0.2, 0.2, 0, 0)


f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.origin[, j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[j] + (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
    f.linear.origin[,j] <- as.matrix(bs.linear[1:n, (((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[j]
    f.nonlinear.origin[,j] <- (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
  alp.origin[i] <- exp(sum(f.origin[i,]))
  y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i])
}

# lambda <- 0.995#sqrt(2*log(n))
# theta <- 0
# lambda.1 <- lambda * theta
# lambda.2 <-. lambda * (1-theta)
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
            second.term[i] <- exp.prime(sum(f[i,]), thres = 10) * log(y.origin[i])
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
# beta.emp <- c(rep(0, no.theta*p), rep(0, p*psi))
beta.emp <- c(as.vector(theta.origin), as.vector(gamma.origin))
beta.map <- optim(beta.emp, fn = log.posterior, #gr = grad.log.posterior, 
                  y.origin = y.origin,
                  # method = "BFGS",
                  method = "CG",
                  # method = "SANN",
                  control = list(fnscale = -1))
# theta.map <- matrix(beta.map$par[1:(2*p)],nrow=2)
theta.map <- beta.map$par[1:(no.theta*p)]
gamma.map <- beta.map$par[(no.theta*p)+1:(psi*p)]

df.theta <- data.frame("seq" = seq(1, (no.theta*p)),
                  theta.map,
                  "theta.true" = as.vector(theta.origin))
df.theta$covariate <- factor(rep(seq(1, 1 + nrow(df.theta) %/% no.theta), each = no.theta, length.out = nrow(df.theta)))
ggplot(df.theta, aes(x = seq)) + 
  geom_point(aes(y = theta.map, color = covariate), size = 1.5) + 
#   geom_smooth(method="gam") +
  geom_point(data = df.theta, aes(y = theta.true, color = "true"), size = 2.5) +
  labs(title=expression("MAP vs True for"~theta)) + 
  # ggtitle(expression(atop(paste("MAP vs True for ", bold(gamma))))) +
#   annotate("text", x = seq(0, 330, length.out=10), y = -1, label = beta, colour = "red", size = 10) +   
#   scale_color_manual(values = c("true"="black"))+
  theme(plot.title = element_text(hjust = 0.5, size = 20))
    
df <- data.frame("seq" = seq(1, (psi*p)), 
                  gamma.map, 
                  "gamma.true" = as.vector(gamma.origin))
df$covariate <- factor(rep(seq(1, 1 + nrow(df) %/% psi), each = psi, length.out = nrow(df)))


ggplot(df, aes(x = seq, y = gamma.map, color = covariate)) + 
  geom_point() + 
  # geom_smooth(method="gam") +
  geom_point(aes(x=seq, y = gamma.true, color = "true")) + 
  # geom_line(aes(x = seq, y = true, color = "true"), linetype = 2) +
  labs(title=expression("MAP vs True for"~gamma)) + 
  # ggtitle(expression(atop(paste("MAP vs True for ", bold(gamma))))) +
#   annotate("text", x = seq(0, 330, length.out=10), y = -1, label = beta, colour = "red", size = 10) +   
  theme(plot.title = element_text(hjust = 0.5, size = 20))

f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
newalpha <- NULL
for (j in 1:p){
  f.linear.new[,j] <- as.matrix(bs.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% matrix(theta.map, nrow=no.theta)[,j]
  f.nonlinear.new[,j] <- bs.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.new[1:n, j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}
new.y <- NULL

# set.seed(100)
for(i in 1:n){  
  temp <- sum(f.new[i,])
  newalpha[i] <- exp(temp)
  new.y[i] <- rPareto(1, 1, alpha = newalpha[i])
}

func.linear.new <- func.nonlinear.new <- func.linear.origin <- func.nonlinear.origin <- func.new <- func.origin <- matrix(, nrow=n, ncol=0)
for (j in 1:p){
  func.origin <- cbind(func.origin, sort(f.origin[,j]))
  func.linear.origin <- cbind(func.linear.origin, sort(f.linear.origin[,j]))
  func.nonlinear.origin <- cbind(func.nonlinear.origin, sort(f.nonlinear.origin[,j]))
  func.new <- cbind(func.new, sort(f.new[,j]))
  func.linear.new <- cbind(func.linear.new, sort(f.linear.new[,j]))
  func.nonlinear.new <- cbind(func.nonlinear.new, sort(f.nonlinear.new[,j]))  
}
# for (j in 1:p){
#   func.new <- cbind(func.new, sort(f.new[,j]))
# }
covariates <- gl(p, n, (p*n))
replicate <- gl(2, n, (p*n))
func.df <- data.frame(seq = seq(0,1,length.out = n), 
                        origin=as.vector(func.origin),
                        origin.linear=as.vector(func.linear.origin),
                        origin.nonlinear=as.vector(func.nonlinear.origin), 
                        new=as.vector(func.new),
                        new.linear=as.vector(func.linear.new),
                        new.nonlinear=as.vector(func.nonlinear.new),
                        covariates=covariates, 
                        replicate=replicate)

ggplot(func.df, aes(x=seq, group=interaction(covariates, replicate))) + 
  geom_line(aes(y=origin, colour = covariates, linetype = "true")) + 
  geom_line(aes(y=new, colour = covariates, linetype = "MAP")) + ylab ("") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  facet_grid(covariates ~ .) + ggtitle("MAP vs True for Smooth Functions") + 
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

ggplot(func.df, aes(x=seq, group=interaction(covariates, replicate))) + 
  geom_line(aes(y=origin.linear, colour = covariates, linetype = "true")) + 
  geom_line(aes(y=new.linear, colour = covariates, linetype = "MAP")) + ylab ("") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  facet_grid(covariates ~ .) + ggtitle("MAP vs True for Linear Functions") + 
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

ggplot(func.df, aes(x=seq, group=interaction(covariates, replicate))) + 
  geom_line(aes(y=origin.nonlinear, colour = covariates, linetype = "true")) + 
  geom_line(aes(y=new.nonlinear, colour = covariates, linetype = "MAP")) + ylab ("") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  facet_grid(covariates ~ .) + ggtitle("MAP vs True for Nonlinear Functions") + 
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20))


data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "trueAlp" = sort(alp.origin),
                            "mapAlp" = sort(newalpha))

plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x)))
print(plt.samp + 
      geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
      geom_line(aes(y = mapAlp, col = paste0("MAP Alpha:",lambda.1,"/",lambda.2,"/",lambda.3)), linewidth = 2.5, linetype = 2) +
      labs(col = "") +
      scale_color_manual(values = c("#e0b430", "red"))+
      theme(text = element_text(size = 15)) + 
      theme(legend.position="top", legend.key.size = unit(1, 'cm')))

# plot(sort(alp.origin))
# plot(sort(new.y), sort(y.origin))
# abline(a=0, b=1, col = "red", lty = 2)
rbind(matrix(theta.map, nrow = no.theta, ncol = p), matrix(gamma.map, nrow = psi, ncol = p))

#### Test Set
f.nonlinear.origin <- f.linear.origin <- f.origin <- f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
true.alpha <- new.alpha <- NULL
for (j in 1:p){
  # f.new[1:n, j] <- as.matrix(xholder.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% matrix(theta.map,nrow=no.theta)[,j] + xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.linear.new[,j] <- as.matrix(xholder.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% matrix(theta.map,nrow=no.theta)[,j]
  f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
  # f.origin[1:n, j] <- as.matrix(xholder.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j] + xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
  f.linear.origin[,j] <- as.matrix(xholder.linear[,(((j-1)*no.theta)+1):(((j-1)*no.theta)+no.theta)]) %*% theta.origin[,j]
  f.nonlinear.origin[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
  f.origin[,j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
}
# set.seed(100)
for(i in 1:n){
  true.alpha[i] <- exp(sum(f.origin[i,]))
  new.alpha[i] <- exp(sum(f.new[i,]))
}

func.linear.new <- func.nonlinear.new <- func.linear.origin <- func.nonlinear.origin <- func.new <- func.origin <- matrix(, nrow=n, ncol=0)
for (j in 1:p){
  func.origin <- cbind(func.origin, sort(f.origin[,j]))
  func.linear.origin <- cbind(func.linear.origin, sort(f.linear.origin))
  func.nonlinear.origin <- cbind(func.nonlinear.origin, sort(f.nonlinear.origin))
  func.new <- cbind(func.new, sort(f.new[,j]))
  func.linear.new <- cbind(func.linear.new, sort(f.linear.new))
  func.nonlinear.new <- cbind(func.nonlinear.new, sort(f.nonlinear.new))
}

covariates <- gl(p, n, (p*n))
replicate <- gl(2, n, (p*n))
func.df <- data.frame(seq = seq(0,1,length.out = n), 
                        origin=as.vector(func.origin), 
                        new=as.vector(func.new), 
                        covariates=covariates, 
                        replicate=replicate)

ggplot(func.df, aes(x=seq, group=interaction(covariates, replicate))) + 
  geom_line(aes(y=origin, colour = covariates, linetype = "true")) + 
  geom_line(aes(y=new, colour = covariates, linetype = "MAP")) + ylab ("") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  facet_grid(covariates ~ .) + ggtitle("MAP vs True for Smooth Functions") + 
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "trueAlp" = sort(true.alpha),
                            "mapAlp" = sort(new.alpha))

plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x)))
print(plt.samp + 
      geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
      geom_line(aes(y = mapAlp, col = paste0("MAP Alpha:",lambda.1,"/",lambda.2,"/",lambda.3)), linewidth = 2.5, linetype = 2) +
      labs(col = "") +
      scale_color_manual(values = c("#e0b430", "red"))+
      theme(text = element_text(size = 15)) + 
      theme(legend.position="top", legend.key.size = unit(1, 'cm')))

# Randomized quantile residuals
r <- matrix(NA, nrow = n, ncol = 20)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 20
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y.origin[i], 1, alpha = newalpha[i]))
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
  # theme(text = element_text(size = 30)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3)) 
# ggsave("../../Figures/residuals_application.pdf")
# ggsave("../Results/residuals_application_IWSM.pdf", width=6, height=6)

#--------------------------------------------------------------------------------------
# plotting bsplines basis function of covariates


plts <- list()
for(j in 1:p){
  x.x <- sort(scale(x.origin[,j]))
  # # Make a matrix containing the thin plate basis
  ndx <- 20
  basis.x <- data.frame(x=x.x, y=y.origin, id = as.factor(ndx))
  # deg <- 3
  # B <- bbase(x.x, min(x.x), max(x.x), nseg = ndx, bdeg = deg)
  knots <- seq(min(x.x), max(x.x), length.out = (ndx))
  B <- basis.tps(x.x, knots, m = 2, rk = FALSE)
  nb1 <- ncol(B)

  # A basis for plotting the fit on the grid xg
  ng <- length(x.x)
  xg <- seq(min(x.x), max(x.x), length.out = ng)
  Bg <- basis.tps(xg, seq(min(x.x), max(x.x), length.out = ndx), m = 2, rk = FALSE)
  # Bg <- bbase(xg, min(x.x), max(x.x), nseg = ndx, bdeg = deg)

  a <- solve(t(B) %*% B, t(B) %*% y.origin, tol = 1e-20)
  z <- Bg %*% a

  # Make a matrix with B-splines scaled by coefficients
  Bsc1 <- Bg %*% diag(c(a))
  # Create data frames for ggplot
  Zf1 <- data.frame(x = xg, y = z, id = as.factor(1))
  Bf1 <- data.frame(x = rep(xg, nb1), y = as.vector(Bsc1), id = as.factor(rep(1:nb1, each = ng)))
  # Bf1$y[abs(Bf1$y) < 0.0001] = NA
  # Bf1 <- na.omit(Bf1)

  # Build the graphs
  plt <- ggplot(Bf1 , aes(x = x, y = y, group = id, colour = id)) +
    ggtitle(paste0("Basis Function of X_" ,j)) +
    # geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line(data = Zf1, linewidth = 1, colour = "grey") +
    # geom_point(data = basis.x, colour = "grey60", size = 0.8, shape=1) +
    geom_line(linewidth = 0.7) +
    # geom_point(data = basis.x, colour = "red", size = 2, shape = 1) +
    xlab("") + ylab ("") +
    JOPS_theme() +
    theme(legend.position = "none") +
    scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

  plts[[j]] <- plt
  # ggsave(paste0("./Results/vert_bsplines_",j,".png"))
}
# plts
grid.arrange(grobs = plts, ncol = 2, nrow = 5)

# Bspline basis
plts <- list()
for(j in 1:p){
  x.x <- sort(scale(x.origin[,j]))
  # # Make a matrix containing the B-spline basis
  ndx <- 13
  basis.x <- data.frame(x=x.x, y=y.origin, id = as.factor(ndx))
  deg <- 3
  B <- bbase(x.x, min(x.x), max(x.x), nseg = ndx, bdeg = deg)
  nb1 <- ncol(B)

  # A basis for plotting the fit on the grid xg
  ng <- length(x.x)
  xg <- seq(min(x.x), max(x.x), length.out = ng)
  # Bg <- basis.tps(xg, seq(min(x.x), max(x.x), length.out = ndx), m = 2, rk = FALSE)
  Bg <- bbase(xg, min(x.x), max(x.x), nseg = ndx, bdeg = deg)

  a <- solve(t(B) %*% B, t(B) %*% y.origin, tol = 1e-20)
  z <- Bg %*% a

  # Make a matrix with B-splines scaled by coefficients
  Bsc1 <- Bg %*% diag(c(a))
  # Create data frames for ggplot
  Zf1 <- data.frame(x = xg, y = z, id = as.factor(1))
  Bf1 <- data.frame(x = rep(xg, nb1), y = as.vector(Bsc1), id = as.factor(rep(1:nb1, each = ng)))
  # Bf1$y[abs(Bf1$y) < 0.0001] = NA
  # Bf1 <- na.omit(Bf1)

  # Build the graphs
  plt <- ggplot(Bf1 , aes(x = x, y = y, group = id, colour = id)) +
    ggtitle(paste("Basis Function of ", colnames(x.origin)[j])) +
    # geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line(data = Zf1, linewidth = 1, colour = "grey") +
    # geom_point(data = basis.x, colour = "grey60", size = 0.8, shape=1) +
    geom_line(linewidth = 0.7) +
    # geom_point(data = basis.x, colour = "red", size = 2, shape = 1) +
    xlab("") + ylab ("") +
    JOPS_theme() +
    theme(legend.position = "none") +
    scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

  plts[[j]] <- plt
  # ggsave(paste0("./Results/vert_bsplines_",j,".png"))
}
# plts
# grid.arrange(grobs = plts, ncol = 2, nrow = 5)

# log.posterior <- function(gamma){
#   log.lik <- function(gamma){
#     f <- matrix(, nrow=n, ncol=p)
#     term <- first.term <- second.term <- NULL
#     for(j in 1:p){
#       coef <- as.matrix(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])
#       f[, j] <- bs.x[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% coef
#     }
    
#     for(i in 1:n){
#       first.term[i] <- sum(f[i,]) - log(y.origin[i])
#       second.term[i] <- exp(sum(f[i,])) * log(y.origin[i])
#       term[i] <- first.term[i] - second.term[i]
#     }
#     return(sum(term))
#   }

#   log.prior <- function(gamma){
#     w <- matrix(0, nrow=(psi-1), ncol = p)
#     term <- first.term <- second.term <- NULL
#     for(j in 1:p){
#       first.term[j] <- -1 * lambda.1 * sum(abs(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)]))
#       for(ps in 1:(psi-1)){
#         w[ps, j] <- abs(gamma[((j-1)*psi)+ps+1] - gamma[((j-1)*psi)+ps])
#       }
#       second.term[j] <- -1 * lambda.2 * sum(w[,j])
#       term[j] <- first.term[j] + second.term[j]
#       # term[j] <- first.term[j] + second.term[j] -(0.9*log(lambda.1^2)-0.1*lambda.1^2-0.9*log(lambda.2^2) -0.1*lambda.2^2)
#     }
#     return(sum(term))
#   }
#   return(log.lik(gamma) + log.prior(gamma))
# }

# gamma.map <- optim(rep(0, (psi*p)), log.posterior, 
#         method = "BFGS",
#         control = list(fnscale = -1))$par

# df <- data.frame("seq"=seq(1,(psi*p)), gamma.map)
# df$jterm <- factor(rep(seq(1, 1 + nrow(df) %/% psi), each = psi, length.out = nrow(df)))

# ggplot(df, aes(x=seq, y = gamma.map, color = jterm)) + 
#   geom_point() + 
#   geom_smooth(method=glm) +
#   ggtitle("Bayesian Fused Lasso") + 
#   annotate("text", x = seq(0, 330, length.out=10), y = -1, label = beta, colour = "red", size = 10) +   
#   theme(plot.title = element_text(hjust = 0.5, size = 20))

# me.log.posterior <- function(gamma){
#   moreau.envelope <- function(w){
#     if(w < -1)  {
#       ans <- -0.5 - w
#     }
#     else if (1 < w){
#       ans <- w - 0.5
#     }
#     else {
#       ans <- w^2 / 2
#     }
#     return(ans)
#   }

#   log.lik <- function(gamma){
#     f <- matrix(, nrow=n, ncol=p)
#     term <- first.term <- second.term <- NULL
#     for(j in 1:p){
#       coef <- as.matrix(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])
#       f[, j] <- bs.x[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% coef
#     }

#     for(i in 1:n){
#       first.term[i] <- sum(f[i,]) - log(y.origin[i])
#       second.term[i] <- exp(sum(f[i,])) * log(y.origin[i])
#       term[i] <- first.term[i] - second.term[i]
#     }
#     return(sum(term))
#   }

#   log.prior <- function(gamma){
#     w <- matrix(0, nrow = (psi-1), ncol = p)
#     term <- first.term <- second.term <- NULL
#     for(j in 1:p){
#       first.term[j] <- -1 * lambda.1 * 
#                         sum(unlist(lapply(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)], moreau.envelope)))
#       for(ps in 1:(psi-1)){
#         w[ps, j] <- moreau.envelope(gamma[((j-1)*psi)+ps+1] - gamma[((j-1)*psi)+ps])
#       }
#       second.term[j] <- -1 * lambda.2 * sum(w[,j])
#       term[j] <- first.term[j] + second.term[j]
#     }
#     return(sum(term))
#   }
#   return(log.lik(gamma) + log.prior(gamma))
# }


# me.gamma.map <- optim(rep(0, psi*p), me.log.posterior, 
#         method = "BFGS",
#         control = list(fnscale = -1))$par
# me.df <- data.frame("seq"=seq(1,(psi*p)), me.gamma.map)
# me.df$jterm <- factor(rep(seq(1, 1 + nrow(df) %/% psi), each = psi, length.out = nrow(df)))

# ggplot(me.df, aes(x=seq, y = me.gamma.map, color = jterm)) + 
#   geom_point() + 
#   geom_smooth(method=glm) +
#   annotate("text", x = seq(0, 330, length.out=10), y = -2.5, label = beta, colour = "red", size = 10) + 
#   ggtitle("Bayesian Fused Lasso via Approximation with Moreau Envelope") + 
#   theme(plot.title = element_text(hjust = 0.5, size = 20))

# new.f.moreau <- new.f <- matrix(, nrow = n, ncol=p)
# newalpha.moreau <- newalpha <- NULL
# for (j in 1:p){
#   new.f[1:n, j] <- new.bs.x[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, nrow = psi, ncol = p)[1:psi, j]
#   new.f.moreau[1:n, j] <- new.bs.x[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(me.gamma.map, nrow = psi, ncol = p)[1:psi, j]
# }
# for(i in 1:n){
#   temp <- sum(new.f[i, 1:p])
#   newalpha[i] <- exp(temp)
#   temp <- sum(new.f.moreau[i, 1:p])
#   newalpha.moreau[i] <- exp(temp)
# }

# trueAlp <- NULL
# for (i in 1:n){
#   trueAlp[i] <- exp(sum(xholder[i, 1:p] %*% beta[1:p]))
# }

# data.scenario <- data.frame("x" = c(1:n),
#                             "constant" = newx,
#                             "trueAlp" = sort(trueAlp),
#                             "mapAlp" = sort(newalpha),
#                             "mapAlp.moreau" = sort(newalpha.moreau))

# plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x)))
# print(plt.samp + ylim(0, 2.5) +
#       geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
#       geom_line(aes(y = mapAlp, col = "Optimised MAP Alpha"), linewidth = 2.5) +
#       geom_line(aes(y = mapAlp.moreau, col = "Optimised MAP Alpha (Moreau)"), linewidth = 2.5, linetype = 2) +
#       labs(col = "") +
#       scale_color_manual(values = c("#e0b430", "red", "darkgreen"))+
#       theme(text = element_text(size = 15)) + 
#       theme(legend.position="top", legend.key.size = unit(1, 'cm')))