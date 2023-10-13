library(JOPS)
library(Pareto)
library(npreg)
library(gridExtra)
library(colorspace)
library(reshape2)
library(splines2)
library(scales)

suppressMessages(library(tidyverse))
# library(ggplotify)

#Scenario 1
set.seed(2)

n <- 5000
psi <- 20
threshold <- 0.90
p <- 10
no.theta <- 1
simul.no <- 50

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
    f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- (bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j])
    f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp(theta.origin[1] + sum(f.origin[i,]))
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
    bs.linear <- cbind(bs.linear, tps[,1:no.theta])
    bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
}

f.nonlinear.new <- f.linear.new <- f.new <- f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
    f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
    f.nonlinear.origin[,j] <- bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
    f.linear.new[,j] <- xholder.linear[, j] * theta.origin[j+1]
    f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
    f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}

true.alpha <- alp.new <- alp.origin <- NULL
for(i in 1:n){
    alp.origin[i] <- exp(theta.origin[1] + sum(f.origin[i,]))
    alp.new[i] <- exp(theta.origin[1] + sum(f.new[i,]))
}

# theta <- 0
# lambda.1 <- lambda * theta
# lambda.2 <-. lambda * (1-theta)


log.posterior <- function(beta, y.origin){
  theta <- beta[1:(p+1)]
  gamma <- matrix(beta[(p+2):(p+1+(psi*p))], ncol=p)
  g <- matrix(, nrow=n, ncol=p)
  lik <- first.lik <- second.lik <- NULL
  for(j in 1:p){
      # coef <- as.matrix(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])
      linear.term <- bs.linear[,j] * theta[j+1]
      nonlinear.term <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[,j]
      g[, j] <- linear.term + nonlinear.term
  }
  for(i in 1:n){
      first.lik[i] <- theta[1] + sum(g[i,]) - log(y.origin[i])
      second.lik[i] <- exp(theta[1] + sum(g[i,])) * log(y.origin[i]/u)
      lik[i] <- first.lik[i] - second.lik[i]
  }
  sum.lik <- sum(lik)

  # lambda.1 <- beta[length(beta)-1]
  # lambda.2 <- beta[length(beta)]
  lambda.1 <- 0.00001
  lambda.2 <- 0.00001
  prior <- first.prior <- second.prior <- NULL
  for(j in 1:p){
      # print(sum(abs(theta[j+1])))
      first.prior[j] <- -1 * lambda.1 * abs(theta[j+1])
      second.prior[j] <- -1 * lambda.2 * sqrt(sum((gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])^2))
      prior[j] <- first.prior[j] + second.prior[j]
  }
  sum.prior <- sum(prior) - (lambda.1 * abs(theta[1])) #+
                # ((p+1) * log(lambda.1)) + (p * psi * log(lambda.2)) +
                # ((0.1-1)*log(lambda.1) - (5 * lambda.1)) + 
                # ((0.01-1)*log(lambda.2) - (0.01 * lambda.2))
  # print(first.prior)
  return(sum.lik + sum.prior)
}

# log.posterior <- function(beta, y.origin){
#     log.lik <- function(beta){
#         exp.prime <- function(x, thres){
#             if(x > thres){ans <- exp(thres) + exp(thres)*(x-thres)}
#             else{ans <- exp(x)}
#             return(ans)
#         }
#         theta.0 <- beta[1]
#         theta <- beta[2:(p+1)]
#         gamma <- matrix(beta[-(1:(p+1))], ncol=p)
#         g <- matrix(, nrow=n, ncol=p)
#         term <- first.term <- second.term <- NULL
#         for(j in 1:p){
#             # coef <- as.matrix(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])
#             linear.term <- bs.linear[,j] * theta[j]
#             nonlinear.term <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma[,j]
#             g[, j] <- rep(theta.0, n) + linear.term + nonlinear.term
#         }
#         for(i in 1:n){
#             first.term[i] <- sum(g[i,]) - log(y.origin[i])
#             second.term[i] <- exp.prime(sum(g[i,]), thres = 10) * log(y.origin[i]/u)
#             term[i] <- first.term[i] - second.term[i]
#         }
#         return(sum(term))
#     }
#     log.prior <- function(beta){
#         moreau.envelope <- function(w){
#             if(w < -1){ans <- -0.5 - w}
#             else if (1 < w){ans <- w - 0.5}
#             else {ans <- (w^2) / 2}
#             return(ans)
#         }
#         # theta.0 <- beta[1]
#         theta <- beta[1:(p+1)]
#         gamma <- matrix(beta[-(1:(p+1))], ncol=p)
#         g.1 <- g.2 <- term <- third.term <- first.term <- second.term <- NULL
#         for(j in 1:p){
#             first.term[j] <- -1 * lambda.1 * sqrt(sum((gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])^2))
#             # first.term[j] <- -1 * lambda.1 * abs(sum(gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)]))
#             second.term[j] <- -1 * lambda.2 * sum((gamma[(((j-1)*psi)+1):(((j-1)*psi)+psi)])^2)
#             third.term[j] <- -1 * lambda.3 * sum(abs(theta[j]))
#             term[j] <- first.term[j] + second.term[j] + third.term[j]
#             # term[j] <- first.term[j] + second.term[j] - lambda.3 * sum(abs(beta))
#         }
#         return(sum(term))
#     }
#     return(log.lik(beta) + log.prior(beta))
# }
# beta.emp <- c(rep(0, no.theta*p), rep(0, p*psi))
beta.emp <- c(as.vector(theta.origin), as.vector(gamma.origin))
beta.map <- optim(beta.emp, fn = log.posterior, #gr = grad.log.posterior, 
                  y.origin = y.origin,
                  # method = "BFGS",
                  method = "CG",
                  # method = "SANN",
                  control = list(trace=2, fnscale = -1, maxit = 1000))
# theta.map <- matrix(beta.map$par[1:(2*p)],nrow=2)
theta.map <- beta.map$par[1:(p+1)]
gamma.map <- beta.map$par[(p+1+1):(p+1+(psi*p))]
lambda.map <- beta.map$par[-c(1:(p+1+(psi*p)))]

df.theta <- data.frame("seq" = seq(1, p+1),
                  theta.map,
                  "theta.true" = theta.origin)
# df.theta$covariate <- factor(rep(seq(1, 1 + nrow(df.theta) %/% no.theta), each = no.theta, length.out = nrow(df.theta)))
df.theta$covariate <- factor(0:p)
df.theta$labels <- factor(0:p)
ggplot(df.theta, aes(x = labels)) + ylab("") + xlab("") +
  geom_point(aes(y = theta.map, col = covariate), size = 6) + 
  geom_point(aes(y = theta.true), color="red", size = 4) +
  # labs(title=expression("MAP vs True for"~theta)) + 
  scale_x_discrete(labels = c(expression(bold(theta[0])),
                              expression(bold(theta[1])),
                              expression(bold(theta[2])),
                              expression(bold(theta[3])),
                              expression(bold(theta[4])),
                              expression(bold(theta[5])),
                              expression(bold(theta[6])),
                              expression(bold(theta[7])),
                              expression(bold(theta[8])),
                              expression(bold(theta[9])),
                              expression(bold(theta[10]))))+
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size=30),
          # axis.ticks.x = element_blank(),
          axis.text.x = element_text(hjust=0.3),
          axis.text = element_text(size = 30),
          panel.grid.minor.x = element_blank())
    
df <- data.frame("seq" = seq(1, (psi*p)), 
                  gamma.map, 
                  "gamma.true" = as.vector(gamma.origin))
df$covariate <- factor(rep(seq(1, 1 + nrow(df) %/% psi), each = psi, length.out = nrow(df)))
df$labels <- factor(1:(psi*p))


ggplot(df, aes(x =labels , y = gamma.map, col = covariate)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("")+
  geom_point(aes(y = gamma.true), color = "red", size=3) + 
  geom_point(size = 4) + ylab("") +
  # geom_smooth(method="gam") +
  # geom_line(aes(x = seq, y = true, color = "true"), linetype = 2) +
  # labs(title=expression("MAP vs True for"~gamma)) + 
  # ggtitle(expression(atop(paste("MAP vs True for ", bold(gamma))))) +
#   annotate("text", x = seq(0, 330, length.out=10), y = -1, label = beta, colour = "red", size = 10) +
  scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+15), 
                    label = c(expression(bold(gamma[1])), 
                              expression(bold(gamma[2])), 
                              expression(bold(gamma[3])), 
                              expression(bold(gamma[4])), 
                              expression(bold(gamma[5])), 
                              expression(bold(gamma[6])),
                              expression(bold(gamma[7])),
                              expression(bold(gamma[8])),
                              expression(bold(gamma[9])),
                              expression(bold(gamma[10]))), 
                    expand=c(0,5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size=30),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(hjust=1),
          axis.text = element_text(size = 30),
          panel.grid.major.x = element_blank())

f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
newalpha <- NULL
for (j in 1:p){
  f.linear.new[,j] <- bs.linear[,j] * theta.map[j+1]
  f.nonlinear.new[,j] <- bs.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.new[1:n, j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}
new.y <- NULL

# set.seed(100)
for(i in 1:n){
  newalpha[i] <- exp(theta.map[1] + sum(f.new[i,]))
  # new.y[i] <- rPareto(1, 1, alpha = newalpha[i])
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
                        x = as.vector(apply(x.origin, 2, sort, method = "quick")),
                        origin=as.vector(func.origin),
                        origin.linear=as.vector(func.linear.origin),
                        origin.nonlinear=as.vector(func.nonlinear.origin), 
                        new=as.vector(func.new),
                        new.linear=as.vector(func.linear.new),
                        new.nonlinear=as.vector(func.nonlinear.new),
                        covariates=covariates, 
                        replicate=replicate)

equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_line(aes(y=origin, colour = covariates, linetype = "true"), linewidth=2) + 
  geom_line(aes(y=new, colour = covariates, linetype = "MAP"), linewidth=2) + 
  ylab("") + xlab ("Smooth Functions") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  facet_grid(covariates ~ .) + ggtitle("MAP vs True for Smooth Functions") +
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_line(aes(y=origin.linear, colour = covariates, linetype = "true"), linewidth=2) + 
  geom_line(aes(y=new.linear, colour = covariates, linetype = "MAP"), linewidth=2) + 
  ylab ("") + facet_grid(covariates ~ .) + xlab("Linear Components") + 
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,-20,0,-30),
        strip.text = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_line(aes(y=origin.nonlinear, colour = covariates, linetype = "true"), linewidth=2) + 
  geom_line(aes(y=new.nonlinear, colour = covariates, linetype = "MAP"), linewidth=2) +
  ylab ("") + facet_grid(covariates ~ .) + xlab("Nonlinear Components") +
  # geom_point(aes(y=origin, shape = replicate)) + geom_point(aes(y=new, shape = replicate)) +
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=45),
        legend.margin=margin(0,0,0,-10),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))


data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "trueAlp" = sort(alp.origin),
                            "mapAlp" = sort(newalpha))

plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x))) + xlab("")
print(plt.samp + 
      geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
      geom_line(aes(y = mapAlp, col = "MAP Alpha"), linewidth = 2.5, linetype = 2) +
      labs(col = "") +
        theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        scale_color_manual(values = c("#e0b430", "red"))+
        theme(text = element_text(size = 15),
                legend.position="bottom", legend.key.size = unit(1, 'cm'),
                axis.text = element_text(size = 20),
                legend.margin=margin(-15,-15,-15,-15),
                legend.box.margin=margin(-25,0,20,0)))

# ggsave(paste0("../Laboratory/Simulation/BayesianFusedLasso/results/map_alpha_n1_train.pdf"), width=10)

# plot(sort(alp.origin))
# plot(sort(new.y), sort(y.origin))
# abline(a=0, b=1, col = "red", lty = 2)
# rbind(matrix(theta.map, nrow = no.theta, ncol = p), matrix(gamma.map, nrow = psi, ncol = p))

#### Test Set
f.nonlinear.origin <- f.linear.origin <- f.origin <- f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
true.alpha <- new.alpha <- NULL
for (j in 1:p){
  f.linear.new[,j] <- xholder.linear[,j] * theta.map[j+1]
  f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.new[,j] <- f.linear.new[,j] + f.nonlinear.new[,j]
  f.linear.origin[,j] <- xholder.linear[,j] * theta.origin[j+1]
  f.nonlinear.origin[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
  f.origin[,j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
}
# set.seed(100)
for(i in 1:n){
  true.alpha[i] <- exp(theta.map[1] + sum(f.origin[i,]))
  new.alpha[i] <- exp(theta.origin[1] + sum(f.new[i,]))
}

func.linear.new <- func.nonlinear.new <- func.linear.origin <- func.nonlinear.origin <- func.new <- func.origin <- matrix(, nrow=n, ncol=0)
for (j in 1:p){
  func.origin <- cbind(func.origin, f.origin[,j])
  func.linear.origin <- cbind(func.linear.origin, f.linear.origin[,j])
  func.nonlinear.origin <- cbind(func.nonlinear.origin, f.nonlinear.origin[,j])
  func.new <- cbind(func.new, f.new[,j])
  func.linear.new <- cbind(func.linear.new, f.linear.new[,j])
  func.nonlinear.new <- cbind(func.nonlinear.new, f.nonlinear.new[,j])
}

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
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=origin, colour = covariates, linetype = "true"), linewidth = 2) + 
  geom_line(aes(y=new, colour = covariates, linetype = "MAP"), linewidth = 2) + 
  ylab ("") + facet_grid(covariates ~ .) + xlab("Smooth Functions") +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave("../Laboratory/Simulation/BayesianFusedLasso/results/map_smooth_n1.pdf", width=10.5, height = 15)

ggplot(func.df, aes(x=seq, group=interaction(covariates, replicate))) +  
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=origin.linear, colour = covariates, linetype = "true"), linewidth = 2) + 
  geom_line(aes(y=new.linear, colour = covariates, linetype = "MAP"), linewidth = 2) + 
  facet_grid(covariates ~ .) + xlab("Linear Component") + ylab ("") +
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,-20,0,-30),
        strip.text = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave("../Laboratory/Simulation/BayesianFusedLasso/results/map_linear_n1.pdf", width=10, height = 15)

ggplot(func.df, aes(x=seq, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=origin.nonlinear, colour = covariates, linetype = "true"), linewidth = 2) + 
  geom_line(aes(y=new.nonlinear, colour = covariates, linetype = "MAP"), linewidth=2) + ylab ("") + xlab("Nonlinear Component") + facet_grid(covariates ~ .) + 
  scale_linetype_manual("functions",values=c("MAP"=3,"true"=1)) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=45),
        legend.margin=margin(0,0,0,-10),
        legend.box.margin=margin(-10,0,-10,0),
        plot.margin = margin(0,0,0,-20),
        strip.text = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave("../Laboratory/Simulation/BayesianFusedLasso/results/map_nonlinear_n1.pdf", width=12, height = 15)


data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "trueAlp" = sort(true.alpha),
                            "mapAlp" = sort(new.alpha))

plt.samp <- ggplot(data = data.scenario, aes(x = constant)) + ylab(expression(alpha(x))) + xlab("")
print(plt.samp + 
      geom_line(aes(y = trueAlp, col = paste0("True Alpha:",n,"/",psi,"/",threshold)), linewidth = 2.5) + 
      geom_line(aes(y = mapAlp, col = "MAP Alpha"), linewidth = 2.5, linetype = 2) +
      labs(col = "") +
        theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
        theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
        scale_color_manual(values = c("#e0b430", "red"))+
        theme(text = element_text(size = 15),
                legend.position="bottom", legend.key.size = unit(1, 'cm'),
                axis.text = element_text(size = 20),
                legend.margin=margin(-15,-15,-15,-15),
                legend.box.margin=margin(-25,0,20,0)))

# ggsave(paste0("../Laboratory/Simulation/BayesianFusedLasso/results/map_alpha_n1_test.pdf"), width=10)

# Randomized quantile residuals
r <- matrix(NA, nrow = n, ncol = 20)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 20
for(i in 1:n){
  for(t in 1:T){
    r[i, t] <- qnorm(pPareto(y.origin[i], u, alpha = newalpha[i]))
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
  theme(text = element_text(size = 30)) + 
  coord_fixed(xlim = c(-3, 3),  
              ylim = c(-3, 3))
# ggsave("../../Figures/residuals_application.pdf")
# ggsave("../Results/residuals_application_IWSM.pdf", width=6, height=6)

#--------------------------------------------------------------------------------------
# plotting bsplines basis function of covariates


# plts <- list()
# for(j in 1:p){
#   x.x <- sort(scale(x.origin[,j]))
#   # # Make a matrix containing the thin plate basis
#   ndx <- 20
#   basis.x <- data.frame(x=x.x, y=y.origin, id = as.factor(ndx))
#   # deg <- 3
#   # B <- bbase(x.x, min(x.x), max(x.x), nseg = ndx, bdeg = deg)
#   knots <- seq(min(x.x), max(x.x), length.out = (ndx))
#   B <- basis.tps(x.x, knots, m = 2, rk = FALSE)
#   nb1 <- ncol(B)

#   # A basis for plotting the fit on the grid xg
#   ng <- length(x.x)
#   xg <- seq(min(x.x), max(x.x), length.out = ng)
#   Bg <- basis.tps(xg, seq(min(x.x), max(x.x), length.out = ndx), m = 2, rk = FALSE)
#   # Bg <- bbase(xg, min(x.x), max(x.x), nseg = ndx, bdeg = deg)

#   a <- solve(t(B) %*% B, t(B) %*% y.origin, tol = 1e-20)
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
#     ggtitle(paste0("Basis Function of X_" ,j)) +
#     # geom_hline(yintercept = 0, linewidth = 0.3) +
#     geom_line(data = Zf1, linewidth = 1, colour = "grey") +
#     # geom_point(data = basis.x, colour = "grey60", size = 0.8, shape=1) +
#     geom_line(linewidth = 0.7) +
#     # geom_point(data = basis.x, colour = "red", size = 2, shape = 1) +
#     xlab("") + ylab ("") +
#     JOPS_theme() +
#     theme(legend.position = "none") +
#     scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

#   plts[[j]] <- plt
#   # ggsave(paste0("./Results/vert_bsplines_",j,".png"))
# }
# # plts
# grid.arrange(grobs = plts, ncol = 2, nrow = 5)

# # Bspline basis
# plts <- list()
# for(j in 1:p){
#   x.x <- sort(scale(x.origin[,j]))
#   # # Make a matrix containing the B-spline basis
#   ndx <- 13
#   basis.x <- data.frame(x=x.x, y=y.origin, id = as.factor(ndx))
#   deg <- 3
#   B <- bbase(x.x, min(x.x), max(x.x), nseg = ndx, bdeg = deg)
#   nb1 <- ncol(B)

#   # A basis for plotting the fit on the grid xg
#   ng <- length(x.x)
#   xg <- seq(min(x.x), max(x.x), length.out = ng)
#   # Bg <- basis.tps(xg, seq(min(x.x), max(x.x), length.out = ndx), m = 2, rk = FALSE)
#   Bg <- bbase(xg, min(x.x), max(x.x), nseg = ndx, bdeg = deg)

#   a <- solve(t(B) %*% B, t(B) %*% y.origin, tol = 1e-20)
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
#     ggtitle(paste("Basis Function of ", colnames(x.origin)[j])) +
#     # geom_hline(yintercept = 0, linewidth = 0.3) +
#     geom_line(data = Zf1, linewidth = 1, colour = "grey") +
#     # geom_point(data = basis.x, colour = "grey60", size = 0.8, shape=1) +
#     geom_line(linewidth = 0.7) +
#     # geom_point(data = basis.x, colour = "red", size = 2, shape = 1) +
#     xlab("") + ylab ("") +
#     JOPS_theme() +
#     theme(legend.position = "none") +
#     scale_color_manual(values = rainbow_hcl(nb1 + 1, start = 10, end = 350))

#   plts[[j]] <- plt
#   # ggsave(paste0("./Results/vert_bsplines_",j,".png"))
# }
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
