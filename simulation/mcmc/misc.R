library(npreg)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(parallel)
library(qqboxplot)

#Scenario 1
set.seed(10)

n <- 15000
psi <- 10
threshold <- 0.95
p <- 5
no.theta <- 1
simul.no <- 50

# Function to generate Gaussian copula
C <- diag(p)
xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
end.holder <- basis.holder <- matrix(, nrow = 2, ncol =0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                        matrix(c(which.min(x.origin[,i]),
                                 which.max(x.origin[,i])), ncol=2))
}
for(i in 1:p){
  knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
  tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  bs.linear <- cbind(bs.linear, tps[,1:no.theta])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
  basis.holder <- cbind(basis.holder, 
                        solve(t(matrix(c(tps[index.holder[i,1], no.theta+1],
                                         tps[index.holder[i,1], no.theta+psi],
                                         tps[index.holder[i,2], no.theta+1],
                                         tps[index.holder[i,2], no.theta+psi]), 
                                       nrow = 2, ncol = 2))))
  end.holder <- cbind(end.holder, 
                      matrix(c(tps[index.holder[i,1], no.theta+1],
                               tps[index.holder[i,1], no.theta+psi],
                               tps[index.holder[i,2], no.theta+1],
                               tps[index.holder[i,2], no.theta+psi]), 
                             nrow = 2, ncol = 2))
}

## Generate sample
gamma.origin <- matrix(, nrow=psi, ncol=p)
for(j in 1:p){
  for (ps in 1:psi){
    if(j %in% c(1,4,5)){gamma.origin[ps,j] <- 0}
    else {
      if(ps==1 || ps==psi){gamma.origin[ps,j] <- 0}
      else{gamma.origin[ps,j] <- -25}
    }
  }
}

# for(j in 1:p){
#   for (ps in 1:psi){
#     if(j %in% c(1,4,5)){gamma.origin[ps,j] <- 0}
#     else if(j==2){
#       if(ps==1 || ps==psi){gamma.origin[ps,j] <- 0}
#       else{gamma.origin[ps,j] <- -25}
#     }
#     else{
#       if(ps==1){gamma.origin[ps,j] <- 0}
#       else if(ps==psi){gamma.origin[ps,j] <- -.5}
#       else if(ps>1 & ps<=(psi/2)){gamma.origin[ps,j] <- -1.5}
#       else if(ps<psi & ps>(psi/2)){gamma.origin[ps,j] <- 1.5}
#     }
#   }
# }
theta.origin <- c(-0.5, 0, -0.5, -0.5, 0, 0)


f.sub.origin <- matrix(, nrow = 2, ncol = p)
for(j in 1:p){
  f.sub.origin[,j] <- as.matrix(bs.nonlinear[index.holder[j,], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))], nrow = 2) %*% gamma.origin[(2:(psi-1)), j]
  gamma.origin[c(1,psi),j] <- -1 * basis.holder[,(((j-1)*2)+1):(((j-1)*2)+2)] %*% as.matrix(f.sub.origin[,j], nrow=2)
}

f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)
for(j in 1:p){
  f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
  f.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
  f.origin[, j] <- f.linear.origin[,j] + f.nonlinear.origin[,j]
}

alp.origin <- y.origin <- NULL
for(i in 1:n){
  alp.origin[i] <- exp(theta.origin[1] + sum(f.origin[i,]))
  y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i]) 
}

u <- quantile(y.origin, threshold)
excess.index <- which(y.origin>u)
x.origin <- as.matrix(x.origin[excess.index,])

y.origin <- y.origin[y.origin > u]
n <- length(y.origin)

xholder.nonlinear <- xholder.linear <-  matrix(,nrow=n, ncol=0)
bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
newx <- seq(0, 1, length.out=n)
xholder <- bs.x <- matrix(, nrow = n, ncol = p)
end.holder <- basis.holder <- matrix(, nrow = 2, ncol =0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                        matrix(c(which.min(x.origin[,i]),
                                 which.max(x.origin[,i])), ncol=2))
}

for(i in 1:p){
  xholder[,i] <- seq(0, 1, length.out = n)  
  test.knot <- seq(0, 1, length.out = psi)
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)  
  tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  basis.holder <- cbind(basis.holder, 
                        solve(t(matrix(c(tps[index.holder[i,1], no.theta+1],
                                         tps[index.holder[i,1], no.theta+psi],
                                         tps[index.holder[i,2], no.theta+1],
                                         tps[index.holder[i,2], no.theta+psi]), 
                                       nrow = 2, ncol = 2))))
  bs.linear <- cbind(bs.linear, tps[,1:no.theta])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)]) 
}

f.nonlinear.new <- f.linear.new <- f.new <- f.nonlinear.origin <- f.linear.origin <- f.origin <- matrix(, nrow = n, ncol = p)

for(j in 1:p){
  f.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
  f.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
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

# as.vector(f.new[,3]*1.8) + (as.vector(f.linear.new[,3]))
data.constraint <- data.frame("x"=newx,
                            "true.smooth" = as.vector(f.new[,2]),
                            # "true.linear" = as.vector(f.linear.new[,2]),
                            # "true.nonlinear" = as.vector(f.nonlinear.new[,2]))
                            # "true.linear" =  as.vector(f.linear.new[,3])*0.1,
                            # "true.nonlinear" = as.vector(f.nonlinear.new[,3]))
                            "true.linear" = as.vector(f.linear.new[,3] * -1) - 0.3,
                            "true.nonlinear" = as.vector(f.new[,3]*1.5) + (as.vector(f.linear.new[,3])) + 0.27)

ggplot(data.constraint, aes(x=x)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=true.smooth), color = "steelblue", linewidth=2.5) + 
  ylab("g(x)") + xlab("x") + ylim(-0.8, 0.8) + 
  theme_minimal(base_size = 30) +
  theme(legend.position="none",
        axis.title = element_text(size = 45),
        axis.text = element_text(size=40))
ggsave(paste0("./simulation/results/illust_smooth1_new.pdf"), width=11.5, height = 8)

ggplot(data.constraint, aes(x=x)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=true.linear), color = "steelblue", linewidth=2.5) + 
  ylab("") + xlab("x") + ylim(-0.8, 0.8) + 
  theme_minimal(base_size = 30) +
  theme(legend.position="none",
        plot.margin = margin(0,0,0,-20),
        axis.title.x = element_text(size = 45),
        axis.text = element_text(size=40))
ggsave(paste0("./simulation/results/illust_linear_nc2_new.pdf"), width=10, height = 7.78)

ggplot(data.constraint, aes(x=x)) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_line(aes(y=true.nonlinear), color = "steelblue", linewidth=2.5) + 
  ylab("") + xlab("x") + ylim(-0.8, 0.8) + 
  theme_minimal(base_size = 30) +
  theme(legend.position="none",
        plot.margin = margin(0,0,0,-20),
        axis.title.x = element_text(size = 45),
        axis.text = element_text(size=40))
ggsave(paste0("./simulation/results/illust_nonlinear_nc2_new.pdf"), width=10, height = 7.78)
