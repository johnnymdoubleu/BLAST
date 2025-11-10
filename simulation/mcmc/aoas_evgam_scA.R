library(npreg)
library(Pareto)
suppressMessages(library(tidyverse))
library(rstan)
library(MESS)
library(evgam)
library(VGAM)
# Scenario A

total.iter <- 2

n <- n.origin <- 15000
psi <- 10
threshold <- 0.95
p <- 5
newp <- p+1
no.theta <- 1

C <- diag(p)
## Generate sample
gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
  for (ps in 1:psi){
    if(j %in% c(1,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
    else {
      if(ps == 1 || ps == psi){gamma.origin[ps, j] <- 0}
      else{gamma.origin[ps, j] <- -25}
    }
  }
}
theta.origin <- c(-0.5, 0, -0.5, -0.5, 0, 0)

model.stan <- "// Stan model for BRSTIR Pareto Uncorrelated Samples
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    array[n] real <lower=1> y; // extreme response
    real <lower=0> atau;
    matrix[2, (2*p)] basisFL;
    array[(p*2)] int indexFL;
    array[n] real <lower=0> trueAlpha; // True Alpha
}
parameters {
    vector[(p+1)] theta; // linear predictor
    array[p] vector[(psi-2)] gammaTemp; // constraint splines coefficient from 2 to psi-1
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    array[p] real <lower=0> tau;
}

transformed parameters {
    array[n] real <lower=0> alpha; // covariate-adjusted tail index
    array[n] real <lower=0> gridalpha; // new tail index
    array[n] real <lower=0> se; // squared error
    real <lower=0> mse; // mean squared error
    array[p] vector[psi] gamma; // splines coefficient
    matrix[n, p] gridgsmooth; // linear component
    
    
    {
      array[p] vector[2] gammaFL;
      matrix[2, p] subgnl;
      matrix[n, p] gnl; // nonlinear component
      matrix[n, p] gl; // linear component
      matrix[n, p] gsmooth; // linear component
      matrix[n, p] gridgnl; // nonlinear component
      matrix[n, p] gridgl; // linear component

      for(j in 1:p){
          gamma[j][2:(psi-1)] = gammaTemp[j][1:(psi-2)];
          subgnl[,j] = bsNonlinear[indexFL[(((j-1)*2)+1):(((j-1)*2)+2)], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j];
          gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * subgnl[,j] * (-1);
          gamma[j][1] = gammaFL[j][1];
          gamma[j][psi] = gammaFL[j][2];  
      };
      
      for (j in 1:p){
          gnl[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
          gridgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
          gl[,j] = bsLinear[,j] * theta[j+1];
          gridgl[,j] = xholderLinear[,j] * theta[j+1];
          gsmooth[,j] = gl[,j] + gnl[,j];
          gridgsmooth[,j] = gridgl[,j] + gridgnl[,j];
      };

      for (i in 1:n){
          alpha[i] = exp(theta[1] + sum(gsmooth[i,])); 
          gridalpha[i] = exp(theta[1] + sum(gridgsmooth[i,]));
          se[i] = pow((gridalpha[i]-trueAlpha[i]), 2);
      };
      mse = mean(se);
    }
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += normal_lpdf(theta[1] | 0, 100);
    target += gamma_lpdf(lambda1 | 1, 1e-3);
    target += gamma_lpdf(lambda2 | 1, 1e-3);
    target += (2*p*log(lambda2));
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, lambda2^2*0.5);
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * (1/tau[j]));
    }
}
"

newgsmooth.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
smooth.scale.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
smooth.1.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
alpha.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
evgam.1.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
evgam.scale.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
alpha.lower.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
alpha.upper.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
mise.scale.container <- mise.1.container <- mise.container <- c()

for(iter in 1:total.iter){
  n <- n.origin
  x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
  
  end.holder <- basis.holder <- matrix(, nrow = 2, ncol =0)
  index.holder <- matrix(, nrow = 0, ncol = 2)
  for(i in 1:p){
    index.holder <- rbind(index.holder, 
                          matrix(c(which.min(x.origin[,i]),
                                    which.max(x.origin[,i])), ncol=2))
  }
  xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
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
  }
  g.sub.origin <- matrix(, nrow = 2, ncol = p)
  for(j in 1:p){
    g.sub.origin[,j] <- as.matrix(bs.nonlinear[index.holder[j,], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))], nrow = 2) %*% gamma.origin[(2:(psi-1)), j]
    gamma.origin[c(1,psi),j] <- -1 * basis.holder[,(((j-1)*2)+1):(((j-1)*2)+2)] %*% as.matrix(g.sub.origin[,j], nrow=2)
  }
  
  g.nonlinear.origin <- g.linear.origin <- g.origin <- matrix(, nrow = n, ncol = p)
  for(j in 1:p){
      g.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
      g.nonlinear.origin[,j] <- bs.nonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
      g.origin[, j] <- g.linear.origin[,j] + g.nonlinear.origin[,j]
  }

  alp.origin <- y.origin <- NULL
  for(i in 1:n){
      alp.origin[i] <- exp(theta.origin[1] + sum(g.origin[i,]))
      y.origin[i] <- rPareto(1, 1, alpha = alp.origin[i])
  }

  u <- quantile(y.origin, threshold)
  x.origin <- x.origin[which(y.origin>u),]
  y.origin <- y.origin[y.origin > u]
  n <- length(y.origin)

  xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
  newx <- seq(0, 1, length.out=n)
  xholder <- bs.x <- matrix(, nrow = n, ncol = p)
  end.holder <- basis.holder <- matrix(, nrow = 2, ncol =0)
  index.holder <- matrix(, nrow = 0, ncol = 2)
  for(i in 1:p){
    index.holder <- rbind(index.holder, 
                          matrix(c(which.min(x.origin[,i]),
                                    which.max(x.origin[,i])), ncol=2))
    xholder[,i] <- seq(0, 1, length.out = n)
  }
  for(i in 1:p){
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

  g.nonlinear.new <- g.linear.new <- g.new <- g.nonlinear.origin <- g.linear.origin <- g.origin <- matrix(, nrow = n, ncol = p)
  for(j in 1:p){
      g.linear.origin[,j] <- bs.linear[, j] * theta.origin[j+1]
      g.nonlinear.origin[,j] <- bs.nonlinear[1:n,(((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
      g.origin[, j] <- g.linear.origin[,j] + g.nonlinear.origin[,j]
      g.linear.new[,j] <- xholder.linear[, j] * theta.origin[j+1]
      g.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% gamma.origin[,j]
      g.new[,j] <- g.linear.new[,j] + g.nonlinear.new[,j]
  }

  true.alpha <- alp.new <- alp.origin <- NULL
  for(i in 1:n){
    alp.origin[i] <- exp(theta.origin[1] + sum(g.origin[i,]))
    alp.new[i] <- exp(theta.origin[1] + sum(g.new[i,]))
  }
  
  data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), basisFL = basis.holder,
                    indexFL = as.vector(t(index.holder)), trueAlpha = alp.new,
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear)
  
  init.alpha <- list(list(gammaTemp = array(rep(2, ((psi-2)*p)), dim=c(p,(psi-2))),
                          theta = rep(0, (p+1)), tau1 = rep(0.1, p),tau2 = rep(0.1, p),
                          lambda1 = 0.1, lambda2 = 1),
                      list(gammaTemp = array(rep(-1, ((psi-2)*p)), dim=c(p,(psi-2))),
                          theta = rep(0, (p+1)), tau1 = rep(0.001, p), tau2=rep(0.001, p),
                          lambda1 = 100, lambda2 = 100),
                      list(gammaTemp = array(rep(-3, ((psi-2)*p)), dim=c(p,(psi-2))),
                          theta = rep(0.1, (p+1)), tau1 = rep(0.5, p),tau2 = rep(0.5, p),
                          lambda1 = 5, lambda2 = 55))
    
  fit1 <- stan(
      model_code = model.stan,  # Stan program
      data = data.stan,    # named list of data
      init = init.alpha,      # initial value
      chains = 3,             # number of Markov chains
      # warmup = 1000,          # number of warmup iterations per chain
      iter = 2000,            # total number of iterations per chain
      cores = parallel::detectCores(), # number of cores (could use one per chain)
      refresh = 500             # no progress shown
  )

  newgsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
  newalpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary
  se.samples <- summary(fit1, par=c("se"), probs = c(0.05,0.5, 0.95))$summary

  
  simul.data <- data.frame(y = y.origin, x.origin)
  vgam.fit.scale <- VGAM::vgam(y ~ s(X1, bs = "tp", k = 10) + s(X2, bs = "tp", k = 10) + s(X3, bs = "tp", k = 10) + s(X4, bs = "tp", k = 10) + s(X5, bs = "tp", k = 10),
                        data = simul.data,
                        family = gpd(threshold= 0,
                                      # lscale="loglink", 
                                      lshape="loglink",
                                      zero = NULL),
                        trace = TRUE,
                        control = vgam.control(maxit = 200))
  fitted.terms <- predict(vgam.fit.scale,newdata = data.frame(xholder), type = "terms")
  vgam.xi.scale <- exp(rowSums(fitted.terms[,c(2, 4, 6, 8, 10)]))
  vgam.fit.1 <- VGAM::vgam(y ~ s(X1, bs = "tp", k = 10) + s(X2, bs = "tp", k = 10) + s(X3, bs = "tp", k = 10) + s(X4, bs = "tp", k = 10) + s(X5, bs = "tp", k = 10),
                        data = simul.data,
                        family = gpd(threshold= 0,
                                      lscale="loglink", 
                                      lshape="loglink",
                                      zero = 1),
                        trace = TRUE,
                        control = vgam.control(maxit = 200))
  fitted.terms <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "terms")
  vgam.xi.1 <- exp(rowSums(fitted.terms[,c(2, 4, 6, 8, 10)]))  
  # gam.scale <- list(y ~ s(X1, bs = "tp", k = 10) + 
  #                       s(X2, bs = "tp", k = 10) + 
  #                       s(X3, bs = "tp", k = 10) + 
  #                       s(X4, bs = "tp", k = 10) + 
  #                       s(X5, bs = "tp", k = 10),
  #                     ~ s(X1, bs = "tp", k = 10) + 
  #                       s(X2, bs = "tp", k = 10) + 
  #                       s(X3, bs = "tp", k = 10) + 
  #                       s(X4, bs = "tp", k = 10) + 
  #                       s(X5, bs = "tp", k = 10))
  # evgam.fit.scale <- evgam::evgam(gam.scale, data = simul.data, family = "gpd")
  # xi.pred.scale <-predict(evgam.fit.scale, newdata = data.frame(xholder), type="response")$shape
  # alpha.pred.scale <- 1/xi.pred.scale

  # xholder.basis.scale <- predict(evgam.fit.scale, newdata = data.frame(xholder), type= "lpmatrix")$shape
  # xi.coef.scale <- tail(evgam.fit.scale$coefficients, (psi-1)*p)
  # gamma.xi.scale <- matrix(xi.coef.scale, ncol = p)
  # alpha.nonlinear.scale <- xi.nonlinear.scale <- matrix(, nrow = n, ncol = p)
  # bs.nonlinear.scale <- xholder.basis.scale[,c(2:((psi-1)*p+1))]
  # for(j in 1:p){
  #   xi.nonlinear.scale[,j] <- bs.nonlinear.scale[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.scale[,j]
  #   alpha.nonlinear.scale[,j] <- 1/(xi.nonlinear.scale[,j])
  # }

  # gam.1 <- list(y ~ 1,
  #                 ~ s(X1, bs = "tp", k = 10) + 
  #                   s(X2, bs = "tp", k = 10) + 
  #                   s(X3, bs = "tp", k = 10) + 
  #                   s(X4, bs = "tp", k = 10) + 
  #                   s(X5, bs = "tp", k = 10))
  # evgam.fit.1 <- evgam::evgam(gam.1, data = simul.data, family = "gpd")
  # xi.pred.1 <-predict(evgam.fit.1, newdata = data.frame(xholder), type="response")$shape
  # alpha.pred.1 <- 1/xi.pred.1

  # xholder.basis.1 <- predict(evgam.fit.1, newdata = data.frame(xholder), type= "lpmatrix")$shape
  # xi.coef.1 <- tail(evgam.fit.1$coefficients, (psi-1)*p)
  # gamma.xi.1 <- matrix(xi.coef.1, ncol = p)
  # alpha.nonlinear.1 <- xi.nonlinear.1 <- matrix(, nrow = n, ncol = p)
  # bs.nonlinear.1 <- xholder.basis.1[,c(2:((psi-1)*p+1))]
  # for(j in 1:p){
  #   xi.nonlinear.1[,j] <- bs.nonlinear.1[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.1[,j]
  #   alpha.nonlinear.1[,j] <- 1/(xi.nonlinear.1[,j])
  # }


  evgam.1.container[,iter] <- vgam.xi.1
  evgam.scale.container[,iter] <- vgam.xi.scale
  alpha.lower.container[,iter] <- 1/(newalpha.samples[,4])
  alpha.container[,iter] <- 1/(newalpha.samples[,5])
  alpha.upper.container[,iter] <- 1/(newalpha.samples[,6])
  newgsmooth.container[,iter] <- as.vector(matrix(newgsmooth.samples[,5], nrow = n, byrow=TRUE))
  # smooth.1.container[,iter] <- as.vector(xi.nonlinear.1)
  # smooth.scale.container <- as.vector(xi.nonlinear.scale)
  
  mise.container[iter] <- auc(newx, se.samples[,5], type="spline")
  mise.1.container[iter] <- auc(newx, ((1/alp.new)-vgam.xi.1)  ,type="spline")
  mise.scale.container[iter] <- auc(newx, ((1/alp.new)-vgam.xi.scale)  ,type="spline")
}



alpha.container$x <- seq(0,1, length.out = n)
evgam.1.container$x <- seq(0,1, length.out = n)
alpha.container$true <- 1/alp.new
alpha.container <- cbind(alpha.container, t(apply(alpha.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
colnames(alpha.container)[(dim(alpha.container)[2]-2):(dim(alpha.container)[2])] <- c("q1","q2","q3")
alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
# alpha.container$q1 <- apply(alpha.lower.container[,1:total.iter], 1, quantile, c(.5))
# alpha.container$q3 <- apply(alpha.upper.container[,1:total.iter], 1, quantile, c(.5))
alpha.container$evgam.1 <- rowMeans(evgam.1.container[,1:total.iter])
alpha.container$evgam.scale <- rowMeans(evgam.scale.container[,1:total.iter])
alpha.container <- as.data.frame(alpha.container)

save(newgsmooth.container, alpha.container, evgam.1.container, evgam.scale.container, mise.container, mise.evgam.container, file="./simulation/results/vgam_mc_scA.Rdata")
# load("./simulation/results/evgam_mc_scA.Rdata")

plt <- ggplot(data = alpha.container, aes(x = x)) + xlab(expression(c)) + labs(col = "") + ylab(expression(xi(c,ldots,c))) #+ ylab("")
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.075, linewidth = 0.7)
  }
} else{
  for(i in 1:total.iter){
  # for(i in 50:100){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.075, linewidth = 0.7)
    plt <- plt + geom_line(data=evgam.1.container, aes(x=x, y = .data[[names(evgam.1.container)[i]]]), alpha = 0.075, linewidth = 0.7)
  }
}

print(plt +
        geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.8) +
        geom_line(aes(y=evgam.1), colour="purple", linewidth = 1.8) +
        geom_line(aes(y=evgam.scale), colour="orange", linewidth = 1.8) +
        scale_fill_manual(values=c("steelblue"), name = "") +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) + #ylim(-1, 2.4) +
        theme(legend.position = "none",
                strip.text = element_blank(),
                axis.text = element_text(size = 18)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_evgam_scA_.pdf"), width=10, height = 7.78)





# equal_breaks <- function(n = 3, s = 0.1,...){
#   function(x){
#     d <- s * diff(range(x)) / (1+2*s)
#     seq = seq(min(x)+d, max(x)-d, length=n)
#     round(seq, -floor(log10(abs(seq[2]-seq[1]))))
#   }
# }

# newgsmooth.container$x <- seq(0,1, length.out = n)
# newgsmooth.container$true <- as.vector(g.new)
# newgsmooth.container <- cbind(newgsmooth.container, t(apply(newgsmooth.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
# colnames(newgsmooth.container)[(dim(newgsmooth.container)[2]-2):(dim(newgsmooth.container)[2])] <- c("q1","q2","q3")
# newgsmooth.container$mean <- rowMeans(newgsmooth.container[,1:total.iter])
# newgsmooth.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]"))
# newgsmooth.container <- as.data.frame(newgsmooth.container)

# plt <- ggplot(data = newgsmooth.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(c))
# if(total.iter <= 50){
#   for(i in 1:total.iter){
#     plt <- plt + geom_line(aes(y = .data[[names(newgsmooth.container)[i]]]), alpha = 0.075, linewidth = 0.7)
#   }
# } else{
#   for(i in 1:total.iter){
#     plt <- plt + geom_line(aes(y = .data[[names(newgsmooth.container)[i]]]), alpha = 0.075, linewidth = 0.7)
#   }
# }
# print(plt + 
#         geom_line(aes(y=true, col = "True"), linewidth = 2, linetype = 2) + 
#         geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5) + 
#         ylim(-1, 1) +
#         facet_grid(covariate ~ ., scales = "free_x", switch = "y",
#                     labeller = label_parsed) +        
#         scale_color_manual(values = c("steelblue", "red"))+
#         guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) +
#         theme_minimal(base_size = 30) +
#         theme(legend.position = "none",
#                 plot.margin = margin(0,0,0,-20),
#                 strip.text.y = element_text(size = 25, colour = "black", angle = 0, face = "bold.italic"),
#                 strip.placement = "outside",
#                 axis.title.x = element_text(size = 35),                
#                 axis.text = element_text(size = 18)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_smooth_scA_",n.origin,".pdf"), width=12.5, height = 15)


print(mean(mise.container))
print(mean(mise.1.container))
print(mean(mise.scale.container))
