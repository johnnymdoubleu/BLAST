library(npreg)
suppressMessages(library(ggplot2))
library(rstan)
library(Pareto)
library(evgam)
library(gridExtra)
library(VGAM)

#Scenario 1
# set.seed(10)
set.seed(111)

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

# x.origin <- apply(x.origin, 2, sort, decreasing=F)
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
# bs.nonlinear <- bs.nonlinear[excess.index,]
# bs.linear <- bs.linear[excess.index,]
y.pre <- y.origin
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
# # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# # x.origin <- as.data.frame(sapply(as.data.frame(x.origin), FUN = range01))
for(i in 1:p){
  # xholder[,i] <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = n)  
  # test.knot <- seq(min(xholder[,i]), max(xholder[,i]), length.out = psi)
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

stan.code <- "// Stan model for BRSTIR Pareto Uncorrelated Samples
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
    array[p] vector[psi] gamma; // splines coefficient
    
    {
      array[p] vector[2] gammaFL;
      matrix[2, p] subgnl;
      matrix[n, p] gnl; // nonlinear component
      matrix[n, p] gl; // linear component
      matrix[n, p] gsmooth; // linear component


      for(j in 1:p){
          gamma[j][2:(psi-1)] = gammaTemp[j][1:(psi-2)];
          subgnl[,j] = bsNonlinear[indexFL[(((j-1)*2)+1):(((j-1)*2)+2)], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j];
          gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * subgnl[,j] * (-1);
          gamma[j][1] = gammaFL[j][1];
          gamma[j][psi] = gammaFL[j][2];  
      };
      
      for (j in 1:p){
          gnl[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
          gl[,j] = bsLinear[,j] * theta[j+1];
          gsmooth[,j] = gl[,j] + gnl[,j];
      };

      for (i in 1:n){
          alpha[i] = exp(theta[1] + sum(gsmooth[i,])); 
      };
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

generated quantities {
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgsmooth; // linear component
    {
      matrix[n, p] newgnl; // nonlinear component
      matrix[n, p] newgl; // linear component
      for (j in 1:p){
            newgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
            newgl[,j] = xholderLinear[,j] * theta[j+1];
            newgsmooth[,j] = newgl[,j] + newgnl[,j];
        };
    }
    for (i in 1:n){ 
      newalpha[i] = exp(theta[1] + sum(newgsmooth[i,]));
    };
}
"

data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                  atau = ((psi+1)/2), basisFL = basis.holder,
                  indexFL = as.vector(t(index.holder)),
                  bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                  xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear)

init.alpha <- list(list(gammaTemp = array(rep(2, ((psi-2)*p)), dim=c(p,(psi-2))),
                        theta = rep(0, (p+1)), tau = rep(0.1, p),
                        lambda1 = 0.1, lambda2 = 1),
                   list(gammaTemp = array(rep(-1, ((psi-2)*p)), dim=c(p,(psi-2))),
                        theta = rep(0, (p+1)), tau = rep(0.001, p),
                        lambda1 = 100, lambda2 = 100),
                   list(gammaTemp = array(rep(-3, ((psi-2)*p)), dim=c(p,(psi-2))),
                        theta = rep(0.1, (p+1)), tau = rep(0.5, p),
                        lambda1 = 5, lambda2 = 55))

system.time(fit1 <- stan(
  model_code = stan.code,  # Stan program
  data = data.stan,    # named list of data
  init = init.alpha,      # initial value
  chains = 3,             # number of Markov chains
  # warmup = 1000,          # number of warmup iterations per chain
  iter = 4000,            # total number of iterations per chain
  cores = parallel::detectCores(), # number of cores (could use one per chain)
  refresh = 500             # no progress shown
))

posterior <- extract(fit1)

tau.samples <- summary(fit1, par=c("tau"), probs = c(0.05,0.5, 0.95))$summary
theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
newgsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05, 0.5, 0.95))$summary
newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary

# save(theta.samples, gamma.samples, lambda.samples, newgsmooth.samples, newalpha.samples, file = "./simulation/results/blast_scA.Rdata")

gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q2 <- gamma.samples[,1]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q2 <- theta.samples[,5]
theta.q3 <- theta.samples[,6]


g.smooth.mean <- as.vector(matrix(newgsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(newgsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(newgsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(newgsmooth.samples[,6], nrow = n, byrow=TRUE))

# equal_breaks <- function(n = 3, s = 0.1,...){
#   function(x){
#     d <- s * diff(range(x)) / (1+2*s)
#     seq = seq(min(x)+d, max(x)-d, length=n)
#     round(seq, -floor(log10(abs(seq[2]-seq[1]))))
#   }
# }

simul.data <- data.frame(y = y.origin, x.origin)

psi <- 10
vgam.fit.scale <- vgam(y ~ s(X1) + s(X2) + s(X3) + s(X4) + s(X5),
# vgam.fit.scale <- vgam(y ~ sm.ps(X1) + sm.ps(X2) + sm.ps(X3) + sm.ps(X4) + sm.ps(X5),
# vgam.fit.scale <- vgam(y ~ bs(X1, knots=knots) + bs(X2, knots=knots) + bs(X3, knots=knots) + bs(X4, knots=knots) + bs(X5, knots=knots),
                        data = simul.data,
                        family = gpd(threshold = u,
                                      lshape="loglink",
                                      zero = NULL),
                        trace = TRUE,
                        control = vgam.control(maxit = 200))
par(mfrow = c(5, 2), mar=c(1.5,1.5,1.5,1.5))
plot(vgam.fit.scale, se = TRUE, shade = TRUE, shcol = "steelblue")
par(mfrow = c(1, 1))
fitted.linear <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "link")
fitted.terms <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "terms")
fitted.response <- predict(vgam.fit.scale,newdata = data.frame(xholder), type = "response")
# summary(vgam.fit)
vgam.xi.scale <- exp(fitted.linear[,2])#exp(rowSums(fitted.terms[,c(2, 4, 6, 8, 10)]))
vgam.fit.1 <- vgam(y ~ s(X1) + s(X2) + s(X3) + s(X4) + s(X5),
# vgam.fit.1 <- vgam(y ~ sm.ps(X1) + sm.ps(X2) + sm.ps(X3) + sm.ps(X4) + sm.ps(X5),
                    data = simul.data,
                    family = gpd(threshold = u, 
                                  lshape="loglink",
                                  zero = 1),
                    trace = TRUE,
                    control = vgam.control(maxit = 200))
par(mfrow = c(2, 3), mar=c(4,4,4,4))
plot(vgam.fit.1, se = TRUE, shade = TRUE, shcol = "steelblue")
par(mfrow = c(1, 1))
fitted.linear <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "link")
fitted.terms <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "terms")
fitted.response <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "response")
# summary(vgam.fit)
vgam.xi.1 <- exp(fitted.linear[,2])#exp(rowSums(fitted.terms[,c(2, 4, 6, 8, 10)]))

simul.evgam <- data.frame(y = y.origin-u, x.origin)
gam.scale <- list(y ~ s(X1, bs = "tp", k = 10) + 
                      s(X2, bs = "tp", k = 10) + 
                      s(X3, bs = "tp", k = 10) + 
                      s(X4, bs = "tp", k = 10) + 
                      s(X5, bs = "tp", k = 10),
                    ~ s(X1, bs = "tp", k = 10) + 
                      s(X2, bs = "tp", k = 10) + 
                      s(X3, bs = "tp", k = 10) + 
                      s(X4, bs = "tp", k = 10) + 
                      s(X5, bs = "tp", k = 10))
evgam.fit.scale <- evgam::evgam(gam.scale, data = simul.evgam, family = "gpd", outer = "Newton")
plot.data <- as.data.frame(evgam.fit.scale$plotdata)
xi.pred.scale <-predict(evgam.fit.scale, newdata = data.frame(xholder), type="response")$shape
# xi.pred <-predict(evgam.fit, newdata = data.frame(xholder))$shape
alpha.pred.scale <- 1/xi.pred.scale

xholder.basis.scale <- predict(evgam.fit.scale, newdata = data.frame(xholder), type= "lpmatrix")$shape
xi.coef.scale <- tail(evgam.fit.scale$coefficients, (psi-1)*p)
gamma.xi.scale <- matrix(xi.coef.scale, ncol = p)
alpha.nonlinear.scale <- xi.nonlinear.scale <- matrix(, nrow = n, ncol = p)
bs.nonlinear.scale <- xholder.basis.scale[,c(2:((psi-1)*p+1))]
for(j in 1:p){
  xi.nonlinear.scale[,j] <- bs.nonlinear.scale[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.scale[,j]
  alpha.nonlinear.scale[,j] <- 1/(xi.nonlinear.scale[,j])
}

gam.1 <- list(y ~ 1,
                ~ s(X1, bs = "tp", k = 10) + 
                  s(X2, bs = "tp", k = 10) + 
                  s(X3, bs = "tp", k = 10) + 
                  s(X4, bs = "tp", k = 10) + 
                  s(X5, bs = "tp", k = 10))
evgam.fit.1 <- evgam::evgam(gam.1, data = simul.evgam, family = "gpd")

# save(evgam.fit.scale, evgam.fit.1, file= "./simulation/results/evgam_scA")

plot.data <- as.data.frame(evgam.fit.1$plotdata)
xi.pred.1 <-predict(evgam.fit.1, newdata = data.frame(xholder), type="response")$shape
alpha.pred.1 <- 1/xi.pred.1

xholder.basis.1 <- predict(evgam.fit.1, newdata = data.frame(xholder), type= "lpmatrix")$shape
xi.coef.1 <- tail(evgam.fit.1$coefficients, (psi-1)*p)
gamma.xi.1 <- matrix(xi.coef.1, ncol = p)
alpha.nonlinear.1 <- xi.nonlinear.1 <- matrix(, nrow = n, ncol = p)
bs.nonlinear.1 <- xholder.basis.1[,c(2:((psi-1)*p+1))]
for(j in 1:p){
  xi.nonlinear.1[,j] <- bs.nonlinear.1[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.1[,j]
  alpha.nonlinear.1[,j] <- 1/(xi.nonlinear.1[,j])
}


alpha.smooth <- data.frame("x" = as.vector(xholder),
                          "true" = as.vector(f.new),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "evgam.scale" = as.vector(alpha.nonlinear.scale),
                          "evgam.1" = as.vector(alpha.nonlinear.1),
                          "covariates" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))


grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(alpha.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=true, colour = "True"), linewidth=1, linetype=2) + 
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1) + 
                  geom_line(aes(y=evgam.1), colour = "orange", linewidth=1) + 
                  ylab("") + xlab(paste0("X[",i,"]")) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue", "red")) +
                  ylim(-1.3, 1.3) + xlim(0,1) +
                  theme_minimal(base_size = 10) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 15),
                          axis.title.x = element_text(size = 15))
  grid.plts[[i]] <- grid.plt
}

grid.arrange(grobs = grid.plts, ncol = 3, nrow = 2)

xi.smooth <- data.frame("x" = as.vector(xholder),
                          "true" = as.vector(1/exp(f.new)),
                          "post.mean" = as.vector(1/exp(g.smooth.mean)),
                          "q1" = as.vector(1/exp(g.smooth.q1)),
                          "q2" = as.vector(1/exp(g.smooth.q2)),
                          "q3" = as.vector(1/exp(g.smooth.q3)),
                          "evgam.scale" = as.vector(xi.nonlinear.scale),
                          "evgam.1" = as.vector(xi.nonlinear.1),
                          "covariates" = gl(p, n, (p*n), labels = c("x[1]", "x[2]", "x[3]", "x[4]", "x[5]")))

grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(xi.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=true, colour = "True"), linewidth=1, linetype = 2) + 
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1) + 
                  geom_line(aes(y=evgam.1), colour = "orange", linewidth=1) + 
                  ylab("") + xlab(paste0("X[",i,"]")) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue", "red")) +
                  ylim(-1.8, 1.8) + xlim(0,1) +
                  theme_minimal(base_size = 10) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 15),
                          axis.title.x = element_text(size = 15))
  grid.plts[[i]] <- grid.plt
}

grid.arrange(grobs = grid.plts, ncol = 3, nrow = 2)

alpha.scenario <- data.frame("x" = newx,
                            "constant" = newx,
                            "true" = (alp.new),
                            "post.mean" = (newalpha.samples[,1]),
                            "post.median" = (newalpha.samples[,5]),
                            "evgam.scale" = 1/xi.pred.scale,
                            "evgam.1" = 1/xi.pred.1,
                            "vgam.scale" = 1/as.vector(vgam.xi.scale),
                            "vgam.1" = 1/as.vector(vgam.xi.1),                            
                            "q1" = (newalpha.samples[,4]),
                            "q3" = (newalpha.samples[,6]))
# "post.mean" = sort(alpha.smooth.new),
# "post.median" = sort(newalpha.samples[,5]),
# "q1" = sort(alpha.smooth.q1),
# "q3" = sort(alpha.smooth.q3))

ggplot(alpha.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,ldots,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y = true, col = "True"), linewidth = 2, linetype=2) + 
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1.5) +
  # geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1.5, linetype=2) +
  # geom_line(aes(y=evgam.1), colour = "orange", linewidth=1.5, linetype=2) +
  geom_line(aes(y=vgam.scale), colour = "purple", linewidth=1.5) +
  geom_line(aes(y=vgam.1), colour = "orange", linewidth=1.5) +  
  scale_color_manual(values=c("steelblue", "red")) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  theme_minimal(base_size = 30) + ylim(0,2.5)+
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",n,"_mcmc_alpha_test_sc1-wi.pdf"), width=10, height = 7.78)

xi.scenario <- data.frame("x" = newx,
                            "constant" = newx,
                            "true" = 1/(alp.new),
                            "post.mean" = 1/(newalpha.samples[,1]),
                            "post.median" = 1/(newalpha.samples[,5]),
                            "evgam.scale" = xi.pred.scale,
                            "evgam.1" = xi.pred.1,
                            "vgam.scale" = as.vector(vgam.xi.scale),
                            "vgam.1" = as.vector(vgam.xi.1),
                            "q1" = 1/(newalpha.samples[,4]),
                            "q3" = 1/(newalpha.samples[,6]))

ggplot(xi.scenario, aes(x=x)) + 
  ylab(expression(xi(c,ldots,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y = true, col = "True"), linewidth = 2, linetype=2) + 
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1.5) +
  # geom_line(aes(y=evgam.scale), colour = "purple", linewidth=1.5, linetype=2) +
  # geom_line(aes(y=evgam.1), colour = "orange", linewidth=1.5, linetype=2) +
  geom_line(aes(y=vgam.scale), colour = "purple", linewidth=1.5) +
  geom_line(aes(y=vgam.1), colour = "orange", linewidth=1.5) +
  scale_color_manual(values=c("steelblue", "red")) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  theme_minimal(base_size = 30) + #ylim(0, 6.1)+
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))

AIC(vgam.fit.1)
AIC(vgam.fit.scale)


# orderred <- rev(sort(y.pre)[14325:length(y.pre)])
orderred <- sort(y.pre, decreasing = TRUE)[1:(n+1)]
n.hill <- length(orderred)
k <- 1:n.hill
loggs <- log(orderred/u)
avesumlog <- cumsum(loggs)/k
xihat <- c(NA, (avesumlog)[2:n.hill])
ses.xi <- xihat * sqrt(k)
xx <- seq(from = n.hill, to = 2)
y.xi <- xihat[xx]
qq <- qnorm(1-(1-threshold)/2)
uu.xi <- y.xi + ses.xi[xx] * qq
ll.xi <- y.xi - ses.xi[xx] * qq

alphahat <- 1/xihat
ses.alpha <- alphahat / sqrt(k)
y.alpha <- alphahat[xx]
uu.alpha <- y.alpha + ses.alpha[xx] * qq
ll.alpha <- y.alpha - ses.alpha[xx] * qq


gamma.alpha <- (xx * y.alpha) / (xx - 1)
gamma.q1.alpha <- 1 / qgamma(0.975, shape = xx, rate = xx * y.alpha)
gamma.q3.alpha <- 1 / qgamma(0.025, shape = xx, rate = xx * y.alpha)
gamma.xi <- (xx * y.xi) / (xx - 1)
gamma.q1.xi <- 1 / qgamma(0.975, shape = xx, rate = xx * y.xi)
gamma.q3.xi <- 1 / qgamma(0.025, shape = xx, rate = xx * y.xi)

data.hill.alpha <- data.frame("k" = c(1:n),
                        "u" = sort(uu.alpha),
                        "l" = sort(ll.alpha),
                        "alpha" = sort(y.alpha),
                        "order" = xx,
                        "blast.mean" = sort(alpha.samples[,1]),
                        "blast.q1" = sort(alpha.samples[,4]),
                        "blast.q3" = sort(alpha.samples[,6]),
                        "gamma.mean" = sort(gamma.alpha),
                        "gamma.q1" = sort(gamma.q1.alpha), 
                        "gamma.q3" = sort(gamma.q3.alpha))
ggplot(data = data.hill.alpha) + 
  geom_ribbon(aes(x = order, ymin = l, ymax = u),
              alpha = 0.2, linetype = "dashed", fill = "black") + 
  geom_line(aes(x = order, y = alpha), linewidth = 1.2, colour = "black") +
  geom_line(aes(x = order, y = blast.mean), linewidth = 1.2, colour = "steelblue") + 
  geom_ribbon(aes(x = order, ymin = blast.q1, ymax = blast.q3),
              alpha = 0.2, linetype = "dashed", fill = "steelblue") + 
  geom_line(aes(x=order, y= gamma.mean), linewidth = 1.2, colour = "brown") + 
  geom_ribbon(aes(x = order, ymin = gamma.q1, ymax = gamma.q3),
              alpha = 0.2, linetype = "dashed", fill = "brown") + 
  labs(x = "Order Statistics", y = expression(alpha)) +
  theme_minimal(base_size = 30) +
  theme(text = element_text(size = 30), 
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")

data.hill.xi <- data.frame("k" = c(1:n),
                        "u" = sort(uu.xi),
                        "l" = sort(ll.xi),
                        "alpha" = sort(y.xi),
                        "order" = xx,
                        "evgam.1" = sort(predict(evgam.fit.1, type="response")$shape),
                        "evgam.scale" = sort(predict(evgam.fit.scale, type="response")$shape),
                        "blast.mean" = sort(1/alpha.samples[,1]),
                        "blast.q1" = sort(1/alpha.samples[,4]),
                        "blast.q3" = sort(1/alpha.samples[,6]),
                        "gamma.mean" = sort(gamma.xi),
                        "gamma.q1" = sort(gamma.q1.xi), 
                        "gamma.q3" = sort(gamma.q3.xi))
ggplot(data = data.hill.xi) + 
  geom_ribbon(aes(x = order, ymin = l, ymax = u),
              alpha = 0.2, linetype = "dashed", fill = "black") + 
  geom_line(aes(x = order, y = alpha), linewidth = 1.2, colour = "black") +
  geom_line(aes(x = order, y = blast.mean), linewidth = 1.2, colour = "steelblue") + 
  geom_ribbon(aes(x = order, ymin = blast.q1, ymax = blast.q3),
              alpha = 0.2, linetype = "dashed", fill = "steelblue") + 
  geom_line(aes(x=order, y= evgam.1), linewidth = 1.2, linetype=3, colour = "purple") + 
  geom_line(aes(x=order, y= gamma.mean), linewidth = 1.2, colour = "brown") + 
  geom_ribbon(aes(x = order, ymin = gamma.q1, ymax = gamma.q3),
              alpha = 0.2, linetype = "dashed", fill = "brown") + 
  labs(x = "Order Statistics", y = expression(xi)) + ylim(0, 32) +
  theme_minimal(base_size = 30) +
  theme(text = element_text(size = 30), 
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")




hill.est <- function(data, k) {
  n <- length(data)
  if (k <= 0 | k >= n) stop("k must be between 1 and n-1")
  sorted.data <- sort(data, decreasing = TRUE)
  
  X.k <- sorted.data[k]
  X.top.k <- sorted.data[1:k]
  H.k <- (1/k) * sum(log(X.top.k) - log(X.k))
  
  return(H.k)
}


bayes.hill.est <- function(y, k) {
  # Inverse gamma = 1 / Gamma, so sample from Gamma and invert
  # 1 / rgamma(length(y), shape = k, rate = k * hill.est(y, k))
  k * hill.est(y,k) / (k-1)
}

bayes.hill.est(y.origin, 100)

