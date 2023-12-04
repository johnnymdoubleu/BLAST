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
library(rstan)
library(ggmcmc)
library(MCMCvis)
library(cmdstanr)
library(ggh4x)



# Scenario B

total.iter <- 2

n <- 5000
psi <- 20
threshold <- 0.9
p <- 5
no.theta <- 1
simul.no <- 50

xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)

# Function to generate Gaussian copula
C <- matrix(c(1, 0.3, 0.5, 0.3, 0.3,
            0.3, 1, 0.95, 0.4, 0.4,
            0.5, 0.95, 1, 0.5, 0.1,
            0.3, 0.4, 0.5 , 1, 0.5,
            0.3, 0.4, 0.5, 0.5, 1), nrow = p)
# C <- diag(p)                
## Generate sample


write("// Stan model for simple linear regression

data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp; 
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector[n] y; // extreme response
    real <lower=0> atau;
}

parameters {
    vector[newp] theta; // linear predictor
    array[p] vector[psi] gamma; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real <lower=0> sigma; //
    array[p] real <lower=0> tau;
}

transformed parameters {
    array[n] real <lower=0> alpha; // tail index
    matrix[n, p] gnl; // nonlinear component
    matrix[n, p] gl; // linear component
    matrix[n, p] gsmooth; // linear component
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgnl; // nonlinear component
    matrix[n, p] newgl; // linear component
    matrix[n, p] newgsmooth; // linear component
    for (j in 1:p){
        gnl[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        gl[,j] = bsLinear[,j] * theta[j+1];
        newgl[,j] = xholderLinear[,j] * theta[j+1];
        gsmooth[,j] = gl[,j] + gnl[,j];
        newgsmooth[,j] = newgl[,j] + newgnl[,j];
    };
    for (i in 1:n){
        alpha[i] = exp(theta[1] + sum(gnl[i,]) + sum(gl[i,]));
        newalpha[i] = exp(theta[1] + sum(newgnl[i,]) + sum(newgl[i,]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 1, 10);
    target += gamma_lpdf(lambda2 | 0.1, 100);
    target += normal_lpdf(theta[1] | 0, 100);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01); // target += double_exponential_lpdf(theta[1] | 0, lambda1)
    target += (newp * log(lambda1) + (p * psi * log(lambda2)));
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, lambda2/sqrt(2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * sqrt(tau[j]) * sqrt(sigma));
    }
}
"
, "model_simulation_sc2.stan")

theta.container <- as.data.frame(matrix(, nrow = newp, ncol= total.iter))
gamma.container <- as.data.frame(matrix(, nrow = (p*psi), ncol = total.iter))
alpha.container <- as.data.frame(matrix(, nrow = n, ncol = total.iter))

for(iter in 1:total.iter){
    n <- 5000
    x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
# x.origin <- scale(x.origin)

    for(i in 1:p){
        knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
        tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
        # bs.linear <- cbind(bs.linear, tps[,1:no.theta])
        bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])  
    }

    gamma.origin <- matrix(, nrow = psi, ncol = p)
    for(j in 1:p){
        for (ps in 1:psi){
            if(j %in% c(1,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
            else if(j==7){
                if(ps <= (psi/2)){gamma.origin[ps, j] <- -0.5}
                else{gamma.origin[ps, j] <- -0.5}
            }
            else {
                if(ps <= (psi/2)){gamma.origin[ps, j] <- 1.7}
                else{gamma.origin[ps, j] <- 1.7}
            }
        }
    }
    theta.origin <- -1.6      
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
    x.origin <- x.origin[which(y.origin>u),]
    y.origin <- y.origin[y.origin > u]
    n <- length(y.origin)

    xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
    newx <- seq(0, 1, length.out=n)
    xholder <- bs.x <- matrix(, nrow = n, ncol = p)
    for(i in 1:p){
        xholder[,i] <- seq(0, 1, length.out = n)  
        test.knot <- seq(0, 1, length.out = psi)
        splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
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

    fit1 <- stan(
        file = "model_simulation_sc2.stan",  # Stan program
        data = data.stan,    # named list of data
        init = init.alpha,      # initial value
        chains = 3,             # number of Markov chains
        warmup = 1000,          # number of warmup iterations per chain
        iter = 2000,            # total number of iterations per chain
        cores = 4,              # number of cores (could use one per chain)
        refresh = 500             # no progress shown
    )

    posterior <- extract(fit1)

    beta.samples <- summary(fit1, par=c("beta"), probs = c(0.05,0.5, 0.95))$summary
    lambda.samples <- summary(fit1, par=c("lambda1"), probs = c(0.05,0.5, 0.95))$summary
    alpha.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
    newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary

    alpha.container[,iter] <- sort(newalpha.samples[,1])
    beta.container[,iter] <- beta.samples[,1]
}

alpha.container$x <- seq(0,1, length.out = n)
alpha.container$true <- alpha.new
alpha.container <- cbind(alpha.container, t(apply(alpha.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
colnames(alpha.container)[(dim(alpha.container)[2]-2):(dim(alpha.container)[2])] <- c("q1","q2","q3")
alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
alpha.container <- as.data.frame(alpha.container)
# for(i in 1:total.iter){
#   # print(.data[[names(data.scenario)[i]]])
#   plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.4,linewidth = 0.7)
# }

ggplot(data = alpha.container, aes(x = x)) + 
    ylab(expression(alpha(x))) + xlab(expression(x)) + 
    geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.5) + 
    geom_line(aes(y=true, col = "True"), linewidth = 1.8) + 
    geom_line(aes(y=mean, col = "Mean"), linewidth = 1.8, linetype = 2) +
    geom_line(aes(y=q2, col = "Median"), linewidth = 1.8, linetype = 2) +
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
    labs(col = "") + #ylim(0, 150) +
    scale_color_manual(values = c("blue","#e0b430","red"))+
    theme_minimal(base_size = 30) +
    theme(plot.title = element_text(hjust = 0.5, size = 30),
            legend.position="top", 
            legend.key.size = unit(1, 'cm'),
            legend.text = element_text(size=18),
            plot.margin = margin(0,0,0,-10),
            strip.text = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size=33),
            axis.title.x = element_text(size = 35))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_alpha_sc1-wi.pdf"), width=10, height = 7.78)


resg <- gather(beta.container,
               key = "group",
               names(beta.container),
               value = "values")
resg$group1 <- factor(rep(1:newp, total.iter))

somelines <- data.frame(value=c(as.vector(beta.origin)),boxplot.nr=c(1:(newp)))
ggplot(resg, aes(group=group1, x = group1, y = values, fill=group1)) + ylim(-0.5,1) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) +
  geom_boxplot() + #coord_cartesian(ylim=c(-1,1))+
  geom_segment(data=somelines,aes(x=boxplot.nr-0.5, xend=boxplot.nr+0.5, 
                                  y=value,yend=value),inherit.aes=FALSE,color="red",linewidth=1.5)+
  labs(x = "", y = "") + 
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
                              expression(bold(theta[10])))) + 
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