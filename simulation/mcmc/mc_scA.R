library(npreg)
library(tmvnsim)
library(reshape2)
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
library(cmdstanr)

# Scenario A

total.iter <- 2

n <- 2500
psi <- 10
threshold <- 0.95
p <- 5
newp <- p+1
no.theta <- 1

# Function to generate Gaussian copula
# C <- matrix(c(1, 0.3, 0.5, 0.3, 0.3,
#             0.3, 1, 0.95, 0.4, 0.4,
#             0.5, 0.95, 1, 0.5, 0.1,
#             0.3, 0.4, 0.5 , 1, 0.5,
#             0.3, 0.4, 0.5, 0.5, 1), nrow = p)
C <- diag(p)
## Generate sample
gamma.origin <- matrix(, nrow = psi, ncol = p)
for(j in 1:p){
    for (ps in 1:psi){
        if(j %in% c(1,4,5,6,9,10)){gamma.origin[ps, j] <- 0}
        else if(j==7){
            if(ps <= (psi/2)){gamma.origin[ps, j] <- -0.1}
            else{gamma.origin[ps, j] <- -0.1}
        }
        else {
            if(ps <= (psi/2)){gamma.origin[ps, j] <- -0.1}
            else{gamma.origin[ps, j] <- -0.1}
        }
    }
}

theta.origin <- c(0.5, 0, -0.2, -0.2, 0, 0)

write("// Stan model for BRSTIR Pareto Uncorrelated Samples

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
    target += gamma_lpdf(lambda2 | 0.1, 10);
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
, "model_simulation_sc1.stan")

# theta.container <- as.data.frame(matrix(, nrow = newp, ncol= total.iter))
# gamma.container <- as.data.frame(matrix(, nrow = (p*psi), ncol = total.iter))
newgsmooth.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
newgl.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
newgnl.container <- as.data.frame(matrix(, nrow = (p*(n*(1-threshold))), ncol = total.iter))
alpha.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
alpha.lower.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))
alpha.upper.container <- as.data.frame(matrix(, nrow = (n*(1-threshold)), ncol = total.iter))

for(iter in 1:total.iter){
    n <- 2500
    x.origin <- pnorm(matrix(rnorm(n*p), ncol = p) %*% chol(C))
    xholder.nonlinear <- xholder.linear <- bs.nonlinear <- bs.linear <- matrix(,nrow=n, ncol=0)
    for(i in 1:p){
        knots <- seq(min(x.origin[,i]), max(x.origin[,i]), length.out = psi)
        tps <- basis.tps(x.origin[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
        bs.linear <- cbind(bs.linear, tps[,1:no.theta])
        bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
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

    data.stan <- list(y = as.vector(y.origin), u = u, p = p, n= n, psi = psi, 
                        atau = ((psi+1)/2), newp = (p+1),
                        bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                        xholderLinear = xholder.linear, 
                        xholderNonlinear = xholder.nonlinear)

    init.alpha <- list(list(gamma = array(rep(0, (psi*p)), dim=c(psi, p)),
                            theta = rep(0, (p+1)), 
                            tau = rep(0.1, p), sigma = 0.1, 
                            lambda1 = 0.1, lambda2 = 0.1),
                      list(gamma = array(rep(0.02, (psi*p)), dim=c(psi, p)),
                            theta = rep(0.01, (p+1)), 
                            tau = rep(0.01, p), sigma = 0.001,
                            lambda1 = 0.01, lambda2 = 0.1),
                      list(gamma = array(rep(0.01, (psi*p)), dim=c(psi, p)),
                            theta = rep(0.05, (p+1)), 
                            tau = rep(0.01, p), sigma = 0.01,
                            lambda1 = 0.1, lambda2 = 0.01))

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

    # theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
    # gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
    lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
    newgl.samples <- summary(fit1, par=c("newgl"), probs = c(0.05, 0.5, 0.95))$summary
    newgnl.samples <- summary(fit1, par=c("newgnl"), probs = c(0.05, 0.5, 0.95))$summary
    newgsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
    newalpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary

    alpha.lower.container[,iter] <- (newalpha.samples[,4])
    alpha.container[,iter] <- (newalpha.samples[,5])
    alpha.upper.container[,iter] <- (newalpha.samples[,6])
    # theta.container[,iter] <- theta.samples[,5]
    # gamma.container[,iter] <- gamma.samples[,5]
    newgl.container[,iter] <- as.vector(matrix(newgl.samples[,5], nrow = n, byrow=TRUE))
    newgnl.container[,iter] <- as.vector(matrix(newgnl.samples[,5], nrow = n, byrow=TRUE))
    newgsmooth.container[,iter] <- as.vector(matrix(newgsmooth.samples[,5], nrow = n, byrow=TRUE))
}

alpha.container$x <- seq(0,1, length.out = n)
alpha.container$true <- sort(alp.new)
alpha.container <- cbind(alpha.container, t(apply(alpha.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
colnames(alpha.container)[(dim(alpha.container)[2]-2):(dim(alpha.container)[2])] <- c("q1","q2","q3")
alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
alpha.container$q1 <- apply(alpha.lower.container[,1:total.iter], 1, quantile, c(.5))
alpha.container$q3 <- apply(alpha.upper.container[,1:total.iter], 1, quantile, c(.5))
alpha.container <- as.data.frame(alpha.container)

plt <- ggplot(data = alpha.container, aes(x = x)) + ylab(expression(alpha(c*bold("1")))) + xlab(expression(c)) + labs(col = "")
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.2, linewidth = 0.7)
  }
} else{
  for(i in 50:100){
    plt <- plt + geom_line(aes(y = .data[[names(alpha.container)[i]]]), alpha = 0.2, linewidth = 0.7)
  }
}
print(plt + geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) + ylim(0.5, 2.5) + 
        geom_line(aes(y=true, col = "True"), linewidth = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5, linetype = 2) +
        scale_fill_manual(values=c("steelblue"), name = "") +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                strip.text = element_blank(),
                axis.text = element_text(size = 35)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_alpha_sc1-wi.pdf"), width=10, height = 7.78)


# resg <- gather(theta.container,
#                key = "group",
#                names(theta.container),
#                value = "values")
# resg$group1 <- factor(rep(1:newp, total.iter))

# somelines <- data.frame(value=c(as.vector(theta.origin)),boxplot.nr=c(1:(newp)))
# ggplot(resg, aes(group=group1, x = group1, y = values, fill=group1)) + ylim(-0.5,1) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + geom_boxplot() + 
#   geom_segment(data=somelines,aes(x=boxplot.nr-0.5, xend=boxplot.nr+0.5, 
#                                   y=value,yend=value),inherit.aes=FALSE,color="red",linewidth=1.5)+
#   labs(x = "", y = "") + 
#   scale_x_discrete(labels = c(expression(bold(theta[0])),
#                               expression(bold(theta[1])),
#                               expression(bold(theta[2])),
#                               expression(bold(theta[3])),
#                               expression(bold(theta[4])),
#                               expression(bold(theta[5])),
#                               expression(bold(theta[6])),
#                               expression(bold(theta[7])),
#                               expression(bold(theta[8])),
#                               expression(bold(theta[9])),
#                               expression(bold(theta[10])))) + 
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#           legend.text.align = 0,
#           legend.title = element_blank(),
#           legend.text = element_text(size=25),
#           legend.margin=margin(0,0,0,-10),
#           legend.box.margin=margin(-10,0,-10,0),
#           plot.margin = margin(0,0,0,-20),
#           axis.text.x = element_text(hjust=0.35),
#           axis.text = element_text(size = 28))
# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_theta_sc2-wi.pdf"), width=10, height = 7.78)

# resg <- gather(gamma.container,
#                key = "group",
#                names(gamma.container),
#                value = "values")
# resg$group1 <- factor(rep(1:(psi*p), total.iter))
# resg$group2 <- factor(rep(1:p, each = psi))

# somelines <- data.frame(value=c(as.vector(gamma.origin)),
#                         boxplot.nr=c(1:(psi*p)),
#                         covariate = factor(rep(1:p, each= psi)))
# ggplot(resg, aes(group=group1, x = group1, y = values, fill=group2)) + ylim(-1.1,1.1) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + labs(x = "", y = "") + 
#   geom_boxplot() + #coord_cartesian(ylim=c(-1,1))+
#   geom_segment(data=somelines,aes(x=boxplot.nr-0.5,xend=boxplot.nr+0.5, 
#                                   y=value,yend=value),inherit.aes=FALSE,
#                                   color="red",
#                                   linewidth=1.5)+
#   scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+7), 
#                     label = c(expression(bold(gamma[1])), 
#                               expression(bold(gamma[2])), 
#                               expression(bold(gamma[3])), 
#                               expression(bold(gamma[4])), 
#                               expression(bold(gamma[5])), 
#                               expression(bold(gamma[6])), 
#                               expression(bold(gamma[7])), 
#                               expression(bold(gamma[8])), 
#                               expression(bold(gamma[9])), 
#                               expression(bold(gamma[10]))),
#                     expand=c(0,3)) +
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 20),
#           legend.title = element_blank(),
#           legend.text = element_text(size=25),
#           legend.margin=margin(0,0,0,-10),
#           legend.box.margin=margin(-10,0,-10,0),
#           plot.margin = margin(0,0,0,-20),
#           axis.text.x = element_text(hjust=0.5),
#           axis.text = element_text(size = 28))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_gamma_sc2-wi.pdf"), width=10, height = 7.78)

equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

newgsmooth.container$x <- seq(0,1, length.out = n)
newgsmooth.container$true <- as.vector(f.new)
newgsmooth.container <- cbind(newgsmooth.container, t(apply(newgsmooth.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
colnames(newgsmooth.container)[(dim(newgsmooth.container)[2]-2):(dim(newgsmooth.container)[2])] <- c("q1","q2","q3")
newgsmooth.container$mean <- rowMeans(newgsmooth.container[,1:total.iter])
newgsmooth.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]", "g[6]"))
newgsmooth.container <- as.data.frame(newgsmooth.container)

plt <- ggplot(data = newgsmooth.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(x))
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(newgsmooth.container)[i]]]), alpha = 0.2, linewidth = 0.7)
    # plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
  }
} else{
  for(i in 1:50){
    plt <- plt + geom_line(aes(y = .data[[names(newgsmooth.container)[i]]]), alpha = 0.2, linewidth = 0.7)
    # plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
  }
}
print(plt + #geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
        geom_line(aes(y=true, col = "True"), linewidth = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5, linetype = 2) + 
        ylim(-0.23, 0.2) +
        # facet_wrap(covariate ~ ., scales = "free_x", nrow = 5,
        #             labeller = label_parsed, strip.position = "left") +
        facet_grid(covariate ~ ., scales = "free_x", switch = "y",
                    labeller = label_parsed) +        
        # scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
        #scale_fill_manual(values=c("steelblue"), name = "") +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                plot.margin = margin(0,0,0,-20),
                strip.text = element_blank(),
                axis.text = element_text(size = 20)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_smooth_sc1-wi.pdf"), width=10, height = 7.78)                


newgl.container$x <- seq(0,1, length.out = n)
newgl.container$true <- as.vector(f.linear.new)
newgl.container <- cbind(newgl.container, t(apply(newgl.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
colnames(newgl.container)[(dim(newgl.container)[2]-2):(dim(newgl.container)[2])] <- c("q1","q2","q3")
newgl.container$mean <- rowMeans(newgl.container[,1:total.iter])
newgl.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]", "g[6]"))
newgl.container <- as.data.frame(newgl.container)

plt <- ggplot(data = newgl.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(x))
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(newgl.container)[i]]]), alpha = 0.2, linewidth = 0.7)
    # plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
  }
}else{
  for(i in 1:50){
    plt <- plt + geom_line(aes(y = .data[[names(newgl.container)[i]]]), alpha = 0.2, linewidth = 0.7)
    # plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
  }
}
print(plt + #geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
        geom_line(aes(y=true, col = "True"), linewidth = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5, linetype = 2) + 
        ylim(-0.23, 0.2) +
        # facet_wrap(covariate ~ ., scales = "free_x", nrow = 5,
        #             labeller = label_parsed, strip.position = "left") +
        facet_grid(covariate ~ ., scales = "free_x", switch = "y",
                    labeller = label_parsed) +  
        # scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
        #scale_fill_manual(values=c("steelblue"), name = "") +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                plot.margin = margin(0,0,0,-20),
                strip.text = element_blank(),
                axis.text = element_text(size = 20)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_linear_sc1-wi.pdf"), width=10, height = 7.78)                

newgnl.container$x <- seq(0,1, length.out = n)
newgnl.container$true <- as.vector(f.nonlinear.new)
newgnl.container <- cbind(newgnl.container, t(apply(newgnl.container[,1:total.iter], 1, quantile, c(0.05, .5, .95))))
colnames(newgnl.container)[(dim(newgnl.container)[2]-2):(dim(newgnl.container)[2])] <- c("q1","q2","q3")
newgnl.container$mean <- rowMeans(newgnl.container[,1:total.iter])
newgnl.container$covariate <- gl(p, n, (p*n), labels = c("g[1]", "g[2]", "g[3]", "g[4]", "g[5]", "g[6]"))
newgnl.container <- as.data.frame(newgnl.container)

plt <- ggplot(data = newgnl.container, aes(x = x, group = covariate)) + ylab("") + xlab(expression(x))
if(total.iter <= 50){
  for(i in 1:total.iter){
    plt <- plt + geom_line(aes(y = .data[[names(newgnl.container)[i]]]), alpha = 0.2, linewidth = 0.7)
    # plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
  }
}else{
  for(i in 1:50){
    plt <- plt + geom_line(aes(y = .data[[names(newgnl.container)[i]]]), alpha = 0.2, linewidth = 0.7)
    # plt <- plt + geom_line(aes(y = .data[[names(data.scenario)[i]]]))
  }
}
print(plt + #geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
        geom_line(aes(y=true, col = "True"), linewidth = 2) + 
        geom_line(aes(y=mean, col = "Mean"), linewidth = 1.5, linetype = 2) + 
        ylim(-0.23, 0.2) +
        # facet_wrap(covariate ~ ., scales = "free_x", nrow = 5,
        #             labeller = label_parsed, strip.position = "left") +
        facet_grid(covariate ~ ., scales = "free_x", switch = "y",
                    labeller = label_parsed) +  
        # scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
        #scale_fill_manual(values=c("steelblue"), name = "") +
        scale_color_manual(values = c("steelblue", "red"))+
        guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
        theme_minimal(base_size = 30) +
        theme(legend.position = "none",
                plot.margin = margin(0,0,0,-20),
                strip.text = element_blank(),
                axis.text = element_text(size = 20)))

# ggsave(paste0("./simulation/results/",Sys.Date(),"_",total.iter,"_MC_nonlinear_sc1-wi.pdf"), width=10, height = 7.78)
