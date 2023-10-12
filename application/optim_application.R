sm <- "// Stan model for simple linear regression
data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> newp; 
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    vector[n] y; // extreme response
    real <lower=0> atau;
}

parameters {
    vector[newp] theta; // linear predictor
    vector[psi] gamma[p]; // splines coefficient
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    real sigma; //
    vector[p] tau;
}

transformed parameters {
    vector[n] alpha; // tail index
    matrix[n, p] gsmooth; // nonlinear component
    for (j in 1:p){
        gsmooth[,j] <- bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
    };
    for (i in 1:n){
        alpha[i] <- exp(theta[1] + dot_product(bsLinear[i], theta[2:newp]) + (gsmooth[i,] *rep_vector(1, p)));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += gamma_lpdf(lambda1 | 1, 5);
    target += gamma_lpdf(lambda2 | 0.1, 0.1);
    target += inv_gamma_lpdf(sigma | 0.01, 0.01);
    target += double_exponential_lpdf(theta[1] | 0, lambda1); // target += normal_lpdf(theta[1] | 0, 0.1);
    for (j in 1:p){
        target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1);
        target += gamma_lpdf(tau[j] | atau, (square(lambda2)/2));
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * tau[j] * sigma);
    }
}
generated quantities {
    // Used in Posterior predictive check
    vector[n] log_lik;
    real y_rep[n] = pareto_rng(u, alpha);
    for (i in 1:n) {
        log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}
"

op1 <- optimizing(sm, data = list(y = as.vector(y), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), newp = (p+1),
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear),
              init = list(gamma = array(rep(0,(psi*p)), dim=c(psi, p)),
                        theta = rep(0, (p+1)), 
                        tau = rep(1, p), sigma = 1, 
                        lambda1 = 1, lambda2 = 1),
              verbose = TRUE)

theta.map <- op1$par[1:(p+1)]
gamma.map <- op1$par[(p+1+1):(p+1+(psi*p))]
lambda.map <- op1$par[(p+2+(psi*p)):(p+3+(psi*p))]
systime <- Sys.time()
Sys.time()
systime <- chartr(":","-",systime)
date <- gsub("-","", substr(systime, 1, 10))
time <- substr(systime, 12, 20)

df.theta <- data.frame("seq" = seq(1, (p+1)),
                  theta.map)
df.theta$covariate <- factor(c("\u03b8",colnames(fwi.scaled)), levels = c("\u03b8",colnames(fwi.scaled)))
df.theta$labels <- factor(c("theta0",colnames(fwi.scaled)))

ggplot(df.theta, aes(x = covariate)) + ylab("") + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) +
  geom_point(aes(y = theta.map, color = covariate), size = 5) + 
  scale_x_discrete(labels = c(expression(bold(theta[0])),
                              expression(bold(theta[1])),
                              expression(bold(theta[2])),
                              expression(bold(theta[3])),
                              expression(bold(theta[4])),
                              expression(bold(theta[5])),
                              expression(bold(theta[6])),
                              expression(bold(theta[7])))) + 
  scale_color_discrete(labels = c(expression(theta[0]),colnames(fwi.scaled))) + 
  theme_minimal(base_size = 30) + xlab('') + #ylim(-0.01, 0.01) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.text = element_text(size=25),
          legend.margin=margin(0,0,0,-10),
          legend.box.margin=margin(-10,0,-10,0),
          plot.margin = margin(0,0,0,-20),
          axis.text.x = element_text(hjust=0.35),
          axis.text = element_text(size = 28))
# ggsave(paste0("./BRSTIR/application/figures/",date,"_map_theta.pdf"), width=10, height = 7.78)
df <- data.frame("seq" = seq(1, (psi*p)), 
                  gamma.map)
df$covariate <- factor(rep(colnames(fwi.scaled), each = psi, length.out = nrow(df)), levels = colnames(fwi.scaled))
df$labels <- factor(1:(psi*p))

ggplot(df, aes(x =labels , y = gamma.map, color = covariate)) + 
  geom_point(size = 4) + ylab("") + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("")+
  scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+10), 
                    label = c(expression(bold(gamma[1])), 
                              expression(bold(gamma[2])), 
                              expression(bold(gamma[3])), 
                              expression(bold(gamma[4])), 
                              expression(bold(gamma[5])), 
                              expression(bold(gamma[6])), 
                              expression(bold(gamma[7]))), 
                    expand=c(0,3)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size=25),
          legend.margin=margin(0,0,0,-10),
          legend.box.margin=margin(-10,0,-10,0),
          plot.margin = margin(0,0,0,-20),
          axis.text.x = element_text(hjust=0.5),
          axis.text = element_text(size = 28))
# ggsave(paste0("./BRSTIR/application/figures/",date,"_map_gamma.pdf"), width=10, height = 7.78)


f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
alpha.new <- cdf <- newalpha <- NULL
for (j in 1:p){
  f.linear.new[,j] <- bs.linear[,j] * theta.map[j+1]
  f.nonlinear.new[,j] <- bs.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.new[1:n, j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}

for(i in 1:n){
  newalpha[i] <- exp(theta.map[1] + sum(f.new[i,]))
}


f.nonlinear.new <- f.linear.new <- f.new <- matrix(, nrow = n, ncol=p)
for (j in 1:p){
  f.linear.new[,j] <- xholder.linear[,j] * theta.map[j+1]
  f.nonlinear.new[,j] <- xholder.nonlinear[, (((j-1)*psi)+1):(((j-1)*psi)+psi)] %*% matrix(gamma.map, ncol = p)[, j]
  f.new[1:n, j] <- f.linear.new[,j] + f.nonlinear.new[,j]
}
for(i in 1:n){
  alpha.new[i] <- exp(theta.map[1] + sum(f.new[i,]))
}

alpha.optim <- op1$par[(158+1):(158+n)]

func.linear.new <- func.nonlinear.new <- func.new <- matrix(, nrow=n, ncol=0)
for (j in 1:p){
  func.new <- cbind(func.new, f.new[,j])
  func.linear.new <- cbind(func.linear.new, f.linear.new[,j])
  func.nonlinear.new <- cbind(func.nonlinear.new, f.nonlinear.new[,j])  
}

covariates <- gl(p, n, (p*n), labels = factor(colnames(fwi.scaled)))
replicate <- gl(2, n, (p*n))
func.df <- data.frame(seq = seq(0,1,length.out = n),
                        # x = as.vector(apply(fwi.scaled, 2, sort, method = "quick")),
                        new=as.vector(func.new),
                        new.linear=as.vector(func.linear.new),
                        new.nonlinear=as.vector(func.nonlinear.new),
                        covariates=covariates, 
                        replicate=replicate,
                        x = rep(newx, p))


equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Smooth Functions") +
  geom_line(aes(y=new, colour = covariates), linewidth=2) + ylab ("") +
  facet_grid(covariates ~ ., scales="free_y") +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + theme_minimal(base_size = 30) + 
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "none",
        plot.margin = margin(0,0,0,-10),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=33),
        axis.title.x = element_text(size = 35))
# ggsave(paste0("./BRSTIR/application/figures/",date,"_map_smooth.pdf"), width=10.5, height = 15)
ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) +  
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Linear Component") + 
  geom_line(aes(y=new.linear, colour = covariates), linewidth=2) + ylab ("") +
  facet_grid(covariates ~ ., scales="free_y") +
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
# ggsave(paste0("./BRSTIR/application/figures/",date,"_map_linear.pdf"), width=11, height = 15)
ggplot(func.df, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + xlab("Nonlinear Component") +
  geom_line(aes(y=new.nonlinear, colour = covariates), linewidth=2) + ylab ("") +
  facet_grid(covariates ~ ., scales="free_y") +
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
# ggsave(paste0("./BRSTIR/application/figures/",date,"_map_nonlinear.pdf"), width=12.5, height = 15)

data.scenario <- data.frame("x" = c(1:n),
                            "constant" = newx,
                            "mapAlp" = sort(newalpha),
                            "optimAlp" = sort(alpha.optim))
ggplot(data = data.scenario, aes(x = constant)) + 
  ylab(expression(alpha(x))) + xlab("") +
  geom_line(aes(y = mapAlp, col = "MAP Alpha"), linewidth = 2.5) +
  geom_line(aes(y = optimAlp, col = "Optim Alpha"),linetype=2, linewidth = 2.5) +
  labs(col = "") +
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
    theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
    scale_color_manual(values = c("red", "green"))+
    theme(text = element_text(size = 15),
            legend.position="bottom", legend.key.size = unit(1, 'cm'),
            axis.text = element_text(size = 20),
            legend.margin=margin(-15,-15,-15,-15),
            legend.box.margin=margin(-25,0,20,0))
# ggsave(paste0("./BRSTIR/application/figures/",date,"_map_nonlinear.pdf"), width=12.5, height = 15)