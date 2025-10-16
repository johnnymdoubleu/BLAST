library(npreg)
library(Pareto)
suppressMessages(library(tidyverse))
library(readxl)
library(gridExtra)
library(corrplot)
library(rstan)
library(loo)
library(qqboxplot)
library(ggdensity)
library(ggforce)
library(ggdist)
options(mc.cores = parallel::detectCores())

# Structure of the FWI System
#DSR : Daily Severity Rating
#FWI : Fire Weather Index
#BUI : Buildup Index
#ISI : Initial Spread Index
#FFMC : Fine FUel Moisture Code
#DMC : Duff Moisture Code
#DC : Drought Code


setwd("C:/Users/Johnny Lee/Documents/GitHub")
df <- read_excel("./BLAST/application/AADiarioAnual.xlsx", col_types = c("date", rep("numeric",40)))
df.long <- gather(df, condition, measurement, "1980":"2019", factor_key=TRUE)
missing.values <- which(!is.na(df.long$measurement))

#NAs on Feb 29 most years, and Feb 14, 1999
#considering the case of leap year, the missing values are the 29th of Feb
#Thus, each year consist of 366 data with either 1 or 0 missing value.
Y <- df.long$measurement[!is.na(df.long$measurement)]
Y[Y==0] <- 1e-5
psi <- 30

multiplesheets <- function(fname) {
    setwd("C:/Users/Johnny Lee/Documents/GitHub")
    # getting info about all excel sheets
    sheets <- excel_sheets(fname)
    tibble <- lapply(sheets, function(x) read_excel(fname, sheet = x, col_types = c("date", rep("numeric", 41))))
    data_frame <- lapply(tibble, as.data.frame)
    # assigning names to data frames
    names(data_frame) <- sheets
    return(data_frame)
}
setwd("C:/Users/Johnny Lee/Documents/GitHub")
path <- "./BLAST/application/DadosDiariosPT_FWI.xlsx"
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
for(i in 1:length(cov)){
    cov.long <- gather(cov[[i]][,1:41], condition, measurement, "1980":"2019", factor_key=TRUE)
    fwi.scaled[,i] <- fwi.index[,i] <- cov.long$measurement[missing.values]
}

fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)
fwi.origin <- fwi.scaled

fwi.origin <- data.frame(fwi.origin, BA=Y)
BA.shifted <- ifelse(fwi.origin$BA == 0, 1e-5, fwi.origin$BA)
fwi.origin$log.BA <- log(BA.shifted)
quant.fit <- quantreg::rq(log.BA ~ DSR + FWI + BUI + ISI + FFMC + DMC + DC, 
                          tau=rep(0.975, length(Y)), data = fwi.origin)#,
                          # method = "fnc", # feasible non-crossing
                          # R = matrix(c(1, rep(0, p)), 1, p+1), # constraint matrix
                          # r = 1) # constraint threshold)
qu <- exp(predict(quant.fit))
y <- Y[which(Y>qu)]
u <- qu[which(Y>qu)]
# u[u < 0] <- 0.001

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled <- as.data.frame(sapply(fwi.scaled[which(Y>qu),], FUN = range01))
head(as.data.frame(lapply(fwi.scaled[which(Y>qu),], FUN = scale)))
n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]

save(fwi.scaled, fwi.origin, qu, y, u, file = "./BLAST/application/wildfire_prep.Rdata")

# qu975 <- qu[which(Y>u)]
# qu975[which(qu975<u)] <- u
# load("./BLAST/application/qrresult.Rdata")
# qu <- predict(quantile.mod)
# qu975 <- qu[which(Y>u),1]

no.theta <- 1 #represents the no. of linear predictors for each smooth functions
newx <- seq(0, 1, length.out=n)
xholder.linear <- xholder.nonlinear <- bs.linear <- bs.nonlinear <- matrix(,nrow=n, ncol=0)
xholder <- matrix(nrow=n, ncol=p)
end.holder <- basis.holder <- matrix(, nrow = 2, ncol = 0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                      matrix(c(which.min(fwi.scaled[,i]),
                              which.max(fwi.scaled[,i])), ncol=2))
}

for(i in 1:p){
  xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
  test.knot <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  knots <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  tps <- basis.tps(fwi.scaled[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  basis.holder <- cbind(basis.holder, 
          solve(matrix(c(tps[index.holder[i,1], no.theta+1],
                  tps[index.holder[i,1], no.theta+psi],
                  tps[index.holder[i,2], no.theta+1],
                  tps[index.holder[i,2], no.theta+psi]), 
                  nrow = 2, ncol = 2)))
  end.holder <- cbind(end.holder, 
                matrix(c(tps[index.holder[i,1], no.theta+1],
                  tps[index.holder[i,1], no.theta+psi],
                  tps[index.holder[i,2], no.theta+1],
                  tps[index.holder[i,2], no.theta+psi]), 
                  nrow = 2, ncol = 2))
  # mid.holder <- cbind(mid.holder,)                  
  bs.linear <- cbind(bs.linear, tps[,1:no.theta])
  bs.nonlinear <- cbind(bs.nonlinear, tps[,-c(1:no.theta)])
}


model.stan <- "data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    vector[n] u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector[n] y; // extreme response
    real <lower=0> atau;
    matrix[2, (2*p)] basisFL;
    array[(p*2)] int indexFL;
}
parameters {
    vector[(p+1)] theta; // linear predictor
    vector[(psi-2)] gammaTemp[p]; // constraint splines coefficient from 2 to psi-1
    real <lower=0> lambda1; // lasso penalty
    real <lower=0> lambda2; // group lasso penalty
    array[p] real <lower=0> tau1;
    array[p] real <lower=0> tau2;
}
transformed parameters {
    array[n] real <lower=0> alpha; // covariate-adjusted tail index
    vector[psi] gamma[p]; // splines coefficient 
    real <lower=0> lambda2o;

    lambda2o=lambda2*100;
    {
      vector[2] gammaFL[p];
      matrix[n, p] gsmooth; // linear component
      
      for(j in 1:p){
          gamma[j][2:(psi-1)] = gammaTemp[j][1:(psi-2)]*100;
          gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * (bsNonlinear[indexFL[(((j-1)*2)+1):(((j-1)*2)+2)], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j]*100) * (-1);
          gamma[j][1] = gammaFL[j][1];
          gamma[j][psi] = gammaFL[j][2];
      };
      for (j in 1:p){
          gsmooth[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j] + bsLinear[,j] * theta[j+1];
      };

      for (i in 1:n){
          alpha[i] = exp(theta[1] + sum(gsmooth[i,]));
      };
    }

}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u[i], alpha[i]);
    }
    target += normal_lpdf(theta[1] | 0, 100);
    target += gamma_lpdf(lambda1 | 1, 1e-3);
    target += gamma_lpdf(lambda2o | 1, 1e-3);
    target += (2*p*log(lambda2o));
    for (j in 1:p){
        target += gamma_lpdf(tau1[j] | 1, lambda1^2*0.5);
        target += normal_lpdf(theta[(j+1)] | 0, sqrt(1/tau1[j]));
        target += gamma_lpdf(tau2[j] | atau, lambda2o^2*0.5);
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * (1/tau2[j]));
    }
}
generated quantities {
    // Used in Posterior predictive check
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgsmooth; // linear component

    for (j in 1:p){
        newgsmooth[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j] + xholderLinear[,j] * theta[j+1];
    };    

    for (i in 1:n){ 
        newalpha[i] = exp(theta[1] + sum(newgsmooth[i,]));
    };
}
"


data.stan <- list(y = as.vector(y), u = as.vector(u), p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), indexFL = as.vector(t(index.holder)),
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear, basisFL = basis.holder)

init.alpha <- list(list(gammaTemp = array(rep(0, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0, (p+1)), 
                        tau1 = rep(0.1, p),tau2 = rep(0.1, p),
                        lambda1 = 0.1, lambda2 = 0.01),
                   list(gammaTemp = array(rep(0, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0, (p+1)), 
                        tau1 = rep(0.001, p),tau2 = rep(0.001, p),
                        lambda1 = 100, lambda2 = 1),
                   list(gammaTemp = array(rep(0, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0.1, (p+1)), 
                        tau1 = rep(0.5, p),tau2 = rep(0.5, p),
                        lambda1 = 5, lambda2 = 5.5))

# stanc("C:/Users/Johnny Lee/Documents/GitHub/BLAST/application/model1.stan")
fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    chains = 3,             # number of Markov chains
    iter = 4000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 2000           # no progress shown
)

posterior <- rstan::extract(fit1)

theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary
# intcpt.samples <- summary(fit1, par=c("intcpt"), probs = c(0.05, 0.5, 0.95))

g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))

data.scenario <- data.frame("x" = xholder[,1],
                            "post.mean" = (alpha.samples[,1]),
                            "post.median" = (alpha.samples[,5]),
                            "q1" = (alpha.samples[,4]),
                            "q3" = (alpha.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,...,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_alpha.pdf"), width=10, height = 7.78)

data.smooth <- data.frame("x" = as.vector(xholder),
                          "true" = as.vector(as.matrix(fwi.scaled)),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  ylim(-1.2, max(data.smooth$q3[((((i-1)*n)+1):(i*n))])) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=-1.2, color = "red", size = 7)
}

grid.arrange(grobs = grid.plts, ncol = 2, nrow = 4)

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_DC.pdf"), grid.plts[[7]], width=10, height = 7.78)