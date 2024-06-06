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
#DSR : Dail Severity Rating
#FWI : Fire Weather Index
#BUI : Buildup Index
#ISI : Initial Spread Index
#FFMC : Fine FUel Moisture Code
#DMC : Duff Moisture Code
#DC : Drought Code


setwd("C:/Users/Johnny Lee/Documents/GitHub")
df <- read_excel("./BRSTIR/application/AADiarioAnual.xlsx", col_types = c("date", rep("numeric",40)))
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
path <- "./BRSTIR/application/DadosDiariosPT_FWI.xlsx"
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
    fwi.index[,i] <- cov.long$measurement[missing.values]
    fwi.scaled[,i] <- cov.long$measurement[missing.values]
}

fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)

psi <- 10
threshold <- 0.975
Y.temp <- Y[1:13801]
u <- quantile(Y, threshold)
y <- Y.temp[Y.temp>u]
fwi.scaled <- fwi.scaled[1:13801,]
fwi.origin <- fwi.scaled <-fwi.scaled[which(Y.temp>u),]
fwi.minmax <- fwi.index[which(Y>quantile(Y,0.975)),]
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled$DSR <- (fwi.scaled$DSR - min(fwi.minmax$DSR))/(max(fwi.minmax$DSR)-min(fwi.minmax$DSR))
fwi.scaled$FWI <- (fwi.scaled$FWI - min(fwi.minmax$FWI))/(max(fwi.minmax$FWI)-min(fwi.minmax$FWI))
fwi.scaled$BUI <- (fwi.scaled$BUI - min(fwi.minmax$BUI))/(max(fwi.minmax$BUI)-min(fwi.minmax$BUI))
fwi.scaled$ISI <- (fwi.scaled$ISI - min(fwi.minmax$ISI))/(max(fwi.minmax$ISI)-min(fwi.minmax$ISI))
fwi.scaled$FFMC <- (fwi.scaled$FFMC - min(fwi.minmax$FFMC))/(max(fwi.minmax$FFMC)-min(fwi.minmax$FFMC))
fwi.scaled$DMC <- (fwi.scaled$DMC - min(fwi.minmax$DMC))/(max(fwi.minmax$DMC)-min(fwi.minmax$DMC))
fwi.scaled$DC <- (fwi.scaled$DC - min(fwi.minmax$DC))/(max(fwi.minmax$DC)-min(fwi.minmax$DC))
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# fwi.scaled <- as.data.frame(sapply(fwi.origin, FUN = range01))
n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]

fwi.origin <- data.frame(fwi.index[which(Y.temp>u),], BA=y)
max.fwi <- fwi.origin[which.max(y),]

ggplot(fwi.origin, aes(x=DSR, y=FFMC)) + 
  geom_point(aes(colour = BA), size= 2.5) + 
  scale_colour_stepsn(colours = c("slategray1", "red"), labels=function(x) format(x, big.mark = ",", scientific = TRUE)) +
  # scale_colour_stepsn(colours = heat.colors(2, rev=TRUE), labels=function(x) format(x, big.mark = ",", scientific = TRUE)) +
  # guides(colour = guide_coloursteps(show.limits = TRUE)) +
  # scale_color_gradientn(colours = heat.colors(2)) +
  geom_density2d(aes(x=DSR, y=FFMC), colour="steelblue", linewidth = 1.3) + 
  # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="steelblue")+ 
  geom_mark_circle(aes(x = max.fwi$DSR, y = max.fwi$FFMC, label = "15th Oct 2017"), con.type = "straight",
                   radius = unit(2.5, "mm"), color = "steelblue", size = 1, 
                   con.colour = "steelblue", con.cap = unit(0, "mm"),
                   label.colour = "steelblue", label.buffer = unit(5, "mm"),
                   label.fill = "transparent")  +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        # plot.margin = margin(0,0,0,-1),
        strip.text = element_blank(),
        axis.title = element_text(size = 30))
# ggsave("./BRSTIR/application/figures/extremeviz.pdf", width = 10, height = 7.78)
# ------------- Explanatory Analaysis
# first.extreme <- which(Y==max(y))
# second.extreme <- which(Y==max(y[-which.max(y)]))
# tenth.extreme <- which(Y==sort(y, decreasing = TRUE)[10])
# ggplot(fwi.index[((first.extreme):(first.extreme+12)),], aes(x=date)) +
#   geom_line(aes(y=DSR, color = "DSR"), linetype = 1) + 
#   geom_line(aes(y=FWI, color = "FWI"), linetype = 2) +
#   geom_line(aes(y=BUI, color = "BUI"), linetype = 3) +
#   geom_line(aes(y=ISI, color = "ISI"), linetype = 4) +
#   geom_line(aes(y=FFMC, color = "FFMC"), linetype = 5) + 
#   geom_line(aes(y=DMC, color = "DMC"), linetype = 6) +
#   geom_line(aes(y=DC, color = "DC"), linetype = 7)  + 
#   ylab("indices") + xlab("dates after extreme fire (sorted by burnt area)") + 
#   scale_color_manual(name = "Indices", values = c(
#     "DSR" = "darkblue", 
#     "FWI" = "red",
#     "BUI" = "green",
#     "ISI" = "yellow",
#     "FFMC" = "orange",
#     "DMC" = "purple",
#     "DC" = "skyblue")) +
#   theme(legend.position="right", 
#       legend.key.size = unit(1, 'cm'),
#       legend.text = element_text(size=20),
#       # plot.margin = margin(0,0,0,-1),
#       axis.title = element_text(size = 20))

# fwi.index[((first.extreme):(first.extreme+12)),]
# fwi.index[13682:13694,]
# fwi.index[second.extreme:(second.extreme+12),]

# df.extreme <- cbind(y, fwi.scaled)
# df.extreme <- as.data.frame(cbind(month = fwi.index$month[which(Y>u)], df.extreme))
# ggplot(df.extreme, aes(x=month, y=y, color=month)) + geom_point(size=6) + theme_minimal() +
#     theme(plot.title = element_text(hjust = 0.5, size = 20),
#         legend.title = element_blank(),
#         legend.text = element_text(size=20),
#         # axis.ticks.x = element_blank(),
#         axis.text.x = element_text(hjust=0.35),
#         axis.text = element_text(size = 20),
#         axis.title = element_text(size = 15))
# ggsave("./BRSTIR/application/figures/datavis.pdf", width=15)



no.theta <- 1
newx <- seq(0, 1, length.out=n)
xholder.linear <- xholder.nonlinear <- bs.linear <- bs.nonlinear <- matrix(,nrow=n, ncol=0)
xholder <- matrix(nrow=n, ncol=p)
end.holder <- basis.holder <- matrix(, nrow = 2, ncol =0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                      matrix(c(which.min(fwi.scaled[,i]),
                              which.max(fwi.scaled[,i])), ncol=2))
}
# for(i in 1:n){
#   xholder[i,] <- centre.fwi + seq(min(PC1), max(PC1), length.out = n)[i] %*% phi[,1]
# }

for(i in 1:p){
  xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
  test.knot <- quantile(fwi.scaled[,i], probs=seq(1/(psi+1),psi/(psi+1),length.out = psi))
  # test.knot <- seq(min(xholder[,i]), max(xholder[,i]), length.out = psi)  
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  knots <- quantile(fwi.scaled[,i], probs=seq(1/(psi+1),psi/(psi+1),length.out = psi))
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


write("data {
    int <lower=1> n; // Sample size
    int <lower=1> p; // regression coefficient size
    int <lower=1> psi; // splines coefficient size
    real <lower=0> u; // large threshold value
    matrix[n,p] bsLinear; // fwi dataset
    matrix[n, (psi*p)] bsNonlinear; // thin plate splines basis
    matrix[n,p] xholderLinear; // fwi dataset
    matrix[n, (psi*p)] xholderNonlinear; // thin plate splines basis    
    vector[n] y; // extreme response
    real <lower=0> atau;
    matrix[2, (2*p)] basisFL;
    array[(p*2)] int indexFL;
    vector[n] newy;
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
    vector[2] gammaFL[p];
    real <lower=0> lambda2o;
    matrix[2, p] subgnl;
    matrix[n, p] gnl; // nonlinear component
    matrix[n, p] gl; // linear component
    matrix[n, p] gsmooth; // linear component
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgnl; // nonlinear component
    matrix[n, p] newgl; // linear component
    matrix[n, p] newgsmooth; // linear component

    lambda2o=lambda2*100;
    for(j in 1:p){
        gamma[j][2:(psi-1)] = gammaTemp[j][1:(psi-2)]*100;
        subgnl[,j] = bsNonlinear[indexFL[(((j-1)*2)+1):(((j-1)*2)+2)], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j]*100;
        gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * subgnl[,j] * (-1);
        gamma[j][1] = gammaFL[j][1];
        gamma[j][psi] = gammaFL[j][2];  
    };

    for (j in 1:p){
        gnl[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        newgnl[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        gl[,j] = bsLinear[,j] * theta[j+1];
        newgl[,j] = xholderLinear[,j] * theta[j+1];
        gsmooth[,j] = gl[,j] + gnl[,j];
        newgsmooth[,j] = newgl[,j] + newgnl[,j];
    };

    for (i in 1:n){
        alpha[i] = exp(theta[1] + sum(gsmooth[i,])); 
        newalpha[i] = exp(theta[1] + sum(newgsmooth[i,]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
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
    vector[n] log_lik;
    vector[n] yrep;
    vector[n] f;
    
    for(i in 1:n){
      yrep[i] = pareto_rng(u, alpha[i]);
      f[i] = (alpha[215]/exp(newy[i]))*(exp(newy[i])/u)^(-alpha[215]); //pareto_rng(u, alpha[i])
      log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}
"
, "model_BRSTIR_doomsday.stan")

data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), indexFL = as.vector(t(index.holder)),
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear, basisFL = basis.holder,
                    newy = seq(log(u), 30, length.out = n))

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

# stanc("C:/Users/Johnny Lee/Documents/GitHub/BRSTIR/application/model1.stan")
fit1 <- stan(
    file = "model_BRSTIR_doomsday.stan",  # Stan program
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    # init_r = 1,
    chains = 3,             # number of Markov chains
    # warmup = 1000,          # number of warmup iterations per chain
    iter = 5000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 500           # no progress shown
)

# saveRDS(fit1, file=paste0("./BRSTIR/application/",Sys.Date(),"_stanfit.rds"))
posterior <- extract(fit1)
# str(posterior)

theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
gl.samples <- summary(fit1, par=c("newgl"), probs = c(0.05, 0.5, 0.95))$summary
gnl.samples <- summary(fit1, par=c("newgnl"), probs = c(0.05, 0.5, 0.95))$summary
# gnlfl.samples <- summary(fit1, par=c("newgfnl", "newglnl"), probs = c(0.05, 0.5, 0.95))$summary
# y.samples <- summary(fit1, par=c("y"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("newgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
# smooth.samples <- summary(fit1,par=c("gsmooth"), probs = c(0.05, 0.5, 0.95))$summary
alp.x.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary
yrep.samples <- summary(fit1, par=c("yrep"), probs = c(0.05,0.5, 0.95))$summary
f.samples <- summary(fit1, par=c("f"), probs = c(0.05,0.5, 0.95))$summary

# summary(fit1, par=c("sigma"), probs = c(0.05,0.5, 0.95))$summary
# summary(fit1, par=c("tau"), probs = c(0.05,0.5, 0.95))$summary
# gammafl.samples <- summary(fit1, par=c("gammaFL"), probs = c(0.05,0.5, 0.95))$summary

gamma.post.mean <- gamma.samples[,1]
gamma.q1 <- gamma.samples[,4]
gamma.q2 <- gamma.samples[,5]
gamma.q3 <- gamma.samples[,6]
theta.post.mean <- theta.samples[,1]
theta.q1 <- theta.samples[,4]
theta.q2 <- theta.samples[,5]
theta.q3 <- theta.samples[,6]


df.theta <- data.frame("seq" = seq(1, (p+1)),
                        "m" = c(theta.q2),
                        "l" = c(theta.q1),
                        "u" = c(theta.q3))
df.theta$covariate <- factor(c("\u03b8",names(fwi.scaled)), levels = c("\u03b8",colnames(fwi.scaled)))
df.theta$labels <- factor(c("\u03b8",colnames(fwi.scaled)))

ggplot(df.theta, aes(x = covariate, y=m, color = covariate)) + ylab("") + xlab('') +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_point(size = 5) + 
  geom_errorbar(aes(ymin = l, ymax = u), width = 0.3, linewidth =1.2) + 
  scale_x_discrete(labels = c(expression(bold(theta[0])),
                              expression(bold(theta[1])),
                              expression(bold(theta[2])),
                              expression(bold(theta[3])),
                              expression(bold(theta[4])),
                              expression(bold(theta[5])),
                              expression(bold(theta[6])),
                              expression(bold(theta[7])))) + 
  scale_color_discrete(labels = c(expression(theta[0]),colnames(fwi.scaled))) + 
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
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_theta.pdf"), width=10, height = 7.78)

df.gamma <- data.frame("seq" = seq(1, ((psi)*p)), 
                  "m" = as.vector(gamma.q2),
                  "l" = as.vector(gamma.q1),
                  "u" = as.vector(gamma.q3))
df.gamma$covariate <- factor(rep(names(fwi.scaled), each = (psi), length.out = nrow(df.gamma)), levels = colnames(fwi.scaled))
df.gamma$labels <- factor(1:((psi)*p))
ggplot(df.gamma, aes(x =labels, y = m, color = covariate)) + 
  geom_errorbar(aes(ymin = l, ymax = u), alpha = 0.4, width = 4, linewidth = 1.2) +
  geom_point(size = 4) + ylab("") + xlab("" ) + #xlim(1,(psi*p)) +
  # geom_ribbon(aes(ymin = l, ymax = u)) +
  # geom_point(size = 4, color = "black") + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  scale_x_discrete(breaks=c(seq(0, (psi*p), psi)+10), 
                    label = c(expression(bold(gamma[1])), 
                              expression(bold(gamma[2])), 
                              expression(bold(gamma[3])), 
                              expression(bold(gamma[4])), 
                              expression(bold(gamma[5])), 
                              expression(bold(gamma[6])), 
                              expression(bold(gamma[7]))),
                    expand=c(0,10)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size=25),
          legend.margin=margin(0,0,0,-10),
          legend.box.margin=margin(-10,0,-10,0),
          plot.margin = margin(0,0,0,-20),
          axis.text.x = element_text(hjust=0.5),
          axis.text = element_text(size = 28))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_gamma.pdf"), width=10, height = 7.78)


# g.linear.mean <- as.vector(matrix(gl.samples[,1], nrow = n, byrow=TRUE))
# g.linear.q1 <- as.vector(matrix(gl.samples[,4], nrow = n, byrow=TRUE))
# g.linear.q2 <- as.vector(matrix(gl.samples[,5], nrow = n, byrow=TRUE))
# g.linear.q3 <- as.vector(matrix(gl.samples[,6], nrow = n, byrow=TRUE))
g.nonlinear.mean <- as.vector(matrix(gnl.samples[,1], nrow = n, byrow=TRUE))
g.nonlinear.q1 <- as.vector(matrix(gnl.samples[,4], nrow = n, byrow=TRUE))
g.nonlinear.q2 <- as.vector(matrix(gnl.samples[,5], nrow = n, byrow=TRUE))
g.nonlinear.q3 <- as.vector(matrix(gnl.samples[,6], nrow = n, byrow=TRUE))
g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))

equal_breaks <- function(n = 3, s = 0.1,...){
  function(x){
    d <- s * diff(range(x)) / (1+2*s)
    seq = seq(min(x)+d, max(x)-d, length=n)
    round(seq, -floor(log10(abs(seq[2]-seq[1]))))
  }
}

data.smooth <- data.frame("x"= as.vector(xholder),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
                          # "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab("") +
  facet_grid(covariates ~ ., scales = "free",
              labeller = label_parsed) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) + 
  # scale_y_continuous(breaks=c(-3,-2,-1,1,2)) +        
          ylim(-3.5, 3.5) +
  # scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
          plot.margin = margin(0,0,0,-20),
          # strip.text = element_blank(),
          axis.text = element_text(size = 20))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_smooth.pdf"), width=12.5, height = 15)

# data.linear <- data.frame("x"= as.vector(xholder),
#                           "post.mean" = as.vector(g.linear.mean),
#                           "q1" = as.vector(g.linear.q1),
#                           "q2" = as.vector(g.linear.q2),
#                           "q3" = as.vector(g.linear.q3),
#                           "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
#                           "fakelab" = rep(1, (p*n)),
#                           "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

# ggplot(data.linear, aes(x=x, group=interaction(covariates, replicate))) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#   # geom_line(aes(y=true, colour = "True"), linewidth=2) + 
#   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#   ylab("") + xlab("") +
#   facet_grid(covariates ~ ., scales = "free",
#               labeller = label_parsed) + 
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values=c("steelblue")) + 
#   guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) + #ylim(-0.65, 0.3) +
#   scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#           plot.margin = margin(0,0,0,-20),
#           # strip.text = element_blank(),
#           axis.text = element_text(size = 20))
# # #ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_linear.pdf"), width=12.5, height = 15)


data.nonlinear <- data.frame("x"=as.vector(xholder),
                          "post.mean" = as.vector(g.nonlinear.mean),
                          "q1" = as.vector(g.nonlinear.q1),
                          "q2" = as.vector(g.nonlinear.q2),
                          "q3" = as.vector(g.nonlinear.q3),
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
                          "fakelab" = rep(1, (p*n)),
                          "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

ggplot(data.nonlinear, aes(x=x, group=interaction(covariates, replicate))) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
  # geom_line(aes(y=true, colour = "True"), linewidth=2) + 
  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
  ylab("") + xlab("") +
  facet_grid(covariates ~ ., scales = "free", #switch = "y",
              labeller = label_parsed) + 
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values=c("steelblue")) + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) + #ylim(-0.65, 0.3) +
  scale_y_continuous(breaks=equal_breaks(n=3, s=0.1)) + 
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
          plot.margin = margin(0,0,0,-20),
          # strip.text = element_blank(),
          axis.text = element_text(size = 20))
# #ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_nonlinear.pdf"), width=12.5, height = 15)

data.scenario <- data.frame("x" = newx,
                            #"x" = seq(min(PC1), max(PC1), length.out = n),
                            "post.mean" = (alpha.samples[,1]),
                            "post.median" = (alpha.samples[,5]),
                            # "post.median" = sort(exp(rowSums(g.q2) + rep(theta.samples[1,5], n)), decreasing = TRUE),
                            "q1" = (alpha.samples[,4]),
                            "q3" = (alpha.samples[,6]))
# data.scenario <- data.frame("x" = seq(0, 1, length.out = n),
#                             "post.mean" = (alp.x.samples[,1]),
#                             "post.median" = (alp.x.samples[,5]),
#                             "q1" = (alp.x.samples[,4]),
#                             "q3" = (alp.x.samples[,6]))

ggplot(data.scenario, aes(x=x)) + 
  ylab(expression(alpha(c,...,c))) + xlab(expression(c)) + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y = true, col = "True"), linewidth = 2) +
  # xlim(-1,1) + #ylim(0, 6.2) + 
  geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue")) + 
  # scale_y_log10() + 
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))

# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_alpha.pdf"), width=10, height = 7.78)

# data.scenario <- data.frame("x" = seq(0, 1, length.out = n),
#                             "post.mean" = (alp.x.samples[,3]),
#                             "post.median" = (alpha.samples[,5]),
#                             # "post.median" = sort(exp(rowSums(g.q2) + rep(theta.samples[1,5], n)), decreasing = TRUE),
#                             "q1" = (alpha.samples[,4]),
#                             "q3" = (alpha.samples[,6]))

# ggplot(data.scenario,aes(x=x)) + 
#   xlab(expression(alpha(bold(x)))) + #ylab(expression(c)) + labs(col = "") +
#   geom_histogram(aes(post.median), fill = "purple", binwidth = 0.3, , alpha = 0.5) + 
#   # geom_histogram(aes(post.mean), fill = "red", binwidth = 0.3, alpha = 0.4) +
#   # geom_histogram(aes(q1), fill = "green", binwidth = 0.3, alpha = 0.1) + 
#   # geom_histogram(aes(q3), fill = "black", binwidth = 0.3, alpha = 0.1) + 
#   # scale_y_log10() + 
#   # guides(color = guide_legend(order = 2), 
#   #         fill = guide_legend(order = 1)) +
#   scale_fill_brewer(palette="Pastel1") +
#   theme_minimal(base_size = 30) +
#   theme(legend.position="bottom", 
#         legend.key.size = unit(1, 'cm'),
#         legend.text = element_text(size=20),
#         # plot.margin = margin(0,0,0,-1),
#         strip.text = element_blank(),
#         axis.title.x = element_text(size = 35))

len <- dim(posterior$alpha)[1]
r <- matrix(, nrow = n, ncol = 100)
# beta <- as.matrix(mcmc[[1]])[, 1:7] 
T <- 100
for(i in 1:n){
  for(t in 1:T){
    # r[i, t] <- qnorm(pPareto(y[i], u, alpha = alp.x.samples[i,5]))
    r[i, t] <- qnorm(pPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i]))
    # r[i, t] <- qnorm(pPareto(y[i], u, alpha = posterior$newalpha[round(runif(1,1,len)),i]))
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
  geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
              fill = "steelblue",
              alpha = 0.4, linetype = "dashed") + 
  geom_line(aes(x = grid, y = trajhat), colour = "steelblue", linetype = "dashed", linewidth = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.2) + 
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 30)) + 
  coord_fixed(xlim = c(-3, 3),
              ylim = c(-3, 3))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_qqplot.pdf"), width=10, height = 7.78)
rp <-c()
for(i in 1:n){
  # rp[i] <- rPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i])
  rp[i] <- qnorm(pPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i]))
}
rp <- data.frame(rp, group = rep("residuals", n))

ggplot(data = rp) + 
  # geom_qqboxplot(aes(factor(group, levels=c("residuals")), y=rp), notch=FALSE, varwidth=TRUE, reference_dist="norm")+ 
  geom_qqboxplot(aes(y=rp), notch=FALSE, varwidth=FALSE, reference_dist="norm")+   
  labs(x = "", y = "Residuals") + ylim(-4,4) + 
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_qqboxplot.pdf"), width = 10, height = 7.78)
             
# saveRDS(data.scenario, file=paste0("Simulation/BayesianPsplines/results/",date,"-",time, "_sc1_data_samp1.rds"))

cat("Finished Running")

# relative_eff(exp(fit.log.lik))
#https://discourse.mc-staqan.org/t/four-questions-about-information-criteria-cross-validation-and-hmc-in-relation-to-a-manuscript-review/13841/3
# y.rep <- as.matrix(fit1, pars = "y_rep")
# ppc_loo_pit_overlay(
#   y = y,
#   yrep = y.rep,
#   lw = weights(fwi.loo$psis_object)
# )

grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),], origin = fwi.scaled[,i] 
                  # c = seq(min(PC1), max(PC1), length.out = n)
                  ), aes(x=x)) + 
  # grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),], origin = fwi.index[which(Y>u),i]), aes(x=x)) +   
                  # geom_point(aes(x= origin, y=q2), alpha = 0.3) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x= origin, y=q2), sides = "b") +
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  ylim(-3, 3) +
                  # scale_y_continuous(breaks=equal_breaks(n=5, s=0.1)) + 
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt
}

grid.arrange(grobs = grid.plts, ncol = 2, nrow = 4)
# grid.plts[[7]]
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_mcmc_smooth.pdf"), width=10, height = 7.78)

data.yrep <- data.frame("y" = seq(log(u), 30, length.out=n),
                          "post.mean" = f.samples[,1],
                          "q1" = f.samples[,4],
                          "q2" = f.samples[,5],
                          "q3" = f.samples[,6])

ggplot(data.yrep, aes(x=y)) + 
  ylab("Density") + xlab("log(Burned Area)") + labs(col = "") +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
  # geom_line(aes(y = true, col = "True"), linewidth = 2) +
  xlim(7.5, 30) + #ylim(0, 6.2) + 
  geom_line(aes(y=q2, col = "Posterior Median"), linewidth=2) +
  scale_fill_manual(values=c("steelblue"), name = "") +
  scale_color_manual(values = c("steelblue")) + 
  annotate("point", x= log(max(Y)), y=0, color = "red", size = 7) +
  guides(color = guide_legend(order = 2), 
          fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 20))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_pareto_post_generative_dooms.pdf"), width = 10, height = 7.78)

#Predictive Distribution check
y.container <- as.data.frame(matrix(, nrow = n, ncol = 0))  
random.alpha.idx <- floor(runif(100, 1, ncol(t(posterior$f))))
for(i in random.alpha.idx){
  y.container <- cbind(y.container, log(rPareto(n, u, t(posterior$alpha)[i])))
}
colnames(y.container) <- paste("col", 1:100, sep="")
y.container$x <- seq(1,n)
y.container$logy <- log(y)
plt <- ggplot(data = y.container, aes(x = logy)) + ylab("density") + xlab("log(Burnt Area)") + labs(col = "")

for(i in names(y.container)){
  plt <- plt + geom_density(aes(x=.data[[i]]), color = "slategray1", alpha = 0.1, linewidht = 0.7)
}

print(plt + geom_density(aes(x=logy), color = "steelblue", linewidth = 2) +
        theme_minimal(base_size = 30) + ylim(0, 2) + #xlim(7.5,25) +
        theme(legend.position = "none",
                axis.text = element_text(size = 35)))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_BRSTIR_predictive_distribution.pdf"), width=10, height = 7.78)

extreme.container <- as.data.frame(matrix(, nrow = n, ncol = 3000))
for(i in 1:3000){
  extreme.container[,i] <- density(log(posterior$yrep[i,]), n=n)$y
}
extreme.container <- cbind(extreme.container, t(apply(extreme.container[,1:3000], 1, quantile, c(0.05, .5, .95))))
colnames(extreme.container)[(dim(extreme.container)[2]-2):(dim(extreme.container)[2])] <- c("q1","q2","q3")
colnames(extreme.container)[1:3000] <- as.character(1:3000)
extreme.container$mean <- rowMeans(extreme.container[,1:3000])
extreme.container$y <- seq(0, 30, length.out = n)
extreme.container <- as.data.frame(extreme.container)


plt <- ggplot(data = extreme.container, aes(x = y)) + xlab("log(Burned Area)") + ylab("Density")+
        # geom_line(aes(y=q2), colour = "steelblue", linewidth = 1.5) +
        geom_line(aes(y=mean), colour = "steelblue", linewidth = 1.5) +
        geom_ribbon(aes(ymin = q1, ymax = q3), fill = "steelblue", alpha = 0.2) + 
        theme_minimal(base_size = 30) + 
        theme(legend.position = "none",
              axis.title = element_text(size = 30))
d <- ggplot_build(plt)$data[[1]]
print(plt + 
        # geom_area(data = subset(d, x>12.44009), aes(x=x,y=y), fill = "slategray1", alpha = 0.5) +
        geom_segment(x=12.44009, xend=12.44009, 
              y=0, yend=approx(x = d$x, y = d$y, xout = 12.4409)$y,
              colour="red", linewidth=1.2, linetype = "dotted"))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_doomsday_generative.pdf"), width = 10, height = 7.78)              

data.extreme <- data.frame("x" = log(y),
                            "post.mean" = log(yrep.samples[,1]),
                            "post.median" = log(yrep.samples[,5]),
                            "q1" = log(yrep.samples[,4]),
                            "q3" = log(yrep.samples[,6]))
# data.extreme <- data.frame("x" = y,
#                             "post.mean" = (f.samples[,1]),
#                             "post.median" = (f.samples[,5]),
#                             "q1" = (f.samples[,4]),
#                             "q3" = (f.samples[,6]))



ggplot(data.extreme, aes(x=x)) + 
  ylab("Density") + xlab("Burned Area") + #labs(col = "") +
  geom_density(aes(x=post.mean), colour = "steelblue") + 
  # geom_density(aes(x=post.median), colour = "steelblue") + 
  # geom_density(aes(x=q1), colour = "purple") + 
  geom_density(aes(x=q3), colour = "red") + 
  # geom_line(aes(x, y = post.mean), colour="steelblue", linewidth = 1) +
  # geom_ribbon(aes(ymin = q1, ymax = q3), fill = "steelblue", alpha = 0.2) +
  # scale_fill_manual(values=c("steelblue"), name = "") +
  # scale_color_manual(values = c("steelblue")) + 
  # scale_y_log10() + 
  # guides(color = guide_legend(order = 2), 
          # fill = guide_legend(order = 1)) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position="none",
        strip.text = element_blank(),
        axis.title.x = element_text(size = 35))


temp.data <- curve_interval(log(posterior$yrep), .width=c(0.05,0.5,0.95))
data.extreme <- data.frame(temp.data[(n+1):(2*n),], logy = log(y), y.origin = y)
# data.extreme <- data.frame(median_qi(log(posterior$yrep), .width=c(0.05,0.5,0.95)),logy = log(y), y.origin = y)
ggplot(data.extreme, aes(x=logy)) + xlab("Burned Area") + ylab("Density")+
  geom_density(aes(x=log(.value)), colour = "steelblue") + 
  # geom_density(aes(x=post.median), colour = "steelblue") + 
  # geom_density(aes(x=log(.lower)), colour = "purple") + 
  geom_density(aes(x=log(.upper)), colour = "red") + xlim(1, 30) + 
  # geom_lineribbon(aes(ymin=sort(.lower), ymax=sort(.upper))) + 
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        legend.position="none",
        plot.margin = margin(0,0,0,-1),
        strip.text = element_blank(),
        axis.title.x = element_text(size = 35))

random.alpha.idx <- floor(runif(1000, 1, ncol(t(posterior$alpha))))
ev.y1 <- ev.y2 <- as.data.frame(matrix(, nrow = 1, ncol = 0))
# for(i ?in 1:ncol(t(posterior$theta))){
ev.alpha.single <- c()  
for(i in random.alpha.idx){
  ev.y1 <- rbind(ev.y1, as.numeric(posterior$yrep[i]))
  # ev.y2 <- rbind(ev.y2, as.numeric(posterior$yrep2[i]))
}
ev.y1 <- as.data.frame(log(ev.y1))
ev.y1$logy <- log(max(y))
colnames(ev.y1) <- c("yrep", "logy")
ev.y1$group <- rep("15th Oct 2017",1000)
# ggplot(data=ev.y, aes(x=yrep, y = group)) +
#   ylab("") + 
#   xlab("log(Burnt Area)") + labs(col = "") +
#   stat_slab(scale = 0.6, colour = "steelblue", fill=NA, slab_linewidth = 1.5, trim = FALSE, expand = TRUE, density = "unbounded", subguide="outside", justification = -0.01) +
#   # stat_spike(aes(linetype = after_stat(at)), at = c("median"), scale=0.7)+
#   stat_dotsinterval(subguide = 'integer', side = "bottom", scale = 0.6, slab_linewidth = NA, position = "dodge") +
#   # geom_point(position = position_jitter(seed = 1, height = 0.05), alpha = 0.1) +  
#   # geom_boxplot(width = 0.2, notch = TRUE, alpha = 0.25, outlier.color = NA) +
#   geom_vline(xintercept = log(max(y)), linetype="dashed", color = "red",) +
#   # geom_label(aes(log(max(y)), 1), label = "Target Length", show.legend = FALSE)+
#   geom_vline(xintercept = log(y[133]), linetype="dashed", color = "black",) +
#   # geom_label(aes(log(y[133]), 1), label = "Target Length", show.legend = FALSE)+
#   theme_minimal(base_size = 30) +  
#   theme(legend.position = "none",
#         plot.margin = margin(0,0,0,25),
#         axis.text.y = element_text(angle = 90, size = 15, vjust = 15, hjust = 0.5),
#         axis.title = element_text(size = 30)) +
#         annotate(x=(log(max(y))+2), y= 0.1, label = "15th Oct 2017", geom="label") +
#         annotate(x=(log(y[133])-2), y= 0.1, label = "18th Jun 2017", geom="label")
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_BRSTIR_two_generative.pdf"), width = 10, height = 7.78)

plt <- ggplot(data = ev.y1, aes(x = yrep)) + ylab("Density") + xlab("log(Burned Area)") + labs(col = "") +
  geom_density(color = "steelblue", linewidth = 1.2) + 
  geom_rug(alpha = 0.1) + 
  # geom_point(aes(x=yrep,y=-Inf),color="steelblue", size = 3.5, alpha = 0.2) +
  xlim(5.5, 40) +
  theme_minimal(base_size = 30) +  
  theme(legend.position = "none",
        axis.title = element_text(size = 30))
  # geom_vline(xintercept = log(max(y)), linetype="dotted", color = "red",) +

d <- ggplot_build(plt)$data[[1]]
print(plt + geom_area(data = subset(d, x>11.48771), aes(x=x,y=y), fill = "slategray1", alpha = 0.5) +
        geom_segment(x=11.48771, xend=11.48771, 
              y=0, yend=approx(x = d$x, y = d$y, xout = 11.48771)$y,
              colour="red", linewidth=1.2, linetype = "dashed"))
# ggsave(paste0("./BRSTIR/application/figures/",Sys.Date(),"_BRSTIR_generative.pdf"), width = 10, height = 7.78)

# density.y <- density(ev.y1$yrep[1:1000]) # see ?density for parameters
# plot(density.y$x,density.y$y, type="l") #can use ggplot for this too
# # set an Avg.position value
# Avg.pos <- 12.44009
# xt <- diff(density.y$x[density.y$x>Avg.pos])
# library(zoo)
# yt <- rollmean(density.y$y[density.y$x>Avg.pos],2)
# # This gives you the area
# sum(xt*yt)

# library(ismev)
# gpd.fit(y, u)

fit.log.lik <- extract_log_lik(fit1)
loo(fit.log.lik, is_method = "sis", cores = 2)
fwi.loo <- loo(fit.log.lik, cores = 2)
plot(fwi.loo, label_points = TRUE)

constraint.elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 2)
constraint.waic <- waic(fit.log.lik, cores = 2)
# save(constraint.elpd.loo, constraint.waic, file = (paste0("./BRSTIR/application/BRSTIR_constraint_",Sys.Date(),"_",psi,"_",floor(threshold*100),"quantile_IC.Rdata")))




