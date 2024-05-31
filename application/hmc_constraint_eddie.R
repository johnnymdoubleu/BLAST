library(npreg)
library(Pareto)
library(tidyr)
library(readxl)
library(rstan)
library(loo)

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

psi <- 30
threshold <- 0.975


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


# fit <- bltir(Y, fwi.scaled, threshold = 0.975, psi = 10, T = 2000)

u <- quantile(Y, threshold)
y <- Y[Y>u]
fwi.origin <- fwi.scaled <-fwi.scaled[which(Y>u),]
# fwi.scaled <- as.data.frame(scale(fwi.scaled))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled <- as.data.frame(sapply(fwi.origin, FUN = range01))

n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]

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

for(i in 1:p){
  xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
  test.knot <- quantile(fwi.scaled[,i], probs=seq(1/(psi+1),psi/(psi+1),length.out = psi))
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  # knots <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
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


model.stan <- "data {
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
    matrix[n, p] gsmooth; // linear component


    lambda2o=lambda2*100;
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
    vector[n] f;
    real yrep;
    array[n] real <lower=0> newalpha; // new tail index
    matrix[n, p] newgsmooth; // linear component    

    for (j in 1:p){
        newgsmooth[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j] + xholderLinear[,j] * theta[j+1];
    };    

    for (i in 1:n){ 
        newalpha[i] = exp(theta[1] + sum(newgsmooth[i,]));
    };    

    yrep = pareto_rng(u, alpha[362]);
    for(i in 1:n){
      f[i] = alpha[362]*(y[i]/u)^(-alpha[362])/y[i]; //pareto_rng(u, alpha[i])
      log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}
"

data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), indexFL = as.vector(t(index.holder)),
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear, basisFL = basis.holder)

init.alpha <- list(list(gammaTemp = array(rep(0, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0, (p+1)), 
                        tau1 = rep(0.1, p), tau2 = rep(0.1, p),
                        lambda1 = 0.1, lambda2 = 0.01),
                   list(gammaTemp = array(rep(0, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0, (p+1)), 
                        tau1 = rep(0.001, p), tau2 = rep(0.001, p),
                        lambda1 = 100, lambda2 = 1),
                   list(gammaTemp = array(rep(0, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0.1, (p+1)), 
                        tau1 = rep(0.5, p), tau2 = rep(0.5, p),
                        lambda1 = 5, lambda2 = 5.5))

# stanc("C:/Users/Johnny Lee/Documents/GitHub/BRSTIR/application/model1.stan")
fit1 <- stan(
    file = model.stan,  # Stan program
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value
    # init_r = 1,
    chains = 3,             # number of Markov chains
    # warmup = 1000,          # number of warmup iterations per chain
    iter = 2000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 500           # no progress shown
)

# saveRDS(fit1, file=paste0("./BRSTIR/application/",Sys.Date(),"_stanfit.rds"))
posterior <- extract(fit1)
save(fit1, posterior, file = "eddie.Rdata")

