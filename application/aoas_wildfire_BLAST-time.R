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
library(pbs)
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
psi <- 30
threshold <- 0.975
# u <- quantile(Y[Y>1], threshold)

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
    fwi.index[,i] <- cov.long$measurement[missing.values]
    fwi.scaled[,i] <- cov.long$measurement[missing.values]
}

fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)
fwi.origin <- fwi.scaled
fwi.origin <- data.frame(fwi.origin, time = c(1:length(Y)), BA=Y)
# fwi.origin <- data.frame(fwi.origin[which(Y>1),], BA=Y[Y>1])
# BA.shifted <- ifelse(fwi.origin$BA == 0, 1e-5, fwi.origin$BA)
# fwi.origin$log.BA <- log(fwi.origin$BA+1)
vars <- colnames(fwi.origin)[1:7]
seq_list <- lapply(vars, function(var) seq(min(fwi.origin[[var]]), max(fwi.origin[[var]]), length.out = length(Y)))
grid.df <- as.data.frame(setNames(seq_list, vars))

fwi.index$BA <- Y
# 1. Prepare the Data
# We create a new dataframe 'df_seasonal' to avoid overwriting your original data
df_seasonal <- fwi.index %>%
  mutate(
    year = as.numeric(as.character(year)),
    Month_Num = match(month, month.abb),
    Season = case_when(
      Month_Num %in% c(12, 1, 2) ~ "Winter",
      Month_Num %in% c(3, 4, 5) ~ "Spring",
      Month_Num %in% c(6, 7, 8) ~ "Summer",
      Month_Num %in% c(9, 10, 11) ~ "Autumn"
    ), SeasonYear = ifelse(Month_Num == 12, year + 1, year)
  )
plot_data <- df_seasonal %>%
  group_by(SeasonYear, Season) %>%
  summarise(
    Q975 = quantile(BA, 0.975, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  arrange(SeasonYear, Season) %>%
  mutate(Label = paste(SeasonYear, Season)) %>%
  mutate(Label = factor(Label, levels = unique(Label)))

df_seasonal <- df_seasonal %>% 
  mutate(Label = paste(SeasonYear, Season)) %>%
  mutate(Label = factor(Label, levels = unique(Label)))
# 3. Plot
ggplot(plot_data, aes(x = Label, y = Q975, group = 1)) +
# ggplot(df_seasonal, aes(x = Label, y = BA, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(aes(color = Season), size = 3) +
  scale_color_manual(values = c(
    "Summer" = "orange", 
    "Winter" = "red", 
    "Spring" = "red", 
    "Autumn" = "red"
  )) +
  labs(y = "BA", x = "Season", color = "Season") +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "none"
  )

u <- quantile(Y, threshold)
y <- Y[Y>u]
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled <- as.data.frame(sapply(fwi.origin[which(Y>u), c(1:7)], FUN = range01))


fwi.origin <- data.frame(fwi.index[which(Y>u),], BA=y)
max.fwi <- fwi.origin[which.max(y),]

# ggplot(fwi.origin, aes(x=DSR, y=FFMC)) + 
#   geom_point(aes(colour = BA), size= 2.5) + 
#   scale_colour_stepsn(colours = c("slategray1", "red"), labels=function(x) format(x, big.mark = ",", scientific = TRUE), breaks=c(0.1e5, 0.5e5, 1e5, 2e5)) +
#   geom_density2d(colour="steelblue", linewidth = 1.3) + 
#   geom_mark_circle(aes(x = max.fwi$DSR, y = max.fwi$FFMC, label = "15th Oct 2017"), con.type = "straight",
#                    radius = unit(2.5, "mm"), color = "steelblue", size = 1, 
#                    con.colour = "steelblue", con.cap = unit(0, "mm"),
#                    label.colour = "steelblue", label.buffer = unit(5, "mm"),
#                    label.fill = "transparent")  +
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 30),
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         strip.text = element_blank(),
#         axis.title = element_text(size = 30))
# # ggsave("./BLAST/application/figures/extremeviz.pdf", width = 10, height = 7.78)

# ggplot(fwi.origin, aes(x=as.numeric(year), y=log(BA), color = BA)) + 
#   ylab("Hectares (log)") + xlab("Time (years)") + 
#   geom_point(size= 2.5, alpha = 0.5) + 
#   scale_colour_stepsn(colours = c("slategray1", "red"), labels=function(y) format(y, big.mark = ",", scientific = TRUE), 
#   breaks = quantile(fwi.origin$BA, probs = seq(0,1,length.out = 20))) + 
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 30),
#         legend.position = "none",
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         strip.text = element_blank(),
#         axis.title = element_text(size = 30))

# # ggsave("./BLAST/application/figures/hectareslog.pdf", width = 10, height = 7.78)

# ggplot(fwi.origin, aes(x=as.numeric(year))) + 
#   ylab("Density") + xlab("Time (years)") + 
#   geom_histogram(aes(y = after_stat(density)), fill = "steelblue", color = "gray", alpha = .2) +
#   geom_rug() +
#   theme_minimal(base_size = 30) +
#   theme(plot.title = element_text(hjust = 0.5, size = 30),
#         legend.position = "none",
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         strip.text = element_blank(),
#         axis.title = element_text(size = 30))



# ggsave("./BLAST/application/figures/intensityfn.pdf", width = 10, height = 7.78)
df_spline <- fwi.index %>%
  mutate(
    # Create a clean Date column (e.g., "1980-01-01")
    Full_Date = dmy(paste(date, month, year)),
    
    # Extract Day of Year (1 to 366)
    DOY = yday(Full_Date)
  )

df_spline <- df_spline[which(Y>u),]

boundary_knots <- c(1, 366) # Covers leap years
internal_knots <- seq(from = boundary_knots[1], 
                      to = boundary_knots[2], 
                      length.out = psi + 2)[-c(1, psi + 2)]

# 3. Generate the Periodic B-Spline Basis
# This creates a matrix where each column is one basis function
pbs_basis <- pbs(
  df_spline$DOY, 
  knots = seq(1, 366, length.out = psi)
)

time.basis <- pbs(
  seq(0, 1,length.out = nrow(fwi.origin)),
  knots = seq(0, 1, length.out = psi)
  # Boundary.knots = c(0,1)
)

grid.basis <- pbs(
  seq(1, 366, length.out = nrow(fwi.origin)),
  knots = seq(1, 366, length.out = psi)
  # Boundary.knots = c(1, 366)
)

n <- nrow(fwi.scaled)
p <- ncol(fwi.scaled)
no.theta <- 1 #represents the no. of linear predictors for each smooth functions
newx <- seq(0, 1, length.out=n)
grid.linear <- grid.nonlinear <- xholder.linear <- xholder.nonlinear <- bs.linear <- bs.nonlinear <- matrix(,nrow=n, ncol=0)
grid.x <- xholder <- matrix(nrow=n, ncol=p)
end.holder <- basis.holder <- matrix(, nrow = 2, ncol =0)
index.holder <- matrix(, nrow = 0, ncol = 2)
for(i in 1:p){
  index.holder <- rbind(index.holder, 
                      matrix(c(which.min(fwi.scaled[,i]),
                              which.max(fwi.scaled[,i])), ncol=2))
}

for(i in 1:p){
  xholder[,i] <- seq(0, 1, length.out = n)
  test.knot <- seq(0, 1, length.out = psi)
  splines <- basis.tps(xholder[,i], test.knot, m=2, rk=FALSE, intercept = FALSE)
  xholder.linear <- cbind(xholder.linear, splines[,1:no.theta])
  xholder.nonlinear <- cbind(xholder.nonlinear, splines[,-c(1:no.theta)])
  grid.x[,i] <- seq(min(df_spline[,i]), max(df_spline[,i]), length.out = n)
  grid.knot <- seq(min(df_spline[,i]), max(df_spline[,i]), length.out = psi)
  splines <- basis.tps(grid.x[,i], grid.knot, m=2, rk=FALSE, intercept = FALSE)
  grid.linear <- cbind(grid.linear, splines[,1:no.theta])
  grid.nonlinear <- cbind(grid.nonlinear, splines[,-c(1:no.theta)])  
  knots <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = psi)
  tps <- basis.tps(fwi.scaled[,i], knots, m = 2, rk = FALSE, intercept = FALSE)
  basis.holder <- cbind(basis.holder, 
          solve(matrix(c(tps[index.holder[i,1], no.theta+1],
                  tps[index.holder[i,1], no.theta+psi],
                  tps[index.holder[i,2], no.theta+1],
                  tps[index.holder[i,2], no.theta+psi]), 
                  nrow = 2, ncol = 2)))
  # end.holder <- cbind(end.holder, 
  #               matrix(c(tps[index.holder[i,1], no.theta+1],
  #                 tps[index.holder[i,1], no.theta+psi],
  #                 tps[index.holder[i,2], no.theta+1],
  #                 tps[index.holder[i,2], no.theta+psi]), 
  #                 nrow = 2, ncol = 2))
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
    matrix[n,p] gridLinear; // fwi dataset
    matrix[n, (psi*p)] gridNonlinear; // thin plate splines basis    
    matrix[n, psi] pbs_basis;
    matrix[n, psi] pbs_time;
    matrix[n, psi] grid_time;
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
    vector[psi] gamma_time;
    real <lower=0> tau_time;
}
transformed parameters {
    array[n] real <lower=0> alpha; // covariate-adjusted tail index
    vector[psi] gamma[p]; // splines coefficient 
    matrix[n, (p+1)] gsmooth; // linear component
    real <lower=0> lambda2o;

    lambda2o=lambda2*100;
    {
      vector[2] gammaFL[p];
      
      for(j in 1:p){
          gamma[j][2:(psi-1)] = gammaTemp[j][1:(psi-2)]*100;
          gammaFL[j] = basisFL[, (((j-1)*2)+1):(((j-1)*2)+2)] * (bsNonlinear[indexFL[(((j-1)*2)+1):(((j-1)*2)+2)], (((j-1)*psi)+2):(((j-1)*psi)+(psi-1))] * gammaTemp[j]*100) * (-1);
          gamma[j][1] = gammaFL[j][1];
          gamma[j][psi] = gammaFL[j][2];
      };
    }
    for (j in 1:p){
        gsmooth[,j] = bsNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j] + bsLinear[,j] * theta[j+1];
    };
    gsmooth[,(p+1)] = pbs_basis * gamma_time;

    for (i in 1:n){
        alpha[i] = exp(theta[1] + sum(gsmooth[i,]));
    };
}

model {
    // likelihood
    for (i in 1:n){
        target += pareto_lpdf(y[i] | u, alpha[i]);
    }
    target += normal_lpdf(theta[1] | 0, 10);
    target += gamma_lpdf(lambda1 | 1, 1e-3);
    target += gamma_lpdf(lambda2o | 1, 0.1);
    target += (2*p*log(lambda2o));
    for (j in 1:p){
        target += gamma_lpdf(tau1[j] | 1, lambda1^2*0.5); // target += exponential_lpdf(tau1[j] | 0.5 * lambda1^2)
        target += normal_lpdf(theta[(j+1)] | 0, sqrt(1/tau1[j])); // target += normal_lpdf(theta[(j+1)] | 0, tau1[j]) target += double_exponential_lpdf(theta[(j+1)] | 0, lambda1)
        target += gamma_lpdf(tau2[j] | atau, lambda2o^2*0.5);
        target += multi_normal_lpdf(gamma[j] | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * (1/tau2[j]));
    }
    target += gamma_lpdf(tau_time | atau, lambda2o^2*0.5);
    target += multi_normal_lpdf(gamma_time | rep_vector(0, psi), diag_matrix(rep_vector(1, psi)) * (1/tau_time));
}

generated quantities {
    // Used in Posterior predictive check    
    vector[n] log_lik;
    real yrep;
    array[n] real <lower=0> newalpha; // new tail index
    vector[n] f;
    matrix[n, (p+1)] gridgsmooth; // smooth component
    matrix[n, p] gridlinear; // linear component
    matrix[n, (p+1)] gridnonlinear; // nonlinear component
    yrep = pareto_rng(u, alpha[1]); 

    {
    matrix[n, (p+1)] newgsmooth; // linear component
      for (j in 1:p){
        newgsmooth[,j] = xholderNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j] + xholderLinear[,j] * theta[j+1];
        gridnonlinear[,j] = gridNonlinear[,(((j-1)*psi)+1):(((j-1)*psi)+psi)] * gamma[j];
        gridlinear[,j] = gridLinear[,j] * theta[j+1];
        gridgsmooth[,j] = gridnonlinear[,j] + gridlinear[,j]; 
      };
      newgsmooth[,(p+1)] = pbs_time * gamma_time;
      gridnonlinear[,(p+1)] = grid_time * gamma_time;
      gridgsmooth[,(p+1)] = gridnonlinear[,(p+1)];

      for (i in 1:n){ 
        newalpha[i] = exp(theta[1] + sum(newgsmooth[i,]));
      };    
    }
    for(i in 1:n){
      f[i] = pareto_rng(u, alpha[i]);
      log_lik[i] = pareto_lpdf(y[i] | u, alpha[i]);
    }
}
"


data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, 
                    atau = ((psi+1)/2), indexFL = as.vector(t(index.holder)),
                    bsLinear = bs.linear, bsNonlinear = bs.nonlinear,
                    xholderLinear = xholder.linear, xholderNonlinear = xholder.nonlinear, basisFL = basis.holder, pbs_basis = pbs_basis, pbs_time = time.basis, grid_time = grid.basis,
                    gridLinear = grid.linear, gridNonlinear = grid.nonlinear)

init.alpha <- list(list(gammaTemp = array(rep(0, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0, (p+1)), gamma_time = rep(0, psi),
                        tau1 = rep(0.1, p), tau2 = rep(0.1, p), tau_time = 0.1,
                        lambda1 = 0.1, lambda2 = 0.01),
                   list(gammaTemp = array(rep(0, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0, (p+1)), gamma_time = rep(0, psi),
                        tau1 = rep(0.001, p), tau2 = rep(0.001, p), tau_time = 0.1,
                        lambda1 = 20, lambda2 = 1),
                   list(gammaTemp = array(rep(0, ((psi-2)*p)), dim=c(p, (psi-2))),
                        theta = rep(0.1, (p+1)), gamma_time = rep(0, psi),
                        tau1 = rep(0.5, p), tau2 = rep(0.5, p), tau_time = 0.1,
                        lambda1 = 5, lambda2 = 5.5))

# stanc("C:/Users/Johnny Lee/Documents/GitHub/BLAST/application/model1.stan")
fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value 
    chains = 3,             # number of Markov chains
    iter = 30000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 5000           # no progress shown
)

# saveRDS(fit1, file=paste0("./BLAST/application/",Sys.Date(),"_stanfit.rds"))
# readRDS(file=paste0("./BLAST/application/2024-11-27","_stanfit.rds"))
posterior <- rstan::extract(fit1)

theta.samples <- summary(fit1, par=c("theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
linear.samples <- summary(fit1, par=c("gridlinear"), probs = c(0.05, 0.5, 0.95))$summary
nonlinear.samples <- summary(fit1, par=c("gridnonlinear"), probs = c(0.05, 0.5, 0.95))$summary
origin.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
alpha.samples <- summary(fit1, par=c("newalpha"), probs = c(0.05,0.5, 0.95))$summary
yrep <- summary(fit1, par=c("yrep"), probs = c(0.05,0.5, 0.95))$summary
f.samples <- summary(fit1, par=c("f"), probs = c(0.05,0.5, 0.95))$summary
loglik.samples <- summary(fit1, par=c("log_lik"), probs = c(0.05,0.5, 0.95))$summary

MCMCvis::MCMCplot(fit1, params = 'theta')
# MCMCvis::MCMCsummary(fit1, params = "gamma")
# MCMCvis::MCMCsummary(fit1, params = "gamma_time")
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))
# posterior::neff_tail(fit1, params = "gamma")

g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))
g.linear.mean <- as.vector(matrix(linear.samples[,1], nrow = n, byrow=TRUE))
g.linear.q1 <- as.vector(matrix(linear.samples[,4], nrow = n, byrow=TRUE))
g.linear.q2 <- as.vector(matrix(linear.samples[,5], nrow = n, byrow=TRUE))
g.linear.q3 <- as.vector(matrix(linear.samples[,6], nrow = n, byrow=TRUE))
g.nonlinear.mean <- as.vector(matrix(nonlinear.samples[,1], nrow = n, byrow=TRUE))
g.nonlinear.q1 <- as.vector(matrix(nonlinear.samples[,4], nrow = n, byrow=TRUE))
g.nonlinear.q2 <- as.vector(matrix(nonlinear.samples[,5], nrow = n, byrow=TRUE))
g.nonlinear.q3 <- as.vector(matrix(nonlinear.samples[,6], nrow = n, byrow=TRUE))

data.smooth <- data.frame("x"= as.vector(cbind(xholder, seq(1, 366, length.out=nrow(fwi.origin)))),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "covariates" = gl((p+1), n, ((p+1)*n), labels = c(names(fwi.scaled),"time")),
                          "replicate" = gl((p+1), n, ((p+1)*n), labels = c(names(fwi.scaled),"time")))

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
          xlim(0, 1) + # ylim(-3.5, 3.5) + 
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
          plot.margin = margin(0,0,0,-20),
          # strip.text = element_blank(),
          axis.text = element_text(size = 20))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_smooth.pdf"), width=12.5, height = 15)

data.scenario <- data.frame("x" = xholder[,1],
                            "post.mean" = (alpha.samples[,1]),
                            "post.median" = (alpha.samples[,5]),
                            "q1" = (alpha.samples[,4]),
                            "q3" = (alpha.samples[,6]),
                            "post.mean.org" = (origin.samples[,1]),
                            "post.median.org" = (origin.samples[,5]),
                            "q1.org" = (origin.samples[,4]),
                            "q3.org" = (origin.samples[,6]))

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



T <- 100
len <- dim(posterior$alpha)[1]
posterior_idx <- sample(1:len, T, replace = TRUE)
r <- sapply(1:T, function(t) {
  alpha_t <- posterior$alpha[posterior_idx[t], ]  # Extract vector of length n
  qnorm(pPareto(y, u, alpha = alpha_t))           # Vectorized across all y
})
quantile.prob <- ppoints(n)
grid <- qnorm(quantile.prob)
traj <- t(apply(r, 2, quantile, probs = quantile.prob, type = 2))
# r <- matrix(, nrow = n, ncol = T)
# for(i in 1:n){
#   for(t in 1:T){
#     r[i, t] <- qnorm(pPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i]))
#   }
# }

# traj <- matrix(NA, nrow = T, ncol = lgrid)
# for (t in 1:T){
#   traj[t, ] <- quantile(r[, t], quantile.prob, type = 2)
# }

l.band <- apply(traj, 2, quantile, prob = 0.025)
trajhat <- apply(traj, 2, quantile, prob = 0.5)
u.band <- apply(traj, 2, quantile, prob = 0.975)
qqplot.df <- data.frame(grid = grid, 
                        l.band = l.band, 
                        trajhat = trajhat, 
                        u.band = u.band)
ggplot(data = ) + 
  geom_ribbon(aes(x = grid, ymin = l.band, ymax = u.band), 
              fill = "steelblue",
              alpha = 0.4, linetype = "dashed") + 
  geom_line(aes(x = grid, y = trajhat), colour = "steelblue", linetype = "dashed", linewidth = 1.2) + 
  geom_abline(intercept = 0, slope = 1, linewidth = 1.2) + 
  labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
  theme_minimal(base_size = 30) +
  theme(axis.text = element_text(size = 20)) + 
  coord_fixed(xlim = c(-3, 3),
              ylim = c(-3, 3))

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_qqplot.pdf"), width=10, height = 7.78)
# save(loglik.samples, data.smooth, data.scenario, qqplot.df, file="./BLAST/application/blast_1.Rdata")

rp <-c()
for(i in 1:n){
  # rp[i] <- rPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i])
  rp[i] <- qnorm(pPareto(y[i], u, alpha = posterior$alpha[round(runif(1,1,len)),i]))
}
rp <- data.frame(rp, group = rep("residuals", n))

ggplot(data = rp) + 
  # geom_qqboxplot(aes(factor(group, levels=c("residuals")), y=rp), notch=FALSE, varwidth=TRUE, reference_dist="norm")+ 
  geom_qqboxplot(aes(y=rp), notch=FALSE, varwidth=FALSE, reference_dist="norm", width = 0.15, qq.colour = "steelblue")+
  labs(x = "", y = "Residuals") + ylim(-4,4) + xlim(-.2,.2)+
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_qqboxplot.pdf"), width = 10, height = 7.78)
             
cat("Finished Running")

# relative_eff(exp(fit.log.lik))
#https://discourse.mc-staqan.org/t/four-questions-about-information-criteria-cross-validation-and-hmc-in-relation-to-a-manuscript-review/13841/3


data.smooth <- data.frame("x"= as.vector(cbind(grid.x, seq(1, 366, length.out=nrow(df_spline)))),
                          "true" = as.vector(as.matrix(df_spline[,c(1:7,13)])),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3),
                          "covariates" = gl((p+1), n, ((p+1)*n), labels = c(names(fwi.scaled),"T")))


grid.plts <- list()
for(i in 1:(p+1)){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab("") + xlab(c(names(fwi.scaled),"T")[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  # ylim(-5, 3) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt #+ annotate("point", x= cbind(fwi.scaled, df_spline$DOY)[which.max(y),i], y=-5, color = "red", size = 7)
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

data.smooth <- data.frame("x"= as.vector(grid.x),
                          "true" = as.vector(as.matrix(df_spline[,c(1:7)])),
                          "post.mean" = as.vector(g.linear.mean),
                          "q1" = as.vector(g.linear.q1),
                          "q2" = as.vector(g.linear.q2),
                          "q3" = as.vector(g.linear.q3),
                          "covariates" = gl((p), n, ((p)*n), labels = names(fwi.scaled)))


grid.plts <- list()
for(i in 1:(p)){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  # ylim(-5, 3) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt #+ annotate("point", x= cbind(fwi.scaled, df_spline$DOY)[which.max(y),i], y=-5, color = "red", size = 7)
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

data.smooth <- data.frame("x"= as.vector(cbind(grid.x, seq(1, 366, length.out=nrow(df_spline)))),
                          "true" = as.vector(as.matrix(df_spline[,c(1:7,13)])),
                          "post.mean" = as.vector(g.nonlinear.mean),
                          "q1" = as.vector(g.nonlinear.q1),
                          "q2" = as.vector(g.nonlinear.q2),
                          "q3" = as.vector(g.nonlinear.q3),
                          "covariates" = gl((p+1), n, ((p+1)*n), labels = c(names(fwi.scaled),"T")))


grid.plts <- list()
for(i in 1:(p+1)){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  # geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab("") + xlab(c(names(fwi.scaled),"T")[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  # ylim(-5, 3) +
                  theme_minimal(base_size = 30) +
                  theme(legend.position = "none",
                          plot.margin = margin(0,0,0,-20),
                          axis.text = element_text(size = 35),
                          axis.title.x = element_text(size = 45))
  grid.plts[[i]] <- grid.plt #+ annotate("point", x= cbind(fwi.scaled, df_spline$DOY)[which.max(y),i], y=-5, color = "red", size = 7)
}

grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_DC.pdf"), grid.plts[[7]], width=10, height = 7.78)

# saveRDS(data.smooth, file="./BLAST/application/figures/comparison/full_stanfit.rds")

#Predictive Distribution check
# y.container <- as.data.frame(matrix(, nrow = n, ncol = 0))  
# random.alpha.idx <- floor(runif(100, 1, ncol(t(posterior$f))))
# for(i in random.alpha.idx){
#   y.container <- cbind(y.container, log(t(posterior$f)[,i]))
# }
# colnames(y.container) <- paste("col", 1:100, sep="")
# y.container$x <- seq(1,n)
# y.container$logy <- log(y)
# plt <- ggplot(data = y.container, aes(x = logy)) + ylab("Density") + xlab("log(Burned Area)") + labs(col = "")

# for(i in names(y.container)){
#   plt <- plt + geom_density(aes(x=.data[[i]]), color = "slategray1", alpha = 0.1, linewidht = 0.7)
# }

# print(plt + geom_density(aes(x=logy), color = "steelblue", linewidth = 2) +
#         theme_minimal(base_size = 30) + ylim(0, 1.25) + xlim(7.5,30) +
#         theme(legend.position = "none",
#                 axis.text = element_text(size = 35)))
# # ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_BLAST_predictive_distribution.pdf"), width=10, height = 7.78)


# extreme.container <- as.data.frame(matrix(, nrow = n, ncol = 3000))
# for(i in 1:3000){
#   extreme.container[,i] <- density(log(posterior$f[i,]), n=n)$y
# }
# extreme.container <- cbind(extreme.container, t(apply(extreme.container[,1:3000], 1, quantile, c(0.05, .5, .95))))
# colnames(extreme.container)[(dim(extreme.container)[2]-2):(dim(extreme.container)[2])] <- c("q1","q2","q3")
# colnames(extreme.container)[1:3000] <- as.character(1:3000)
# extreme.container$mean <- rowMeans(extreme.container[,1:3000])
# extreme.container$y <- seq(0, 30, length.out = n)
# extreme.container <- as.data.frame(extreme.container)


# plt <- ggplot(data = extreme.container, aes(x = y)) + xlab("log(Burned Area)") + ylab("Density")+
#         geom_line(aes(y=mean), colour = "steelblue", linewidth = 1.5) +
#         geom_ribbon(aes(ymin = q1, ymax = q3), fill = "steelblue", alpha = 0.2) + 
#         theme_minimal(base_size = 30) + 
#         theme(legend.position = "none",
#               axis.title = element_text(size = 30))
# d <- ggplot_build(plt)$data[[1]]
# print(plt + 
#         geom_segment(x=12.44009, xend=12.44009, 
#               y=0, yend=approx(x = d$x, y = d$y, xout = 12.4409)$y,
#               colour="red", linewidth=1.2, linetype = "dotted"))
# # ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_post_generative.pdf"), width = 10, height = 7.78)

# random.alpha.idx <- floor(runif(1000, 1, ncol(t(posterior$alpha))))
# ev.y1 <- ev.y2 <- as.data.frame(matrix(, nrow = 1, ncol = 0))
# ev.alpha.single <- c()  
# for(i in random.alpha.idx){
#   ev.y1 <- rbind(ev.y1, as.numeric(posterior$yrep[i]))
# }
# ev.y1 <- as.data.frame(log(ev.y1))
# ev.y1$logy <- max(log(y))
# colnames(ev.y1) <- c("yrep", "logy")
# ev.y1$group <- rep("15th Oct 2017",1000)
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
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_BLAST_two_generative.pdf"), width = 10, height = 7.78)

# plt <- ggplot(data = ev.y1, aes(x = yrep)) + ylab("Density") + xlab("log(Burned Area)") + labs(col = "") +
#   geom_density(color = "steelblue", linewidth = 1.2) + 
#   geom_rug(alpha = 0.1) + 
#   xlim(5.5, 40) +
#   theme_minimal(base_size = 30) +  
#   theme(legend.position = "none",
#         axis.title = element_text(size = 30))

# d <- ggplot_build(plt)$data[[1]]
# print(plt + geom_area(data = subset(d, x>12.44009), aes(x=x,y=y), fill = "slategray1", alpha = 0.5) +
#         geom_segment(x=12.44009, xend=12.44009, 
#               y=0, yend=approx(x = d$x, y = d$y, xout = 12.4409)$y,
#               colour="red", linewidth=1.2, linetype = "dotted"))
# # ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_BLAST_generative.pdf"), width = 10, height = 7.78)

# library(ismev)
# gpd.fit(y-u, u)

# fit.log.lik <- extract_log_lik(fit1)
# constraint.elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 2)
# # save(constraint.elpd.loo, constraint.waic, file = (paste0("./BLAST/application/BLAST_constraint_",Sys.Date(),"_",psi,"_",floor(threshold*100),"quantile_IC.Rdata")))




# library(qrmtools)
# # library(Dowd)
# library(ReIns)
# hill_qrm <- Hill_estimator(Y[Y>1], k = c(5, 1000))
# Hill_plot(Y[Y>1], k = c(5, 1000), conf.level = 0.95, log="")

# tail_size <- length(Y) / 4 - 5 # The tail.size must be < n/4
# # pickands_result <- PickandsEstimator(Y, tail.size = tail_size)
# # PickandsPlot(pickands_result)
# # print(paste("Pickands Estimator:", pickands_result))
# moment_result <- Moment(Y[Y>1], plot = TRUE)
# moment_k <- moment_result$k
# moment_gamma <- moment_result$gamma
# stable_estimate <- median(moment_result$gamma[moment_result$k > 50])
# print(paste("DEdH Tail Index Estimate:", stable_estimate))

# H <- Hill(Y[Y>1], plot=FALSE)
# M <- Moment(Y[Y>1])
# gH <- genHill(Y[Y>1], gamma=H$gamma)
# # Plot estimates
# plot(H$k[5:500], M$gamma[5:500], xlab="k", ylab=expression(gamma), type="l")
# # lines(H$k[5:500], gH$gamma[5:500], lty=2)
# lines(H$k[5:500], gH$gamma[5:500], lty=3)
# legend("bottomright", c("Moment", "Generalised Hill"), lty=1:2)
