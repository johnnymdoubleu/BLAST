library(VGAM)
library(mgcv)
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
library(mboost)
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
psi.origin <- psi <- 30
threshold <- 0.95

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

# era5 <- read_excel("./BLAST/application/ERA_5.xlsx")
# era5 <- era5[era5$year>1979,]
# era5 <- era5[!(era5$year == 1999 & era5$month == 2 & era5$day == 14), ]
# fwi.index$ERA5 <- fwi.scaled$ERA5 <- as.numeric(era5$ERA_5)
fwi.scaled$time <- fwi.index$time <- c(1:length(Y))
fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)



above.0 <- which(Y>0)
Y <- Y[above.0]
fwi.scaled <- fwi.scaled[above.0,]
fwi.index <- fwi.index[above.0,]

fwi.season <- fwi.index %>%
  mutate(
    year = as.numeric(as.character(year)), 
    Month_Num = match(month, month.abb), 
    
    season = case_when(
      Month_Num %in% c(12, 1, 2) ~ "Winter",
      Month_Num %in% c(3, 4, 5) ~ "Spring",
      Month_Num %in% c(6, 7, 8) ~ "Summer",
      Month_Num %in% c(9, 10, 11) ~ "Autumn"
    ),

    SeasonYear = ifelse(Month_Num == 12, year + 1, year)
  )
fwi.season <- fwi.season %>% 
  mutate(Label = paste(SeasonYear, season)) %>%
  mutate(Label = factor(Label, levels = unique(Label)))

fwi.season$code <- factor(fwi.season$season, levels=c("Winter", "Spring", "Summer", "Autumn"), labels = c(1,2,3,4))
fwi.season$winter <- ifelse(fwi.season$season == "Winter", 1, 0)
fwi.season$spring <- ifelse(fwi.season$season == "Spring", 1, 0)
fwi.season$summer <- ifelse(fwi.season$season == "Summer", 1, 0)
fwi.season$autumn <- ifelse(fwi.season$season == "Autumn", 1, 0)
fwi.season$BA <- Y
fwi.season$log.BA <- log(fwi.season$BA)
fwi.season$cos.time <- cos(2*pi*fwi.season$time / 365)
fwi.season$sin.time <- sin(2*pi*fwi.season$time / 365)

fwi.scaled <- fwi.df <- fwi.season[,c(3,4,5,6,7,16)]
# save(fwi.season, fwi.scaled, file="./BLAST/application/wildfire_season.Rdata")
load(paste0("./BLAST/application/qr-",threshold*1000,"-t.Rdata"))

# u <- quantile(Y, threshold)
# excess <- which(Y>u)
preds <- preds.qgam
excess <- which(Y>preds)
u <- preds[excess]

y <- Y[excess]
# u <- rep(quantile(Y, threshold),length(y))
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
n <- dim(fwi.df[excess,])[[1]]
p <- dim(fwi.df[excess,])[[2]]-1

fwi.scaled <- fwi.scaled[excess,]
fwi.scaled[,c(1:p)] <- as.data.frame(sapply(fwi.scaled[,c(1:p)], FUN = range01))

fwi.origin <- data.frame(fwi.df[excess,], BA=y)
max.fwi <- fwi.origin[which.max(y),]
fwi.grid <- data.frame(lapply(fwi.origin[,c(1:p)], function(x) seq(min(x), max(x), length.out = nrow(fwi.scaled))))
fwi.minmax <- sapply(fwi.origin[,c(1:p)], function(x) max(x)-min(x))
fwi.min <- sapply(fwi.origin[,c(1:p)], function(x) min(x))

psi <- psi - 2
model_data <- list()       # Stores parameters (U, Lambda, Coefs) for later
train_basis_flat <- list() # Stores actual Training Matrices for cbind
grid_basis_flat <- list()  # Stores actual Grid Matrices for cbind
season_code_full <- as.factor(fwi.scaled[, 6]) 
covariates <- colnames(fwi.scaled[, c(1:p)])

# Initialize storage
fwi_basis_list <- list()
linear_basis_list <- list()
nonlinear_basis_list <- list()
model_data <- list() # <--- NEW: Stores params for the grid

for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  x_vec <- fwi.scaled[, i]
  x_true <- fwi.origin[,i]
  # Initialize sub-list for this variable
  model_data[[var_name]] <- list() 
  
  sm_list <- smoothCon(mgcv::s(x_vec, by=season_code_full, bs = "tp", k = psi + 2), 
                       data = data.frame(x_vec = x_vec, season_code_full = season_code_full), 
                       knots = NULL)
  
  for (j in 1:length(sm_list)) {
    season_label <- levels(season_code_full)[j]
    sm_spec <- sm_list[[j]] 
    
    mask <- as.numeric(season_code_full == season_label)
    
    # 1. Linear Basis Construction
    vec_intercept <- mask
    vec_slope <- x_vec * mask
    
    lin_col_name <- paste(var_name, "S", season_label, sep="_")
    linear_basis_list[[lin_col_name]] <- matrix(vec_slope, ncol=1, dimnames=list(NULL, lin_col_name))

    fwi_basis_list[[lin_col_name]] <- matrix(x_true*mask, ncol=1, dimnames=list(NULL, lin_col_name))
    
    # Define Null Space
    X_lin_local <- cbind(vec_intercept, vec_slope)
    
    qr_lin <- qr(X_lin_local)
    Q_lin <- qr.Q(qr_lin)
    R_lin <- qr.R(qr_lin)

    # 2. Spectral Decomposition
    X_raw <- sm_spec$X
    S     <- sm_spec$S[[1]] 
    
    eig <- eigen(S, symmetric = TRUE)
    pos_idx <- which(eig$values > max(eig$values) * 1e-8)
    
    if(length(pos_idx) == 0) next 
    
    U_pen <- eig$vectors[, pos_idx]       
    Lambda_pen <- diag(eig$values[pos_idx]) 
    Lambda_sqrt_inv <- solve(sqrt(Lambda_pen))
    
    Z_spectral <- X_raw %*% U_pen %*% Lambda_sqrt_inv
    
    # 3. Orthogonalization
    Gamma_Q <- t(Q_lin) %*% Z_spectral
    Gamma_Original <- backsolve(R_lin, Gamma_Q)
    Z_orth <- Z_spectral - X_lin_local %*% Gamma_Original
    
    # 4. Cleanup & Store Training Basis
    keep_cols <- colSums(Z_orth^2) > 1e-9
    Z_final <- Z_orth[, keep_cols, drop = FALSE]
    
    train_scale <- apply(Z_final, 2, sd)
    train_scale[train_scale < 1e-12] <- 1 
    Z_final <- scale(Z_final, center = FALSE, scale = train_scale)
    Z_final <- Z_final * mask
    
    col_names <- paste(var_name, "S", season_label, 1:ncol(Z_final), sep="_")
    colnames(Z_final) <- col_names
    nonlinear_basis_list[[paste(var_name, season_label, sep="_")]] <- Z_final
    
    # --- NEW: SAVE PARAMETERS FOR GRID ---
    model_data[[var_name]][[season_label]] <- list(
      sm_spec = sm_spec, 
      U_pen = U_pen, 
      Lambda_sqrt_inv = Lambda_sqrt_inv, 
      projection_coefs = Gamma_Original, 
      keep_cols = keep_cols, 
      scale_stats = train_scale
    )
  }
}

# Create Final Training Matrices
true.linear <- do.call(cbind, fwi_basis_list)
fwi.linear <- bs.linear <- do.call(cbind, linear_basis_list)
bs.linear <- cbind(rep(1,n), bs.linear)
bs.nonlinear <- do.call(cbind, nonlinear_basis_list)

cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,2,3,4,5)]) %*% bs.nonlinear[,c((1):(4*psi))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,6,7,8,9)]) %*% bs.nonlinear[,c((4*(psi+1)):(4*(psi*2)))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,10,11,12,13)]) %*% bs.nonlinear[,c((4*(psi*2+1)):(4*(psi*3)))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,14,15,16,17)]) %*% bs.nonlinear[,c((4*(psi*3+1)):(4*(psi*4)))]), "\n")
cat("Orthogonality Check (Linear vs Nonlinear):", sum(t(bs.linear[,c(1,18,19,20,21)]) %*% bs.nonlinear[,c((4*(psi*4+1)):(4*(psi*5)))]), "\n")


newx <- seq(0, 1, length.out = n)
xholder <- do.call(cbind, lapply(1:p, function(j) {newx}))
colnames(xholder) <- covariates

# Note: 'n' must match the row dimension expected by your Stan data block.
# We no longer define a global 'newx' here because every season gets its own grid.

grid_linear_list <- list()
grid_nonlinear_list <- list()

for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  
  # 1. Retrieve Global Constants for this variable
  # (Required to scale the seasonal grid back to the model's 0-1 scale)
  global_min <- fwi.min[i]
  global_range <- fwi.minmax[i]
  
  if (is.null(model_data[[var_name]])) next
  
  for (season_label in names(model_data[[var_name]])) {
    params <- model_data[[var_name]][[season_label]]
    lin_col_name <- paste(var_name, "S", season_label, sep="_")
    
    # --- NEW LOGIC START ---
    
    # A. Find the indices for this specific season
    season_idx <- which(season_code_full == season_label)
    
    # B. Extract the RAW data for this season (to find observed bounds)
    raw_season_vals <- fwi.origin[season_idx, i]
    season_min <- min(raw_season_vals, na.rm = TRUE)
    season_max <- max(raw_season_vals, na.rm = TRUE)
    
    # C. Create a grid covering ONLY the observed seasonal range
    # We create a sequence of length 'n' between the seasonal min and max
    raw_grid <- seq(season_min, season_max, length.out = n)
    
    # D. Scale this grid to the Global [0, 1] domain
    # This is crucial: The model coefficients 'theta' and 'gamma' were trained 
    # on data scaled globally. We must provide inputs on that same scale.
    grid_vals <- (raw_grid - global_min) / global_range
    
    # --- NEW LOGIC END ---

    # Store the grid values (on 0-1 scale) for plotting x-axis later
    grid_linear_list[[lin_col_name]] <- matrix(grid_vals, ncol=1, dimnames=list(NULL, lin_col_name))
    
    pred_df <- data.frame(
      x_vec = grid_vals,
      season_code_full = factor(season_label, levels = levels(season_code_full))
    )
    
    # Projection (unchanged)
    X_raw_grid <- PredictMat(params$sm_spec, pred_df)
    Z_spectral_grid <- X_raw_grid %*% params$U_pen %*% params$Lambda_sqrt_inv
    intercept_grid <- rep(1, length(grid_vals))
    slope_grid     <- grid_vals
    X_lin_grid_local <- cbind(intercept_grid, slope_grid)
    Z_orth_grid <- Z_spectral_grid - X_lin_grid_local %*% params$projection_coefs
    Z_final_grid <- Z_orth_grid[, params$keep_cols, drop = FALSE]
    Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)
    col_names <- paste(var_name, "S", season_label, 1:ncol(Z_final_grid), sep="_")
    colnames(Z_final_grid) <- col_names
    grid_nonlinear_list[[paste(var_name, season_label, sep="_")]] <- Z_final_grid
  }
}

xholder.linear <- do.call(cbind, grid_linear_list)
xholder.nonlinear <- do.call(cbind, grid_nonlinear_list)

# grid_linear_list <- list()
# grid_nonlinear_list <- list()

# for (i in seq_along(covariates)) {
#   var_name <- covariates[i]
#   grid_vals <- xholder[, i]
  
#   if (is.null(model_data[[var_name]])) next
#   for (season_label in names(model_data[[var_name]])) {
#     params <- model_data[[var_name]][[season_label]]
#     lin_col_name <- paste(var_name, "S", season_label, sep="_")
#     grid_linear_list[[lin_col_name]] <- matrix(grid_vals, ncol=1, dimnames=list(NULL, lin_col_name))
#     pred_df <- data.frame(
#       x_vec = grid_vals,
#       season_code_full = factor(season_label, levels = levels(season_code_full))
#     )
    
#     X_raw_grid <- PredictMat(params$sm_spec, pred_df)
#     Z_spectral_grid <- X_raw_grid %*% params$U_pen %*% params$Lambda_sqrt_inv
#     intercept_grid <- rep(1, length(grid_vals))
#     slope_grid     <- grid_vals
#     X_lin_grid_local <- cbind(intercept_grid, slope_grid)
#     Z_orth_grid <- Z_spectral_grid - X_lin_grid_local %*% params$projection_coefs
#     Z_final_grid <- Z_orth_grid[, params$keep_cols, drop = FALSE]
#     Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)
#     col_names <- paste(var_name, "S", season_label, 1:ncol(Z_final_grid), sep="_")
#     colnames(Z_final_grid) <- col_names
#     grid_nonlinear_list[[paste(var_name, season_label, sep="_")]] <- Z_final_grid
#   }
# }

# xholder.linear <- do.call(cbind, grid_linear_list)
# xholder.nonlinear <- do.call(cbind, grid_nonlinear_list)

scale_stats_list <- list()
for (var_name in names(model_data)) {
  for (season_label in names(model_data[[var_name]])) {
    stats <- model_data[[var_name]][[season_label]]$scale_stats
    scale_stats_list[[paste(var_name, "S", season_label, sep="_")]] <- stats
  }
}

Z_scales <- unlist(scale_stats_list)
grid_nonlinear_flat <- list()
grid_linear_flat <- list()

for (i in seq_along(covariates)) {
  var_name <- covariates[i]
  
  if (is.null(model_data[[var_name]])) next
  
  # 1. Global Constants (for scaling back to model units)
  global_min <- fwi.min[i]
  global_range <- fwi.minmax[i]
  
  # Determine how many points you want in your grid (matching original fwi.grid size)
  grid_len <- nrow(fwi.grid)
  
  # Note: We do NOT calculate x_vals_raw here anymore because it must be season-specific.

  for (season_label in names(model_data[[var_name]])) {
    params <- model_data[[var_name]][[season_label]]
    lin_col_name <- paste(var_name, "S", season_label, sep="_")
    
    # --- NEW LOGIC START ---
    
    # A. Identify the season indices in the ORIGINAL data
    season_idx <- which(season_code_full == season_label)
    
    # B. Get the observed range for this specific season
    raw_season_vals <- fwi.origin[season_idx, i]
    season_min <- min(raw_season_vals, na.rm = TRUE)
    season_max <- max(raw_season_vals, na.rm = TRUE)
    
    # C. Generate a sequence specifically for this season
    x_vals_raw <- seq(season_min, season_max, length.out = grid_len)
    
    # D. Scale to Global [0, 1] (Model Space)
    x_vals_scaled <- (x_vals_raw - global_min) / global_range
    
    # --- NEW LOGIC END ---
    
    # 2. Re-create the Linear Basis Matrix for Orthogonalization 
    # (Must be done INSIDE the loop now, as x_vals_scaled changes)
    X_lin_grid_local <- cbind(rep(1, length(x_vals_scaled)), x_vals_scaled)

    # 3. Store Linear Part (Raw Value - Global Min)
    # This aligns with your 'sweep' logic later
    grid_linear_flat[[lin_col_name]] <- matrix((x_vals_raw - global_min), ncol=1, dimnames=list(NULL, lin_col_name))
    
    grid_df <- data.frame(
      x_vec = x_vals_scaled,
      season_code_full = factor(season_label, levels = levels(season_code_full))
    )
    
    # 4. Projection & Orthogonalization
    X_raw_grid <- PredictMat(params$sm_spec, grid_df)
    Z_spectral_grid <- X_raw_grid %*% params$U_pen %*% params$Lambda_sqrt_inv
    
    # Use the local linear matrix for projection
    Z_orth_grid <- Z_spectral_grid - X_lin_grid_local %*% params$projection_coefs
    
    Z_final_grid <- Z_orth_grid[, params$keep_cols, drop = FALSE]
    Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)
    
    col_names <- paste(var_name, "S", season_label, 1:ncol(Z_final_grid), sep="_")
    colnames(Z_final_grid) <- col_names
    grid_nonlinear_flat[[paste(var_name, season_label, sep="_")]] <- Z_final_grid
  }
}

xgrid.linear <- do.call(cbind, grid_linear_flat)
fwi.true <- sweep(xgrid.linear, 2, rep(fwi.min, each=4), "+")
xgrid.nonlinear <- do.call(cbind, grid_nonlinear_flat)
# grid_nonlinear_flat <- list()
# grid_linear_flat <- list()
# for (i in seq_along(covariates)) {
#   var_name <- covariates[i]
  
#   if (is.null(model_data[[var_name]])) next
#   x_vals_raw <- fwi.grid[, i]
#   x_vals_scaled <- (x_vals_raw - fwi.min[i]) / fwi.minmax[i]
#   X_lin_grid_local <- cbind(rep(1, length(x_vals_scaled)), x_vals_scaled)
#   for (season_label in names(model_data[[var_name]])) {
#     params <- model_data[[var_name]][[season_label]]
#     lin_col_name <- paste(var_name, "S", season_label, sep="_")
#     grid_linear_flat[[lin_col_name]] <- matrix((x_vals_raw - fwi.min[i]), ncol=1, dimnames=list(NULL, lin_col_name))
    
#     grid_df <- data.frame(
#       x_vec = x_vals_scaled,
#       season_code_full = factor(season_label, levels = levels(season_code_full))
#     )
    
#     X_raw_grid <- PredictMat(params$sm_spec, grid_df)
#     Z_spectral_grid <- X_raw_grid %*% params$U_pen %*% params$Lambda_sqrt_inv
#     Z_orth_grid <- Z_spectral_grid - X_lin_grid_local %*% params$projection_coefs
#     Z_final_grid <- Z_orth_grid[, params$keep_cols, drop = FALSE]
#     Z_final_grid <- scale(Z_final_grid, center = FALSE, scale = params$scale_stats)
#     col_names <- paste(var_name, "S", season_label, 1:ncol(Z_final_grid), sep="_")
#     colnames(Z_final_grid) <- col_names
#     grid_nonlinear_flat[[paste(var_name, season_label, sep="_")]] <- Z_final_grid
#   }
# }

# xgrid.linear <- do.call(cbind, grid_linear_flat)
# fwi.true <- sweep(xgrid.linear, 2, rep(fwi.min, each=4), "+")
# xgrid.nonlinear <- do.call(cbind, grid_nonlinear_flat)

X_sd   <- apply(bs.linear[,-1], 2, sd)
bs.linear[,-1] <- scale(bs.linear[,-1], center = FALSE, scale = X_sd)

model.stan <- "// Stan model for BLAST Pareto Samples
data {
  int <lower=1> n; // Sample size
  int <lower=1> p; // regression coefficient size
  int <lower=1> psi; // splines coefficient size
  int <lower=1> n_seasons;
  matrix[n, (p*n_seasons)] bsLinear; // fwi dataset
  matrix[n, (psi*p*n_seasons)] bsNonlinear; // thin plate splines basis
  matrix[n, (p*n_seasons)] xholderLinear; // fwi dataset
  matrix[n, (psi*p*n_seasons)] xholderNonlinear; // thin plate splines basis    
  matrix[n, (p*n_seasons)] gridL; // fwi dataset
  matrix[n, (psi*p*n_seasons)] gridNL; // thin plate splines basis
  vector[n] u; // large threshold value
  vector[n] y; // extreme response
  real <lower=0> atau;

  vector[(psi*p*n_seasons)] Z_scales;
  vector[p * n_seasons] X_minmax;
  vector[p * n_seasons] X_sd;
}

parameters {
  real theta0;
  vector[p * n_seasons] theta; // linear predictor
  array[p, n_seasons] vector[psi] gamma_raw;
  array[p, n_seasons] real <lower=0> lambda1; // lasso penalty //array[p] real <lower=0> 
  array[p, n_seasons] real <lower=0> lambda2; // lambda2 group lasso penalty
  array[p, n_seasons] real <lower=0> tau;
}

transformed parameters {
  vector[n] alpha;        // Tail index
  array[p, n_seasons] vector[psi] gamma; // Scaled spline coefficients
  
  {
    vector[n] lin_predictor;
    lin_predictor = rep_vector(theta0, n); 
    for (j in 1:p){
      for (s in 1:n_seasons){
        int lin_idx = (j-1) * n_seasons + s;
        int nl_start = (lin_idx-1) * psi + 1; 
        int nl_end   = lin_idx * psi;
        
        for (k in 1:psi){
           int flat_idx = nl_start + k - 1;
           gamma[j,s][k] = gamma_raw[j,s][k] * sqrt(tau[j,s]) * Z_scales[flat_idx];
        }
        
        lin_predictor += col(bsLinear, lin_idx) * theta[lin_idx] 
                       + block(bsNonlinear, 1, nl_start, n, psi) * gamma[j,s];
      }
    }
    alpha = exp(lin_predictor);
  }
}

model {
  // likelihood
  target += pareto_lpdf(y | u, alpha);
  target += normal_lpdf(theta0 | 0, 10);
  for (j in 1:p){
    for (s in 1:n_seasons){
      int idx = (j-1) * n_seasons + s;
      target += gamma_lpdf(lambda1[j,s] | 1, 1); 
      target += gamma_lpdf(lambda2[j,s] | 1e-2, 1e-2);
      target += double_exponential_lpdf(theta[idx] | 0, 1 / lambda1[j,s]);
      target += gamma_lpdf(tau[j,s] | atau, square(lambda2[j,s]) * 0.5);
      target += std_normal_lpdf(gamma_raw[j,s]);
    }
  }
}

generated quantities {
  matrix[n, n_seasons] alphaseason;
  vector[n] log_lik;
  matrix[n, p * n_seasons] gridgl;  // Linear component
  matrix[n, p * n_seasons] gridgnl; // Non-linear component
  matrix[n, p * n_seasons] gridgsmooth; // Non-linear component
  matrix[n, p * n_seasons] fwismooth; // Total smooth effect

  vector[p * n_seasons] theta_origin = theta ./ X_sd;
  vector[p * n_seasons] theta_fwi = theta_origin ./ X_minmax;

  for (j in 1:p){
    for (s in 1:n_seasons){
      int lin_idx = (j-1) * n_seasons + s;
      int nl_start = (lin_idx-1) * psi + 1;
      
      gridgl[, lin_idx] = col(xholderLinear, lin_idx) * theta_origin[lin_idx];
      gridgnl[, lin_idx] = block(xholderNonlinear, 1, nl_start, n, psi) * gamma[j,s];

      fwismooth[, lin_idx] = block(gridNL, 1, nl_start, n, psi) * gamma[j,s] 
                            + col(gridL, lin_idx) * theta_fwi[lin_idx]; 
    }
  }
  gridgsmooth = gridgl + gridgnl;
  
  {
    for (s in 1:n_seasons) {
      vector[n] current_season_pred = rep_vector(theta0, n);     
      for (j in 1:p) {
        int lin_idx = (j-1) * n_seasons + s;
        current_season_pred += gridgsmooth[, lin_idx];
      }
      alphaseason[, s] = exp(current_season_pred);
    }
  }
  
  for (i in 1:n){
    log_lik[i] = pareto_lpdf(y[i] | u[i], alpha[i]);
  }
}
"

# {
#   vector[n] grid_lin_pred = rep_vector(theta0, n);
#   for(k in 1:(p * n_seasons)) {
#     gridgsmooth[,k] = gridgl[, k] + gridgnl[, k];
#     grid_lin_pred += gridgsmooth[,k];
#   }
#   gridalpha = exp(grid_lin_pred);
# }

# data.stan <- list(y = as.vector(y), u = u, p = p, n= n, psi = psi, Z_scales= Z_scales,
#                   atau = ((psi+1)/2), xholderNonlinear = xholder.nonlinear, 
#                   bsLinear = bs.linear[,-1], bsNonlinear = bs.nonlinear, X_min=fwi.min,
#                   xholderLinear = xholder.linear[,-1], X_minmax = fwi.minmax, 
#                   X_means = X_means, X_sd = X_sd,
#                   gridL = xgrid.linear[,-1], gridNL = xgrid.nonlinear)

# init.alpha <- list(list(gamma_raw= array(rep(0.2, (psi*p)), dim=c(p, psi)), #rho = 1, 
#                         theta = rep(-0.1, (p+1)), tau = rep(0.1, p), #sigma_ridge = rep(0.1, p),
#                         lambda1 = rep(0.1,p), lambda2 = rep(1, p)),
#                    list(gamma_raw = array(rep(-0.15, (psi*p)), dim=c(p, psi)),# rho = 1,
#                         theta = rep(0.05, (p+1)), tau = rep(0.2, p), #sigma_ridge = rep(0.2, p), 
#                         lambda1 = rep(2,p), lambda2 = rep(5, p)),
#                    list(gamma_raw= array(rep(0.1, (psi*p)), dim=c(p, psi)), #rho = 1,
#                         theta = rep(0.1, (p+1)), tau = rep(0.1, p), #sigma_ridge = rep(0.1, p),
#                         lambda1 = rep(0.1,p), lambda2 = rep(0.1, p)))

data.stan <- list(
  y = as.vector(y), 
  u = u, 
  n = n, 
  p = p, 
  n_seasons = 4,               # <--- NEW
  psi = psi, 
  bsLinear = bs.linear[,-1], 
  bsNonlinear = bs.nonlinear,
  xholderLinear = xholder.linear, 
  xholderNonlinear = xholder.nonlinear, 
  gridL = xgrid.linear,
  gridNL = xgrid.nonlinear,
  atau = ((psi+1)/2), 
  Z_scales = Z_scales,
  X_minmax = rep(fwi.minmax, each = 4),
  X_sd = X_sd
)

init.alpha <- list(
  list(
    theta0 = 0.1,                                  
    theta = rep(-0.1, p * 4),
    gamma_raw = array(rep(0.1, p * 4 * psi), dim = c(p, 4, psi)), 
    tau = array(rep(0.1, p * 4), dim = c(p, 4)),
    lambda1 = array(rep(0.1, p * 4), dim = c(p, 4)), 
    lambda2 = array(rep(1, p * 4), dim = c(p, 4))
  ),
  list(
    theta0 = 0.05,
    theta = rep(0.05, p * 4),
    gamma_raw = array(rep(-0.1, p * 4 * psi), dim = c(p, 4, psi)),
    tau = array(rep(0.2, p * 4), dim = c(p, 4)),
    lambda1 = array(rep(2, p * 4), dim = c(p, 4)), 
    lambda2 = array(rep(5, p * 4), dim = c(p, 4))
  ),
  list(
    theta0 = -0.1,
    theta = rep(0.1, p * 4),
    gamma_raw = array(rep(0.05, p * 4 * psi), dim = c(p, 4, psi)),
    tau = array(rep(0.1, p * 4), dim = c(p, 4)),
    lambda1 = array(rep(0.1, p * 4), dim = c(p, 4)), 
    lambda2 = array(rep(0.1, p * 4), dim = c(p, 4))
  )
)

fit1 <- stan(
    model_code = model.stan,
    model_name = "BLAST",
    data = data.stan,    # named list of data
    init = init.alpha,      # initial value 
    chains = 3,             # number of Markov chains
    iter = 3000,            # total number of iterations per chain
    cores = parallel::detectCores(), # number of cores (could use one per chain)
    refresh = 1500           # no progress shown
)

posterior <- rstan::extract(fit1)
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit1, pars="lp__") + ylab("") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text = element_text(size = 18))

theta.samples <- summary(fit1, par=c("theta0", "theta"), probs = c(0.05,0.5, 0.95))$summary
gamma.samples <- summary(fit1, par=c("gamma"), probs = c(0.05,0.5, 0.95))$summary
lambda.samples <- summary(fit1, par=c("lambda1", "lambda2"), probs = c(0.05,0.5, 0.95))$summary
gsmooth.samples <- summary(fit1, par=c("gridgsmooth"), probs = c(0.05, 0.5, 0.95))$summary
fwismooth.samples <- summary(fit1, par=c("fwismooth"), probs = c(0.05, 0.5, 0.95))$summary
gridgl.samples <- summary(fit1, par=c("gridgl"), probs = c(0.05, 0.5, 0.95))$summary
gridgnl.samples <- summary(fit1, par=c("gridgnl"), probs = c(0.05, 0.5, 0.95))$summary
origin.samples <- summary(fit1, par=c("alpha"), probs = c(0.05,0.5, 0.95))$summary
# alpha.samples <- summary(fit1, par=c("gridalpha"), probs = c(0.05,0.5, 0.95))$summary
season.samples <- summary(fit1, par=c("alphaseason"), probs = c(0.05,0.5, 0.95))$summary
loglik.samples <- summary(fit1, par=c("log_lik"), probs = c(0.05,0.5, 0.95))$summary

# MCMCvis::MCMCplot(fit1, params = 'theta')
# MCMCvis::MCMCplot(fit1, params = "gamma")

g.smooth.mean <- as.vector(matrix(gsmooth.samples[,1], nrow = n, byrow=TRUE))
g.smooth.q1 <- as.vector(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.smooth.q2 <- as.vector(matrix(gsmooth.samples[,5], nrow = n, byrow=TRUE))
g.smooth.q3 <- as.vector(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))
g.linear.mean <- as.vector(matrix(gridgl.samples[,1], nrow = n, byrow=TRUE))
g.linear.q1 <- as.vector(matrix(gridgl.samples[,4], nrow = n, byrow=TRUE))
g.linear.q2 <- as.vector(matrix(gridgl.samples[,5], nrow = n, byrow=TRUE))
g.linear.q3 <- as.vector(matrix(gridgl.samples[,6], nrow = n, byrow=TRUE))
g.nonlinear.mean <- as.vector(matrix(gridgnl.samples[,1], nrow = n, byrow=TRUE))
g.nonlinear.q1 <- as.vector(matrix(gridgnl.samples[,4], nrow = n, byrow=TRUE))
g.nonlinear.q2 <- as.vector(matrix(gridgnl.samples[,5], nrow = n, byrow=TRUE))
g.nonlinear.q3 <- as.vector(matrix(gridgnl.samples[,6], nrow = n, byrow=TRUE))
fwi.smooth.mean <- as.vector(matrix(fwismooth.samples[,1], nrow = n, byrow=TRUE))
fwi.smooth.q1 <- as.vector(matrix(fwismooth.samples[,4], nrow = n, byrow=TRUE))
fwi.smooth.q2 <- as.vector(matrix(fwismooth.samples[,5], nrow = n, byrow=TRUE))
fwi.smooth.q3 <- as.vector(matrix(fwismooth.samples[,6], nrow = n, byrow=TRUE))

g.min <- as.data.frame(matrix(gsmooth.samples[,4], nrow = n, byrow=TRUE))
g.min.samples <- sapply(g.min, min)
g.max <- as.data.frame(matrix(gsmooth.samples[,6], nrow = n, byrow=TRUE))
g.max.samples <- sapply(g.max, max)

fwi.smooth.min <- as.data.frame(matrix(fwismooth.samples[,4], nrow = n, byrow=TRUE))
fwi.min.samples <- sapply(fwi.smooth.min, min)
fwi.smooth.max <- as.data.frame(matrix(fwismooth.samples[,6], nrow = n, byrow=TRUE))
fwi.max.samples <- sapply(fwi.smooth.max, max)

# data.smooth <- data.frame("x"= as.vector(xholder),
#                           "post.mean" = as.vector(g.smooth.mean),
#                           "q1" = as.vector(g.smooth.q1),
#                           "q2" = as.vector(g.smooth.q2),
#                           "q3" = as.vector(g.smooth.q3),
#                           "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)),
#                           "replicate" = gl(p, n, (p*n), labels = names(fwi.scaled)))

# ggplot(data.smooth, aes(x=x, group=interaction(covariates, replicate))) + 
#   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#   ylab("") + xlab("") +
#   facet_grid(covariates ~ ., scales = "free",
#               labeller = label_parsed) + 
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values=c("steelblue")) + 
#   guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) + 
#           # ylim(-3.5, 3.5) +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#           plot.margin = margin(0,0,0,-20),
#           # strip.text = element_blank(),
#           axis.text = element_text(size = 20))
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_smooth.pdf"), width=12.5, height = 15)

# xholder.evgam <- data.frame(xholder, fwi.season[excess,c(17,18,19,20)])
# simul.data <- data.frame(BA = y-u, fwi.scaled[,c(1:p)], 
                          # fwi.season[excess,c(17,18,19,20)])
xholder.evgam <- data.frame(xholder, code = fwi.season[excess,16])
simul.data <- data.frame(BA = y-u, fwi.scaled[,c(1:p)], 
                          code = fwi.season[excess,16])

gam.scale <- list(BA ~ s(BUI, by=code, bs = "ts", k = 30) + 
                      s(ISI, by=code, bs = "ts", k = 30) + 
                      s(FFMC, by=code, bs = "ts", k = 30) +
                      s(DMC, by=code, bs = "ts", k = 30) + 
                      s(DC, by=code, bs = "ts", k = 30),
                    ~ s(BUI, by=code, bs = "ts", k = 30) + 
                      s(ISI, by=code, bs = "ts", k = 30) + 
                      s(FFMC, by=code, bs = "ts", k = 30) +
                      s(DMC, by=code, bs = "ts", k = 30) +
                      s(DC, by=code, bs = "ts", k = 30))
# evgam.fit.scale <- evgam::evgam(gam.scale, data = simul.data, family = "gpd")
# pred.scale <- predict(evgam.fit.scale, newdata = xholder.evgam, type="response", se.fit = TRUE)
# xi.pred.scale <-pred.scale$fitted$shape
# xi.se.scale <- pred.scale$se.fit$shape
# xi.low.scale <- xi.pred.scale - (1.96 * xi.se.scale)
# xi.high.scale <- xi.pred.scale + (1.96 * xi.se.scale)
# alpha.pred.scale <- 1/xi.pred.scale
# xholder.basis.scale <- predict(evgam.fit.scale, newdata = xholder, type= "lpmatrix")$shape
# psi <- 15
# xi.coef.scale <- tail(evgam.fit.scale$coefficients, (psi-1)*p)
# gamma.xi.scale <- matrix(xi.coef.scale, ncol = p)
# alpha.nonlinear.scale <- xi.nonlinear.scale <- matrix(, nrow = n, ncol = p)
# bs.nonlinear.scale <- xholder.basis.scale[,c(2:((psi-1)*p+1))]
# for(j in 1:p){
#   xi.nonlinear.scale[,j] <- bs.nonlinear.scale[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.scale[,j]
#   alpha.nonlinear.scale[,j] <- 1/(xi.nonlinear.scale[,j])
# }

n_grid <- nrow(bs.linear)
season_levels <- levels(simul.data$code) 
covariate_names <- c("BUI", "ISI", "FFMC", "DMC", "DC")
plot_data_list <- list()

for(var in covariate_names) {
  var_seq <- seq(min(simul.data[[var]]), max(simul.data[[var]]), length.out = n_grid)
  base_df <- simul.data %>% 
    summarise(across(all_of(covariate_names), median)) %>% 
    slice(rep(1, n_grid))
  base_df[[var]] <- var_seq
  expanded_df <- expand_grid(
    base_df, 
    code = factor(season_levels, levels = season_levels)
  )
  preds <- predict(evgam.fit.scale, newdata = expanded_df, type = "response")
  expanded_df$fitted_xi <- preds$fitted$shape
  expanded_df$covariate_name <- var
  expanded_df$x_value <- expanded_df[[var]] # Generic x-column for faceting
  
  plot_data_list[[var]] <- expanded_df
}

# Combine all into one big plotting dataframe
final_plot_data <- bind_rows(plot_data_list)


gam.1 <- list(BA ~ 1,
                ~ s(BUI, by=code, bs = "ts", k = 30) + 
                  s(ISI, by=code, bs = "ts", k = 30) + 
                  s(FFMC, by=code, bs = "ts", k = 30) +
                  s(DMC, by=code, bs = "ts", k = 30) +
                  s(DC, by=code, bs = "ts", k = 30) + 
                  winter + spring + summer + autumn)
evgam.fit.1 <- evgam::evgam(gam.1, data = simul.data, family = "gpd")
pred.1 <- predict(evgam.fit.1, newdata = xholder.evgam, type="response", se.fit = TRUE)
xi.pred.1 <-pred.1$fitted$shape
# xi.se.1 <- pred.1$se.fit$shape
# xi.low.1 <- xi.pred.1 - (1.96 * xi.se.1)
# xi.high.1 <- xi.pred.1 + (1.96 * xi.se.1)
# alpha.pred.1 <- 1/xi.pred.1



# xholder.basis.1 <- predict(evgam.fit.1, newdata = xholder, type= "lpmatrix")$shape
# xi.coef.1 <- tail(evgam.fit.1$coefficients, (psi-1)*p)
# gamma.xi.1 <- matrix(xi.coef.1, ncol = p)
# alpha.nonlinear.1 <- xi.nonlinear.1 <- matrix(, nrow = n, ncol = p)
# bs.nonlinear.1 <- xholder.basis.1[,c(2:((psi-1)*p+1))]
# for(j in 1:p){
#   xi.nonlinear.1[,j] <- bs.nonlinear.1[,(((j-1)*(psi-1))+1):(((j-1)*(psi-1))+(psi-1))] %*% gamma.xi.1[,j]
#   alpha.nonlinear.1[,j] <- 1/(xi.nonlinear.1[,j])
# }
# simul.data <- data.frame(BA = y, fwi.scaled[,c(1:p)])
# vgam.fit.scale <- vgam(BA ~ sm.ps(BUI, ps.int = 28) + 
#                            sm.ps(ISI, ps.int = 28) + 
#                            sm.ps(FFMC, ps.int = 28) + 
#                            sm.ps(DMC, ps.int = 28) + 
#                            sm.ps(DC, ps.int = 28),
#                         data = simul.data,
#                         family = gpd(threshold= u,
#                                       lshape="loglink",
#                                       zero = NULL),
#                         trace = TRUE,
#                         control = vgam.control(maxit = 200))
# fitted.linear <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "link")
# fitted.terms <- predict(vgam.fit.scale, newdata = data.frame(xholder), type = "terms")
# vgam.xi.scale <- exp(fitted.linear[,2])
# vgam.sigma.scale <- exp(fitted.linear[,1])

# vgam.fit.1 <- vgam(BA ~ sm.ps(BUI, ps.int = 28) + 
#                        sm.ps(ISI, ps.int = 28) + 
#                        sm.ps(FFMC, ps.int = 28) + 
#                        sm.ps(DMC, ps.int = 28) + 
#                        sm.ps(DC, ps.int = 28),
#                       data = simul.data,
#                       family = gpd(threshold= u,
#                                     lshape="loglink",
#                                     zero = 1),
#                       trace = TRUE,
#                       control = vgam.control(maxit = 200))
# fitted.linear <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "link")
# fitted.terms <- predict(vgam.fit.1, newdata = data.frame(xholder), type = "terms")
# vgam.xi.1 <- exp(fitted.linear[,2])
# vgam.sigma.1 <- exp(fitted.linear[,1])


alpha.mean <- as.vector(matrix(season.samples[,1], nrow = n, byrow=TRUE))
alpha.q1 <- as.vector(matrix(season.samples[,4], nrow = n, byrow=TRUE))
alpha.q2 <- as.vector(matrix(season.samples[,5], nrow = n, byrow=TRUE))
alpha.q3 <- as.vector(matrix(season.samples[,6], nrow = n, byrow=TRUE))

data.alpha <- data.frame("x" = newx,
                          "post.mean" = alpha.mean,
                          "q1" = alpha.q1,
                          "q2" = alpha.q2,
                          "q3" = alpha.q3)
seasons <- c("Winter", "Spring", "Summer", "Autumn")

grid.plts <- list()
for(i in 1:4){
  grid.plt <- ggplot(data = data.frame(data.alpha[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
                  geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
                  geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
                  # geom_rug(aes(x=true, y=q2), sides = "b") +
                  ylab(expression(alpha(c,...,c))) + xlab(seasons[i]) +
                  scale_fill_manual(values=c("steelblue"), name = "") + 
                  scale_color_manual(values=c("steelblue")) +
                  theme_minimal(base_size = 20) +
                  theme(legend.position = "none",
                        plot.margin = margin(5, 5, 5, 5),
                        plot.title = element_text(hjust = 0.5, face = "bold"),
                        axis.text = element_text(size = 18),
                        axis.title.x = element_text(size = 22))
  grid.plts[[i]] <- grid.plt
}

marrangeGrob(grobs = grid.plts, nrow = 2, ncol = 2)

summary(fit1, par=c("theta_fwi"), probs = c(0.05,0.5, 0.95))$summary

# data.scenario <- data.frame("x" = newx,
#                             "post.mean" = (alpha.samples[,1]),
#                             "post.median" = (alpha.samples[,5]),
#                             "q1" = (alpha.samples[,4]),
#                             "q3" = (alpha.samples[,6]))

# ggplot(data.scenario, aes(x=x)) + 
#   ylab(expression(alpha(c,...,c))) + xlab(expression(c)) + labs(col = "") +
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
#   geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values = c("steelblue")) + 
#   guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(size = 20))

# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_alpha.pdf"), width=10, height = 7.78)

# xi.scenario <- data.frame("x" = newx,
#                             "post.mean" = (1/alpha.samples[,1]),
#                             "post.median" = (1/alpha.samples[,5]),
#                             "q1" = (1/alpha.samples[,4]),
#                             "q3" = (1/alpha.samples[,6]),
#                             "evgam.1" = xi.pred.1,
#                             # "evgam.1.q1" = xi.low.1,
#                             # "evgam.1.q3" = xi.high.1,
#                             # "vgam.1" = vgam.xi.1,
#                             # "vgam.scale" = vgam.xi.scale,
#                             # "evgam.scale.q1" = xi.low.scale,
#                             # "evgam.scale.q3" = xi.high.scale,
#                             "evgam.scale" = xi.pred.scale)

# ggplot(xi.scenario, aes(x=x)) + 
#   ylab(expression(xi(c,...,c))) + xlab(expression(c)) + labs(col = "") +
#   geom_ribbon(aes(ymin = q1, ymax = q3, fill="Credible Band"), alpha = 0.2) +
#   geom_line(aes(y=post.median, col = "Posterior Median"), linewidth=1) +
#   # geom_ribbon(aes(ymin = evgam.scale.q1, ymax = evgam.scale.q3), fill= "orange", alpha = 0.2) +
#   geom_line(aes(y=evgam.scale), colour = "orange", linewidth=1, linetype=3) +
#   # geom_ribbon(aes(ymin = evgam.1.q1, ymax = evgam.1.q3), fill= "purple", alpha = 0.2) +
#   geom_line(aes(y=evgam.1), colour = "purple", linewidth=1, linetype=3) +
#   # geom_line(aes(y=vgam.scale), colour = "orange", linewidth=1, linetype=4) +
#   # geom_line(aes(y=vgam.1), colour = "purple", linewidth=1, linetype=4) +  
#   scale_fill_manual(values=c("steelblue"), name = "") +
#   scale_color_manual(values = c("steelblue")) + 
#   guides(color = guide_legend(order = 2), 
#           fill = guide_legend(order = 1)) +
#   theme_minimal(base_size = 30) +
#   theme(legend.position = "none",
#         strip.text = element_blank(),
#         axis.text = element_text(size = 20))

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
                        # evgam.traj = evgam.traj, 
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

# scaled_fwi_grid <- sweep(bs.linear[,-1]*X_sd, 2, fwi.min, "-")
# scaled_fwi_grid <- sweep(scaled_fwi_grid, 2, fwi.minmax, "/")

data.smooth <- data.frame("x" = newx,
                          "true" = as.vector(as.matrix(fwi.linear)),
                          "post.mean" = as.vector(g.smooth.mean),
                          "q1" = as.vector(g.smooth.q1),
                          "q2" = as.vector(g.smooth.q2),
                          "q3" = as.vector(g.smooth.q3))


# grid.plts <- list()
# for(i in 1:(p*4)){
#   grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
#                   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#                   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#                   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#                   geom_rug(data = subset(data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), true>0), aes(x=true, y=q2), sides = "b") +
#                   geom_rug(aes(x=true, y=q2), sides = "b") +
#                   ylab("") + xlab(colnames(bs.linear[,-1])[i]) +
#                   scale_fill_manual(values=c("steelblue"), name = "") + 
#                   scale_color_manual(values=c("steelblue")) +
#                   ylim(g.min.samples[i], g.max.samples[i]) +
#                   theme_minimal(base_size = 30) +
#                   theme(legend.position = "none",
#                           plot.margin = margin(0,0,0,-20),
#                           axis.text = element_text(size = 35),
#                           axis.title.x = element_text(size = 45))
#   grid.plts[[i]] <- grid.plt #+ annotate("point", x= fwi.linear[which.max(y),i], y=g.min.samples[i], color = "red", size = 7)
# }

# marrangeGrob(grobs = grid.plts, nrow = 1, ncol = 4)

# # ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_mcmc_DC.pdf"), grid.plts[[7]], width=10, height = 7.78)


# data.linear <- data.frame("x" = newx,
#                           "true" = as.vector(as.matrix(fwi.linear)),
#                           "post.mean" = as.vector(g.linear.mean),
#                           "q1" = as.vector(g.linear.q1),
#                           "q2" = as.vector(g.linear.q2),
#                           "q3" = as.vector(g.linear.q3),
#                           "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


# grid.plts <- list()
# for(i in 1:(p*4)){
#   grid.plt <- ggplot(data = data.frame(data.linear[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
#                   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#                   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#                   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#                   geom_rug(data = subset(data.frame(data.linear[((((i-1)*n)+1):(i*n)),]), true>0), aes(x=true, y=q2), sides = "b") +
#                   # geom_rug(aes(x=true, y=q2), sides = "b") +
#                   ylab("") + xlab(colnames(bs.linear[,-1])[i]) +
#                   scale_fill_manual(values=c("steelblue"), name = "") + 
#                   scale_color_manual(values=c("steelblue")) +
#                   ylim(g.min.samples[i], g.max.samples[i]) +
#                   theme_minimal(base_size = 30) +
#                   theme(legend.position = "none",
#                           plot.margin = margin(0,0,0,-20),
#                           axis.text = element_text(size = 35),
#                           axis.title.x = element_text(size = 45))
#   grid.plts[[i]] <- grid.plt #+ annotate("point", x= fwi.linear[which.max(y),i], y=g.min.samples[i], color = "red", size = 7)
# }

# # grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)
# marrangeGrob(grobs = grid.plts, nrow = 1, ncol = 4)


# data.nonlinear <- data.frame("x" = newx,
#                             "true" = as.vector(as.matrix(fwi.linear)),
#                             "post.mean" = as.vector(g.nonlinear.mean),
#                             "q1" = as.vector(g.nonlinear.q1),
#                             "q2" = as.vector(g.nonlinear.q2),
#                             "q3" = as.vector(g.nonlinear.q3),
#                             "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


# grid.plts <- list()
# for(i in 1:(p*4)){
#   grid.plt <- ggplot(data = data.frame(data.nonlinear[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
#                   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#                   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#                   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#                   geom_rug(data = subset(data.frame(data.nonlinear[((((i-1)*n)+1):(i*n)),]), true>0), aes(x=true, y=q2), sides = "b") +
#                   ylab("") + xlab(colnames(bs.linear[,-1])[i]) +
#                   scale_fill_manual(values=c("steelblue"), name = "") + 
#                   scale_color_manual(values=c("steelblue")) +
#                   ylim(g.min.samples[i], g.max.samples[i]) +
#                   theme_minimal(base_size = 30) +
#                   theme(legend.position = "none",
#                           plot.margin = margin(0,0,0,-20),
#                           axis.text = element_text(size = 35),
#                           axis.title.x = element_text(size = 45))
#   grid.plts[[i]] <- grid.plt #+ annotate("point", x= fwi.linear[which.max(y),i], y=g.min.samples[i], color = "red", size = 7)
# }
# marrangeGrob(grobs = grid.plts, nrow = 1, ncol = 4)

data.fwi <- data.frame("x" = as.vector(as.matrix(fwi.true)),
                        "true" = as.vector(as.matrix(true.linear)),
                        "post.mean" = as.vector(fwi.smooth.mean),
                        "q1" = as.vector(fwi.smooth.q1),
                        "q2" = as.vector(fwi.smooth.q2),
                        "q3" = as.vector(fwi.smooth.q3))


# grid.plts <- list()
# for(i in 1:(p*4)){
#   grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), aes(x=x)) + 
#                   geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
#                   geom_ribbon(aes(ymin = q1, ymax = q3, fill = "Credible Band"), alpha = 0.2) +
#                   geom_line(aes(y=q2, colour = "Posterior Median"), linewidth=1) + 
#                   geom_rug(data = subset(data.frame(data.smooth[((((i-1)*n)+1):(i*n)),]), true>0), aes(x=true, y=q2), sides = "b") +
#                   ylab("") + xlab(colnames(bs.linear[,-1])[i]) +
#                   scale_fill_manual(values=c("steelblue"), name = "") + 
#                   scale_color_manual(values=c("steelblue")) +
#                   ylim(fwi.min.samples[i], fwi.max.samples[i]) +
#                   theme_minimal(base_size = 30) +
#                   theme(legend.position = "none",
#                           plot.margin = margin(0,0,0,-20),
#                           axis.text = element_text(size = 35),
#                           axis.title.x = element_text(size = 45))
#   grid.plts[[i]] <- grid.plt #+ annotate("point", x= fwi.true[which.max(y),i], y=fwi.min.samples[i], color = "red", size = 7)
# }

# marrangeGrob(grobs = grid.plts, nrow = 1, ncol = 4)

# summary(fit1, par=c("theta_fwi"), probs = c(0.05,0.5, 0.95))$summary

# 1. Configuration
global_y_min <- min(data.fwi$q1, na.rm = TRUE)
global_y_max <- max(data.fwi$q3, na.rm = TRUE)
seasons <- c("Winter", "Spring", "Summer", "Autumn")
n_seasons <- 4
sorted_plots <- list()
plot_counter <- 1

for(s in 1:n_seasons) {
  for(j in 1:p) {
    stan_idx <- (j - 1) * n_seasons + s
    start_row <- (stan_idx - 1) * n + 1
    end_row   <- stan_idx * n
    df_slice <- data.smooth[start_row:end_row, ]
    current_var <- covariates[j]
    current_season <- seasons[s]
    x_min <- min(df_slice$x, na.rm = TRUE)
    x_max <- max(df_slice$x, na.rm = TRUE)    
    my_breaks <- pretty(c(x_min, x_max), n = 5)
    p_plot <- ggplot(data = df_slice, aes(x = x)) + 
      geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 1.2) + 
      geom_ribbon(aes(ymin = q1, ymax = q3), fill = "steelblue", alpha = 0.2) +
      geom_line(aes(y = q2), color="steelblue", linewidth = 1) + 
      geom_rug(data = subset(df_slice, true > 1e-6), 
               aes(x = true, y = q2), sides = "b", alpha = 0.6) +
      ylim(global_y_min, global_y_max) +
      ylab("") + xlab(paste(current_var, "-", current_season)) +
      # scale_fill_manual(values = c("steelblue"), name = "") + 
      # scale_color_manual(values = c("steelblue")) +
      scale_x_continuous(breaks = my_breaks, limits = range(my_breaks)) +
      # ylim(fwi.min.samples[stan_idx], fwi.max.samples[stan_idx]) +
      
      theme_minimal(base_size = 20) +
      theme(legend.position = "none",
            plot.margin = margin(5, 5, 5, 5),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text = element_text(size = 18),
            axis.title.x = element_text(size = 22))
    sorted_plots[[plot_counter]] <- p_plot
    plot_counter <- plot_counter + 1
  }
}
marrangeGrob(grobs = sorted_plots, nrow = 1, ncol = 5, top = NULL)


sorted_plots <- list()
plot_counter <- 1
for(s in 1:n_seasons) {
  for(j in 1:p) {
    stan_idx <- (j - 1) * n_seasons + s
    start_row <- (stan_idx - 1) * n + 1
    end_row   <- stan_idx * n
    df_slice <- data.fwi[start_row:end_row, ]
    current_var <- covariates[j]
    current_season <- seasons[s]
    x_min <- min(df_slice$x, na.rm = TRUE)
    x_max <- max(df_slice$x, na.rm = TRUE)
    my_breaks <- floor(seq(from = x_min, to = x_max, length.out = 5))
    p_plot <- ggplot(data = df_slice, aes(x = x)) + 
      geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 1.2) + 
      geom_ribbon(aes(ymin = q1, ymax = q3), fill = "steelblue", alpha = 0.2) +
      geom_line(aes(y = q2), color="steelblue", linewidth = 1) + 
      geom_rug(data = subset(df_slice, true > 0), 
               aes(x = true, y = q2), sides = "b", alpha = 0.6) +
      ylim(global_y_min, global_y_max) +
      ylab("") + xlab(paste(current_var, "-", current_season)) +
      scale_x_continuous(breaks = my_breaks, limits = range(my_breaks)) +
      # ylim(fwi.min.samples[stan_idx], fwi.max.samples[stan_idx]) +
      
      theme_minimal(base_size = 20) +
      theme(legend.position = "none",
            plot.margin = margin(5, 5, 5, 5),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text = element_text(size = 18),
            axis.title.x = element_text(size = 22))
    sorted_plots[[plot_counter]] <- p_plot
    plot_counter <- plot_counter + 1
  }
}

marrangeGrob(grobs = sorted_plots, nrow = 1, ncol = 5, top = NULL)

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
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_BLAST_predictive_distribution.pdf"), width=10, height = 7.78)


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
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_pareto_post_generative.pdf"), width = 10, height = 7.78)

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
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_BLAST_generative.pdf"), width = 10, height = 7.78)

# library(ismev)
# gpd.fit(y-u, u)

# fit.log.lik <- extract_log_lik(fit1)
# constraint.elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 2)
# save(constraint.elpd.loo, constraint.waic, file = (paste0("./BLAST/application/BLAST_constraint_",Sys.Date(),"_",psi,"_",floor(threshold*100),"quantile_IC.Rdata")))


# x_data <- sort(Y[Y > 0], decreasing = TRUE)
# n <- length(x_data)
# k_range_pick <- 15:500
# pick_res <- data.frame(k = k_range_pick, xi = NA, lower = NA, upper = NA)

# for(i in seq_along(k_range_pick)) {
#   k <- k_range_pick[i]
  
#   # Quantiles: X_{n-k+1}, X_{n-2k+1}, X_{n-4k+1}
#   q1 <- x_data[k]
#   q2 <- x_data[2*k]
#   q3 <- x_data[4*k]
  
#   # 1. Point Estimate
#   # Formula: (1/ln2) * ln( (q1-q2) / (q2-q3) )
#   xi_hat <- (1/log(2)) * log((q1 - q2) / (q2 - q3))
  
#   # 2. Asymptotic Variance Formula for Pickands
#   # Var = [ xi^2 * (2^(2xi+1) + 1) ] / [ (2 * (2^xi - 1) * ln2)^2 ]
#   # Handle the xi=0 case to prevent division by zero (limit is approx 3.24)
#   if(abs(xi_hat) < 1e-6) {
#     asy_var <- 3.24 # Approx limit
#   } else {
#     num <- (xi_hat^2) * ((2^(2*xi_hat + 1)) + 1)
#     den <- (2 * (2^xi_hat - 1) * log(2))^2
#     asy_var <- num / den
#   }
  
#   # Standard Error = sqrt(Var / k)
#   se <- sqrt(asy_var / k)
  
#   pick_res$xi[i] <- xi_hat
#   pick_res$lower[i] <- xi_hat - 1.96 * se
#   pick_res$upper[i] <- xi_hat + 1.96 * se
# }

# # --- 3. MOMENT (DedH) ESTIMATOR ---
# # Constraint: Uses all k up to n-1
# k_range_mom <- 15:500
# mom_res <- data.frame(k = k_range_mom, xi = NA, lower = NA, upper = NA)

# for(i in seq_along(k_range_mom)) {
#   k <- k_range_mom[i]
#   y_k <- x_data[1:k]
  
#   # 1. Point Estimate
#   log_excess <- log(y_k) - log(x_data[k+1])
#   M1 <- mean(log_excess)
#   M2 <- mean(log_excess^2)
  
#   xi_hat <- M1 + 1 - 0.5 * (1 - (M1^2 / M2))^(-1)
  
#   # 2. Asymptotic Variance (for xi > 0)
#   # Var = 1 + xi^2
#   asy_var <- 1 + xi_hat^2
#   se <- sqrt(asy_var / k)
  
#   mom_res$xi[i] <- xi_hat
#   mom_res$lower[i] <- xi_hat - 1.96 * se
#   mom_res$upper[i] <- xi_hat + 1.96 * se
# }


# p1 <- ggplot(pick_res, aes(x=k, y=xi)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), fill="darkgreen", alpha=0.2) +
#   geom_line(color="darkgreen", linewidth=1) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_cartesian(ylim=c(-0.15, 1.5)) + # Zoom to relevant range
#   labs(title="Pickands Estimator", 
#        y="Extreme Value Index") +
#   theme_minimal(base_size = 30) + ylim(-10,10) +
#   theme(legend.position = "none",
#           axis.text = element_text(size = 35),
#           axis.title.x = element_text(size = 45))

# p2 <- ggplot(mom_res, aes(x=k, y=xi)) +
#   geom_ribbon(aes(ymin=lower, ymax=upper), fill="purple", alpha=0.2) +
#   geom_line(color="purple", linewidth=1) +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_cartesian(ylim=c(-0.15, 1.5)) + 
#   labs(title="Moment Estimator", y="") +
#   theme_minimal(base_size = 30) + ylim(0, 2) +
#   theme(legend.position = "none",
#           axis.text = element_text(size = 35),
#           axis.text.y = element_blank(),
#           axis.title.x = element_text(size = 45))

# grid.plt <- grid.arrange(p1, p2, nrow=1)
# ggsave(paste0("./BLAST/application/figures/",Sys.Date(),"_heavytail.pdf"), grid.plt, width=22, height = 7.78)
fit.log.lik <- extract_log_lik(fit1)
elpd.loo <- loo(fit.log.lik, is_method = "sis", cores = 4)
elpd.loo
# save(elpd.loo, file = (paste0("./BLAST/application/BLAST_full_",Sys.Date(),"_",psi+2,"_",floor(threshold*100),"quantile_IC.Rdata")))
