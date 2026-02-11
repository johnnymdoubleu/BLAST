setwd("../BLAST/simulation/results")
iter <- 2
n <- 15000
EV <- 1
threshold <- 0.95
if(EV==TRUE){
  file_pattern <- paste0("evgam_mc_scC2_",n,"_.*.Rdata")
  # file_pattern <- paste0("evgam_mc_scA_.*.Rdata")
}else if(EV==FALSE){
  file_pattern <- paste0("2026-02-09_",iter,"_MC_scC_",n,"_.*.Rdata")
}

file_list <- list.files(pattern = file_pattern)

# Check if files were found
if (length(file_list) == 0) {
  stop("No files found matching the pattern.")
}

# 2. Initialize the 4 final containers
# We use lists initially because they are faster to append to than dataframes
alpha.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)
newgsmooth.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)
gridgsmooth.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)
gridgl.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)
gridgnl.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)
evgam.1.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)
evgam.scale.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)
vgam.1.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)
vgam.scale.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)
mise.container <- mise.evgam.1.container <- mise.evgam.scale.container <- mise.vgam.1.container <- mise.vgam.scale.container <- c()
qqplot.container <- matrix(, nrow=ceiling(n*(1-threshold)), ncol=0)

# 3. Iterate through the files
for (f in file_list) {
  message(paste("Loading file:", f))
  
  temp_env <- new.env()
  
  if(EV == TRUE){
    tryCatch({
      load(f, envir = temp_env)
      alpha.container <- cbind(alpha.container, temp_env$alpha.container[,c(1:iter)])
      gridgsmooth.container <- cbind(gridgsmooth.container, temp_env$newgsmooth.container[,c(1:iter)])
      evgam.1.container <- cbind(evgam.1.container, temp_env$evgam.1.container[,c(1:iter)])
      evgam.scale.container <- cbind(evgam.scale.container, temp_env$evgam.scale.container[,c(1:iter)])
      vgam.1.container <- cbind(vgam.1.container, temp_env$vgam.1.container[,c(1:iter)])
      vgam.scale.container <- cbind(vgam.scale.container, temp_env$vgam.scale.container[,c(1:iter)])      
      mise.container <- c(mise.container, temp_env$mise.container)
      mise.vgam.1.container <- c(mise.vgam.1.container, temp_env$mise.vgam.1.container)
      mise.vgam.scale.container <- c(mise.vgam.scale.container, temp_env$mise.vgam.scale.container)      
      mise.evgam.1.container <- c(mise.evgam.1.container, temp_env$mise.evgam.1.container)
      mise.evgam.scale.container <- c(mise.evgam.scale.container, temp_env$mise.evgam.scale.container)
    }, error = function(e) {
      warning(paste("Error loading file", f, ":", e$message))
    })
  }
  else {
    tryCatch({
      load(f, envir = temp_env)
      alpha.container <- cbind(alpha.container, temp_env$alpha.container[,c(1:iter)])
      gridgsmooth.container <- cbind(gridgsmooth.container, temp_env$gridgsmooth.container[,c(1:iter)])
      gridgl.container <- cbind(gridgl.container, temp_env$gridgl.container[,c(1:iter)])
      gridgnl.container <- cbind(gridgnl.container, temp_env$gridgnl.container[,c(1:iter)])
      mise.container <- c(mise.container, temp_env$mise.container)
      qqplot.container <- cbind(qqplot.container, temp_env$qqplot.container[,c(1:iter)])
    }, error = function(e) {
      warning(paste("Error loading file", f, ":", e$message))
    })    
  }
}

if(EV==FALSE){  
  total.iter <- length(file_list) * iter
  colnames(alpha.container) <- paste0("V", 1:total.iter)
  colnames(gridgsmooth.container) <- paste0("V", 1:total.iter)
  colnames(gridgl.container) <- paste0("V", 1:total.iter)
  colnames(gridgnl.container) <- paste0("V", 1:total.iter)
  colnames(qqplot.container) <- paste0("V", 1:total.iter)
  alpha.container$x <- temp_env$alpha.container$x
  alpha.container$true <- temp_env$alpha.container$true
  alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
  alpha.container <- as.data.frame(alpha.container)
  gridgsmooth.container$x <- temp_env$gridgsmooth.container$x
  gridgsmooth.container$true <- temp_env$gridgsmooth.container$true
  gridgsmooth.container$mean <- rowMeans(gridgsmooth.container[,1:total.iter])
  gridgsmooth.container$covariate <- temp_env$gridgsmooth.container$covariate
  gridgsmooth.container <- as.data.frame(gridgsmooth.container)
  gridgl.container$x <- temp_env$gridgl.container$x
  gridgl.container$true <- temp_env$gridgl.container$true
  gridgl.container$mean <- rowMeans(gridgl.container[,1:total.iter])
  gridgl.container$covariate <- temp_env$gridgl.container$covariate
  gridgl.container <- as.data.frame(gridgl.container)
  gridgnl.container$x <- temp_env$gridgnl.container$x
  gridgnl.container$true <- temp_env$gridgnl.container$true
  gridgnl.container$mean <- rowMeans(gridgnl.container[,1:total.iter])
  gridgnl.container$covariate <- temp_env$gridgnl.container$covariate
  gridgnl.container <- as.data.frame(gridgnl.container)
  qqplot.container$grid <- temp_env$qqplot.container$grid
  qqplot.container$mean <- rowMeans(qqplot.container[,1:total.iter])

  
  # save(alpha.container, gridgnl.container, gridgl.container, gridgsmooth.container, mise.container, qqplot.container, file = paste0(Sys.Date(),"_",total.iter,"_MC_scC_",n,".Rdata"))
} else {
  total.iter <- length(file_list) * iter
  colnames(alpha.container) <- paste0("V", 1:total.iter)
  colnames(evgam.1.container) <- paste0("V", 1:total.iter)
  colnames(evgam.scale.container) <- paste0("V", 1:total.iter)
  colnames(vgam.1.container) <- paste0("V", 1:total.iter)
  colnames(vgam.scale.container) <- paste0("V", 1:total.iter)  
  colnames(gridgsmooth.container) <- paste0("V", 1:total.iter)
  alpha.container$x <- temp_env$alpha.container$x
  alpha.container$true <- temp_env$alpha.container$true
  alpha.container$mean <- rowMeans(alpha.container[,1:total.iter])
  alpha.container$evgam.1 <- rowMeans(evgam.1.container[,1:total.iter])
  alpha.container$evgam.scale <- rowMeans(evgam.scale.container[,1:total.iter])  
  alpha.container$vgam.1 <- rowMeans(vgam.1.container[,1:total.iter])
  alpha.container$vgam.scale <- rowMeans(vgam.scale.container[,1:total.iter])   
  alpha.container <- as.data.frame(alpha.container)
  gridgsmooth.container$x <- temp_env$gridgsmooth.container$x
  gridgsmooth.container$true <- temp_env$gridgsmooth.container$true
  gridgsmooth.container$mean <- rowMeans(gridgsmooth.container[,1:total.iter])
  gridgsmooth.container$covariate <- temp_env$gridgsmooth.container$covariate
  gridgsmooth.container <- as.data.frame(gridgsmooth.container)
  
  cat("BLAST:   ", mean(mise.container, na.rm=TRUE), "±", sd(mise.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.container))), "\n", 
    "EVGAM:  ", mean(mise.evgam.1.container, na.rm=TRUE), "±", sd(mise.evgam.1.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.evgam.1.container))), "\n",
    "EVGAM-σ:", mean(mise.evgam.scale.container, na.rm=TRUE), "±", sd(mise.evgam.scale.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.evgam.scale.container))), "\n",
    "VGAM:   ", mean(mise.vgam.1.container, na.rm=TRUE), "±", sd(mise.vgam.1.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.vgam.1.container))), "\n",
    "VGAM-σ: ", mean(mise.vgam.scale.container, na.rm=TRUE), "±", sd(mise.vgam.scale.container, na.rm=TRUE)/sqrt(sum(!is.na(mise.vgam.scale.container))), "\n")

  # save(alpha.container, gridgsmooth.container, vgam.1.container, vgam.scale.container, evgam.1.container, evgam.scale.container, mise.container, mise.evgam.1.container, mise.evgam.scale.container, mise.vgam.1.container, mise.vgam.scale.container, file = paste0(Sys.Date(),"_evgam_mc_scC2_",(n*0.05),".Rdata"))  
}
