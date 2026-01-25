setwd("../BLAST/simulation/results")
iter <- 2
n <- 15000
EVGAM <- FALSE
file_pattern <- paste0("2026-01-24_",iter,"_MC_scD_15000_.*.Rdata")
file_list <- list.files(pattern = file_pattern)

# Check if files were found
if (length(file_list) == 0) {
  stop("No files found matching the pattern.")
}

# 2. Initialize the 4 final containers
# We use lists initially because they are faster to append to than dataframes
alpha.container <- matrix(, nrow=n*0.05, ncol=0)
newgsmooth.container <- matrix(, nrow=n*0.05, ncol=0)
gridgsmooth.container <- matrix(, nrow=n*0.05, ncol=0)
gridgl.container <- matrix(, nrow=n*0.05, ncol=0)
gridgnl.container <- matrix(, nrow=n*0.05, ncol=0)
mise.container <- c()
qqplot.container <- matrix(, nrow=n*0.05, ncol=0)

# 3. Iterate through the files
for (f in file_list) {
  message(paste("Loading file:", f))
  
  temp_env <- new.env()
  
  if(EVGAM == TRUE){
    tryCatch({
      load(f, envir = temp_env)
      alpha.container <- cbind(alpha.container, temp_env$alpha.container[,c(1:iter)])
      gridgsmooth.container <- cbind(gridgsmooth.container, temp_env$gridgsmooth.container[,c(1:iter)])
      evgam.1.container <- cbind(gridgl.container, temp_env$gridgl.container[,c(1:iter)])
      evgam.scale.container <- cbind(gridgnl.container, temp_env$gridgnl.container[,c(1:iter)])
      mise.container <- c(mise.container, temp_env$mise.container)
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

if(EVGAM==FALSE){  
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
  # newgsmooth.container$x <- temp_env$newgsmooth.container$x
  # newgsmooth.container$true <- temp_env$newgsmooth.container$true
  # newgsmooth.container$covarite <- temp_env$newgsmooth.container$covariate
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

  
  # save(alpha.container, gridgnl.container, gridgl.container, gridgsmooth.container, mise.container, qqplot.container, file = paste0(Sys.Date(),"_",total.iter,"_MC_scD_",n,".Rdata"))
}