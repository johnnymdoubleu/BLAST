setwd("../BLAST/simulation/mcmc")
iter <- 2
file_pattern <- paste0("2026-01-23_",iter,"_MC_scO_15000_.*.Rdata")
file_list <- list.files(pattern = file_pattern)

# Check if files were found
if (length(file_list) == 0) {
  stop("No files found matching the pattern.")
}

# 2. Initialize the 4 final containers
# We use lists initially because they are faster to append to than dataframes
alpha.container <- matrix(, nrow=750, ncol=0)
newgsmooth.container <- matrix(, nrow=750, ncol=0)
gridgsmooth.container <- matrix(, nrow=750, ncol=0)
gridgl.container <- matrix(, nrow=750, ncol=0)
gridgnl.container <- matrix(, nrow=750, ncol=0)
mise.container <- c()
qqplot.container <- matrix(, nrow=750, ncol=0)

# 3. Iterate through the files
for (f in file_list) {
  message(paste("Loading file:", f))
  
  temp_env <- new.env()
  
  tryCatch({
    load(f, envir = temp_env)
    
    # 4. Append contents to the final containers
    # We wrap the loaded object in list() to maintain structure before combining
    alpha.container <- cbind(alpha.container, temp_env$alpha.container[,c(1:iter)])
    # newgsmooth.container <- cbind(newgsmooth.container, temp_env$newgsmooth.container[,c(1:5)])
    gridgsmooth.container <- cbind(gridgsmooth.container, temp_env$gridgsmooth.container[,c(1:iter)])
    gridgl.container <- cbind(gridgl.container, temp_env$gridgl.container[,c(1:iter)])
    gridgnl.container <- cbind(gridgnl.container, temp_env$gridgnl.container[,c(1:iter)])
    mise.container <- c(mise.container, temp_env$mise.container)
    qqplot.container <- cbind(qqplot.container, temp_env$qqplot.container[,c(1:iter)])
    
  }, error = function(e) {
    warning(paste("Error loading file", f, ":", e$message))
  })
}

colnames(alpha.container) <- paste0("V", 1:(length(file_list)*iter))
colnames(gridgsmooth.container) <- paste0("V", 1:(length(file_list)*iter))
colnames(gridgl.container) <- paste0("V", 1:(length(file_list)*iter))
colnames(gridgnl.container) <- paste0("V", 1:(length(file_list)*iter))
colnames(qqplot.container) <- paste0("V", 1:(length(file_list)*iter))
alpha.container$x <- temp_env$alpha.container$x
alpha.container$true <- temp_env$alpha.container$true
# newgsmooth.container$x <- temp_env$newgsmooth.container$x
# newgsmooth.container$true <- temp_env$newgsmooth.container$true
# newgsmooth.container$covarite <- temp_env$newgsmooth.container$covariate
gridgsmooth.container$x <- temp_env$gridgsmooth.container$x
gridgsmooth.container$true <- temp_env$gridgsmooth.container$true
gridgsmooth.container$covarite <- temp_env$gridgsmooth.container$covariate
gridgl.container$x <- temp_env$gridgl.container$x
gridgl.container$true <- temp_env$gridgl.container$true
gridgl.container$covarite <- temp_env$gridgl.container$covariate
gridgnl.container$x <- temp_env$gridgnl.container$x
gridgnl.container$true <- temp_env$gridgnl.container$true
gridgnl.container$covarite <- temp_env$gridgnl.container$covariate
qqplot.container$grid <- temp_env$qqplot.container$grid

save(alpha.container, gridgnl.container, gridgl.container, gridgsmooth.container, mise.container, qqplot.container, file = (paste0(Sys.Date(),"_",length(file_list)*iter,"_MC_sc0_",nrow(alpha.container),".Rdata")))
