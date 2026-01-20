# 1. Get list of files matching the pattern
# This regex looks for files starting with the date, having a number, and ending with your specific suffix
setwd("../BLAST/simulation")
file_pattern <- "2026-01-19_2_MC_scA-2_15000_.*_.Rdata" 
file_list <- list.files(pattern = file_pattern)

# Check if files were found
if (length(file_list) == 0) {
  stop("No files found matching the pattern.")
}

# 2. Initialize the 4 final containers
# We use lists initially because they are faster to append to than dataframes
alpha.container <- matrix(, nrow=750, ncol=0)
newgsmooth.container <- matrix(, nrow=750, ncol=0)
mise.container <- c()
qqplot.container <- matrix(, nrow=750, ncol=0)

# 3. Iterate through the files
for (f in file_list) {
  message(paste("Loading file:", f))
  
  # Create a temporary environment to load the data
  # This prevents the loaded objects from overwriting variables in your global environment
  temp_env <- new.env()
  
  tryCatch({
    load(f, envir = temp_env)
    
    # 4. Append contents to the final containers
    # We wrap the loaded object in list() to maintain structure before combining
    alpha.container <- cbind(alpha.container, temp_env$alpha.container[,c(1,2)])
    newgsmooth.container <- cbind(newgsmooth.container, temp_env$newgsmooth.container[,c(1,2)])
    mise.container <- c(mise.container, temp_env$mise.container)
    qqplot.container <- cbind(qqplot.container, temp_env$qqplot.container[,c(1,2)])
    
  }, error = function(e) {
    warning(paste("Error loading file", f, ":", e$message))
  })
}

colnames(alpha.container) <- paste0("V", 1:50)
colnames(newgsmooth.container) <- paste0("V", 1:50)
colnames(qqplot.container) <- paste0("V", 1:50)
alpha.container$x <- temp_env$alpha.container$x
alpha.container$true <- temp_env$alpha.container$true
newgsmooth.container$x <- temp_env$newgsmooth.container$x
newgsmooth.container$true <- temp_env$newgsmooth.container$true
newgsmooth.container$mise.container <- temp_env$newgsmooth.container$covariate
qqplot.container$grid <- temp_env$qqplot.container$grid
