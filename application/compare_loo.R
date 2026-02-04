library(loo)
setwd("C:/Users/Johnny Lee/Documents/GitHub")
p <- 7
threshold <- 0.97
date <- "2026-02-04"
method <- "time"

load(paste0("./BLAST/application/BLAST_full_",date,"_",30,"_",floor(threshold*100),"quantile_IC_",p,"_", method ,".Rdata"))
elpd.loo.full.30 <- elpd.loo
load(paste0("./BLAST/application/BLAST_linear_",date,"_",30,"_",floor(threshold*100),"quantile_IC_",p,"_", method ,".Rdata"))
elpd.loo.linear <- elpd.loo
load(paste0("./BLAST/application/BHST_full_",date,"_",30,"_",floor(threshold*100),"quantile_IC_",p,"_", method ,".Rdata"))
elpd.loo.hs <- elpd.loo
load(paste0("./BLAST/application/BRIT_full_",date,"_",30,"_",floor(threshold*100),"quantile_IC_",p,"_", method ,".Rdata"))
elpd.loo.ridge.full <- elpd.loo
load(paste0("./BLAST/application/BRIT_linear_",date,"_",30,"_",floor(threshold*100),"quantile_IC_",p,"_", method ,".Rdata"))
elpd.loo.ridge.linear <- elpd.loo


compare <- loo_compare(list(BLAST_full = elpd.loo.full.30, 
                            BLAST_linear=elpd.loo.linear, 
                            BHST_full = elpd.loo.hs,
                            BRIT_full = elpd.loo.ridge.full,
                            BRIT_linear = elpd.loo.ridge.linear))
print(compare, simplify = FALSE)

