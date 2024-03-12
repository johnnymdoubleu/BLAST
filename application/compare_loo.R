library(loo)
setwd("C:/Users/Johnny Lee/Documents/GitHub")
threshold = 0.95
load(paste0("./BRSTIR/application/BRTIR_",Sys.Date(),"_",floor(threshold*100),"quantile_IC.Rdata"))
load(paste0("./BRSTIR/application/BRSTIR_",Sys.Date(),"_",floor(threshold*100),"quantile_IC.Rdata"))
load(paste0("./BRSTIR/application/BRSTIR_constraint_",Sys.Date(),"_",floor(threshold*100),"quantile_IC.Rdata"))

loo_compare(brtir.elpd.loo, brstir.elpd.loo)
loo_compare(constraint.elpd.loo, brtir.elpd.loo)
loo_compare(brtir.waic, brstir.waic)
loo_compare(constraint.waic, brstir.waic)
