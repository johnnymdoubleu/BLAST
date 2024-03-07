library(loo)

threshold = 0.95
load(paste0("./BRSTIR/application/BRTIR_",Sys.Date(),"_",floor(threshold*100),"quantile_IC.Rdata"))
load(paste0("./BRSTIR/application/BRSTIR_",Sys.Date(),"_",floor(threshold*100),"quantile_IC.Rdata"))

loo_compare(brtir.elpd.loo, brstir.elpd.loo)
loo_compare(brtir.waic, brstir.waic)
