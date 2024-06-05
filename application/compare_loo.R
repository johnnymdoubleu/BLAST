library(loo)
setwd("C:/Users/Johnny Lee/Documents/GitHub")
threshold <- 0.97
date <- "2024-06-04"
load(paste0("./BRSTIR/application/BRTIR_",date,"_",floor(threshold*100),"quantile_IC.Rdata"))
# load(paste0("./BRSTIR/application/BRSTIR_",Sys.Date(),"_",floor(threshold*100),"quantile_IC.Rdata"))
# load(paste0("./BRSTIR/application/BRSTIR_constraint_",date,"_",5,"_",floor(threshold*100),"quantile_IC.Rdata"))
# constraint.elpd.loo.5 <- constraint.elpd.loo
load(paste0("./BRSTIR/application/BRSTIR_constraint_",date,"_",10,"_",floor(threshold*100),"quantile_IC.Rdata"))
constraint.elpd.loo.10 <- constraint.elpd.loo
load(paste0("./BRSTIR/application/BRSTIR_constraint_",date,"_",20,"_",floor(threshold*100),"quantile_IC.Rdata"))
constraint.elpd.loo.20 <- constraint.elpd.loo
load(paste0("./BRSTIR/application/BRSTIR_constraint_",date,"_",30,"_",floor(threshold*100),"quantile_IC.Rdata"))
constraint.elpd.loo.30 <- constraint.elpd.loo
# load(paste0("./BRSTIR/application/BRSTIR_constraint_",date,"_",40,"_",floor(threshold*100),"quantile_IC.Rdata"))
# constraint.elpd.loo.40 <- constraint.elpd.loo


compare <- loo_compare(constraint.elpd.loo.10, constraint.elpd.loo.20,constraint.elpd.loo.30, brtir.elpd.loo)
print(compare, simplify = FALSE)

# loo_compare(brtir.elpd.loo, brstir.elpd.loo)
# loo_compare(brtir.waic, brstir.waic)
# loo_compare(constraint.waic, brstir.waic)
# compare <- loo_compare(constraint.elpd.loo.10, constraint.elpd.loo.20)
# print(compare, simplify = FALSE)
