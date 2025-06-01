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
library(patchwork)
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
u <- quantile(Y, threshold)
y <- Y[Y>u]

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
fwi.origin <- fwi.scaled <-fwi.scaled[which(Y>u),]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
fwi.scaled <- as.data.frame(sapply(fwi.origin, FUN = range01))

n <- dim(fwi.scaled)[[1]]
p <- dim(fwi.scaled)[[2]]

fwi.origin <- data.frame(fwi.index[which(Y>u),], BA=y)
max.fwi <- fwi.origin[which.max(y),]
xholder <- matrix(nrow=n, ncol=p)
for(i in 1:p){
  xholder[,i] <- seq(min(fwi.scaled[,i]), max(fwi.scaled[,i]), length.out = n)
}

data.full <- readRDS("./BLAST/application/figures/comparison/full_stanfit.rds")
data.hs <- readRDS("./BLAST/application/figures/comparison/hs_stanfit.rds")
data.linear <- readRDS("./BLAST/application/figures/comparison/linear_stanfit.rds")

data.smooth <- data.frame("x" = as.vector(xholder),
                          "true" = data.full$true,
                          "full.mean" = data.full$post.mean,
                          "full.q1" = data.full$q1,
                          "full.q2" = data.full$q2,
                          "full.q3" = data.full$q3,
                          "hs.mean" = data.hs$post.mean,
                          "hs.q1" = data.hs$q1,
                          "hs.q2" = data.hs$q2,
                          "hs.q3" = data.hs$q3,
                          "linear.mean" = data.linear$post.mean,
                          "linear.q1" = data.linear$q1,
                          "linear.q2" = data.linear$q2,
                          "linear.q3" = data.linear$q3,
                          "covariates" = gl(p, n, (p*n), labels = names(fwi.scaled)))


grid.plts <- list()
for(i in 1:p){
  grid.plt <- ggplot(data = data.frame(data.smooth[((((i-1)*n)+1):(i*n)),])
                  , aes(x=x)) + 
                  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", linewidth = 2) + 
                  geom_ribbon(aes(ymin = linear.q1, ymax = linear.q3, fill = "Linear"), alpha = 0.2) +
                  geom_line(aes(y=linear.q2, colour = "Linear"), linewidth=1, linetype = 3) + 
                  geom_rug(aes(x=true, y=full.q2), sides = "b") +
                  geom_ribbon(aes(ymin = hs.q1, ymax = hs.q3, fill = "Horseshoe"), alpha = 0.2) +
                  geom_line(aes(y=hs.q2, colour = "Horseshoe"), linewidth=1, linetype = 2) + 
                  geom_ribbon(aes(ymin = full.q1, ymax = full.q3, fill = "Full"), alpha = 0.2) +
                  geom_line(aes(y=full.mean, colour = "Full"), linewidth=1.4) + 
                  ylab("") + xlab(names(fwi.scaled)[i]) +
                  scale_fill_manual(values=c("steelblue", "red", "darkgreen"), name = "") + 
                  scale_color_manual(values=c("steelblue", "red", "darkgreen")) +
                  ylim(-5.05, 2.6) +
                  guides(color = guide_legend("Versions"),
                            fill = guide_legend("Versions")) +
                  theme_minimal(base_size = 30)
                  if (i == 1) {
                    grid.plt <- grid.plt + theme(plot.margin = margin(0,0,0,-20),
                                  legend.position = "bottom",
                                  legend.title = element_blank(),
                                  axis.text = element_text(size = 22),
                                  axis.title.x = element_text(size = 35))
                  }
                  else if(i == 4) {
                    grid.plt <- grid.plt + theme(plot.margin = margin(0,25,0,-20),
                                  legend.position = "bottom",
                                  legend.title = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.text.x = element_text(size = 22),
                                  axis.title.x = element_text(size = 35))
                  }
                  else if(i == 5) {
                    grid.plt <- grid.plt + theme(plot.margin = margin(0,0,0,-20),
                                  legend.position = "none",
                                  legend.title = element_blank(),
                                  axis.text = element_text(size = 22),
                                  axis.title.x = element_text(size = 35))
                  }
                  else if(i == 6 || i == 7){
                     grid.plt <- grid.plt + theme(plot.margin = margin(0,0,0,-20),
                                  legend.position = "none",
                                  legend.title = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.text = element_text(size = 22),
                                  axis.title.x = element_text(size = 35))                   
                  }                  
                  else {
                    grid.plt <- grid.plt + theme(plot.margin = margin(0,0,0,-20),
                                  legend.position = "bottom",
                                  legend.title = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.text.x = element_text(size = 22),
                                  axis.title.x = element_text(size = 35))
                  }
  grid.plts[[i]] <- grid.plt + annotate("point", x= fwi.scaled[which.max(y),i], y=-5.05, color = "red", size = 4)
}
# grid.arrange(grobs = grid.plts, ncol = 4, nrow = 2)


pw <- grid.plts[[1]] + grid.plts[[2]] + grid.plts[[3]] + grid.plts[[4]] + grid.plts[[5]] + grid.plts[[6]] + grid.plts[[7]] + guide_area()
pw + plot_layout(ncol=4, guides = "collect") & theme(legend.position = "right", legend.title = element_blank())

pw <- (grid.plts[[1]] | grid.plts[[2]] | grid.plts[[3]] | grid.plts[[4]]) / (grid.plts[[5]] | grid.plts[[6]] | grid.plts[[7]])
pw <- (grid.plts[[1]] | grid.plts[[2]] | grid.plts[[3]] | grid.plts[[4]]) /
      (wrap_elements(grid.plts[[5]]) | wrap_elements(grid.plts[[6]]) | wrap_elements(grid.plts[[7]]))
pw + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.title = element_blank())
# ggsave(paste0("./BLAST/application/figures/comparison/pareto_mcmc_comparisons.pdf"), width=18, height = 10)



