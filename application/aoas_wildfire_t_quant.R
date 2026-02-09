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
library(splines)
library(quantreg)
library(qgam)
library(mgcViz)
library(brms)
library(evgam)
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
# Y[Y==0] <- 1e-5
psi <- 30

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
    fwi.scaled[,i] <- fwi.index[,i] <- cov.long$measurement[missing.values]
}

fwi.index$date <- substr(cov.long$...1[missing.values],9,10)
fwi.index$month <- factor(format(as.Date(substr(cov.long$...1[missing.values],1,10), "%Y-%m-%d"),"%b"),
                            levels = c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
fwi.index$date <- as.numeric(fwi.index$date)
fwi.index$year <- substr(as.Date(cov.long$condition[missing.values], "%Y"),1,4)
fwi.origin <- fwi.scaled

fwi.origin <- data.frame(fwi.origin, time = c(1:length(Y)), BA=Y)


# fit.qr <- brm(bf(BA ~ s(time), quantile = 0.975),
#                 data = fwi.origin,
#                 cores = 3,
#                 chain = 3,
#                 family = asym_laplace())

# evgam.time <- BA ~ s(time, k = 30)
# ald.time.fit <- evgam(evgam.time, data = fwi.origin, family = "ald", ald.args=list(tau = 0.975))
# evgam.cov <- BA ~ s(BUI, k = 30) + s(ISI, k = 30) + s(FFMC, k = 30) + s(DMC, k = 30) + s(DC, k = 30)
# ald.cov.fit <- evgam(evgam.cov, data = fwi.origin, family = "ald", ald.args=list(tau = 0.975))
# save(ald.time.fit, ald.cov.fit, file="quant-evgam.Rdata")
# load("./BLAST/application/quant-evgam-scaled.Rdata")
# evgam.cov.pred <- predict(ald.cov.fit, type = "response")$location
# evgam.time.pred <- predict(ald.time.fit, type = "response")$location
# save(ald.time.fit, ald.cov.fit, evgam.cov.pred, evgam.time.pred, file="./BLAST/application/quant-evgam-scaled.Rdata")

# quant.fit <- qgam(BA ~ s(time,k=30), data = fwi.origin, qu = 0.975)
# quant.fit <- qgam(BA ~ s(BUI,bs="ts", k=10) +
#                       s(ISI,bs="ts", k=10) +
#                       s(FFMC,bs="ts", k=10) +
#                       s(DMC,bs="ts", k=10) +
#                       s(DC,bs="ts", k=10), data = fwi.origin, qu = 0.95)
# save(quant.fit, file = './BLAST/application/qgam_95_10_ts.Rdata')
# load("./BLAST/application/quant-time.Rdata")
# load("./BLAST/application/qgam_5_none.Rdata")
load("./BLAST/application/quant-t_30.Rdata")
# load("./BLAST/application/qgam_96_30_ts.Rdata")
# print(plot(quant.fit, allTerms = TRUE), pages = 1)
check(quant.fit)
quant.viz <- getViz(quant.fit, nsim = 20)
print(plot(quant.viz, allTerms = TRUE), pages = 1)
# check1D(quant.viz, fwi.origin[,4]) + l_gridQCheck1D(qu = 0.975)

preds <- predict(quant.fit)
# save(preds, quant.fit, file="./BLAST/application/quant-time.Rdata")
empirical_qu <- mean(fwi.origin$BA < preds)
print(paste("Empirical Quantile:", empirical_qu))

# fwi.origin <- data.frame(fwi.origin[which(Y>1),], BA=Y[Y>1])
# BA.shifted <- ifelse(fwi.origin$BA == 0, 1e-5, fwi.origin$BA)
# fwi.origin$log.BA <- log(fwi.origin$BA+1)
vars <- colnames(fwi.origin)[1:7]
seq_list <- lapply(vars, function(var) seq(min(fwi.origin[[var]]), max(fwi.origin[[var]]), length.out = length(Y)))
grid.df <- as.data.frame(setNames(seq_list, vars))

library(tidyverse)
fwi.index$BA <- fwi.origin$BA
df_seasonal <- fwi.index %>%
  mutate(
    year = as.numeric(as.character(year)), 
    Month_Num = match(month, month.abb), 
    
    Season = case_when(
      Month_Num %in% c(12, 1, 2) ~ "Winter",
      Month_Num %in% c(3, 4, 5) ~ "Spring",
      Month_Num %in% c(6, 7, 8) ~ "Summer",
      Month_Num %in% c(9, 10, 11) ~ "Autumn"
    ),

    SeasonYear = ifelse(Month_Num == 12, year + 1, year)
  )
  
plot_data <- df_seasonal %>%
  group_by(SeasonYear, Season) %>%
  summarise(
    Q975 = quantile(BA, 0.975, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
  arrange(SeasonYear, Season) %>%
  mutate(Label = paste(SeasonYear, Season)) %>%
  mutate(Label = factor(Label, levels = unique(Label)))

df_seasonal <- df_seasonal %>% 
  mutate(Label = paste(SeasonYear, Season)) %>%
  mutate(Label = factor(Label, levels = unique(Label)))
# 3. Plot
ggplot(plot_data, aes(x = Label, y = Q975, group = 1)) +
# ggplot(df_seasonal, aes(x = Label, y = BA, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_hline(aes(yintercept=quantile(Y,0.975)), linetype="dashed") +
  geom_point(aes(color = Season), size = 3) +
  scale_color_manual(values = c(
    "Summer" = "orange", 
    "Winter" = "red", 
    "Spring" = "red", 
    "Autumn" = "red"
  )) +
  labs(y = "BA", x = "Season", color = "Season") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


df_seasonal$fitted_975 <- predict(quant.fit)
load("./BLAST/application/qgam_96_30_ts.Rdata")
df_seasonal$fitted_cov <- predict(quant.fit)

block_summary <- df_seasonal %>%
  group_by(Label, Season) %>%
  summarise(
    Empirical = (quantile(BA, probs = 0.975, na.rm = TRUE)),
    Model_Smooth = median(fitted_975),#(quantile(fitted_975, probs = 0.975, na.rm = TRUE)), # Average smooth value for that block
    Model_cov = median(fitted_cov), #(quantile(fitted_cov, probs = 0.95, na.rm = TRUE)),
    Exceedance_Count = sum(BA > Model_cov, na.rm = TRUE),
    .groups = 'drop'
)%>%
mutate(across(c(Empirical, Model_Smooth, Model_cov), ~ifelse(. <= 0, 1, .)))

max_log_val <- log10(max(block_summary$Empirical, na.rm = TRUE))
max_count <- max(block_summary$Exceedance_Count, na.rm = TRUE)
count_scale <- ifelse(max_count == 0, 1, max_log_val / max_count)

ggplot(block_summary, aes(x = Label, group = 1)) +
  geom_col(aes(y = 10^(Exceedance_Count * count_scale), fill = Season), 
           alpha = 0.3, width = 0.7) +
  geom_hline(aes(yintercept=(quantile(Y, 0.975))), linetype="dashed") +
  geom_line(aes(y = Model_Smooth), color = "steelblue", linewidth = 1.2) +
  geom_line(aes(y = Model_cov), color = "darkgreen", linewidth = 1, linetype = 2) +
  geom_point(aes(y = Empirical, color = Season), size = 3) +
  scale_y_log10(
    name = "Burnt Area (Log10)",
    sec.axis = sec_axis(~ log10(.) / count_scale, 
                        name = "Exceedance Count",
                        breaks = seq(0, max_count, by = 1))
  ) +
  scale_color_manual(values = c("Summer" = "orange", "Winter" = "red", "Spring" = "red", "Autumn" = "red")) +
  scale_fill_manual(values = c("Summer" = "orange", "Winter" = "red", "Spring" = "red", "Autumn" = "red")) +
  labs(x = "Season Year") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


ggplot(block_summary, aes(x = Label, group = 1)) +
  geom_point(aes(y = Empirical, color = Season), size = 3) +
  # geom_line(aes(y = Empirical), color = "grey80", linetype = "dashed") +
  geom_hline(aes(yintercept=(quantile(Y, 0.975))), linetype="dashed") +
  geom_line(aes(y = Model_Smooth), color = "steelblue", linewidth = 1.2) +
  geom_line(aes(y = Model_cov), color = "darkgreen", linewidth = 0.7, linetype=2) + scale_y_log10() +
  scale_color_manual(values = c(
    "Summer" = "orange", 
    "Winter" = "red", 
    "Spring" = "red", 
    "Autumn" = "red"
  )) +
  labs(
    y = "Burnt Area", 
    x = "Time (blockwise)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )



coverage_summary <- df_seasonal %>%
  left_join(block_summary %>% select(Label, Model_Smooth), by = "Label") %>%
  group_by(Label, Season) %>%
  summarise(
    # Using the block's 97.5th percentile prediction as the boundary
    n_obs = n(),
    n_below = sum(BA <= Model_Smooth),
    empirical_coverage = n_below / n_obs,
    # empirical_coverage = mean(BA <= Model_Predicted),
    .groups = 'drop'
  )

# 3. View the summary
print(coverage_summary)

tail(coverage_summary)
block_summary$origin_Empirical <- block_summary$Empirical
block_summary$origin_Model_Smooth <- block_summary$Model_Smooth
block_summary$origin_Model_cov <- block_summary$Model_cov

df_final <- df_seasonal %>%
  left_join(block_summary %>% 
              select(Label, origin_Empirical, origin_Model_Smooth, origin_Model_cov), 
            by = "Label")

fwi.dd <- df_final %>% mutate(excess = BA > origin_Model_Smooth)
fwi.dd <- df_final %>% mutate(excess = BA > origin_Model_cov)
tail(fwi.dd[which(fwi.dd$excess==TRUE),])
# save(preds, quant.fit, fwi.dd, file="./BLAST/application/quant-975_30.Rdata")

length(which(fwi.dd$excess==TRUE))
