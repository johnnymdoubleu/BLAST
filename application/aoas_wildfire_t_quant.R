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

quant.fit <- qgam(BA ~ s(time), data = fwi.origin, qu = 0.975)
# load("./BLAST/application/quant-time.Rdata")
print(plot(quant.fit, allTerms = TRUE), pages = 1)
quant.viz <- getViz(quant.fit, nsim = 20)
print(plot(quant.viz))
check1D(quant.viz, fwi.origin[,8]) + l_gridQCheck1D(qu = 0.95)

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

# 2. Compute the Blockwise (Seasonal) Empirical Quantile
# This calculates the actual 97.5th percentile for each specific season/year
block_summary <- df_seasonal %>%
  group_by(Label, Season) %>%
  summarise(
    Empirical_975 = log(quantile(BA, probs = 0.975, na.rm = TRUE)),
    Model_Smooth_975 = log(quantile(fitted_975, probs = 0.975, na.rm = TRUE)), # Average smooth value for that block
    .groups = 'drop'
)

ggplot(block_summary, aes(x = Label, group = 1)) +
  geom_point(aes(y = Empirical_975, color = Season), size = 3) +
  geom_line(aes(y = Empirical_975), color = "grey80", linetype = "dashed") +
  geom_hline(aes(yintercept=log(quantile(Y,0.975))), linetype="dashed") +
  geom_line(aes(y = Model_Smooth_975), color = "steelblue", linewidth = 1.2) +
  
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
  left_join(block_summary %>% select(Label, Model_Smooth_975), by = "Label") %>%
  group_by(Label, Season) %>%
  summarise(
    # Using the block's 97.5th percentile prediction as the boundary
    n_obs = n(),
    n_below = sum(log(BA) <= Model_Smooth_975),
    empirical_coverage = n_below / n_obs,
    # empirical_coverage = mean(BA <= Model_Predicted_975),
    .groups = 'drop'
  )

# 3. View the summary
print(coverage_summary)

tail(coverage_summary)
block_summary$origin_Empirical_975 <- exp(block_summary$Empirical_975)
block_summary$origin_Model_Smooth_975 <- exp(block_summary$Model_Smooth_975)

df_final <- df_seasonal %>%
  left_join(block_summary %>% 
              select(Label, origin_Empirical_975, origin_Model_Smooth_975), 
            by = "Label")

fwi.dd <- df_final %>% mutate(excess = BA > origin_Model_Smooth_975)
tail(fwi.dd[which(fwi.dd$excess==TRUE),])
save(preds, quant.fit,fwi.dd, file="./BLAST/application/quant-time.Rdata")

length(which(fwi.dd$excess==TRUE))
