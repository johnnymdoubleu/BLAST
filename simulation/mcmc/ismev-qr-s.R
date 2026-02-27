library(ismev)
library(evgam)
library(ggplot2)

data(wooster)
negated.wooster <- -1 * wooster

time.seq <- 1:length(wooster)
evgam.df <- data.frame(
  y = negated.wooster,
  sin.time = sin(2*pi*time.seq / 365),
  cos.time = cos(2*pi*time.seq / 365)
)
evgam.wooster <- y ~ 1 + cos.time + sin.time
ald.wooster.80 <- evgam(evgam.wooster, data = evgam.df, family = "ald", ald.args = list(tau = 0.80))
u.vec.80 <- (predict(ald.wooster.80)$location)
ald.wooster.85 <- evgam(evgam.wooster, data = evgam.df, family = "ald", ald.args = list(tau = 0.85))
u.vec.85 <- (predict(ald.wooster.85)$location)
ald.wooster.90 <- evgam(evgam.wooster, data = evgam.df, family = "ald", ald.args = list(tau = 0.90))
u.vec.90 <- (predict(ald.wooster.90)$location)
ald.wooster.91 <- evgam(evgam.wooster, data = evgam.df, family = "ald", ald.args = list(tau = 0.91))
u.vec.91 <- (predict(ald.wooster.91)$location)
ald.wooster.95 <- evgam(evgam.wooster, data = evgam.df, family = "ald", ald.args = list(tau = 0.95))
u.vec.95 <- (predict(ald.wooster.95)$location)
ald.wooster.99 <- evgam(evgam.wooster, data = evgam.df, family = "ald", ald.args = list(tau = 0.99))
u.vec.99 <- (predict(ald.wooster.99)$location)

wooster.df <- data.frame(y = negated.wooster,
                         u.80 = u.vec.80,
                         u.85 = u.vec.85,
                         u.90 = u.vec.90,
                         u.91 = u.vec.91,
                         u.95 = u.vec.95,
                         u.99 = u.vec.99,
                         time = as.Date("1983-01-01") + 0:(length(wooster)-1))

ggplot(data = wooster.df) +
  geom_point(aes(y=y, x=time), size= 1.2) +
  # geom_line(aes(y=u.80, x = time), linewidth = 0.8, color = "red") +
  # geom_line(aes(y=u.85, x = time), linewidth = 0.8, color = "orange") +
  # geom_line(aes(y=u.90, x = time), linewidth = 0.8, color = "darkgreen") +
  geom_line(aes(y=u.91, x = time), linewidth = 1) +
  # geom_line(aes(y=u.95, x = time), linewidth = 0.8, color = "steelblue") +
  # geom_line(aes(y=u.99, x = time), linewidth = 0.8, color = "purple") +
  ylab("Daily Minimum Temperature (Degree Below 0F)") + xlab("Year") +
  scale_x_date(
    breaks = as.Date(c("1983-01-01", "1984-01-01", "1985-01-01",
                       "1986-01-01", "1987-01-01", "1988-01-01")),
    date_labels = "%Y"
  ) +
  scale_y_continuous(
    breaks = seq(-60, 20, by = 20)
  ) +  
  theme_minimal(base_size = 15) + 
  theme(legend.position = "none",
          axis.title = element_text(size = 27), 
          axis.text = element_text(size = 30))
