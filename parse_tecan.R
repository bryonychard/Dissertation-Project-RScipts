library(cowplot)
library(tidyverse)

df = read_csv('PAO1-Tobra-P98-re.csv')
df$wellf= factor(df$well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))
ggplot(df, aes(x=time_hours, y=mean_od)) +
  geom_line(size=2, color='#0072B2') +
  scale_x_continuous('Time (hours)') +
  scale_y_continuous('OD600', limits=c(0, NA)) +
  facet_wrap(~wellf, nrow=8, ncol=12) +
  theme_bw()


#try this

ggplot(df, aes(x=time_hours, y=mean_od)) +
  geom_line(color='#0072B2') +
  scale_x_continuous('Time (hours)') +
  scale_y_continuous('OD600', limits=c(0, NA)) +
  facet_wrap(~wellf, nrow=8, ncol=12) +
  theme_bw()
