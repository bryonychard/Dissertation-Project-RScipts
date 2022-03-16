library(cowplot)
library(tidyverse)
library(DescTools)

#load in data
df = read_csv('PAO1-Ceftazidime-P103-abc-re.csv')
df$well= factor(df$well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))

# to create factors for rows and columns (effectively splitting well into its letter and number part)
df = df %>% mutate(row=gsub('^([A-H])(.*)', '\\1', well), 
                   column=factor(gsub('^([A-H])(.*)', '\\2', well), levels=1:12))


df <- df %>% mutate(replicate=factor(replicate))


# Add phage and antibiotic concentration columns

scientific_10 <- function(x) {
  parse(text=gsub("\\de", " 10^", scales::scientific_format()(x)))
}

phage_concentrations=c(0, 0, 10^c(2,3,4,5,6,7,8))
antibiotic_values=c(0, 0, 1 * 2^c(0:10))

df = df %>% mutate(antibiotic_conc=antibiotic_values[as.integer(column)]) %>%
  mutate(phage_conc=phage_concentrations[match(row, LETTERS)])


# Plot G11 well

well_G11_df = df %>% filter(well=='G11')
ggplot(well_G11_df, aes(x=time_hours, y=mean_od, color=replicate)) +
  geom_line() +
  scale_x_continuous('Time (hours)') +
  scale_y_continuous('OD600') +
  theme_bw()

#Plot all wells

ggplot(df, aes(x=time_hours, y=mean_od, color=replicate)) +
  geom_line() +
  scale_x_continuous('Time (hours)') +
  scale_y_continuous('OD600') +
  facet_wrap(~well, nrow=8) +
  theme_bw()

library(plotrix)

# Summary statistics (min OD, max OD, mean, stdev)

data_df_AUC = df %>% group_by(well) %>% 
  summarise(AUC_trap=AUC(x=time_hours, 
                         y=mean_od,
                         method='trapezoid'))

data_df_AUC = df %>%
  group_by(well, row, column, phage_conc, antibiotic_conc) %>%
  summarise(AUC_trap=AUC(x=time_hours, y=mean_od, method='trapezoid'),
            min_od=min(mean_od),
            max_od=max(mean_od),
            mean_mean_od=mean(mean_od),
            stdev_od=sd(mean_od),
            std.error=std.error(mean_od))

stdev = data_df_AUC %>% select(stdev_od)

print((stdev)/sqrt(3))

# replicate 1 differs from others - human error (not adding phage) or something wrong with phage stock 

#linear model of log(phage_conc) against std.error 
# Show that error is greater at lower concentrations - more evolutionary capacity for resistance
#at higher concentrations everything is killed so theres less variation 
# Plot may also show you which dose optimal for reducing variance (greater concentrations have less effect because phage self-replicating)

#Change x axis to start at 0, strange theme, need to +1 for phage?

phage_concentrations=c(1, 1, 10^c(2,3,4,5,6,7,8))
data_df_AUC = data_df_AUC %>% mutate(antibiotic_conc=antibiotic_values[as.integer(column)]) %>%
  mutate(phage_conc=phage_concentrations[match(row, LETTERS)])

ggplot(data_df_AUC, aes(x=log10(phage_conc), y=std.error)) +
  geom_point() +
  geom_smooth() +
  scale_x_continuous('CPL P103 Concentration (Log10)') +
  scale_y_continuous('Standard Error') +
  theme_bw()



scale_x_continuous('CPL P98 Concentration (Log10)') +

data_df_lm <- lm(std.error ~ phage_conc, data_df_AUC)
summary(data_df_lm)





  
