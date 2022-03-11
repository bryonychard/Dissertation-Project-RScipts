library(cowplot)
library(tidyverse)
library(DescTools)
df = read_csv('PAO1-tobramycin-P98.csv')
df$well= factor(df$well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))

growthcurves_plot = ggplot(df, aes(x=time_hours, y=mean_od)) +
  geom_line(color='#0072B2') +
  scale_x_continuous('Time (hours)') +
  scale_y_continuous('OD600') +
  facet_wrap(~well, nrow=8) +
  theme_bw()

# to create factors for rows and columns (effectively splitting well into its letter and number part).
# remember, we need to use a factor with defined levels, or it'll work alphabetically and put 10 after 1.

df = df %>% mutate(row=gsub('^([A-H])(.*)', '\\1', well), 
                   column=factor(gsub('^([A-H])(.*)', '\\2', well), levels=1:12))

#filter on column==1 for your blanks
# Then we can pull out the other data (i.e. anything NOT in column 1,
# which is !=)

blks = df %>%
  filter(column==1)


data_df = df %>%
  filter(column !=1)

ggplot(blks, aes(x=time_hours, y=mean_od)) +
  geom_line(color='#0072B2') +
  scale_x_continuous('Time (hours)') +
  scale_y_continuous('OD600') +
  facet_wrap(~well, nrow=1) +
  theme_bw()

# Now we have the blanks, we can summarise the data to get a mean blank score for each
# time_hours value

blks_summary = blks %>% group_by(time_hours) %>%
  summarise(mean_blk_od=mean(mean_od))


# Now we can add that to the data_df that contains the non-blank data

data_df = data_df %>% left_join(blks_summary) %>%
  mutate(normalised_mean_od = mean_od - mean_blk_od) %>%
  mutate(normalised_mean_od = ifelse(normalised_mean_od <0, 0, normalised_mean_od))

#This is a utility function to correctly format the x axis exponents.

scientific_10 <- function(x) {
  parse(text=gsub("\\de", " 10^", scales::scientific_format()(x)))
}

phage_concentrations=c(0, 10^c(2,3,4,5,6,7,8))
antibiotic_values=c(0, 0.03125 * 2^c(0:10))

data_df = data_df %>% mutate(antibiotic_conc=antibiotic_values[as.integer(column)-1]) %>%
  mutate(phage_conc=phage_concentrations[match(row, LETTERS)])

make_interaction_growth_plot<-function(target_phage_conc, target_antibiotic_conc){
  test = data_df %>% 
    filter(phage_conc %in% c(0, target_phage_conc) & antibiotic_conc %in% c(0,target_antibiotic_conc) & phage_conc + antibiotic_conc >0) %>%
    mutate(treatment=ifelse(phage_conc >0 & antibiotic_conc >0, 'Combination',
                            ifelse(phage_conc >0, 'Phage', 'Tobramycin'))) %>%
    mutate(treatment=factor(treatment, levels=c('Combination', 'Phage', 'Tobramycin')))
  
  ggplot(test, aes(x=time_hours, y=normalised_mean_od, color=treatment)) +
    geom_point() +
    geom_line() +
    scale_x_continuous('Time (hours)') +
    scale_y_continuous('OD600') +
    theme_bw(16) +
    ggtitle(paste(target_phage_conc, "PFU /", target_antibiotic_conc,"ug/mL"))
}

make_interaction_growth_plot(1e3, 1)

