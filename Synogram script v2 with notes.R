library(cowplot)
library(tidyverse)
library(DescTools)
df = read_csv('PAO1-Tobra-P67S-re.csv')
df$well= factor(df$well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))


growthcurves_plot = ggplot(df, aes(x=time_hours, y=mean_od)) +
  geom_line(color='#0072B2') +
  scale_x_continuous('Time (hours)') +
  scale_y_continuous('OD600') +
  facet_wrap(~well, nrow=8) +
  theme_bw()


#Let's first mimic the heatmaps in the synogram paper. 

#To render a heatmap, we need to know the values for the row and column, so first we need 
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


# In the synogram paper, they picked an arbitrary time point for generating
# their comparisons so we'll do the same. Let's pick the 8 hour timepoint

df_8hrs = data_df %>% filter(time_hours==8)


# Just changed the A1 to A2 here as I think that's your positive control well.
ctrl_mean_od = df_8hrs %>% filter(well=='A2') %>% pull(normalised_mean_od)


# The values they use in the paper in their heatmap are
# calculated as
# (OD_growthcontrol - OD_treatment)/(OD_growthcontrol) X 100
# so let's do the same

df_8hrs = df_8hrs %>% mutate(reduction= ((ctrl_mean_od - normalised_mean_od)/ctrl_mean_od)*100)


# Now we can plot the heatmap using geom_tile
# You'll notice in the paper that they have a three-point scale for their color. 
# negative numbers go from black to orange, positive numbers go from orange to white.
# scale_fill_gradient colors the tiles in, from low color to high color.
# scale_fill_gradient2 allows for a midpoint, so that's what we'll use here.
# you can find a list of named colors here: http://derekogle.com/NCGraphing/resources/colors
# or you can use the hex values.

# I've picked the dark orange from the Okabe and Ito palette here to mimic the synogram paper
# Check out this link for color use https://clauswilke.com/dataviz/color-pitfalls.html

# I've also added some tweaks to make it look more like their plots:
# 1. The border is made thicker (panel.border)
# 2. Added a dashed line to show the just phage and just antibiotic rows.
# 3. The guide for reduction colours is made longer and formatted (the guides bit)


#This is a utility function to correctly format the x axis exponents.

scientific_10 <- function(x) {
  parse(text=gsub("\\de", " 10^", scales::scientific_format()(x)))
}

phage_exponents=c(0, 10^c(2,3,4,5,6,7,8))
antibiotic_values=c(0, 0.03125 * 2^c(0:10))

heatmap_plot = ggplot(df_8hrs, aes(x=row, y=column, fill=reduction)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient2(low='black', midpoint = 0, mid = "#D55E00", high='white') +
  scale_x_discrete('P67_Chicken_PAO1', labels=c(scientific_10(phage_exponents))) +
  scale_y_discrete(expression(paste('Tobramycin (', mu, 'g/mL)')), labels=antibiotic_values) +
  theme_bw(16) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  geom_vline(xintercept =1.5, linetype='dashed') +
  geom_hline(yintercept=1.5, linetype='dashed') +
  guides(fill = guide_colourbar(title="Reduction %", 
                                title.theme =element_text(size=16, face='bold'),
                                barwidth = 1.5, barheight = 24, 
                                frame.colour = 'black', frame.linewidth=2,
                                ticks.colour='black', ticks.linewidth=1))


# for some reason this often causes lots and lots of rows unless you put the things
# you want to keep (like row and column) in the group_by function.

# Let's look at virulence.


# work out the area under the curve
# note, because we want to keep the row and column, we just add them to the 
# summarise method

df_AUC = data_df %>% group_by(well, row, column) %>%
  summarise(AUC_trap=AUC(x=time_hours, y=normalised_mean_od, method='trapezoid'))


# changed A1 to A2 (control curve value)

ctrl_AUC = df_AUC %>% filter(well=='A2') %>% pull(AUC_trap)


#calculate virulence
#If a value is <0, we'll set it to 0.

df_AUC = df_AUC %>%
  mutate(virulence = 1 - (AUC_trap/ctrl_AUC)) %>% 
  mutate(virulence=ifelse(virulence<0, 0, virulence))


# Now see if you can do something similar for virulence on a heatmap.

virulence_plot = ggplot(df_AUC, aes(x=row, y=column, fill=virulence)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient2(low='black', midpoint = 0, mid = "#009E73", high='white') +
  geom_vline(xintercept =1.5, linetype='dashed') +
  geom_hline(yintercept=1.5, linetype='dashed') +
  scale_x_discrete('P67_Chicken_PAO1', labels=c(scientific_10(phage_exponents))) +
  scale_y_discrete(expression(paste('Tobramycin (', mu, 'g/mL)')), labels=antibiotic_values) +
  theme_bw(16) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  guides(fill = guide_colourbar(title="Virulence", 
                                title.theme =element_text(size=16, face='bold'),
                                barwidth = 1.5, barheight = 24, 
                                frame.colour = 'black', frame.linewidth=2,
                                ticks.colour='black', ticks.linewidth=1))


# finally we'll use plot_grid from the cowplot package to plot the two
# plots together (notice I stored the plots as variables above for this purpose.)

plot_grid(heatmap_plot, virulence_plot, labels = "AUTO",
          label_size=16, rel_widths = c(1, 1))






