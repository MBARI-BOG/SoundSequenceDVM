#061523
#kpitz

# Test significance of day - night values at different depths

# Load Libraries -----------------------------------------------------------------
library(readr) #read csv files
library(lubridate) #for date modifications
library(dplyr)
library(ggplot2)
library(ggthemes)
library(magrittr)
library(tidyr)
library(RColorBrewer) #colors for plotting
library(forcats) 
library(stringr)
library(wesanderson)
library(cowplot)
library(tibble)


# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- 'figures/Fig2/'
# Set directory to retrieve data
data_directory = "Data/filtered_seq_data/"


paletteDayNight <- c(wes_palette("Chevalier1", type = "discrete")[2], wes_palette("Darjeeling2", type = "discrete")[2], 'darkgrey')
paletteEcolCat <- c('blue3', 'darkorchid', 'deepskyblue','chartreuse' )

# Import Data -------------------------------------------------------------

#ASV table
print('ASV table')
file = paste("CN19S_",marker,"_otu_filtered.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
otu.c <- read_csv(filepath) %>% rename('ASV' = 1)

#taxa table
print('taxa table')
file = paste("CN19S_",marker,"_taxa_filtered.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
tax.c <- read_csv(filepath) %>% rename('ASV' = 1)

#metadata table
print('metadata table')
#file = "A2W_12S_meta_Filtered.csv"
file = paste("CN19S_",marker,"_meta_filtered.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
samp.c <- read_csv(filepath) %>% rename('SampleID' = 1)

#OTU table long format with percent total reads
potu.c <- otu.c %>%
  tidyr::pivot_longer( -ASV, names_to ='SampleID',values_to = 'reads' ) %>%
  group_by(SampleID) %>%
  mutate(per_tot = reads / sum(reads) *100) %>%
  ungroup() %>%
  arrange(-reads)
head(potu.c)


#get lowest taxonomic annotation level for ASVs
species_label <- tax.c %>%
  mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='s_'| Species =='no_hit' ~as.character(Genus),
                             TRUE ~ as.character(Species))) %>%
  mutate(Species = case_when(Species=='unassigned' | Species =='unknown' | Species =='g_'| Species =='no_hit'~as.character(Family),
                             TRUE ~ as.character(Species))) %>%
  mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Order),
                             TRUE ~ as.character(Species))) %>%
  mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Class),
                             TRUE ~ as.character(Species))) %>%
  mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Phylum),
                             TRUE ~ as.character(Species))) %>%
  mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character('Unknown'),
                             TRUE ~ as.character(Species)))


# Import species ecological categories that have been manually determined

filepath = "data/metadata/CN19S_Taxa_Categories.csv"
# species designations tibble:
sp_desig <- read_csv(filepath)

# Create Seasonal Variables -----------------------------------------------

meta <- samp.c %>% 
  mutate(time = ymd_hms(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         str_detect(SampleID, '_V')==TRUE ~'ROV',
                         TRUE~'CTD')) %>%
  mutate(month =  month(time)) %>%
  mutate(hour = hour(time)) %>%
  mutate(day =  day(time)) %>%
  mutate(year =  year(time)) %>%
  mutate(jday = yday(time)) %>%
  mutate(month_char = as.character(month)) %>%
  mutate(year_char = as.character(year)) %>%
  #more detailed depth bins
  mutate(depth_bin = case_when(depth <=25 ~ "00_0-25m",
                               depth >25 & depth <=50 ~ "01_25-50m",
                               depth >50 & depth <=75 ~ "02_50-75m",
                               depth >75 & depth <=100 ~ "03_75-100m",
                               depth >100 & depth <=150 ~ "04_100-150m",
                               depth >150 & depth <=200 ~ "05_150-200m",
                               depth >200 & depth <=250 ~ "06_200-250m",
                               depth >250 & depth <=300 ~ "07_250-300m",
                               depth >300 & depth <=400 ~ "08_300-400m",
                               depth >400 & depth <=500 ~ "09_400-500m",
                               depth >400 & depth <=600 ~ "10_500-600m",
                               depth >600 & depth <=750 ~ "11_600-750m", TRUE ~ "unknown"
  )) 

# Also include different depth bin

meta %<>% mutate(depth_bin2 = case_when(depth <=100 ~ "0-100m",
                                        depth >100 & depth <=200 ~ "100-200m",
                                        depth >200 & depth <=300 ~ "200-300m",
                                        depth >300 & depth <=400 ~ "300-400m",
                                        depth >400 & depth <=500 ~ "400-500m",
                                        depth >400 & depth <=600 ~ "500-600m",
                                        depth >600 & depth <=750 ~ "600-700m", TRUE ~ "unknown"
)) 


# KS and Wilcox Tests ----------------

#create dataframe merged by ecological category
df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  left_join(meta %>% select(SampleID, FilterID)) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot, FilterID) %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic'))


## Epipelagic and Mesopelagic Groups -------------------------
for (ecol_group in c('epipelagic', 'mesopelagic')){
  #merge with metadata, limit to one ecological group
  stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
    #filter(Ecological_Category == 'mesopelagic') %>%
    filter(Ecological_Category == ecol_group) %>%
    filter(diel %in% c('day', 'night')) %>%
    select(depth, sum_per_tot, Ecological_Category, depth_bin2, diel)
  
  ### Shallow samples 0-100m ------------------
  test_day <- stats %>% filter(depth>=0, depth<=99) %>% filter(diel == 'day')
  test_night <- stats %>% filter(depth>=0, depth<=99) %>% filter(diel == 'night')
  
  ks_stat <- ks.test(test_day$sum_per_tot, test_night$sum_per_tot)
  ks_stat

  #Mann-Whitney test
  wilcox <- wilcox.test(test_day$sum_per_tot, test_night$sum_per_tot, alternative = "two.sided")
  
  night_density_curve <- as_tibble(density(test_night$sum_per_tot, from=0, to=100)[c("x", "y")])
  day_density_curve <- as_tibble(density(test_day$sum_per_tot, from=0, to=100)[c("x", "y")])
  
  p <- full_join(night_density_curve %>% mutate(diel = 'night'), day_density_curve%>% mutate(diel = 'day')) %>%
    ggplot(aes(x=x, y=y, color=diel)) +
    geom_line(size=2) +
    ggtitle(paste(ecol_group,' 0-99m\n',
                  'KS pvalue: ', ks_stat[2],
                  '; statistic: ', ks_stat[1],
                  '\nWilcox pvalue: ', wilcox[3],
                  '; statistic: ', wilcox[1],
                  '\n#day observations:', nrow(test_day),
                  ' #night observations:', nrow(test_night),
                  sep=''))+
    labs(x='Percent Total Reads', y='Kernel Density Estimate') +
    scale_color_manual(values=paletteDayNight )+
    scale_fill_manual(values=paletteDayNight)+
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"))
  p
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_shallow.png', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_shallow.svg', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  
  ###  Deep samples 100-500m ------------------
  test_day <- stats %>% filter(depth>=100, depth<=500) %>% filter(diel == 'day')
  test_night <- stats %>% filter(depth>=100, depth<=500) %>% filter(diel == 'night')
  ks_stat <- ks.test(test_day$sum_per_tot, test_night$sum_per_tot)
  ks_stat

  #Mann-Whitney test
  wilcox <- wilcox.test(test_day$sum_per_tot, test_night$sum_per_tot, alternative = "two.sided")
  
  night_density_curve <- as_tibble(density(test_night$sum_per_tot, from=0, to=100)[c("x", "y")])
  day_density_curve <- as_tibble(density(test_day$sum_per_tot, from=0, to=100)[c("x", "y")])
  
  p <- full_join(night_density_curve %>% mutate(diel = 'night'), day_density_curve%>% mutate(diel = 'day')) %>%
    ggplot(aes(x=x, y=y, color=diel)) +
    geom_line(size=2) +
    ggtitle(paste(ecol_group,' 100-500m\n',
                  'KS pvalue: ', ks_stat[2],
                  '; statistic: ', ks_stat[1],
                  '\nWilcox pvalue: ', wilcox[3],
                  '; statistic: ', wilcox[1],
                  '\n#day observations:', nrow(test_day),
                  ' #night observations:', nrow(test_night),
                  sep=''))+
    labs(x='Percent Total Reads', y='Kernel Density Estimate') +
    scale_color_manual(values=paletteDayNight )+
    scale_fill_manual(values=paletteDayNight)+
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"))
  p
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_deep.png', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_deep.svg', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  
  
}


  
#  Look at stats ------

meta %>% distinct(FilterID, .keep_all = TRUE) %>% filter(diel %in% c('day', 'night')) %>% glimpse()

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  left_join(meta %>% select(SampleID, FilterID)) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot, FilterID) %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic'))


stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
  mutate(count=1) %>%
  group_by(diel, Ecological_Category, depth_bin2) %>%
  mutate(median = median(sum_per_tot)) %>%
  mutate(mean = mean(sum_per_tot)) %>%
  mutate(max = max(sum_per_tot)) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(diel, Ecological_Category, depth_bin2, median, mean, sum_count, max, .keep_all = FALSE) %>%
  # make values for segments, end is current median, start is median from previous depth
  mutate(xend = NaN) %>% # will hold depth bin values
  #mutate(yend = NaN) %>% # will hold per tot values
  mutate(xend = case_when(depth_bin2 == "0-100m" ~ "0-100m",
                          depth_bin2 == "100-200m" ~ "0-100m",
                          depth_bin2 == "200-300m" ~ "100-200m",
                          depth_bin2 == "300-400m" ~ "200-300m",
                          depth_bin2 == "400-500m" ~ "300-400m",
                          depth_bin2 == "500-600m" ~ "400-500m",
                          depth_bin2 == "600-700m" ~ "500-600m", TRUE ~ "unknown"))

# ANOVA

stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  filter(Ecological_Category == 'mesopelagic') %>%
  #filter(Ecological_Category == 'epipelagic') %>%
  # filter(depth_bin2 =='100-200m') %>%
  #mutate(hour = as.integer(hour)) %>%
  #mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  select(depth,sum_per_tot, Ecological_Category, depth_bin2, diel) #%>%
  #mutate(count=1)

library(ggridges)
ggplot(stats, aes(y = depth_bin2, x = sum_per_tot)) + geom_density_ridges()+ facet_wrap(~diel)
ggplot(stats, aes(y = diel, x = sum_per_tot)) + geom_density_ridges()+ facet_wrap(~depth_bin2)


# plot 1: Density 
stats %>%
  filter(depth_bin2 =='0-100m') %>%
  ggplot(aes(x=sum_per_tot, group=diel, fill=diel, alpha=0.5)) + 
  #facet_wrap(~depth_bin2)+
  #geom_density(adjust=1.5) +
  geom_density() +
  ggtitle('0-100m')

stats %>%
  filter(depth_bin2 =='100-200m') %>%
  ggplot(aes(x=sum_per_tot, group=diel, fill=diel, alpha=0.5)) + 
  #facet_wrap(~depth_bin2)+
  #geom_density(adjust=1.5) +
  geom_density() +
  ggtitle('100-200m')

stats %>%
  filter(depth_bin2 =='200-300m') %>%
  ggplot(aes(x=sum_per_tot, group=diel, fill=diel, alpha=0.5)) + 
  #facet_wrap(~depth_bin2)+
  #geom_density(adjust=1.5) +
  geom_density() +
  ggtitle('200-300m')

stats %>%
  #filter(depth_bin2 =='200-300m') %>%
  ggplot(aes(x=sum_per_tot, group=depth_bin2, fill=diel, alpha=0.5)) + 
  #facet_wrap(~depth_bin2)+
  #geom_density(adjust=1.5) +
  geom_density() 


# try hex plot - looks neat.   ###########
library(hexbin)
stats %>% 
  #filter(depth<400) %>%
  ggplot(aes(y=sum_per_tot, x=depth)) +
  geom_hex(bins=10) +
  geom_smooth(aes(color=diel)) +
  #facet_wrap(~diel)+
  #scale_x_log10() +
  scale_y_sqrt() +
  scale_fill_viridis_c() 

x <- rnorm(50)
y <- runif(30)
# Do x and y come from the same distribution?
ks.test(x, y)


# Imitating a ridgeline plot
stats %>%
  filter(depth<400) %>%
  ggplot(aes(sum_per_tot, colour = factor(depth_bin2))) +
  geom_ribbon(
    stat = "density", outline.type = "upper",
    aes(
      fill = after_scale(alpha(colour, 0.3)),
      ymin = after_stat(group),
      ymax = after_stat(group + ndensity)
    )
  )

# calculate density then
#perform Kolmogorov-Smirnov test
ks.test(data1, data2)
#https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test

# Create refernce default `density()` curve
default_density_curve <- as_tibble(density(mtcars$mpg)[c("x", "y")])
# create density curves for day and night at single depth
#depth_val = '100-200m'
#depth_val = '200-300m'
depth_val = '0-100m'
test_day <- stats %>% filter(depth_bin2 ==depth_val) %>% filter(diel == 'day')
test_night <- stats %>% filter(depth_bin2 ==depth_val) %>% filter(diel == 'night')

night_density_curve <- as_tibble(density(test_night$sum_per_tot, from=0, to=100)[c("x", "y")])
ggplot(night_density_curve, aes(x=x, y=y)) +
  geom_line()
day_density_curve <- as_tibble(density(test_day$sum_per_tot, from=0, to=100)[c("x", "y")])
ggplot(day_density_curve, aes(x=x, y=y)) +
  geom_line()
print(depth_val)
ks_stat <- ks.test(day_density_curve$y, night_density_curve$y)
ks_stat[2]
# plot
test <- full_join(night_density_curve %>% mutate(diel = 'night'), day_density_curve%>% mutate(diel = 'day')) %>%
  ggplot(aes(x=x, y=y, color=diel)) +
  geom_line() +
  ggtitle(paste(depth_val,': KS pvalue', ks_stat[2], 'statistic', ks_stat[1]))



#depth_val = '100-200m'
#depth_val = '200-300m'
#depth_val = '0-100m'
depth_val = '300-400m'
#test_day <- stats %>% filter(depth_bin2 ==depth_val) %>% filter(diel == 'day')
#test_night <- stats %>% filter(depth_bin2 ==depth_val) %>% filter(diel == 'night')

## This works  #########
#https://www.statology.org/kolmogorov-smirnov-test-r/#:~:text=The%20Kolmogorov%2DSmirnov%20test%20is,use%20this%20function%20in%20practice.

stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  #filter(Ecological_Category == 'mesopelagic') %>%
  filter(Ecological_Category == 'epipelagic') %>%
  # filter(depth_bin2 =='100-200m') %>%
  #mutate(hour = as.integer(hour)) %>%
  #mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  select(depth,sum_per_tot, Ecological_Category, depth_bin2, diel)

test_day <- stats %>% filter(depth>=100, depth<=500) %>% filter(diel == 'day')
test_night <- stats %>% filter(depth>=100, depth<=500) %>% filter(diel == 'night')
test_day <- stats %>% filter(depth>=0, depth<=99) %>% filter(diel == 'day')
test_night <- stats %>% filter(depth>=0, depth<=99) %>% filter(diel == 'night')
ks_stat <- ks.test(test_day$sum_per_tot, test_night$sum_per_tot)
ks_stat
#https://www.graphpad.com/guides/prism/latest/statistics/interpreting_results_kolmogorov-smirnov_test.htm
#Mann-Whitney test
wilcox.test(test_day$sum_per_tot, test_night$sum_per_tot, alternative = "two.sided")

night_density_curve <- as_tibble(density(test_night$sum_per_tot, from=0, to=100)[c("x", "y")])
ggplot(night_density_curve, aes(x=x, y=y)) +
  geom_line()
day_density_curve <- as_tibble(density(test_day$sum_per_tot, from=0, to=100)[c("x", "y")])
ggplot(day_density_curve, aes(x=x, y=y)) +
  geom_line()
test <- full_join(night_density_curve %>% mutate(diel = 'night'), day_density_curve%>% mutate(diel = 'day')) %>%
  ggplot(aes(x=x, y=y, color=diel)) +
  geom_line() +
  ggtitle(paste('100-500m',': KS pvalue', ks_stat[2], 'statistic', ks_stat[1]))
test


library(hexbin)
stats %>% 
  #filter(depth<400) %>%
  filter(depth>=100, depth<400) %>%
  ggplot(aes(y=sum_per_tot, x=-depth)) +
  geom_hex(bins=10) +
  geom_smooth(aes(color=diel)) +
  #facet_wrap(~diel)+
  #scale_x_log10() +
  scale_y_sqrt() +
  coord_flip()+
  scale_fill_viridis_c() 

stats %>% 
  #filter(depth<400) %>%
  filter(depth>=0, depth<400) %>%
  ggplot(aes(y=sum_per_tot, x=depth)) +
  geom_hex(bins=10) +
  geom_smooth(aes(color=diel)) +
  #facet_wrap(~diel)+
  #scale_x_log10() +
  scale_y_sqrt() +
  scale_fill_viridis_c() 



library(ggridges)
ggplot(stats, aes(y = depth_bin2, x = sum_per_tot)) + geom_density_ridges()+ facet_wrap(~diel)
ggplot(stats, aes(y = diel, x = sum_per_tot)) + geom_density_ridges()+ facet_wrap(~depth_bin2)

ggplot(stats, aes(x = depth_bin2, y = sum_per_tot))+
  #geom_density_ridges()+
  #facet_wrap(~depth_bin2)+
  geom_violin()+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)


ggplot(stats, aes(y = diel, x = sum_per_tot))+
  #geom_density_ridges()+
  #facet_wrap(~depth_bin2)+
  geom_violin()+ 
  geom_dotplot(binaxis='x', stackdir='center', dotsize=0.2)

stats %>%
  filter(depth <400) %>%
  ggplot(aes(x = depth, y = sum_per_tot, color=diel))+
  geom_point() +
  geom_smooth(span=0.7)


#MW U
res <- wilcox.test(sum_per_tot~ diel,
                   data = stats,
                   exact = FALSE)
res

lapply(split(stats, factor(stats$depth_bin2)), function(x)wilcox.test(data=x, sum_per_tot ~ diel, exact = FALSE))


lapply(split(stats, factor(stats$depth_bin2)), function(x)t.test(data=x, sum_per_tot ~ diel, paired=FALSE))


# Pairwise comparisons
pwc <- stats %>%
  pairwise_t_test(sum_per_tot ~ diel, p.adjust.method = "bonferroni")
pwc


ad_aov <- aov(sum_per_tot ~ diel+depth_bin2 , data = stats)
# look at effects and interactions
summary(ad_aov)

library(broom)
# this extracts ANOVA output into a nice tidy dataframe
tidy_ad_aov <- tidy(ad_aov)

# call and save the pair.t.test
ad_pairwise <- pairwise.t.test(stats$sum_per_tot,stats$diel:stats$depth_bin2, 
                               p.adj = "none")
# call and tidy the tukey posthoc
tidy_ad_tukey <- tidy(TukeyHSD(ad_aov, which ='diel:depth_bin2'))

# Create Diel Boxplot of Mesopelagic and Epipelagic EcolCats --------------

meta %>% distinct(FilterID, .keep_all = TRUE) %>% filter(diel %in% c('day', 'night')) %>% glimpse()

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  left_join(meta %>% select(SampleID, FilterID)) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot, FilterID) %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic'))

stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
  mutate(count=1) %>%
  group_by(diel, Ecological_Category, depth_bin2) %>%
  mutate(median = median(sum_per_tot)) %>%
  mutate(max = max(sum_per_tot)) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(diel, Ecological_Category, depth_bin2, median, sum_count, max) %>%
  # make values for segments, end is current median, start is median from previous depth
  mutate(xend = NaN) %>% # will hold depth bin values
  mutate(xend = case_when(depth_bin2 == "0-100m" ~ "0-100m",
                          depth_bin2 == "100-200m" ~ "0-100m",
                          depth_bin2 == "200-300m" ~ "100-200m",
                          depth_bin2 == "300-400m" ~ "200-300m",
                          depth_bin2 == "400-500m" ~ "300-400m",
                          depth_bin2 == "500-600m" ~ "400-500m",
                          depth_bin2 == "600-700m" ~ "500-600m", TRUE ~ "unknown"))
#make yend df
# xend as index, merge back in to get median value
stats2 <- stats %>%
  left_join(stats %>% select(depth_bin2,median, Ecological_Category, diel) %>% rename(yend=median) %>% rename(xend=depth_bin2))

# assign text colour
textcol <- "grey40"
print("Begin plotting...")

p2 <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
  #add in stats df
  full_join(stats2) %>%
  #now plot
  ggplot(aes(x=Ecological_Category, y=sum_per_tot, fill=diel)) +
  geom_boxplot(aes(x=fct_rev(depth_bin2)), alpha=0.5)+
  geom_point(aes(y=median, x=depth_bin2, color=diel)) +
  # add in count of number of samples
  geom_text(data=.%>%filter(Ecological_Category=='mesopelagic') %>% distinct(depth_bin2,diel, .keep_all = TRUE), aes(label=sum_count, y=110, x=depth_bin2, color=diel), position = position_dodge(width = 0.9))+
  geom_segment(aes(y=median,yend=yend, x=depth_bin2,xend=xend, color=diel)) +
  facet_grid(. ~ Ecological_Category) +
  #formatting...
  scale_color_manual(values=paletteDayNight )+
  scale_fill_manual(values=paletteDayNight)+
  coord_flip() +
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  labs(y="Percent Total Reads",x="Depth Bin")+
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=8,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8,colour=textcol),
    axis.text.y=element_text(size=8,colour=textcol),
    plot.background=element_blank(),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())

p2

filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

#  Look at stats ------
stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
  mutate(count=1) %>%
  group_by(diel, Ecological_Category, depth_bin2) %>%
  mutate(median = median(sum_per_tot)) %>%
  mutate(mean = mean(sum_per_tot)) %>%
  mutate(max = max(sum_per_tot)) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(diel, Ecological_Category, depth_bin2, median, mean, sum_count, max, .keep_all = FALSE) %>%
  # make values for segments, end is current median, start is median from previous depth
  mutate(xend = NaN) %>% # will hold depth bin values
  #mutate(yend = NaN) %>% # will hold per tot values
  mutate(xend = case_when(depth_bin2 == "0-100m" ~ "0-100m",
                          depth_bin2 == "100-200m" ~ "0-100m",
                          depth_bin2 == "200-300m" ~ "100-200m",
                          depth_bin2 == "300-400m" ~ "200-300m",
                          depth_bin2 == "400-500m" ~ "300-400m",
                          depth_bin2 == "500-600m" ~ "400-500m",
                          depth_bin2 == "600-700m" ~ "500-600m", TRUE ~ "unknown"))

# Join plots --------------

#Put all together
pf1 <- p1  + theme_minimal() + theme(text = element_text(size=14), 
                                    legend.position = 'none',
                                    axis.line = element_line(color = "black"),
                                    axis.ticks = element_line(color = "black"),
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    axis.title.y=element_blank(),
                                    axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()
)
pf2 <- p2  + theme_minimal() + theme(text = element_text(size=14), 
                                    legend.position = 'none',
                                    axis.line = element_line(color = "black"),
                                    axis.ticks = element_line(color = "black"),
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()
)


plot_grid(pf2, pf1,ncol = 2, align = "h", axis = "bt", rel_widths = c(1, 0.5))

filename = paste(plot_directory, 'Acoustic_EpiMeso_Boxplot.png', sep='')
filename
ggsave(filename, width=8, height=6, units = 'in')
filename = paste(plot_directory, 'Acoustic_EpiMeso_Boxplot.svg', sep='')
filename
ggsave(filename, width=8, height=6, units = 'in')



# scrap -------------------------------------------------------------------

# Import Acoustic data -------------------------------------------

filepath <- 'data/acoustic_data/Spring canon acoustics summaries_bydepth_overall.csv'
ac_sum <- read_csv(filepath) 

# # Plot Acoustic Data -----------------------------------------
# glimpse(ac_sum)
# p1 <- ac_sum %>%
#   ggplot(aes(x = -depth, y=mean, color=diel, group=diel), alpha=0.4)+
#   geom_line(alpha=0.8)+
#   geom_errorbar(aes(ymin=mean-SD, ymax =mean+SD), width=.4,alpha=0.8)+
#   geom_point(aes(shape=diel), alpha=0.7, size=6)+
#   #coord_cartesian(ylim = c(100, -700), expand = FALSE)+
#   scale_color_manual(values=paletteDayNight )+
#   scale_fill_manual(values=paletteDayNight)+
#   theme_minimal() +
#   coord_flip(xlim = c(-660, 60), expand = FALSE)+
#   scale_x_continuous(n.breaks = 6)+
#   #scale_x_continuous(breaks=x_ticks_hr) +
#   xlab('Depth (m)')+
#   #ylab('Mean Scattering (dB)')+
#   labs(y=expression(paste("Mean Scattering "," (dB re 1", m^-1, ')')))+
#   labs(fill='Diel', color='Diel', shape='Diel')
# 
# p1
# filename = paste(plot_directory, 'Acoustic_DayNight_point.png', sep='')
# #print('Plot of top 20 Genus average by month:')
# print(filename)
# ggsave(filename,height = 6, width =5, units = 'in')
# filename = paste(plot_directory, 'Acoustic_DayNight_point.svg', sep='')
# #print('Plot of top 20 Genus average by month:')
# print(filename)
# ggsave(filename,height = 6, width =5, units = 'in')

# Get Acoustic difference -------------
ac_diff <- ac_sum %>%
  pivot_wider(id_cols=depth, names_from = diel, values_from = mean) %>%
  mutate(log_difference = day-night)

ac_diff
