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
#library(cowplot)
library(tibble)


# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- 'figures/Plate_stats/'
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
# species_label <- tax.c %>%
#   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='s_'| Species =='no_hit' ~as.character(Genus),
#                              TRUE ~ as.character(Species))) %>%
#   mutate(Species = case_when(Species=='unassigned' | Species =='unknown' | Species =='g_'| Species =='no_hit'~as.character(Family),
#                              TRUE ~ as.character(Species))) %>%
#   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Order),
#                              TRUE ~ as.character(Species))) %>%
#   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Class),
#                              TRUE ~ as.character(Species))) %>%
#   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character(Phylum),
#                              TRUE ~ as.character(Species))) %>%
#   mutate(Species = case_when(Species=='unassigned' | Species =='unknown'| Species =='no_hit' ~as.character('Unknown'),
#                              TRUE ~ as.character(Species)))


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
# http://www.physics.csbsju.edu/stats/KS-test.html

#create dataframe merged by ecological category

## Plates JJ and RR ----------
print('Plates JJ and RR')
df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  left_join(meta %>% select(SampleID, FilterID, PlateID)) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot, FilterID, PlateID) %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  filter(PlateID %in% c('JJ', 'RR'))


### Epipelagic and Mesopelagic Groups 
for (ecol_group in c('epipelagic', 'mesopelagic')){
  print(ecol_group)
  #merge with metadata, limit to one ecological group
  stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
    #filter(Ecological_Category == 'mesopelagic') %>%
    filter(Ecological_Category == ecol_group) %>%
    filter(diel %in% c('day', 'night')) %>%
    #filter(PlateID %in% c('CE', 'BT')) %>%
    select(depth, sum_per_tot, Ecological_Category, depth_bin2, diel)
  
  #### Shallow samples 0-99m 
  print('Depth 0-<100m')
  test_day <- stats %>% filter(depth>=0, depth<100) %>% filter(diel == 'day')
  test_night <- stats %>% filter(depth>=0, depth<100) %>% filter(diel == 'night')
  print('Number of Day Observations:')
  print(nrow(test_day))
  print('Number of Night Observations:')
  print(nrow(test_night))
  ks_stat <- ks.test(test_day$sum_per_tot, test_night$sum_per_tot)
  print(ks_stat)

  #Mann-Whitney test
  wilcox <- wilcox.test(test_day$sum_per_tot, test_night$sum_per_tot, alternative = "two.sided")
  print(wilcox)
  
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
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_shallow_JJRR.png', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_shallow_JJRR.svg', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  
  #### Deep samples 100-500m 
  print('Depth 100m-500m')
  test_day <- stats %>% filter(depth>=100, depth<=500) %>% filter(diel == 'day')
  test_night <- stats %>% filter(depth>=100, depth<=500) %>% filter(diel == 'night')
  print('Number of Day Observations:')
  print(nrow(test_day))
  print('Number of Night Observations:')
  print(nrow(test_night))
  ks_stat <- ks.test(test_day$sum_per_tot, test_night$sum_per_tot)
  print(ks_stat)
  
  #Mann-Whitney test
  wilcox <- wilcox.test(test_day$sum_per_tot, test_night$sum_per_tot, alternative = "two.sided")
  print(wilcox)
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
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_deep_JJRR.png', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_deep_JJRR.svg', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  
  
}

# Now Limit to just 300m depth to directly compare with samples in CE and BT:


print('Limit to 300m')

### Epipelagic and Mesopelagic Groups 
for (ecol_group in c('epipelagic', 'mesopelagic')){
  print(ecol_group)
  #merge with metadata, limit to one ecological group
  stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
    #filter(Ecological_Category == 'mesopelagic') %>%
    filter(Ecological_Category == ecol_group) %>%
    filter(diel %in% c('day', 'night')) %>%
    #filter(PlateID %in% c('CE', 'BT')) %>%
    select(depth, sum_per_tot, Ecological_Category, depth_bin2, diel)
  
  #### samples 100-300m 
  print('Depth 100-<=300m')
  test_day <- stats %>% filter(depth>=100, depth<=300) %>% filter(diel == 'day')
  test_night <- stats %>% filter(depth>=100, depth<=300) %>% filter(diel == 'night')
  print('Number of Day Observations:')
  print(nrow(test_day))
  print('Number of Night Observations:')
  print(nrow(test_night))
  ks_stat <- ks.test(test_day$sum_per_tot, test_night$sum_per_tot)
  print(ks_stat)
  
  #Mann-Whitney test
  wilcox <- wilcox.test(test_day$sum_per_tot, test_night$sum_per_tot, alternative = "two.sided")
  print(wilcox)
  
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
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_100_300_JJRR.png', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_100_300_JJRR.svg', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
}


## Plates CE and BT ----------
# Plates CE, BT pippin samples
# highest depth is 300m
print('Plates CE and BT')
df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  left_join(meta %>% select(SampleID, FilterID, PlateID)) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot, FilterID, PlateID) %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  filter(PlateID %in% c('CE', 'BT'))


### Epipelagic and Mesopelagic Groups 
for (ecol_group in c('epipelagic', 'mesopelagic')){
  print(ecol_group)
  #merge with metadata, limit to one ecological group
  stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
    #filter(Ecological_Category == 'mesopelagic') %>%
    filter(Ecological_Category == ecol_group) %>%
    filter(diel %in% c('day', 'night')) %>%
    #filter(PlateID %in% c('CE', 'BT')) %>%
    select(depth, sum_per_tot, Ecological_Category, depth_bin2, diel)
  
  #### Shallow samples 0-99m 
  print('Depth 0-<100m')
  test_day <- stats %>% filter(depth>=0, depth<100) %>% filter(diel == 'day')
  test_night <- stats %>% filter(depth>=0, depth<100) %>% filter(diel == 'night')
  print('Number of Day Observations:')
  print(nrow(test_day))
  print('Number of Night Observations:')
  print(nrow(test_night))
  ks_stat <- ks.test(test_day$sum_per_tot, test_night$sum_per_tot)
  print(ks_stat)
  
  #Mann-Whitney test
  wilcox <- wilcox.test(test_day$sum_per_tot, test_night$sum_per_tot, alternative = "two.sided")
  print(wilcox)
  
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
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_shallow_CEBT.png', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_shallow_CEBT.svg', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  
  #### Deep samples 100-500m 
  print('Depth 100m-300m')
  test_day <- stats %>% filter(depth>=100, depth<=500) %>% filter(diel == 'day')
  test_night <- stats %>% filter(depth>=100, depth<=500) %>% filter(diel == 'night')
  print('Number of Day Observations:')
  print(nrow(test_day))
  print('Number of Night Observations:')
  print(nrow(test_night))
  ks_stat <- ks.test(test_day$sum_per_tot, test_night$sum_per_tot)
  print(ks_stat)
  
  #Mann-Whitney test
  wilcox <- wilcox.test(test_day$sum_per_tot, test_night$sum_per_tot, alternative = "two.sided")
  print(wilcox)
  night_density_curve <- as_tibble(density(test_night$sum_per_tot, from=0, to=100)[c("x", "y")])
  day_density_curve <- as_tibble(density(test_day$sum_per_tot, from=0, to=100)[c("x", "y")])
  
  p <- full_join(night_density_curve %>% mutate(diel = 'night'), day_density_curve%>% mutate(diel = 'day')) %>%
    ggplot(aes(x=x, y=y, color=diel)) +
    geom_line(size=2) +
    ggtitle(paste(ecol_group,' 100-300m\n',
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
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_deep_CEBT.png', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  filename = paste(plot_directory, 'KS_test_',ecol_group,'_deep_CEBT.svg', sep='')
  filename
  ggsave(filename, width=8, height=6, units = 'in')
  
  
}


# # Pipe-friendly version of Wilcox test ---------------------------------------------------
# # https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/
# 
# library(rstatix)
# library(ggpubr)
# 
# # ecol_group in c('epipelagic', 'mesopelagic')
# ecol_group = 'mesopelagic'
# stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
#   #filter(Ecological_Category == 'mesopelagic') %>%
#   filter(Ecological_Category == ecol_group) %>%
#   filter(diel %in% c('day', 'night')) %>%
#   select(depth, sum_per_tot, Ecological_Category, depth_bin2, diel)
# 
# # Wilcox test across depth bins by day and night
# for (bin_depth in c('0-100m','100-200m', '300-400m', '400-500m', '500-600m')){
#   stat.test <- stats %>%
#     #filter(depth>=100, depth<=500) %>%
#     filter(depth_bin2==bin_depth) %>%
#     wilcox_test(sum_per_tot ~ diel, paired=FALSE)
#   print(bin_depth)
#   print(stat.test)
# }
# 
# # Wilcox test across "shallow" and "deep" samples by day and night
# #shallow
# stat.test <- stats %>%
#   filter(depth>0, depth<100) %>%
#   #filter(depth_bin2==bin_depth) %>%
#   wilcox_test(sum_per_tot ~ diel, paired=FALSE)
# print(stat.test)
# #deep
# stat.test <- stats %>%
#   filter(depth>=100, depth<=500) %>%
#   #filter(depth_bin2==bin_depth) %>%
#   wilcox_test(sum_per_tot ~ diel, paired=FALSE)
# print(stat.test)
# 
# # epipelagic
# ecol_group = 'epipelagic'
# stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
#   #filter(Ecological_Category == 'mesopelagic') %>%
#   filter(Ecological_Category == ecol_group) %>%
#   filter(diel %in% c('day', 'night')) %>%
#   select(depth, sum_per_tot, Ecological_Category, depth_bin2, diel)
# 
# # Wilcox test across depth bins by day and night
# for (bin_depth in c('0-100m','100-200m', '300-400m', '400-500m', '500-600m')){
#   stat.test <- stats %>%
#     #filter(depth>=100, depth<=500) %>%
#     filter(depth_bin2==bin_depth) %>%
#     wilcox_test(sum_per_tot ~ diel, paired=FALSE)
#   print(bin_depth)
#   print(stat.test)
# }
# 
# # Wilcox test across "shallow" and "deep" samples by day and night
# #shallow
# stat.test <- stats %>%
#   filter(depth>0, depth<100) %>%
#   #filter(depth_bin2==bin_depth) %>%
#   wilcox_test(sum_per_tot ~ diel, paired=FALSE)
# print(stat.test)
# #deep
# stat.test <- stats %>%
#   filter(depth>=100, depth<=500) %>%
#   #filter(depth_bin2==bin_depth) %>%
#   wilcox_test(sum_per_tot ~ diel, paired=FALSE)
# print(stat.test)

