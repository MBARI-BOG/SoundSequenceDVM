#040221
#kpitz

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
library(viridis) #color maps
library(wesanderson) #colors

# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- './figures/Plate_stats/'
# Set directory to retrieve data
#data_directory = "Data/filtered_seq_data/"
data_directory = "./Data/Dada2_seq_data/"

paletteDayNight <- c(wes_palette("Chevalier1", type = "discrete")[2], wes_palette("Darjeeling2", type = "discrete")[2], 'darkgrey')

# Import Data -------------------------------------------------------------

#ASV table
print('ASV table')
#file = paste("CN19S_",marker,"_otu_filtered.csv", sep='')
file = paste("CN19S_",marker,"_Dada2_otu_merged.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
otu.c <- read_csv(filepath) %>% rename('ASV' = 1)

#taxa table
print('taxa table')
#file = paste("CN19S_",marker,"_taxa_filtered.csv", sep='')
file = paste("CN19S_",marker,"_Dada2_taxa_merged.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
tax.c <- read_csv(filepath) %>% rename('ASV' = 1)

#metadata table
print('metadata table')
#file = "A2W_12S_meta_Filtered.csv"
#file = paste("CN19S_",marker,"_meta_filtered.csv", sep='')
file = paste("CN19S_",marker,"_Dada2_meta_merged.csv", sep='')
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

# Import species ecological categories that have been manually determined

filepath = "data/metadata/CN19S_Taxa_Categories.csv"
# species designations tibble:
sp_desig <- read_csv(filepath)

# Merge taxa table and ecological categories
# adjust Species label

species_label <- full_join(tax.c, sp_desig) %>%
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




# Create Seasonal Variables -----------------------------------------------

meta <- samp.c %>% 
  # create date/time format column from local time
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
  mutate(depth_bin = case_when(depth <=25 ~ "00_0-25m",
                               depth >25 & depth <=75 ~ "01_25-75m",
                               #depth >50 & depth <=75 ~ "02_50-75m",
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

# Alternative depth bins

meta %<>% mutate(depth_bin2 = case_when(depth <=100 ~ "0-100m",
                                        depth >100 & depth <=200 ~ "100-200m",
                                        depth >200 & depth <=300 ~ "200-300m",
                                        depth >300 & depth <=400 ~ "300-400m",
                                        depth >400 & depth <=500 ~ "400-500m",
                                        depth >400 & depth <=600 ~ "500-600m",
                                        depth >600 & depth <=750 ~ "600-700m", TRUE ~ "unknown"
)) 

# Create 'cast' IDs for ESP samples:
# Adjust time-label to work? want consistent night bins?
# label with the night of the preceding day

meta %<>% 
  # if hour is <6 and diel == 'night', then new time label should be for preceding day
  mutate(night_label = case_when(hour<6 & diel == 'night' ~ paste(format(as.Date(time - days(1)), "%m-%d"), diel, sep=' '),
                                 hour<9 & diel == 'transition' ~ paste(format(as.Date(time ), "%m-%d"), 'dawn', sep=' '),
                                 hour>15 & diel == 'transition' ~ paste(format(as.Date(time ), "%m-%d"), 'dusk', sep=' '),
                                 TRUE ~ paste(format(as.Date(time), "%m-%d"), diel, sep=' ')))
meta %<>% 
  mutate(SAMPLING_rdepth = case_when(is.na(SAMPLING_rdepth) ~ depth,
                                     TRUE ~ SAMPLING_rdepth))


# Merged ASV table -------------------------------------

# Make merged OTU table by taxonomic ID; can keep ASV column to mark unique taxon
potu.merged <- left_join(potu.c, species_label) %>%
  #left_join(sp_desig) %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, 
           Species, Ecological_Category, SampleID) %>%
  mutate(reads = sum(reads)) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, Ecological_Category, SampleID, .keep_all = TRUE)


# Proportion Mesopelagic across plates -------

# Focus on depth_bin "200-300m" = depth >200 & depth <=300
# Also show difference in 100-300m samples.

### 100-300m ----------

# assign text colour for plots
textcol <- "grey40"

# p1 <- potu.merged %>%
#   left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
#   filter(sample_type == 'environmental') %>%
#   filter(depth <=300) %>%  # try to remove depth effects
#   filter(depth >=100) %>%
#   filter(reads >0) %>%
#   filter(Phylum == 'Chordata') %>%
#   # get percent Mesopelagic in each sample
#   group_by(Ecological_Category, SampleID) %>%
#   mutate(per_tot = sum(per_tot)) %>%
#   mutate(reads = sum(reads)) %>%
#   ungroup() %>%
#   distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
#   filter(Ecological_Category=='mesopelagic') %>%
#   #filter(diel=='night') %>%
#   ggplot(aes(y=per_tot, x=PlateID, color=PlateID, group=PlateID, shape=ESP)) +
#   geom_boxplot() +
#   geom_jitter()+
#   facet_grid(~diel)+
#   ggtitle('Proportion Mesopelagic Species 100-300m') +
#   theme_minimal() +
#   guides(fill=guide_legend(ncol=2)) +
#   theme(
#     #legend
#     legend.position="bottom",legend.direction="vertical",
#     legend.text=element_text(colour=textcol,size=10,face="bold"),
#     legend.key.height=grid::unit(0.3,"cm"),
#     legend.key.width=grid::unit(0.3,"cm"),
#     legend.title=element_text(colour=textcol,size=10,face="bold"),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10,colour=textcol),
#     #axis.text.x=element_text(size=7,colour=textcol),
#     axis.text.y=element_text(size=10,colour=textcol),
#     axis.title.y = element_text(size=10),
#     plot.background=element_blank(),
#     panel.border=element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(size = .25),
#     plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))
# 
# filename = paste(plot_directory, marker,'_propMeso_byplate_100_300m.png', sep='')
# #print('Plot of top 20 Genus average by month:')
# print(filename)
# ggsave(filename,height = 6, width =6, units = 'in')

# look just at transition samples (6 or 18 hour)
p1 <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >=100) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel=='transition') %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) %>%
  ggplot(aes(y=per_tot, x=migration, color=PlateID, group=migration, shape=ESP)) +
  geom_boxplot() +
  geom_jitter()+
  #facet_grid(~diel)+
  ggtitle('Proportion Mesopelagic Species 100-300m') +
  labs(y='Percent Mesopelagic Species')+
  theme_minimal() +
  #guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=10,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=10,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=10,colour=textcol),
    axis.title.y = element_text(size=10),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))

filename = paste(plot_directory, marker,'_propMeso_byplate_100_300m_transition.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 4, width =4, units = 'in')
p1


# with STATS now 
library(ggpubr)

trans_df <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >=100) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel=='transition') %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) 

p1 <- ggboxplot(trans_df, x = "migration", y = "per_tot",
               #color = "PlateID",
               #shape = "ESP",
               palette = "jco",
               add = "jitter") +
  stat_compare_means() +
  labs(y='Percent Mesopelagic Species')

filename = paste(plot_directory, marker,'_propMeso_byplate_100_300m_transition_wilcox.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 4, width =4, units = 'in')
p1
# Change method
#p + stat_compare_means(method = "t.test")


# trans_df <- potu.merged %>%
#   left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
#   filter(sample_type == 'environmental') %>%
#   filter(depth <=300) %>%  # try to remove depth effects
#   filter(depth >=100) %>%
#   filter(reads >0) %>%
#   filter(Phylum == 'Chordata') %>%
#   # get percent Mesopelagic in each sample
#   group_by(Ecological_Category, SampleID) %>%
#   mutate(per_tot = sum(per_tot)) %>%
#   mutate(reads = sum(reads)) %>%
#   ungroup() %>%
#   distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
#   filter(Ecological_Category=='mesopelagic') %>%
#   #filter(diel=='night') %>%
#   filter(diel %in% c('day', 'night')) %>%
#   #filter(PlateID %in% c('CE', 'BT')) %>%
#   filter(PlateID %in% c('JJ', 'RR')) %>%
#   # bin morning and evening migrations
#   mutate(migration = case_when(hour==6 ~ 'morning_downwards',
#                                hour %in% c(18,19,20) ~ 'evening_upwards')) 
# 
# p1 <- ggboxplot(trans_df, x = "diel", y = "per_tot",
#                 color = "diel",
#                 #shape = "ESP",
#                 palette = "jco",
#                 add = "jitter") +
#   stat_compare_means()
# 
# # p1 <- ggboxplot(trans_df, x = "PlateID", y = "per_tot",
# #                 color = "PlateID",
# #                 #shape = "ESP",
# #                 palette = "jco",
# #                 add = "jitter") +
# #   stat_compare_means()
# 
# filename = paste(plot_directory, marker,'_propMeso_byplate_100_300m_diel_TD_wilcox.png', sep='')
# #print('Plot of top 20 Genus average by month:')
# print(filename)
# ggsave(filename,height = 4, width =4, units = 'in')
# p1

# seperate by diel group
# Day
trans_df <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >=100) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel=='day') %>%
  #filter(diel %in% c('day', 'night')) %>%
  #filter(PlateID %in% c('CE', 'BT')) %>%
  #filter(PlateID %in% c('JJ', 'RR')) %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) 

p1 <- ggboxplot(trans_df, x = "PlateID", y = "per_tot",
                color = "PlateID",
                #shape = "ESP",
                palette = "jco",
                add = "jitter") +
  stat_compare_means() +
  labs(y='Percent Mesopelagic Reads') +
  ggtitle('Day Samples')

filename = paste(plot_directory, marker,'_propMeso_byplate_100_300m_day_KW.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 4, width =4, units = 'in')
p1

#Night
trans_df <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >=100) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel=='night') %>%
  #filter(diel %in% c('day', 'night')) %>%
  #filter(PlateID %in% c('CE', 'BT')) %>%
  #filter(PlateID %in% c('JJ', 'RR')) %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) 

p1 <- ggboxplot(trans_df, x = "PlateID", y = "per_tot",
                color = "PlateID",
                #shape = "ESP",
                palette = "jco",
                add = "jitter") +
  stat_compare_means() +
  labs(y='Percent Mesopelagic Reads') +
  ggtitle('Night Samples')

filename = paste(plot_directory, marker,'_propMeso_byplate_100_300m_night_KW.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 4, width =4, units = 'in')
p1

# together
trans_df <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >=100) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel %in% c('day', 'night')) %>%
  # bin PCR methods:
  mutate(PCR = case_when(PlateID %in% c('JJ', 'RR') ~ 'Method1: Plates JJ, RR',
                         PlateID %in% c('CE', 'BT') ~ 'Method2: Plates CE, BT')) %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) 


p1 <- ggboxplot(trans_df, x = "diel", y = "per_tot",
                color = "PCR",
                #shape = "ESP",
                palette = "Dark2",
                add = "jitter") +
  stat_compare_means(aes(group = PCR)) +
  stat_compare_means(label.y = 120) +
  ggtitle('100-300m Mesopelagic Percent Reads')+
  labs(y='Percent Mesopelagic Reads')

filename = paste(plot_directory, marker,'_propMeso_byplate_100_300m_diel_PCR_stat.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 6, width =6, units = 'in')
p1

## 200-300m ------------------
# Focus on depth_bin "200-300m" = depth >200 & depth <=300
# assign text colour for plots
textcol <- "grey40"

# look just at transition samples (6 or 18 hour)
p1 <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >200) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel=='transition') %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) %>%
  ggplot(aes(y=per_tot, x=migration, color=PlateID, group=migration, shape=ESP)) +
  geom_boxplot() +
  geom_jitter()+
  #facet_grid(~diel)+
  ggtitle('Proportion Mesopelagic Species 200-300m') +
  labs(y='Percent Mesopelagic Species')+
  theme_minimal() +
  #guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=10,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=10,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=10,colour=textcol),
    axis.title.y = element_text(size=10),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))

filename = paste(plot_directory, marker,'_propMeso_byplate_200_300m_transition.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 4, width =4, units = 'in')
p1


# with STATS now 
library(ggpubr)

trans_df <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >200) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel=='transition') %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) 

p1 <- ggboxplot(trans_df, x = "migration", y = "per_tot",
                #color = "PlateID",
                #shape = "ESP",
                palette = "jco",
                add = "jitter") +
  stat_compare_means() +
  labs(y='Percent Mesopelagic Species')

filename = paste(plot_directory, marker,'_propMeso_byplate_200_300m_transition_wilcox.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 4, width =4, units = 'in')
p1

# seperate by diel group
# Day
trans_df <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >200) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel=='day') %>%
  #filter(diel %in% c('day', 'night')) %>%
  #filter(PlateID %in% c('CE', 'BT')) %>%
  #filter(PlateID %in% c('JJ', 'RR')) %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) 

p1 <- ggboxplot(trans_df, x = "PlateID", y = "per_tot",
                color = "PlateID",
                #shape = "ESP",
                palette = "jco",
                add = "jitter") +
  stat_compare_means() +
  labs(y='Percent Mesopelagic Reads') +
  ggtitle('Day Samples')

filename = paste(plot_directory, marker,'_propMeso_byplate_200_300m_day_KW.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 4, width =4, units = 'in')
p1

#Night
trans_df <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >200) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel=='night') %>%
  #filter(diel %in% c('day', 'night')) %>%
  #filter(PlateID %in% c('CE', 'BT')) %>%
  #filter(PlateID %in% c('JJ', 'RR')) %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) 

p1 <- ggboxplot(trans_df, x = "PlateID", y = "per_tot",
                color = "PlateID",
                #shape = "ESP",
                palette = "jco",
                add = "jitter") +
  stat_compare_means() +
  labs(y='Percent Mesopelagic Reads') +
  ggtitle('Night Samples')

filename = paste(plot_directory, marker,'_propMeso_byplate_200_300m_night_KW.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 4, width =4, units = 'in')
p1

# together
trans_df <- potu.merged %>%
  left_join(meta %>% select(SampleID, PlateID, depth, ESP, diel, sample_type, hour)) %>%
  filter(sample_type == 'environmental') %>%
  filter(depth <=300) %>%  # try to remove depth effects
  filter(depth >200) %>%
  filter(reads >0) %>%
  filter(Phylum == 'Chordata') %>%
  # get percent Mesopelagic in each sample
  group_by(Ecological_Category, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Ecological_Category, SampleID, per_tot, reads, depth, ESP, diel, PlateID, hour) %>%
  filter(Ecological_Category=='mesopelagic') %>%
  filter(diel %in% c('day', 'night')) %>%
  # bin PCR methods:
  mutate(PCR = case_when(PlateID %in% c('JJ', 'RR') ~ 'Method1: Plates JJ, RR',
                         PlateID %in% c('CE', 'BT') ~ 'Method2: Plates CE, BT')) %>%
  # bin morning and evening migrations
  mutate(migration = case_when(hour==6 ~ 'morning_downwards',
                               hour %in% c(18,19,20) ~ 'evening_upwards')) 


p1 <- ggboxplot(trans_df, x = "diel", y = "per_tot",
                color = "PCR",
                #shape = "ESP",
                palette = "Dark2",
                add = "jitter") +
  stat_compare_means(aes(group = PCR)) +
  stat_compare_means(label.y = 120) +
  ggtitle('200-300m Mesopelagic Percent Reads')+
  labs(y='Percent Mesopelagic Reads')

filename = paste(plot_directory, marker,'_propMeso_byplate_200_300m_diel_PCR_stat.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 6, width =6, units = 'in')
p1

