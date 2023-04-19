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
plot_directory <- 'figures/Fig4/'
# Set directory to retrieve data
data_directory = "Data/filtered_seq_data/"

paletteDayNight <- c(wes_palette("Chevalier1", type = "discrete")[2], wes_palette("Darjeeling2", type = "discrete")[2], 'darkgrey')

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
  left_join(sp_desig) %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species, Ecological_Category, SampleID) %>%
  mutate(reads = sum(reads)) %>%
  mutate(per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, Ecological_Category, SampleID, .keep_all = TRUE)

# Proportion Per Time -------------

avg.time <- full_join(potu.merged, meta,  by = c("SampleID")) %>%
  # average by day first and then by diel pattern? (night_label)
  # average within depth bins by diel
  group_by(Species,night_label, depth_bin2) %>%
  mutate(reads = mean(reads)) %>%
  mutate(per_tot = mean(per_tot)) %>%
  ungroup() %>%
  distinct(Species,night_label, depth_bin2, .keep_all = TRUE)

# top 5 species:
top_taxa <- avg.time %>%
  filter(Ecological_Category %in% c('mesopelagic')) %>%
  filter(Species !='unassigned') %>%
  group_by(Species) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  distinct(Species,.keep_all = TRUE ) %>%
  arrange(-sum_per_tot) %>%
  select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
  ungroup() %>%
  select(Species, sum_per_tot) %>%
  top_n(5)

# assign text colour
textcol <- "grey40"
  
bp_top <- left_join(potu.c, meta,  by = c("SampleID")) %>% #join with metadata
  left_join(species_label,  by = c("ASV")) %>%  #join with taxonomy
  right_join(top_taxa) %>% #limit to top taxa
  # take proportion per time point
  group_by(night_label) %>%
  mutate(prop_per_tot = per_tot/sum(per_tot)) %>%
  ungroup() %>%
  ggplot(aes(x=night_label, y=prop_per_tot))+
  geom_bar(stat='identity', aes(fill = Species))+
  scale_fill_manual(values = wes_palette(5, name = "Darjeeling1", type = "continuous"), name = "") +
  labs(x="",y="Percent Total Reads")+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="bottom",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=8,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=6,colour=textcol),
    axis.title.y = element_text(size=6),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())
bp_top

filename = paste(plot_directory, marker,'_rDepth_Day_Diel_ESPCTD_prop.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 3, width =5, units = 'in')

filename = paste(plot_directory, marker,'_rDepth_Day_Diel_ESPCTD_prop.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 3, width =5, units = 'in')

# Day Night Average -------------------------------------------------------

# Bar plot of average percent total reads of highest relative abundance mesopelagic species
# (top 5)

avg.diel <- full_join(potu.merged, meta,  by = c("SampleID")) %>%
  # average by day first and then by diel pattern (night_label)
  # average within depth bins by diel
  group_by(Species,night_label, depth_bin2) %>%
  mutate(reads = mean(reads)) %>%
  mutate(per_tot = mean(per_tot)) %>%
  ungroup() %>%
  distinct(Species,night_label, depth_bin2, .keep_all = TRUE) %>%
  # average over diel Day or Night
  filter(diel %in% c('day','night')) %>%
  group_by(Species,diel, depth_bin2) %>%
  mutate(reads = mean(reads)) %>%
  mutate(per_tot = mean(per_tot)) %>%
  ungroup() %>%
  distinct(Species,night_label, depth_bin2, .keep_all = TRUE)
  
# Plot
# assign text colour
textcol <- "grey40"
  
# top 10 species:
top_taxa <- avg.diel %>%
  filter(Ecological_Category %in% c('mesopelagic')) %>%
  filter(Species !='unassigned') %>%
  group_by(Species) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  distinct(Species,.keep_all = TRUE ) %>%
  arrange(-sum_per_tot) %>%
  select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
  ungroup() %>%
  select(Species, sum_per_tot) %>%
  top_n(5)

bp_top <- avg.diel %>% 
  right_join(top_taxa) %>% #limit to top taxa
  ggplot(aes(x=depth_bin2, y=per_tot))+
  #geom_bar(aes( y=0.5),stat='identity', fill = "grey",alpha=0.8, width=20)+
  #geom_bar(data = samples, stat='identity', aes(x=SAMPLING_rdepth, y=100), color='lightgrey', fill='lightgrey', alpha=0.2,width=5)+
  geom_bar(stat='identity', aes(fill = Species))+
  coord_flip()+
  scale_x_discrete(limits=rev) +
  #scale_x_reverse()+
  #facet_wrap(~ SAMPLING_station_number)+
  facet_wrap(~diel, nrow=1) +
  #scale_colour_brewer(palette = "Set1") +
  #scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values = wes_palette(5, name = "Darjeeling1", type = "continuous"), name = "") +
  #scale_color_manual(values = wes_palette(5, name = "Darjeeling1", type = "continuous"), name = "") +
  labs(x="",y="Percent Total Reads")+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="bottom",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=8,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=6,colour=textcol),
    axis.title.y = element_text(size=6),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())
bp_top

filename = paste(plot_directory, marker,'_Meso_Diel_avgDayfirst.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 3, width =5, units = 'in')

filename = paste(plot_directory, marker,'_Meso_Diel_avgDayfirst.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 3, width =5, units = 'in')

# Number of Samples by Depth Bin--------

summary_samples <- meta %>%
  mutate(count =1 ) %>%
  group_by(night_label) %>%
  mutate(count = sum(count)) %>%
  #mutate(per_tot = mean(per_tot)) %>%
  ungroup() %>%
  distinct(night_label,count, .keep_all = FALSE)

summary_samples_depth <- meta %>%
  mutate(count =1 ) %>%
  filter(diel %in% c('day', 'night')) %>%
  group_by(diel, depth_bin2) %>%
  mutate(count = sum(count)) %>%
  #mutate(per_tot = mean(per_tot)) %>%
  ungroup() %>%
  distinct(diel, depth_bin2, count, .keep_all = FALSE)


#plots
avg.time.count <- meta %>%
  mutate(count =1 ) %>%
  # full_join(potu.merged, meta,  by = c("SampleID")) %>%
  # average by day first and then by diel pattern? (night_label)
  # average within depth bins by diel
  group_by(night_label, depth_bin2) %>%
  mutate(count = sum(count)) %>%
  #mutate(per_tot = mean(per_tot)) %>%
  ungroup() %>%
  distinct(night_label, depth_bin2, count, .keep_all = FALSE)

# assign text colour
textcol <- "grey40"
  
bp_top <- avg.time.count %>%
  ggplot(aes(x=night_label, y=count))+
  geom_bar(stat='identity', aes(fill = depth_bin2))+
  #scale_fill_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Spectral") +
  #scale_fill_manual(values = wes_palette(5, name = "Darjeeling1", type = "continuous"), name = "") +
  labs(x="",y="sample count")+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="bottom",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=8,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=6,colour=textcol),
    axis.title.y = element_text(size=6),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())
bp_top

filename = paste(plot_directory, marker,'_num_samples_bydepth.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 3, width =5, units = 'in')

filename = paste(plot_directory, marker,'_num_samples_bydepth.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 3, width =5, units = 'in')

# Number of Samples by Type of sampling -----

avg.time.count <- meta %>%
  mutate(count =1 ) %>%
  # full_join(potu.merged, meta,  by = c("SampleID")) %>%
  # average by day first and then by diel pattern? (night_label)
  # average within depth bins by diel
  group_by(night_label, ESP) %>%
  mutate(count = sum(count)) %>%
  #mutate(per_tot = mean(per_tot)) %>%
  ungroup() %>%
  distinct(night_label, ESP, count, .keep_all = FALSE)

# assign text colour
textcol <- "grey40"
  
bp_top <- avg.time.count %>%
  #left_join(potu.c, meta,  by = c("SampleID")) %>% #join with metadata
  #left_join(species_label,  by = c("ASV")) %>%  #join with taxonomy
  #right_join(top_taxa) %>% #limit to top taxa
  # take proportion per time point
  #group_by(night_label) %>%
  #mutate(prop_per_tot = per_tot/sum(per_tot)) %>%
  #ungroup() %>%
  ggplot(aes(x=night_label, y=count))+
  geom_bar(stat='identity', aes(fill = ESP))+
  scale_fill_manual(values = wes_palette(5, name = "Darjeeling2", type = "continuous"), name = "") +
  labs(x="",y="sample count")+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="bottom",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=8,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=6,colour=textcol),
    axis.title.y = element_text(size=6),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())
bp_top

filename = paste(plot_directory, marker,'_num_samples_bytype.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 3, width =5, units = 'in')

filename = paste(plot_directory, marker,'_num_samples_bytype.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 3, width =5, units = 'in')
