#040221
#kpitz

# Set directory to save plots
directory <- 'figures/Mesopelagic_barplots/'


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
marker = sym("12S")

# Import Data -------------------------------------------------------------
data_directory = "Data/filtered_seq_data/"

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


# Import species designations

filepath = "/Users/kpitz/Projects/CN19all/Canon2019Spring_Taxa_Categories_KBB_KP.csv"

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

# Bar by sp designation ------------------------------------------------------------

### ALL SAMPLES ###

#plot each sample by depth by species designation
test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)


# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  labs(x="",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
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

filename = paste(directory, marker,'_EcolCat_bar_ALL.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 6, width =20, units = 'in')

### Bar plot Grouped by binned depth levels ----------

#plot each sample by depth by species designation
test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)


# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  #ggplot(aes(x=fct_reorder(SampleID, desc(depth)), y=sum_per_tot, fill=Ecological_Category)) +
  ggplot(aes(x=fct_reorder(SampleID, desc(time)), y=sum_per_tot, fill=Ecological_Category)) +
  #geom_col(position='dodge', color='black') +
  geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  facet_grid(diel ~ depth_bin, margins=TRUE, scales = "free")+
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  labs(x="",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
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

filename = paste(directory, marker,'_EcolCat_bar_ALL_depthbin.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 10, width =20, units = 'in')

## Boxplot

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  #filter(depth<600) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  #ggplot(aes(x=fct_reorder(SampleID, desc(depth)), y=sum_per_tot, fill=Ecological_Category)) +
  ggplot(aes(x=fct_reorder(SampleID, desc(time)), y=sum_per_tot, fill=Ecological_Category)) +
  #geom_col(position='dodge', color='black') +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  geom_boxplot(aes(x=Ecological_Category))+
  geom_jitter(aes(x=Ecological_Category, color=Ecological_Category), size=0.4, alpha=0.6)+
  facet_grid(diel ~ depth_bin, margins=TRUE, scales = "free")+
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  labs(x="",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="bottom",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=12,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=12,face="bold"),
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

filename = paste(directory, marker,'_EcolCat_bar_ALL_depthbin_boxplot.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 10, width =12, units = 'in')


### Seperate by Day and Night

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  #filter(diel == 'night') %>%
  ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  facet_grid(diel ~ ., margins=TRUE)+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  labs(x="",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
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

filename = paste(directory, marker,'_EcolCat_bar_ALL_DayNight.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 6, width =20, units = 'in')

# LIMIT TO CTD samples?

#plot each sample by depth by species disignation
test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

#limit samples
samp_lim <- samp.c %>%
  select(SampleID,SAMPLING_cruise, depth, ESP, SAMPLING_station, SC, SAMPLING_station_number, diel ) %>%
  #filter(depth<600) %>%
  filter(depth>=-1) %>%
  filter(SAMPLING_station %in% c('MARS', 'OFFMARS_E')) %>%
  #filter(diel == 'day') %>%
  #filter(SC!=57) %>%
  filter(ESP %in% c('KOA', 'MV1') ==FALSE)

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- right_join(test, samp_lim,  by = c("SampleID")) %>% #join with metadata
  ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  labs(x="",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
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

filename = paste(directory, marker,'_EcolCat_bar_CTD.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 6, width =20, units = 'in')

# ROV

#limit samples
samp_lim <- samp.c %>%
  select(SampleID,SAMPLING_cruise, SAMPLING_platform, depth, ESP, SAMPLING_station, SC, SAMPLING_station_number, diel ) %>%
  filter(depth>=-1) %>%
  filter(SAMPLING_platform == 'ROV') 

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- right_join(test, samp_lim,  by = c("SampleID")) %>% #join with metadata
  ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  labs(x="",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
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
filename = paste(directory, marker,'_EcolCat_bar_ROV.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 6, width =8, units = 'in')


# Scatter plot of reads by depth -----------------------------------------

test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  #filter(depth<600) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  ggplot(aes(x = depth, y = sum_per_tot)) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_point(aes(color = Ecological_Category), alpha=0.7) +
  geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  facet_grid(diel ~ ., margins=TRUE)+
  #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  labs(x="Depth (m)",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_scatter_Diel_tab10.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 10, width =6, units = 'in')

### Make Ecol categories smoothed plot with depth on the y-axis ####

test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  #filter(depth<600) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  ggplot(aes(x = -depth, y = sum_per_tot)) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_point(aes(color = Ecological_Category), alpha=0.7) +
  geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  facet_grid(diel ~ ., margins=TRUE)+
  #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  labs(x="Depth (m)",y="Percent Total Reads")+
  coord_flip() +
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_scatter_Diel_tab10_depthy.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 20, width =6, units = 'in')

----- ## Make separate graphs ------

test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  mutate(SAMPLING_platform = replace(SAMPLING_platform,SAMPLING_platform== 'daphne', 'LRAUV'))%>%
  mutate(SAMPLING_platform = replace(SAMPLING_platform,SAMPLING_platform== 'makai', 'LRAUV'))%>%
  filter(diel %in% c('day')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  ggplot(aes(x = -depth, y = sum_per_tot)) +
  geom_point(aes(color = Ecological_Category, shape=SAMPLING_platform), alpha=0.7) +
  geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5, span=0.6, level=0.95) +
  scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_shape_manual(values=c(17, 8, 1)) + 
  labs(x="Depth (m)",y="Percent Total Reads")+
  coord_flip() +
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_scatter_Diel_depthy_DAY_all_markers.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')
filename = paste(directory, marker,'_EcolCat_scatter_Diel_depthy_DAY_all_markers.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

# NIGHT
# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  mutate(SAMPLING_platform = replace(SAMPLING_platform,SAMPLING_platform== 'daphne', 'LRAUV'))%>%
  mutate(SAMPLING_platform = replace(SAMPLING_platform,SAMPLING_platform== 'makai', 'LRAUV'))%>%
  filter(diel %in% c('night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  ggplot(aes(x = -depth, y = sum_per_tot)) +
  geom_point(aes(color = Ecological_Category, shape=SAMPLING_platform), alpha=0.7) +
  geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5, span=0.6, level=0.95) +
  scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_shape_manual(values=c(17, 1)) + 
  labs(x="Depth (m)",y="Percent Total Reads")+
  coord_flip() +
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_scatter_Diel_depthy_NIGHT_all_markers.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')
filename = paste(directory, marker,'_EcolCat_scatter_Diel_depthy_NIGHT_all_markers.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

# TRANSITION
# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  mutate(SAMPLING_platform = replace(SAMPLING_platform,SAMPLING_platform== 'daphne', 'LRAUV'))%>%
  mutate(SAMPLING_platform = replace(SAMPLING_platform,SAMPLING_platform== 'makai', 'LRAUV'))%>%
  filter(diel %in% c('transition')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  ggplot(aes(x = -depth, y = sum_per_tot)) +
  geom_point(aes(color = Ecological_Category, shape=SAMPLING_platform), alpha=0.7) +
  geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5, span=0.6, level=0.95) +
  scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_shape_manual(values=c(17, 1)) + 
  labs(x="Depth (m)",y="Percent Total Reads")+
  coord_flip() +
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_scatter_Diel_depthy_TRANS_all_markers.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')
filename = paste(directory, marker,'_EcolCat_scatter_Diel_depthy_TRANS_all_markers.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

# # NIGHT
# # assign text colour
# textcol <- "grey40"
# print("Begin plotting...")
# bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
#   filter(depth>-1) %>%
#   #filter(SAMPLING_platform=='WESTERN FLYER') %>%
#   #filter(depth<600) %>%
#   #filter(diel %in% c('day', 'night', 'transition')) %>%
#   filter(diel %in% c('night')) %>%
#   filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
#                                     'cosmopolitan')) %>%
#   ggplot(aes(x = -depth, y = sum_per_tot)) +
#   #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
#   geom_point(aes(color = Ecological_Category), alpha=0.7) +
#   geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5, span=0.6, level=0.95) +
#   #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
#   #facet_grid(diel ~ ., margins=TRUE)+
#   #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
#   #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
#   #scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
#   scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
#   scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
#   #scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
#   #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
#   labs(x="Depth (m)",y="Percent Total Reads")+
#   coord_flip() +
#   #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
#   #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
#   theme_minimal() +
#   guides(fill=guide_legend(ncol=2)) +
#   theme(
#     #legend
#     legend.position="right",legend.direction="vertical",
#     legend.text=element_text(colour=textcol,size=8,face="bold"),
#     legend.key.height=grid::unit(0.3,"cm"),
#     legend.key.width=grid::unit(0.3,"cm"),
#     legend.title=element_text(colour=textcol,size=8,face="bold"),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=textcol),
#     #axis.text.x=element_text(size=7,colour=textcol),
#     axis.text.y=element_text(size=6,colour=textcol),
#     axis.title.y = element_text(size=6),
#     plot.background=element_blank(),
#     panel.border=element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(size = .25),
#     plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
#     plot.title=element_blank())
# 
# bp_top
# filename = paste(directory, marker,'_EcolCat_scatter_Diel_depthy_NIGHT_all.png', sep='')
# #print('Plot of top 20 Genus average by month:')
# print(filename)
# ggsave(filename,height = 5, width =5, units = 'in')

# # transition
# textcol <- "grey40"
# print("Begin plotting...")
# bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
#   filter(depth>-1) %>%
#   #filter(depth<600) %>%
#   #filter(diel %in% c('day', 'night', 'transition')) %>%
#   #filter(SAMPLING_platform=='WESTERN FLYER') %>%
#   filter(diel %in% c('transition')) %>%
#   filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
#                                     'cosmopolitan')) %>%
#   ggplot(aes(x = -depth, y = sum_per_tot)) +
#   #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
#   geom_point(aes(color = Ecological_Category), alpha=0.7) +
#   geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5, span=0.6, level=0.95) +
#   #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
#   #facet_grid(diel ~ ., margins=TRUE)+
#   #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
#   #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
#   #scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
#   scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
#   scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
#   #scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
#   #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
#   labs(x="Depth (m)",y="Percent Total Reads")+
#   coord_flip() +
#   #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
#   #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
#   theme_minimal() +
#   guides(fill=guide_legend(ncol=2)) +
#   theme(
#     #legend
#     legend.position="right",legend.direction="vertical",
#     legend.text=element_text(colour=textcol,size=8,face="bold"),
#     legend.key.height=grid::unit(0.3,"cm"),
#     legend.key.width=grid::unit(0.3,"cm"),
#     legend.title=element_text(colour=textcol,size=8,face="bold"),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=5,colour=textcol),
#     #axis.text.x=element_text(size=7,colour=textcol),
#     axis.text.y=element_text(size=6,colour=textcol),
#     axis.title.y = element_text(size=6),
#     plot.background=element_blank(),
#     panel.border=element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(size = .25),
#     plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
#     plot.title=element_blank())
# 
# bp_top
# filename = paste(directory, marker,'_EcolCat_scatter_Diel_depthy_TRNS_all.png', sep='')
# #print('Plot of top 20 Genus average by month:')
# print(filename)
# ggsave(filename,height = 5, width =5, units = 'in')


### Plot sampling intensity ####

textcol <- "grey40"
print("Begin plotting...")

bp_top <- samp.c %>%
  select(depth, diel, SAMPLING_platform, SampleID) %>%
  mutate(count=1) %>%
  mutate( ints = cut(depth ,breaks = 25, labels=FALSE)) %>% 
  group_by(ints, diel, SAMPLING_platform) %>% 
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(ints,diel, SAMPLING_platform,.keep_all = TRUE) %>%
  ggplot(aes(y= ints, x=SAMPLING_platform, fill=sum_count)) +
  facet_grid(. ~ diel, margins=FALSE, scales = "free")+
  geom_tile()+
  scale_y_reverse()+
  scale_fill_viridis()+
  #labs(y="Percent Total Reads")+
  theme_minimal() +
  #guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_sampling_effort_by vehicle.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =3, units = 'in')
filename = paste(directory, marker,'_sampling_effort_byvehicle.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =3, units = 'in')

# merged by diel
bp_top <- samp.c %>%
  select(depth, diel, SAMPLING_platform, SampleID) %>%
  mutate(count=1) %>%
  mutate( ints = cut(depth ,breaks = 25, labels=FALSE)) %>% 
  group_by(ints, diel, SAMPLING_platform) %>% 
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(ints,diel, SAMPLING_platform,.keep_all = TRUE) %>%
  ggplot(aes(y= ints, x=diel, fill=sum_count)) +
  #facet_grid(. ~ diel, margins=FALSE, scales = "free")+
  geom_tile()+
  scale_y_reverse()+
  scale_fill_viridis()+
  #labs(y="Percent Total Reads")+
  theme_minimal() +
  #guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_sampling_effort_all.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =3, units = 'in')
filename = paste(directory, marker,'_sampling_effort_all.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =3, units = 'in')


### Plot EcolCats by Hour and Depth ####

test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

# assign text colour
library(viridis)
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  #filter(depth<600) %>%
  filter(sum_per_tot > 0) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  #filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
  #                                  'cosmopolitan')) %>%
  filter(Ecological_Category %in% c('epipelagic')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour) %>%
  ggplot(aes(x = hour, y = -depth,size=sum_per_tot )) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_point(aes(color = sum_per_tot), alpha=0.7) +
  #geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  facet_grid(Ecological_Category ~ ., margins=FALSE)+
  scale_color_viridis()+
  #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #labs(x="Depth (m)",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_scatter_hour_depth_epi.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

#bin data and average?
bp_top <- full_join(test, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  #filter(depth<600) %>%
  filter(sum_per_tot > 0) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  #filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
  #                                  'cosmopolitan')) %>%
  filter(Ecological_Category %in% c('epipelagic')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour) %>%
  ggplot(aes(x = hour, y = -depth,z=sum_per_tot )) +
  #geom_point(aes(color = sum_per_tot), alpha=0.7) +
  stat_summary_hex(fun = function(x) sd(x), bins=10)+
  #stat_summary_hex(fun = function(x) mean(x), bins=10)+
  scale_fill_viridis(option='turbo')+
  theme_minimal() 

bp_top
filename = paste(directory, marker,'_EcolCat_scatter_hour_depth_epi_std.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

### boxplots of Ecol Cats through depth ####

test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

# assign text colour
library(viridis)
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  #filter(depth<600) %>%
  filter(sum_per_tot > 0) %>%
  filter(diel %in% c('day', 'night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  #filter(Ecological_Category %in% c('epipelagic')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin, diel) %>%
  ggplot(aes(x=Ecological_Category, y=sum_per_tot, fill=Ecological_Category)) +
  #geom_col(position='dodge', color='black') +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  geom_boxplot(aes(x=fct_rev(depth_bin)), alpha=0.5)+
  #geom_jitter(aes(x=fct_rev(depth_bin), color=Ecological_Category), size=0.7, alpha=0.8)+
  facet_grid(Ecological_Category ~ diel, scales = "free")+
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
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
filename = paste(directory, marker,'_EcolCat_boxplot_depthbin.png', sep='')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

bp_top <- full_join(test, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  #filter(depth<600) %>%
  filter(sum_per_tot > 0) %>%
  filter(diel %in% c('day', 'night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  #filter(Ecological_Category %in% c('epipelagic')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin, diel) %>%
  ggplot(aes(x=Ecological_Category, y=sum_per_tot, fill=diel)) +
  #geom_col(position='dodge', color='black') +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  geom_boxplot(aes(x=fct_rev(depth_bin)), alpha=0.5)+
  #geom_jitter(aes(x=fct_rev(depth_bin), color=diel), size=0.7, alpha=0.8)+
  #facet_grid(Ecological_Category ~ ., scales = "free")+
  facet_grid(. ~ Ecological_Category, scales = "free")+
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = -1)+
  scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = -1)+
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
filename = paste(directory, marker,'_EcolCat_boxplot_depthbin_dielside.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

### Plot same ecological category, color by Diel ######

#old way:






# New way? (not working)
smoothed_lines <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  #filter(depth<600) %>%
  #filter(diel %in% c('day', 'night')) %>%
  filter(diel %in% c('day')) %>%
  #filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
  #                                  'cosmopolitan')) %>%
  filter(Ecological_Category %in% c('mesopelagic'))
  
  
# Create model that will do the same thing as under the hood in ggplot2
model <- loess( sum_per_tot ~ -depth, data = smoothed_lines, span = 0.6)

library(broom)
# Add predicted values from model to original dataset using broom library
smoothed_lines2 <- augment(model, smoothed_lines)

# Plot both lines 
ggplot(data=smoothed_lines2, aes(-depth,sum_per_tot)) + 
  geom_line(aes(-depth, .fitted), color = "red") +
  stat_smooth(method = "loess", span = 0.6)


# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, samp.c,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  filter(SAMPLING_platform=='WESTERN FLYER') %>%
  #filter(depth<600) %>%
  filter(diel %in% c('day', 'night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  ggplot(aes(x = -depth, y = sum_per_tot)) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_point(aes(color = diel), alpha=0.4) +
  #geom_smooth(aes(color = diel, fill=diel, linetype=diel),alpha=0.5, span=0.6, level=0.95, se=FALSE) +
  geom_smooth(aes(color = diel, fill=diel, linetype=diel),alpha=0.5, span=0.6, level=0.95) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  facet_grid(Ecological_Category ~ ., margins=FALSE)+
  facet_grid(. ~ Ecological_Category, margins=FALSE)+
  #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  scale_fill_manual(values = c('chocolate1', 'royalblue3', 'darkgreen')) +
  scale_color_manual(values = c('chocolate1', 'royalblue3', 'darkgreen')) +
  labs(x="Depth (m)",y="Percent Total Reads")+
  coord_flip()+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_scatter_Diel_tab10_flipped_CTD.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 4, width =6, units = 'in')



# Ecol Scatter through depth and time -------------------------------------

test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

# assign text colour
library(viridis)
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  #filter(depth<600) %>%
  filter(sum_per_tot > 0) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  select(time, depth, sum_per_tot, Ecological_Category) %>%
  ggplot(aes(x = time, y = -depth,size=sum_per_tot )) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_point(aes(color = sum_per_tot), alpha=0.7) +
  #geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  facet_grid(Ecological_Category ~ ., margins=FALSE)+
  scale_color_viridis()+
  #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #labs(x="Depth (m)",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_scatter_depth_time_nozeros_colorpertot.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 10, width =12, units = 'in')

## color by sample type
test <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

# assign text colour
library(viridis)
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(test, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  #filter(depth<600) %>%
  filter(sum_per_tot >0) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  ggplot(aes(x = time, y = -depth,size=sum_per_tot )) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_point(aes(color = SAMPLING_platform, shape=ESP ), alpha=0.7) +
  #geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  facet_grid(Ecological_Category ~ ., margins=FALSE)+
  #scale_color_viridis()+
  #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #labs(x="Depth (m)",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_scatter_depth_time_samptype.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 10, width =12, units = 'in')



# Look at Daphne samples --------------------------------------------------

# ## color by sample type
# test <- tax.c %>% left_join(sp_desig) %>%
#   left_join(potu.c) %>%
#   group_by(Ecological_Category, SampleID) %>%
#   mutate(sum_reads = sum(reads)) %>%
#   mutate(sum_per_tot = sum(per_tot)) %>%
#   ungroup() %>%
#   distinct(SampleID, Ecological_Category, .keep_all=TRUE)

top_taxa <- inner_join(potu.c, meta,  by = c("SampleID")) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(depth>-1) %>%
  filter(depth<40) %>%
  filter(time > mdy('06/02/2019')) %>%
  filter(SAMPLING_platform == 'daphne') %>%
  group_by(Species) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  arrange(-sum_per_tot) %>%
  distinct(Species,.keep_all = TRUE ) %>%
  select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
  ungroup() %>%
  select(Species, sum_per_tot) %>%
  top_n(20)
print(top_taxa)


# assign text colour
library(viridis)
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(potu.c, meta,  by = c("SampleID")) %>% #join with metadata
  left_join(species_label) %>%
  right_join(top_taxa) %>%
  filter(depth>-1) %>%
  filter(depth<40) %>%
  filter(time > mdy_hm('06/04/2019 12:00')) %>%
  filter(SAMPLING_platform == 'daphne') %>%
  filter(sum_per_tot >0) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  #filter(Ecological_Category %in% c('mesopelagic')) %>%
  ggplot(aes(x = time, y =per_tot )) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  #geom_point(aes(color = SAMPLING_platform, shape=ESP ), alpha=0.7) +
  #geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5) +
  geom_bar(stat = "identity", aes(fill = Species))+
  #facet_grid(Ecological_Category ~ ., margins=FALSE)+
  #scale_color_viridis()+
  #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #labs(x="Depth (m)",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_barplot_daphne.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 10, width =12, units = 'in')

## raw reads

top_taxa <- inner_join(potu.c, meta,  by = c("SampleID")) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(depth>-1) %>%
  filter(depth<40) %>%
  filter(time > mdy('06/02/2019')) %>%
  filter(SAMPLING_platform == 'daphne') %>%
  group_by(Species) %>%
  mutate(sum_reads = sum(reads)) %>%
  arrange(-sum_reads) %>%
  distinct(Species,.keep_all = TRUE ) %>%
  select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_reads) %>%
  ungroup() %>%
  select(Species, sum_reads) %>%
  top_n(20)
print(top_taxa)


# assign text colour
library(viridis)
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(potu.c, meta,  by = c("SampleID")) %>% #join with metadata
  left_join(species_label) %>%
  right_join(top_taxa) %>%
  filter(depth>-1) %>%
  filter(depth<40) %>%
  filter(time > mdy_hm('06/04/2019 12:00')) %>%
  filter(SAMPLING_platform == 'daphne') %>%
  filter(sum_reads >0) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  #filter(Ecological_Category %in% c('mesopelagic')) %>%
  ggplot(aes(x = time, y =reads )) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_reads)) +
  #geom_point(aes(color = SAMPLING_platform, shape=ESP ), alpha=0.7) +
  #geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5) +
  geom_bar(stat = "identity", aes(fill = Species))+
  #facet_grid(Ecological_Category ~ ., margins=FALSE)+
  #scale_color_viridis()+
  #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #labs(x="Depth (m)",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_barplot_daphne_reads.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 10, width =12, units = 'in')

# Steno through time and depth:

# assign text colour
library(viridis)
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(potu.c, meta,  by = c("SampleID")) %>% #join with metadata
  left_join(species_label) %>%
  #right_join(top_taxa) %>%
  filter(Species == 'Stenobrachius leucopsarus') %>%
  filter(depth>-1) %>%
  filter(depth<40) %>%
  filter(time > mdy_hm('06/04/2019 12:00')) %>%
  filter(SAMPLING_platform == 'daphne') %>%
  filter(reads >0) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  #filter(Ecological_Category %in% c('mesopelagic')) %>%
  ggplot(aes(x = time, y =-depth)) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_reads)) +
  geom_point(aes(size=per_tot, color=diel ), alpha=0.7) +
  #geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5) +
  #geom_bar(stat = "identity", aes(fill = Species))+
  #facet_grid(Ecological_Category ~ ., margins=FALSE)+
  #scale_color_viridis()+
  #scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  #scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
  #labs(x="Depth (m)",y="Percent Total Reads")+
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
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
filename = paste(directory, marker,'_EcolCat_barplot_daphne__Steno_percent.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 10, width =12, units = 'in')

# plot lat/lon of LRAUV samples:

map_plot <- meta %>%
  #filter(SAMPLING_platform == 'daphne') %>%
  ggplot(aes(x=decimalLongitude, y=decimalLatitude, color=SAMPLING_platform),alpha=0.5) +
  geom_point()
map_plot


# Just Mesopelagic fishes -------------------------------------------------

meso_taxa <- tax.c %>% left_join(sp_desig) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'meso and bathypelagic'))


taxas <- c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
#### Percent Reads
for (val in taxas) {
  taxa_level = sym(val)
  top_taxa <- potu.c %>%
    right_join(meso_taxa %>% select(ASV)) %>%
    left_join(species_label) %>%
    # filter(!!taxa_level != 'Unknown') %>%
    # filter(!!taxa_level !='no_hit') %>%
    filter(!!taxa_level !='unassigned') %>%
    # filter(!!taxa_level !='unknown') %>%
    # filter(!!taxa_level !='s_') %>%
    # filter(!!taxa_level !='g_') %>%
    group_by(!!taxa_level) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    distinct(!!taxa_level,.keep_all = TRUE ) %>%
    arrange(-sum_per_tot) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    #print(n = Inf) %>%
    ungroup() %>%
    select(!!taxa_level, sum_per_tot) %>%
    top_n(10)
  
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>% #join with metadata
    inner_join(species_label,  by = c("ASV")) %>%  #join with taxonomy
    right_join(top_taxa) %>% #limit to top taxa
    ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = !!taxa_level))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y="Percent Total Reads")+
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
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
  
  filename = paste(directory, marker,'_top10',taxa_level,'_mesotaxa_bar.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 6, width =20, units = 'in')
}

####  Raw Reads
for (val in taxas) {
  taxa_level = sym(val)
  top_taxa <- potu.c %>%
    right_join(meso_taxa %>% select(ASV)) %>%
    left_join(species_label) %>%
    # filter(!!taxa_level != 'Unknown') %>%
    # filter(!!taxa_level !='no_hit') %>%
    filter(!!taxa_level !='unassigned') %>%
    # filter(!!taxa_level !='unknown') %>%
    # filter(!!taxa_level !='s_') %>%
    # filter(!!taxa_level !='g_') %>%
    group_by(!!taxa_level) %>%
    mutate(sum_per_tot = sum(reads)) %>%
    distinct(!!taxa_level,.keep_all = TRUE ) %>%
    arrange(-sum_per_tot) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    #print(n = Inf) %>%
    ungroup() %>%
    select(!!taxa_level, sum_per_tot) %>%
    top_n(10)
  
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>% #join with metadata
    inner_join(species_label,  by = c("ASV")) %>%  #join with taxonomy
    right_join(top_taxa) %>% #limit to top taxa
    ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = reads)) +
    geom_bar(stat = "identity", aes(fill = !!taxa_level))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y="Total Reads")+
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
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
  
  filename = paste(directory, marker,'_top10',taxa_level,'_mesotaxa_bar_rawreads.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 6, width =20, units = 'in')
}

#######   Just CTD samples <600m  #########
samp_lim

#### Percent Reads
for (val in taxas) {
  taxa_level = sym(val)
  top_taxa <- potu.c %>%
    right_join(meso_taxa %>% select(ASV)) %>%
    left_join(species_label) %>%
    right_join(samp_lim %>% select(SampleID)) %>%
    # filter(!!taxa_level != 'Unknown') %>%
    # filter(!!taxa_level !='no_hit') %>%
    filter(!!taxa_level !='unassigned') %>%
    # filter(!!taxa_level !='unknown') %>%
    # filter(!!taxa_level !='s_') %>%
    # filter(!!taxa_level !='g_') %>%
    group_by(!!taxa_level) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    distinct(!!taxa_level,.keep_all = TRUE ) %>%
    arrange(-sum_per_tot) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    #print(n = Inf) %>%
    ungroup() %>%
    select(!!taxa_level, sum_per_tot) %>%
    top_n(10)
  
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- inner_join(potu.c, samp_lim,  by = c("SampleID")) %>% #join with metadata
    inner_join(species_label,  by = c("ASV")) %>%  #join with taxonomy
    right_join(top_taxa) %>% #limit to top taxa
    ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = !!taxa_level))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y="Percent Total Reads")+
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
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
  
  filename = paste(directory, marker,'_top10',taxa_level,'_mesotaxa_bar_limCTD.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 6, width =20, units = 'in')
}


for (val in taxas) {
  taxa_level = sym(val)
  top_taxa <- potu.c %>%
    right_join(meso_taxa %>% select(ASV)) %>%
    left_join(species_label) %>%
    right_join(samp_lim %>% select(SampleID)) %>%
    # filter(!!taxa_level != 'Unknown') %>%
    # filter(!!taxa_level !='no_hit') %>%
    filter(!!taxa_level !='unassigned') %>%
    # filter(!!taxa_level !='unknown') %>%
    # filter(!!taxa_level !='s_') %>%
    # filter(!!taxa_level !='g_') %>%
    group_by(!!taxa_level) %>%
    mutate(sum_per_tot = sum(reads)) %>%
    distinct(!!taxa_level,.keep_all = TRUE ) %>%
    arrange(-sum_per_tot) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    #print(n = Inf) %>%
    ungroup() %>%
    select(!!taxa_level, sum_per_tot) %>%
    top_n(10)
  
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- inner_join(potu.c, samp_lim,  by = c("SampleID")) %>% #join with metadata
    inner_join(species_label,  by = c("ASV")) %>%  #join with taxonomy
    right_join(top_taxa) %>% #limit to top taxa
    ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = reads)) +
    geom_bar(stat = "identity", aes(fill = !!taxa_level))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y="Total Reads")+
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
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
  
  filename = paste(directory, marker,'_top10',taxa_level,'_mesotaxa_bar_rawreads_limCTD.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 6, width =20, units = 'in')
}

# Mesopelagic Species Scatter ---------------------------------------------

meso_taxa <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  filter(Ecological_Category %in% c('mesopelagic')) %>%
  group_by(Genus, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(Genus, SampleID, .keep_all=TRUE)


taxas <- c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
#### Percent Reads

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- inner_join(meso_taxa, samp.c,  by = c("SampleID")) %>% #join with metadata
  #inner_join(species_label,  by = c("ASV")) %>%  #join with taxonomy
  #right_join(top_taxa) %>% #limit to top taxa
  filter(sum_per_tot >0) %>%
  filter(depth>-1) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  filter(Genus %in% c('Diaphus', 'Leuroglossus', 'Stenobrachius')) %>%
  ggplot(aes(x = depth, y = sum_per_tot)) +
  #ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = sum_per_tot)) +
  geom_point(aes(color = Genus), alpha=0.7) +
  geom_smooth(aes(color = Genus, fill=Genus), alpha=0.5) +
  #geom_bar(stat = "identity", aes(fill = Ecological_Category))+
  facet_grid(diel ~ ., margins=TRUE)+
  scale_color_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  scale_fill_tableau(palette = "Classic Color Blind", type = c("regular"), direction = -1)+
  labs(x="Depth (m)",y="Percent Total Reads")+
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
filename = paste(directory, marker,'_Genus_mesotaxa_scatter.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')


# Stacked Bar by Cast -----------------------------------------------------

taxas <- c('Class','Order', 'Family', 'Genus', 'Species')
#### Percent Reads
for (val in taxas) {
  taxa_level = sym(val)
  top_taxa <- potu.c %>%
    right_join(meso_taxa %>% select(ASV)) %>%
    left_join(species_label) %>%
    right_join(samp_lim) %>%
    filter(SAMPLING_station_number >=1) %>%
    filter(SAMPLING_station=='MARS') %>%
    filter(SAMPLING_station_number %in% c('11', '13', '16', '20', '23', '27') ==FALSE) %>%
    # filter(!!taxa_level != 'Unknown') %>%
    # filter(!!taxa_level !='no_hit') %>%
    filter(!!taxa_level !='unassigned') %>%
    # filter(!!taxa_level !='unknown') %>%
    # filter(!!taxa_level !='s_') %>%
    # filter(!!taxa_level !='g_') %>%
    group_by(!!taxa_level) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    distinct(!!taxa_level,.keep_all = TRUE ) %>%
    arrange(-sum_per_tot) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    #print(n = Inf) %>%
    ungroup() %>%
    select(!!taxa_level, sum_per_tot) %>%
    top_n(10)
  
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- left_join(potu.c, samp_lim,  by = c("SampleID")) %>% #join with metadata
    left_join(species_label,  by = c("ASV")) %>%  #join with taxonomy
    right_join(top_taxa) %>% #limit to top taxa
    filter(SAMPLING_station_number >=1) %>%
    filter(SAMPLING_station=='MARS') %>%
    #arrange(SAMPLING_station_number) %>%
    mutate(SAMPLING_station_number = replace(SAMPLING_station_number, SAMPLING_station_number=='3', '03')) %>%
    filter(SAMPLING_station_number %in% c('11', '13', '16', '20', '23', '27') ==FALSE) %>%
    ggplot(aes(x=depth, y=per_tot))+
    #geom_bar(aes( y=0.5),stat='identity', fill = "grey",alpha=0.8, width=20)+
    geom_bar(stat='identity', aes(fill = !!taxa_level), width=20)+
    coord_flip()+
    scale_x_reverse()+
    facet_wrap(~ SAMPLING_station_number)+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
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
  
  filename = paste(directory, marker,'_top10',taxa_level,'_meso_bar_bycast.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 8, width =10, units = 'in')
}





# SCRAP BELOW ---------------

# Plot Top Taxa Compositionally --------------------------------------

#get top 5 most abundant taxa for dataset
#want to limit to certain groups for each marker.
#markers <- c("12S","COI","18S","12SnoCL", "COInoH")
#could plot as percent of total or percent of group

if (marker == 'COI') {
  print(marker)
  phyla <- c('Arthropoda', 'Cnidaria', 'Bacillariophyta')
}
if (marker =='18S') {
  print(marker)
  #phyla <- c('Bacillariophyta', 'unknown')
  #phyla <- c('Bacillariophyta_Dinophyceae')
  #phyla <- c('Bacillariophyta')
  #phyla <- c('Dinophyceae')
  phyla <- c('Dinophyceae','Bacillariophyta','Arthropoda', 'Cnidaria')
}
if (marker =='12S') {
  print(marker)
  phyla <- c('Chordata')
}
if (marker =='12SnoCL') {
  print(marker)
  phyla <- c('Chordata')
}
if (marker == 'COInoH') {
  print(marker)
  phyla <- c('Arthropoda')
}

for (val in phyla) {
  var = sym(val)
  
  #Split by Diel group
  top_taxa <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    filter(Phylum == var) %>%  #limit to a certain group
    group_by(Species) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    arrange(-sum_per_tot) %>%
    distinct(Species,.keep_all = TRUE ) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    ungroup() %>%
    filter(Species != var) %>% #remove general phylum level annot
    filter(Species != 'Dinophyceae') %>%
    filter(Species != 'Dinophyceae sp. CCMP1878') %>% #remove unclassified dino
    filter(Species != 'Dinophyceae sp. UDMS0803') %>% #remove unclassified dino
    #filter(Genus != 'unknown') %>% #remove unclassified genera
    select(Species, sum_per_tot) %>%
    top_n(20)
  print(top_taxa)
  
  
  # assign text colour
  textcol <- "grey40"
  
  bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    right_join(top_taxa) %>% #limit to most abundant taxa
    unite(Label, Class, Order, Family, Genus, Species, sep="_", remove='False') %>%
    mutate(depth_char = as.character(depth)) %>%
    filter(diel %in% c('day', 'night', 'transition')) %>%
    #fct_reorder(name, desc(val))
    ggplot(aes(x = fct_reorder(depth_char, desc(depth)), y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = Label))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    facet_wrap(~diel, nrow =4) +
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    theme_minimal() +
    theme(
      #legend
      legend.position="bottom",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=5,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=5,face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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
  
  
  filename = paste(directory, marker,'_top20sp_bar_',var,'_comp_nightday.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  filename = paste(directory, marker,'_top20sp_bar_',var,'_comp_nightday.svg', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  
  # get compositional data in that group
  #get top 20 species in group
  top_taxa <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    filter(Phylum == var) %>%  #limit to a certain group
    group_by(Species) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    arrange(-sum_per_tot) %>%
    distinct(Species,.keep_all = TRUE ) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    ungroup() %>%
    filter(Species != var) %>% #remove general phylum level annot
    filter(Species != 'Dinophyceae') %>%
    filter(Species != 'Dinophyceae sp. CCMP1878') %>% #remove unclassified dino
    filter(Species != 'Dinophyceae sp. UDMS0803') %>% #remove unclassified dino
    #filter(Genus != 'unknown') %>% #remove unclassified genera
    select(Species, sum_per_tot) %>%
    top_n(10)
  print(top_taxa)
  
  
  # assign text colour
  textcol <- "grey40"
  
  bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    right_join(top_taxa) %>% #limit to most abundant taxa
    unite(Label, Class, Order, Family, Genus, Species, sep="_", remove='False') %>%
    #fct_reorder(name, desc(val))
    ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = Label))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    theme_minimal() +
    theme(
      #legend
      legend.position="bottom",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=5,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=5,face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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
  
  
  filename = paste(directory, marker,'_top10sp_bar_',var,'_comp.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =20, units = 'in')
  filename = paste(directory, marker,'_top20sp_bar_',var,'_comp.svg', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =20, units = 'in')
  
  #Split by Diel group
  top_taxa <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    filter(Phylum == var) %>%  #limit to a certain group
    group_by(Species) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    arrange(-sum_per_tot) %>%
    distinct(Species,.keep_all = TRUE ) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    ungroup() %>%
    filter(Species != var) %>% #remove general phylum level annot
    filter(Species != 'Dinophyceae') %>%
    filter(Species != 'Dinophyceae sp. CCMP1878') %>% #remove unclassified dino
    filter(Species != 'Dinophyceae sp. UDMS0803') %>% #remove unclassified dino
    #filter(Genus != 'unknown') %>% #remove unclassified genera
    select(Species, sum_per_tot) %>%
    top_n(20)
  print(top_taxa)
  
  
  # assign text colour
  textcol <- "grey40"
  
  bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    right_join(top_taxa) %>% #limit to most abundant taxa
    unite(Label, Class, Order, Family, Genus, Species, sep="_", remove='False') %>%
    mutate(depth_char = as.character(depth)) %>%
    filter(diel %in% c('7am-6pm', '8pm-5am', 'night', 'day')) %>%
    #fct_reorder(name, desc(val))
    ggplot(aes(x = fct_reorder(depth_char, desc(depth)), y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = Label))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    facet_wrap(~diel, nrow =4) +
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    theme_minimal() +
    theme(
      #legend
      legend.position="bottom",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=5,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=5,face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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
  
  
  filename = paste(directory, marker,'_top20sp_bar_',var,'_comp_nightday.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  filename = paste(directory, marker,'_top20sp_bar_',var,'_comp_nightday.svg', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  
  #Split by Diel group, just ESP samples
  top_taxa <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    filter(Sampling_method == 'ESP') %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    filter(Phylum == var) %>%  #limit to a certain group
    group_by(Species) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    arrange(-sum_per_tot) %>%
    distinct(Species,.keep_all = TRUE ) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    ungroup() %>%
    filter(Species != var) %>% #remove general phylum level annot
    filter(Species != 'Dinophyceae') %>%
    filter(Species != 'Dinophyceae sp. CCMP1878') %>% #remove unclassified dino
    filter(Species != 'Dinophyceae sp. UDMS0803') %>% #remove unclassified dino
    #filter(Genus != 'unknown') %>% #remove unclassified genera
    select(Species, sum_per_tot) %>%
    top_n(20)
  print(top_taxa)
  
  
  # assign text colour
  textcol <- "grey40"
  
  bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    filter(Sampling_method == 'ESP') %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    right_join(top_taxa) %>% #limit to most abundant taxa
    unite(Label, Class, Order, Family, Genus, Species, sep="_", remove='False') %>%
    mutate(depth_char = as.character(depth)) %>%
    filter(diel %in% c('7am-6pm', '8pm-5am', 'night', 'day')) %>%
    #fct_reorder(name, desc(val))
    ggplot(aes(x = fct_reorder(depth_char, desc(depth)), y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = Label))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    facet_wrap(~diel, nrow =4) +
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    theme_minimal() +
    theme(
      #legend
      legend.position="bottom",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=5,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=5,face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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
  
  
  filename = paste(directory, marker,'_top20sp_bar_',var,'_comp_nightday_ESP.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  
  #Split by Diel group, just ESP, raw read number
  top_taxa <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    filter(Sampling_method == 'ESP') %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    filter(Phylum == var) %>%  #limit to a certain group
    group_by(Species) %>%
    mutate(sum_reads = sum(reads)) %>%
    arrange(-sum_reads) %>%
    distinct(Species,.keep_all = TRUE ) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_reads) %>%
    ungroup() %>%
    filter(Species != var) %>% #remove general phylum level annot
    filter(Species != 'Dinophyceae') %>%
    filter(Species != 'Dinophyceae sp. CCMP1878') %>% #remove unclassified dino
    filter(Species != 'Dinophyceae sp. UDMS0803') %>% #remove unclassified dino
    #filter(Genus != 'unknown') %>% #remove unclassified genera
    select(Species, sum_reads) %>%
    top_n(20)
  print(top_taxa)
  
  
  # assign text colour
  textcol <- "grey40"
  
  bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    filter(Sampling_method == 'ESP') %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    right_join(top_taxa) %>% #limit to most abundant taxa
    unite(Label, Class, Order, Family, Genus, Species, sep="_", remove='False') %>%
    mutate(depth_char = as.character(depth)) %>%
    filter(diel %in% c('7am-6pm', '8pm-5am', 'night', 'day')) %>%
    #fct_reorder(name, desc(val))
    ggplot(aes(x = fct_reorder(depth_char, desc(depth)), y = reads)) +
    geom_bar(stat = "identity", aes(fill = Label))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    facet_wrap(~diel, nrow =4) +
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    theme_minimal() +
    theme(
      #legend
      legend.position="bottom",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=5,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=5,face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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
  
  
  filename = paste(directory, marker,'_top20sp_bar_',var,'_rawreads_nightday_ESP.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  
  #Split by Diel group, just ESP, raw read number, no Anchovy
  top_taxa <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    filter(Sampling_method == 'ESP') %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    filter(Phylum == var) %>%  #limit to a certain group
    group_by(Species) %>%
    mutate(sum_reads = sum(reads)) %>%
    arrange(-sum_reads) %>%
    distinct(Species,.keep_all = TRUE ) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_reads) %>%
    ungroup() %>%
    filter(Species != var) %>% #remove general phylum level annot
    filter(Species != 'Dinophyceae') %>%
    filter(Species != 'Dinophyceae sp. CCMP1878') %>% #remove unclassified dino
    filter(Species != 'Dinophyceae sp. UDMS0803') %>% #remove unclassified dino
    filter(Species != 'Engraulis mordax') %>% #remove anchovy here
    #filter(Genus != 'unknown') %>% #remove unclassified genera
    select(Species, sum_reads) %>%
    top_n(20)
  print(top_taxa)
  
  
  # assign text colour
  textcol <- "grey40"
  
  bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    filter(Sampling_method == 'ESP') %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    right_join(top_taxa) %>% #limit to most abundant taxa
    unite(Label, Class, Order, Family, Genus, Species, sep="_", remove='False') %>%
    mutate(depth_char = as.character(depth)) %>%
    filter(diel %in% c('7am-6pm', '8pm-5am', 'night', 'day')) %>%
    #fct_reorder(name, desc(val))
    ggplot(aes(x = fct_reorder(depth_char, desc(depth)), y = reads)) +
    geom_bar(stat = "identity", aes(fill = Label))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    facet_wrap(~diel, nrow =4) +
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    theme_minimal() +
    theme(
      #legend
      legend.position="bottom",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=5,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=5,face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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
  
  
  filename = paste(directory, marker,'_top20sp_bar_',var,'_rawreads_nightday_ESP_noAnch.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  
  #Split by Diel group, just ESP, raw read number, no Anchovy
  top_taxa <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    filter(Sampling_method == 'ESP') %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    filter(Phylum == var) %>%  #limit to a certain group
    group_by(Species) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    arrange(-sum_per_tot) %>%
    distinct(Species,.keep_all = TRUE ) %>%
    select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
    ungroup() %>%
    filter(Species != var) %>% #remove general phylum level annot
    filter(Species != 'Dinophyceae') %>%
    filter(Species != 'Dinophyceae sp. CCMP1878') %>% #remove unclassified dino
    filter(Species != 'Dinophyceae sp. UDMS0803') %>% #remove unclassified dino
    filter(Species != 'Engraulis mordax') %>% #remove anchovy here
    #filter(Genus != 'unknown') %>% #remove unclassified genera
    select(Species, sum_per_tot) %>%
    top_n(20)
  print(top_taxa)
  
  
  # assign text colour
  textcol <- "grey40"
  
  bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    filter(Sampling_method == 'ESP') %>%
    inner_join(species_label,  by = c("ASV")) %>%
    mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                              TRUE ~ as.character(Phylum))) %>%
    right_join(top_taxa) %>% #limit to most abundant taxa
    unite(Label, Class, Order, Family, Genus, Species, sep="_", remove='False') %>%
    mutate(depth_char = as.character(depth)) %>%
    filter(diel %in% c('7am-6pm', '8pm-5am', 'night', 'day')) %>%
    #fct_reorder(name, desc(val))
    ggplot(aes(x = fct_reorder(depth_char, desc(depth)), y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = Label))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    facet_wrap(~diel, nrow =4) +
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    theme_minimal() +
    theme(
      #legend
      legend.position="bottom",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=5,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=5,face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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
  
  
  filename = paste(directory, marker,'_top20sp_bar_',var,'_comp_nightday_ESP_noAnch.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  
  ## Number of Unique taxa by sample by Diel group - color by Order
  
  bp_top <- inner_join(potu.c, species_label,  by = c("ASV")) %>%
    group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species, SampleID) %>%
    mutate(sum_reads = sum(reads)) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    ungroup() %>%
    distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, SampleID, .keep_all=TRUE) %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species, SampleID, sum_reads,sum_per_tot ) %>%
    filter(sum_reads >= 5) %>% #just in case
    mutate(count = 1) %>%
    group_by(SampleID, Order) %>%
    mutate(count = sum(count)) %>%
    ungroup() %>%
    distinct(SampleID,Order, count, .keep_all=FALSE) %>%
    left_join(samp.c) %>%
    filter(Sampling_method=='ESP') %>%
    mutate(depth_char = as.character(depth)) %>%
    filter(diel %in% c('7am-6pm', '8pm-5am', 'night', 'day')) %>%
    ggplot(aes(x = fct_reorder(depth_char, desc(depth)), y = count)) +
    geom_bar(stat = "identity", aes(fill = Order))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    facet_wrap(~diel, nrow =4) +
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    theme_minimal() +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    theme(
      #legend
      legend.position="bottom",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=5,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=5,face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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
  filename = paste(directory, marker,'_NumTaxa_',var,'_byOrder.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  
  ## Number of Unique taxa by sample by Diel group - color by Order
  
  bp_top <- inner_join(potu.c, species_label,  by = c("ASV")) %>%
    mutate(count=1) %>%
    filter(reads>0) %>% #just in case
    filter(Species!='Engraulis mordax') %>%  #remove anchovy if want
    group_by(Order, SampleID) %>%
    mutate(count = sum(count)) %>%
    ungroup() %>%
    distinct(Order, SampleID,count, .keep_all=FALSE) %>%
    left_join(samp.c) %>%
    filter(Sampling_method=='ESP') %>%
    mutate(depth_char = as.character(depth)) %>%
    filter(diel %in% c('7am-6pm', '8pm-5am', 'night', 'day')) %>%
    ggplot(aes(x = fct_reorder(depth_char, desc(depth)), y = count)) +
    geom_bar(stat = "identity", aes(fill = Order))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    facet_wrap(~diel, nrow =4) +
    #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    theme_minimal() +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    theme(
      #legend
      legend.position="bottom",legend.direction="vertical",
      legend.text=element_text(colour=textcol,size=5,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour=textcol,size=5,face="bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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
  filename = paste(directory, marker,'_NumASVs_',var,'_byOrder_noanch.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  
  
  
}


# SCRAP BELOW ###########

#stat summary and geom smooth showing distribution of organisms:
#percent of total reads

#how does # of ASVs correspond with # of reads? Per_tot? w/depth and time
#merge ASVS with same taxonomy
test <- inner_join(potu.c, meta %>% select(SampleID, depth, local_time, time_label,diel, PlateID, ESP),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Stenobrachius') %>%
  #filter(Genus == 'Diaphus') %>%
  group_by(Genus, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Genus, SampleID, .keep_all=TRUE) %>%
  filter(ESP!='Bongo') 


test %>% 
  filter(diel !='transition') %>%
  #filter(reads>10) %>%
  ggplot(aes(x=depth, y=per_tot, color=diel, fill=diel))+
  geom_point()+
  stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20,  geom='bar', alpha=0.4)+
  stat_summary(geom = "errorbar")+
  #stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 40, aes(color=diel), geom='line')
  #stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20, aes(color=diel))#+
  #geom_smooth(aes(color=diel),method = 'loess') +
  geom_smooth()
#scale_y_reverse()+
#scale_x_discrete(position = "top") 
#scale_x_continuous(position="top")+
#labs(x='Percent Total Reads', y="Depth(m)")#+
#coord_flip()

# look at number of reads:
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Stenobrachius') %>%
  #filter(Genus == 'Diaphus') %>%
  group_by(Genus, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Genus, SampleID, .keep_all=TRUE) %>%
  filter(ESP!='Bongo') 


test %>% 
  filter(diel !='transition') %>%
  #filter(reads>10) %>%
  ggplot(aes(x=depth, y=reads, color=diel, fill=diel))+
  geom_point()+
  stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20,  geom='bar', alpha=0.4)+
  stat_summary(geom = "errorbar")+
  #stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 40, aes(color=diel), geom='line')
  #stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20, aes(color=diel))#+
  #geom_smooth(aes(color=diel),method = 'loess') +
  geom_smooth()
#scale_y_reverse()+
#scale_x_discrete(position = "top") 
#scale_x_continuous(position="top")+
#labs(x='Percent Total Reads', y="Depth(m)")#+
#coord_flip()

# go cast by cast and look at distribution through time.

test <- inner_join(potu.c, samp.c %>% select(SampleID, SAMPLING_station,SAMPLING_station_number,depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  #filter(Genus == 'Stenobrachius') %>%
  #filter(Genus == 'Diaphus') %>%
  filter(Genus == 'Engraulis') %>%
  group_by(Genus, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Genus, SampleID, .keep_all=TRUE) %>%
  filter(ESP!='Bongo') 


test %>% 
  filter(SAMPLING_station_number >=1) %>%
  arrange(SAMPLING_station_number) %>%
  mutate(SAMPLING_station_number = replace(SAMPLING_station_number, SAMPLING_station_number=='3', '03')) %>%
  unite(label, SAMPLING_station_number,local_time, SAMPLING_station,sep='_') %>%
  ggplot(aes(x=depth, y=per_tot, fill=diel, color=diel))+
  geom_bar(stat='identity')+
  geom_point()+
  coord_flip()+
  scale_x_reverse()+
  facet_wrap(.~label, nrow=2)

filename = paste(directory, marker,'_Engraulis_casts_pertot.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =20, units = 'in')


# By reads?
# go cast by cast and look at distribution through time.

test <- inner_join(potu.c, samp.c %>% select(SampleID, SAMPLING_station,SAMPLING_station_number,depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Stenobrachius') %>%
  #filter(Genus == 'Diaphus') %>%
  #filter(Genus == 'Engraulis') %>%
  filter(reads>1) %>%
  mutate(ASVs=1) %>%
  group_by(Genus, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  mutate(ASVs= sum(ASVs)) %>%
  ungroup() %>%
  distinct(Genus, SampleID, .keep_all=TRUE) %>%
  filter(ESP!='Bongo') 


test %>% 
  filter(SAMPLING_station_number >=1) %>%
  arrange(SAMPLING_station_number) %>%
  mutate(SAMPLING_station_number = replace(SAMPLING_station_number, SAMPLING_station_number=='3', '03')) %>%
  unite(label, SAMPLING_station_number,local_time, SAMPLING_station,sep='_') %>%
  ggplot(aes(x=depth, y=ASVs, fill=diel, color=diel))+
  geom_bar(stat='identity')+
  geom_point()+
  coord_flip()+
  scale_x_reverse()+
  facet_wrap(.~label, nrow=2)

filename = paste(directory, marker,'_Stenobrachius_casts_ASVs_2reads.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =20, units = 'in')

# Just MARS Casts, through time, showing depth 'center' of abundance?
test <- inner_join(potu.c, samp.c %>% select(SampleID, SAMPLING_station,SAMPLING_station_number,depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Stenobrachius') %>%
  #filter(Genus == 'Diaphus') %>%
  #filter(Genus == 'Engraulis') %>%
  #filter(reads>1) %>%
  #mutate(ASVs=1) %>%
  group_by(Genus, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  #mutate(ASVs= sum(ASVs)) %>%
  ungroup() %>%
  distinct(Genus, SampleID, .keep_all=TRUE) %>%
  filter(ESP!='Bongo') %>%
  filter(SAMPLING_station_number >=1) %>%
  filter(SAMPLING_station=='MARS') %>%
  arrange(SAMPLING_station_number) %>%
  mutate(SAMPLING_station_number = replace(SAMPLING_station_number, SAMPLING_station_number=='3', '03')) %>%
  mutate(local_time = mdy_hm(local_time))

center <- test %>%
  #group_by(SampleID) %>%
  mutate(value = depth*per_tot) %>%
  mutate(count=1) %>%
  #ungroup() %>%
  group_by(SAMPLING_station_number) %>%
  mutate(sum_value = sum(value)) %>%
  #mutate(count_value = sum(count)) %>%
  mutate(count_value = sum(per_tot)) %>%
  mutate(mean_depth2 = sum_value/count_value) %>%
  #mutate(mean_depth = mean(value)) %>%
  distinct(SAMPLING_station_number, .keep_all=TRUE) %>%
  select(mean_depth2, local_time)


test %>%
  left_join(center) %>%
  filter(SAMPLING_station_number %in% c('11', '13', '16', '20', '23', '27') ==FALSE) %>%
  #unite(label, SAMPLING_station_number,local_time, SAMPLING_station,sep='_') %>%
  ggplot(aes(y=depth, x=local_time))+
  #geom_bar(stat='identity')+
  geom_point(aes(size=per_tot, color=diel))+
  geom_point(aes(y=mean_depth2,size=1),color='black')+
  geom_line(aes(y=mean_depth2), color='black', size=1, linetype = "dashed") +
  #geom_line(data = center, aes(x=local_time, y=mean_depth2)) %>%
  #coord_flip()+
  scale_y_reverse()

filename = paste(directory, marker,'_Stenobrachius_casts_throughtime.png', sep='')
print(var)
filename
ggsave(filename,height = 5, width =10, units = 'in')

## Try Other species:
test <- inner_join(potu.c, samp.c %>% select(SampleID, SAMPLING_station,SAMPLING_station_number,depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  #filter(Genus == 'Stenobrachius') %>%
  #filter(Genus == 'Diaphus') %>%
  filter(Genus == 'Engraulis') %>%
  #filter(reads>1) %>%
  #mutate(ASVs=1) %>%
  group_by(Genus, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  #mutate(ASVs= sum(ASVs)) %>%
  ungroup() %>%
  distinct(Genus, SampleID, .keep_all=TRUE) %>%
  filter(ESP!='Bongo') %>%
  filter(SAMPLING_station_number >=1) %>%
  filter(SAMPLING_station=='MARS') %>%
  arrange(SAMPLING_station_number) %>%
  mutate(SAMPLING_station_number = replace(SAMPLING_station_number, SAMPLING_station_number=='3', '03')) %>%
  mutate(local_time = mdy_hm(local_time))

center <- test %>%
  #group_by(SampleID) %>%
  mutate(value = depth*per_tot) %>%
  mutate(count=1) %>%
  #ungroup() %>%
  group_by(SAMPLING_station_number) %>%
  mutate(sum_value = sum(value)) %>%
  #mutate(count_value = sum(count)) %>%
  mutate(count_value = sum(per_tot)) %>%
  mutate(mean_depth2 = sum_value/count_value) %>%
  #mutate(mean_depth = mean(value)) %>%
  distinct(SAMPLING_station_number, .keep_all=TRUE) %>%
  select(mean_depth2, local_time)


test %>%
  left_join(center) %>%
  filter(SAMPLING_station_number %in% c('11', '13', '16', '20', '23', '27') ==FALSE) %>%
  #unite(label, SAMPLING_station_number,local_time, SAMPLING_station,sep='_') %>%
  ggplot(aes(y=depth, x=local_time))+
  #geom_bar(stat='identity')+
  geom_point(aes(size=per_tot, color=diel))+
  geom_point(aes(y=mean_depth2,size=1),color='black')+
  geom_line(aes(y=mean_depth2), color='black', size=1, linetype = "dashed") +
  #geom_line(data = center, aes(x=local_time, y=mean_depth2)) %>%
  #coord_flip()+
  scale_y_reverse()

filename = paste(directory, marker,'_Engraulis_casts_throughtime.png', sep='')
print(var)
filename
ggsave(filename,height = 5, width =10, units = 'in')

## Plot ESP Samples


test <- inner_join(potu.c, samp.c %>% select(SampleID, SAMPLING_station,SAMPLING_station_number,depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Stenobrachius') %>%
  #filter(Genus == 'Diaphus') %>%
  #filter(Genus == 'Engraulis') %>%
  group_by(Genus, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Genus, SampleID, .keep_all=TRUE) %>%
  filter(ESP!='Bongo') %>%
  filter(ESP =='ESP') %>%
  #filter(SAMPLING_station_number >=1) %>%
  #filter(SAMPLING_station=='MARS') %>%
  #arrange(SAMPLING_station_number) %>%
  #mutate(SAMPLING_station_number = replace(SAMPLING_station_number, SAMPLING_station_number=='3', '03')) %>%
  mutate(local_time = mdy_hm(local_time))

# center <- test %>%
#   #group_by(SampleID) %>%
#   mutate(value = depth*per_tot) %>%
#   mutate(count=1) %>%
#   #ungroup() %>%
#   group_by(SAMPLING_station_number) %>%
#   mutate(sum_value = sum(value)) %>%
#   #mutate(count_value = sum(count)) %>%
#   mutate(count_value = sum(per_tot)) %>%
#   mutate(mean_depth2 = sum_value/count_value) %>%
#   #mutate(mean_depth = mean(value)) %>%
#   distinct(SAMPLING_station_number, .keep_all=TRUE) %>%
#   select(mean_depth2, local_time)


test %>%
  #left_join(center) %>%
  #filter(SAMPLING_station_number %in% c('11', '13', '16', '20', '23', '27') ==FALSE) %>%
  #unite(label, SAMPLING_station_number,local_time, SAMPLING_station,sep='_') %>%
  ggplot(aes(y=depth, x=local_time))+
  #geom_bar(stat='identity')+
  geom_point(aes(size=per_tot, color=diel))+
  #geom_point(aes(y=mean_depth2,size=1),color='black')+
  #geom_line(aes(y=mean_depth2), color='black', size=1, linetype = "dashed") +
  #geom_line(data = center, aes(x=local_time, y=mean_depth2)) %>%
  #coord_flip()+
  scale_y_reverse()

filename = paste(directory, marker,'_Engraulis_ESPcasts_throughtime.png', sep='')
print(var)
filename
ggsave(filename,height = 5, width =10, units = 'in')






### Working on BELOW ####

facet_grid(.~SAMPLING_station_number)
#facet_grid(local_time ~.)



#facet_grid(. ~SAMPLING_station_number)
#facet_grid(SAMPLING_station_number ~ .)

#stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20,  geom='bar', alpha=0.4)+
#stat_summary(geom = "errorbar")+
#stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 40, aes(color=diel), geom='line')
#stat_summary_bin(fun.data = "mean_cl_boot", binwidth = 20, aes(color=diel))#+
#geom_smooth(aes(color=diel),method = 'loess') +
#geom_smooth()





test %>% ggplot(aes(y=depth, x=per_tot))+
  geom_point()



filter(reads>1) %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(sum_ASVs = sum(count)) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID,.keep_all=TRUE)







### Species profiles with depth ------

# Look at top 10 most abundant species and their reads by depth value
# During Day and During Night
# Engraulis mordax               10743.
# 2 Stenobrachius leucopsarus       2340.
# 3 Merluccius productus            1907.
# 4 Leuroglossus schmidti           1415.
# 5 Diaphus theta                   1387.
# 6 Sebastes                        1165.
# 7 Lipolagus ochotensis             400.
# 8 Delphinidae                      207.
# 9 Cyclothone acclinidens           181.
# 10 Hydrolagus colliei               179.

library(stringr)
#how does # of ASVs correspond with # of reads? Per_tot? w/depth and time
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Engraulis') %>%
  filter(reads>1) %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(sum_ASVs = sum(count)) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID,.keep_all=TRUE)

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=diel))
p

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=depth))
p
library(viridis)
p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=time, size=depth, shape=diel))+
  scale_color_viridis()
p

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_per_tot))+ 
  geom_point(aes(color=time, size=depth, shape=ESP))+
  scale_color_viridis()
p
p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=ESP, size=depth, shape=diel))
p

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_per_tot))+ 
  geom_point(aes(color=ESP, size=depth, shape=diel))
p

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=PlateID, size=depth, shape=ESP))
p

#pull out top 10
top_taxa <- potu.c %>%
  full_join(species_label) %>%
  filter(Species !='unassigned') %>%
  filter(!!taxa_level !='s_') %>%
  filter(!!taxa_level !='g_') %>%
  group_by(Species) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  distinct(Species,.keep_all = TRUE ) %>%
  arrange(-sum_per_tot) %>%
  select(Kingdom, Phylum, Class, Order, Family,Genus, Species, sum_per_tot) %>%
  ungroup() %>%
  select(Species, sum_per_tot) %>%
  top_n(10)
top_taxa


# plot read number by depth:
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel),  by = c("SampleID")) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Species == 'Engraulis mordax') %>%
  ggplot(aes(y=depth, x=reads, color=diel))+
  geom_point(alpha=0.8)
test

test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel),  by = c("SampleID")) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Species == 'Engraulis mordax') %>%
  ggplot(aes(y=depth, x=local_time, color=diel))+
  geom_point(alpha=0.8)
test


test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Species == 'Engraulis mordax') %>%
  ggplot(aes(y=depth, x=time, fill=diel))+
  geom_point(aes(color=diel),alpha=0.8)+
  scale_y_reverse()
#geom_contour_filled(aes(z=reads))
test

#binned density plot?
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Engraulis') %>%
  filter(reads>0)
#number of ASVs within each sample that are Engraulis
p <- ggplot(test, aes(x = time, y = depth, z=reads)) + stat_binhex()
p

#how does # of ASVs correspond with # of reads? Per_tot?
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Engraulis') %>%
  filter(reads>0) %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(sum_count = sum(count)) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID,.keep_all=TRUE)

p <- test %>% ggplot(aes(x = sum_count, y = sum_reads))+ 
  geom_point(aes(color=diel))
p

p <- test %>% ggplot(aes(x = sum_count, y = sum_reads))+ 
  geom_point(aes(color=depth))
p


#SCRAP
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Species == 'Engraulis mordax') %>%
  ggplot(aes(y=depth, x=time, z=reads))+
  stat_summary_hex(function(z){log(sum(z))})

geom_point(aes(fill=reads),alpha=0.6, shape=21, size=5)+
  #scale_fill_gradient(name = "reads", trans = "log")+
  scale_colour_gradient(name = "reads", trans = "log")
#geom_contour_filled(aes(z=reads),bins=5)
test


bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
                            TRUE ~ as.character(Phylum))) %>%
  right_join(top_taxa) %>% #limit to most abundant taxa
  unite(Label, Class, Order, Family, Genus, Species, sep="_", remove='False') %>%
  mutate(depth_char = as.character(depth)) %>%
  filter(diel %in% c('day', 'night', 'transition')) %>%
  #fct_reorder(name, desc(val))
  ggplot(aes(x = fct_reorder(depth_char, desc(depth)), y = per_tot)) +
  geom_bar(stat = "identity", aes(fill = Label))+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  labs(x="",y=paste("Percent ",var," Reads", sep=""))+
  facet_wrap(~diel, nrow =4) +
  #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
  theme_minimal() +
  theme(
    #legend
    legend.position="bottom",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=5,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=5,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
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


