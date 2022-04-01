#Species Diversity and relative abundance with depth
#kpitz
#120821


# Set directory to save plots
directory <- 'figures/Depth_sp_patterns/'


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
library(viridis)
library(stringr)

#markers <- c("12S")

marker = sym("12S")
  
# Import Data -------------------------------------------------------------
data_directory = "Data/filtered_seq_data/"

#ASV table
print('ASV table')
file = paste("CN19S_",marker,"_otu_Filtered.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
otu.c <- read_csv(filepath) %>% rename('ASV' = 1)

#taxa table
print('taxa table')
file = paste("CN19S_",marker,"_taxa_Filtered.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
tax.c <- read_csv(filepath) %>% rename('ASV' = 1)

#metadata table
print('metadata table')
#file = "A2W_12S_meta_Filtered.csv"
file = paste("CN19S_",marker,"_meta_Filtered.csv", sep='')
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



# Limit by Taxa and Depth -------------------------------------------------

limit_by_taxonomy <- function(taxa, potu) {
  new_potu <- left_join(taxa, potu) %>%
    select(ASV, SampleID, reads) %>%
    group_by(SampleID) %>%
    mutate(per_tot = reads / sum(reads) *100) %>%
    ungroup()
  return(new_potu)
}

#can't return multiple objects, have to do it as a list
limit_by_meta <- function(meta, taxa, potu) {
  new_potu <- left_join(meta %>% select(SampleID), potu) %>%
    #remove ASVs with 0 reads across entire dataset:
    group_by(ASV) %>%
    mutate(total_reads = sum(reads))%>%
    ungroup() %>%
    filter(total_reads>0)%>%
    select(-total_reads)
  new_taxa <- left_join(new_potu %>% select(ASV), taxa) %>%
    distinct(ASV, .keep_all = TRUE)
  out <- list(new_potu, new_taxa)
  return(out)
}

#limit by taxonomy
taxa_filt <- tax.c %>%
  filter(Family!='Engraulidae') %>%
  filter(Family!='Merlucciidae')

potu_filt <- limit_by_taxonomy(taxa_filt, potu.c)

meta_filt <- samp.c %>%
  filter(depth <600)

out <- limit_by_meta(meta_filt, taxa_filt, potu_filt) 
potu_filt <- out[[1]]
taxa_filt <- out[[2]]

# Need to limit samples with <500 reads

samp_500 <- potu_filt %>%
  group_by(SampleID) %>%
  mutate(total_reads = sum(reads))%>%
  ungroup() %>%
  distinct(SampleID,total_reads, .keep_all=FALSE) %>%
  filter(total_reads >=100)

out <- limit_by_meta(samp_500, taxa_filt, potu_filt) 
potu_filt <- out[[1]]
taxa_filt <- out[[2]]
meta_filt <- left_join(samp_500, meta_filt)

# Plot --------------------------------------------------------------------

#look at number of reads without hake and anchovy:
potu_filt %>%
  group_by(SampleID) %>%
  mutate(total_reads = sum(reads))%>%
  ungroup() %>%
  distinct(SampleID,total_reads) %>%
  arrange(total_reads) %>%
  glimpse()

# Just MARS Casts, through time, showing depth 'center' of abundance?
# Need to show all samples taken - zero is meaningful
# Plus column for amount of specific taxa

taxas = c('Stenobrachius','Lipolagus', 'Diaphus',
          'Sebastes', 'Leuroglossus', 'Nannobrachium', 'Sardinops', 'Scomber')
for (val in taxas) {
  taxa_level = sym(val)
  test <- left_join(potu_filt, meta_filt %>% 
                      select(SampleID, SAMPLING_station,SAMPLING_station_number,
                             depth, local_time, time_label,diel, PlateID),  
                    by = c("SampleID")) %>%
    mutate(time = mdy_hm(local_time)) %>%
    mutate(time_since = as.numeric(time)) %>%
    mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                           str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                           TRUE~'CTD')) %>%
    left_join(species_label,  by = c("ASV")) %>%
    #filter(Genus == 'Stenobrachius') %>%
    #filter(Genus == 'Lipolagus') %>%
    filter(Genus == taxa_level) %>%
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
    #geom_point(aes(y=mean_depth2,size=1),color='black')+
    geom_line(aes(y=mean_depth2), color='black', size=1, linetype = "dashed") +
    #geom_line(data = center, aes(x=local_time, y=mean_depth2)) %>%
    #coord_flip()+
    scale_color_manual(values = c("darkgoldenrod1", "deepskyblue", "chartreuse2"))+
    scale_y_reverse()
  
  filename = paste(directory, marker,'_',taxa_level, '_casts_throughtime_nohake_noanch.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 5, width =10, units = 'in')
  
  #Try to plot bars instead of points through time.
  #By cast?
  
  test %>%
    left_join(center) %>%
    filter(SAMPLING_station_number %in% c('11', '13', '16', '20', '23', '27') ==FALSE) %>%
    #unite(label, SAMPLING_station_number,local_time, SAMPLING_station,sep='_') %>%
    ggplot(aes(x=depth, y=per_tot))+
    geom_bar(aes( y=100),stat='identity', fill = "grey",alpha=0.8, width=10)+
    geom_bar(stat='identity', aes(fill = diel), width=10)+
    coord_flip()+
    scale_x_reverse()+
    facet_wrap(~ SAMPLING_station_number)+
    labs(x='depth (m)', y='Percent Total Reads', title=taxa_level) +
    #day night transition
    scale_fill_manual(values = c("darkgoldenrod1", "deepskyblue", "chartreuse2"))
  
  filename = paste(directory, marker,'_',taxa_level,'_bycast_nohake_noanch.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
}


# Stacked Bar by Cast -----------------------------------------------------


taxas <- c('Class','Order', 'Family', 'Genus', 'Species')
#### Percent Reads
for (val in taxas) {
  taxa_level = sym(val)
  top_taxa <- left_join(potu_filt, meta_filt,  by = c("SampleID")) %>% #join with metadata
    left_join(species_label,  by = c("ASV")) %>%
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
  bp_top <- left_join(potu_filt, meta_filt,  by = c("SampleID")) %>% #join with metadata
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
  
  filename = paste(directory, marker,'_top10',taxa_level,'_bar_bycast.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 8, width =10, units = 'in')
}


# Total Number of Reads by Cast -------------------------------------------


# Number of total reads per sample by cast:
textcol <- "grey40"
bp_top <- left_join(potu_filt, meta_filt,  by = c("SampleID")) %>% #join with metadata
  filter(SAMPLING_station_number >=1) %>%
  filter(SAMPLING_station=='MARS') %>%
  mutate(SAMPLING_station_number = replace(SAMPLING_station_number, SAMPLING_station_number=='3', '03')) %>%
  filter(SAMPLING_station_number %in% c('11', '13', '16', '20', '23', '27') ==FALSE) %>%
  distinct(SampleID, .keep_all=TRUE) %>%
  ggplot(aes(x=depth, y=total_reads))+
  geom_bar(aes( y=30000),stat='identity', fill = "grey",alpha=0.8, width=20)+
  geom_bar(stat='identity', fill ='black', width=20)+
  coord_flip()+
  scale_x_reverse()+
  facet_wrap(~ SAMPLING_station_number)+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  labs(x="",y="Total Reads without Anchovy or Hake")+
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
filename = paste(directory, marker,'_reads_bycast.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 8, width =10, units = 'in')



# Center of Relative Abundance ---------------------------------------------------

#above 200m and >= 200m?
taxas = c('Stenobrachius','Lipolagus', 'Diaphus',
          'Sebastes', 'Leuroglossus', 'Nannobrachium', 'Sardinops', 'Scomber')
#taxas = c('Stenobrachius')

for (val in taxas) {
  taxa_level = sym(val)
  test <- left_join(potu_filt, meta_filt %>% 
                      select(SampleID, SAMPLING_station,SAMPLING_station_number,
                             depth, local_time, time_label,diel, PlateID),  
                    by = c("SampleID")) %>%
    mutate(time = mdy_hm(local_time)) %>%
    mutate(time_since = as.numeric(time)) %>%
    mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                           str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                           TRUE~'CTD')) %>%
    left_join(species_label,  by = c("ASV")) %>%
    filter(Genus == taxa_level) %>%
    group_by(Genus, SampleID) %>%
    mutate(per_tot = sum(per_tot)) %>%
    mutate(reads = sum(reads)) %>%
    ungroup() %>%
    distinct(Genus, SampleID, .keep_all=TRUE) %>%
    filter(ESP!='Bongo') %>%
    filter(SAMPLING_station_number >=1) %>%
    filter(SAMPLING_station=='MARS') %>%
    arrange(SAMPLING_station_number) %>%
    mutate(SAMPLING_station_number = replace(SAMPLING_station_number, SAMPLING_station_number=='3', '03')) %>%
    mutate(local_time = mdy_hm(local_time))
  
  
  test %>%
    #left_join(center) %>%
    filter(SAMPLING_station_number %in% c('11', '13', '16', '20', '23', '27') ==FALSE) %>%
    #unite(label, SAMPLING_station_number,local_time, SAMPLING_station,sep='_') %>%
    mutate(depth_bin = case_when(depth >= 200 ~ '>=200m', TRUE ~ '<200m')) %>%
    select(SampleID, per_tot, depth, depth_bin, SAMPLING_station_number, diel) %>%
    ggplot(aes(x=depth, y=per_tot))+
    geom_boxplot(aes(group=depth_bin, fill=diel))+
    #geom_bar(aes( y=100),stat='identity', fill = "grey",alpha=0.8, width=10)+
    #geom_bar(stat='identity', aes(fill = diel), width=10)+
    geom_point(aes(color=diel))+
    #geom_smooth(span=1)+
    #geom_line()+
    coord_flip()+
    scale_x_reverse()+
    facet_wrap(~ SAMPLING_station_number)+
    labs(x='depth (m)', y='Percent Total Reads', title=taxa_level) +
    #day night transition
    scale_fill_manual(values = c("darkgoldenrod1", "deepskyblue", "chartreuse2"))+
    scale_color_manual(values = c("darkgoldenrod1", "deepskyblue", "chartreuse2"))+
    ylim(0,100)
  
  filename = paste(directory, marker,'_',taxa_level,'_bycast_nohake_noanch_boxplot.png', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
}







test <- inner_join(potu_filt, meta_filt %>% select(SampleID, SAMPLING_station,SAMPLING_station_number,depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  #filter(Genus == 'Stenobrachius') %>%
  filter(Genus == 'Lipolagus') %>%
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

# All Samples -------------------------------------------------------------



# All samples

test <- inner_join(potu_filt, meta_filt %>% select(SampleID, SAMPLING_station,SAMPLING_station_number,depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  #filter(Genus == 'Stenobrachius') %>%
  filter(Genus == 'Leuroglossus') %>%
  group_by(Genus, SampleID) %>%
  mutate(per_tot = sum(per_tot)) %>%
  mutate(reads = sum(reads)) %>%
  ungroup() %>%
  distinct(Genus, SampleID, .keep_all=TRUE) %>%
  filter(ESP!='Bongo') %>%
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
  #filter(ESP=='ESP') %>%
  filter(SAMPLING_station_number %in% c('11', '13', '16', '20', '23', '27') ==FALSE) %>%
  #unite(label, SAMPLING_station_number,local_time, SAMPLING_station,sep='_') %>%
  ggplot(aes(y=depth, x=local_time))+
  #geom_bar(stat='identity')+
  geom_point(aes(size=reads, color=ESP), alpha=0.8)+
  #geom_point(aes(y=mean_depth2,size=1),color='black')+
  #geom_line(aes(y=mean_depth2), color='black', size=1, linetype = "dashed") +
  #geom_line(data = center, aes(x=local_time, y=mean_depth2)) %>%
  #coord_flip()+
  scale_y_reverse() +
  scale_color_manual(values = c("darkgoldenrod1", "deepskyblue", "chartreuse2"))

filename = paste(directory, marker,'_Leuroglossus_casts_throughtime_nohake_noanch_reads.png', sep='')
print(var)
filename
ggsave(filename,height = 5, width =10, units = 'in')




















# Stenobrachius
test <- inner_join(potu_filt, 
                   meta_filt %>% select(SampleID, depth, local_time, time_label,diel, PlateID),  
                   by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(hour=hour(time)) %>%
  mutate(hour_char=as.character(hour)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Stenobrachius') %>%
  filter(reads>1) %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(sum_ASVs = sum(count)) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID,.keep_all=TRUE)

p <- test %>% ggplot(aes(x = time, y = depth))+ 
  geom_point(aes(color=fct_reorder(hour_char, hour), size=per_tot, shape=ESP))+
  scale_y_reverse()+
  scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
p

p <- test %>% ggplot(aes(x = depth, y = sum_per_tot))+ 
  geom_point() %>%
  facet_wrap(~ diel)
#geom_point(aes(color=fct_reorder(hour_char, hour), size=per_tot, shape=ESP))+
facet_wrap(~diel) %>%
  #scale_y_reverse()+
  scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
p


p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_per_tot, shape=ESP, color=PlateID))+ 
  geom_point() +
  facet_wrap(~diel, ncol=1, nrow=3)
p









# Plot Unfiltered ---------------------------------------------------------

# Stenobrachius
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(hour=hour(time)) %>%
  mutate(hour_char=as.character(hour)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Stenobrachius') %>%
  filter(reads>1) %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(sum_ASVs = sum(count)) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID,.keep_all=TRUE)


p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_per_tot, shape=ESP, color=PlateID))+ 
  geom_point() +
  facet_wrap(~diel, ncol=1, nrow=3)
p

p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_reads, shape=ESP, color=PlateID))+ 
  geom_point() +
  facet_wrap(~diel, ncol=1, nrow=3)
p



p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=diel, size=depth, shape=ESP))
p
filename = paste(directory, marker,'_NumASVs_NumReads_byDiel_Stenobrachius.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=PlateID, size=depth, shape=ESP))
p
filename = paste(directory, marker,'_NumASVs_NumReads_byPlateID_Stenobrachius.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')


p <- test %>% ggplot(aes(x = time, y = depth))+ 
  geom_point(aes(color=diel, size=sum_reads, shape=ESP))+
  scale_y_reverse()
  #geom_point(aes(color=diel, size=sum_reads, shape=diel))
p
filename = paste(directory, marker,'_NumReads_throughtime_Stenobrachius.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')

p <- test %>% ggplot(aes(x = time, y = depth))+ 
  geom_point(aes(color=diel, size=sum_per_tot, shape=ESP))+
  scale_y_reverse()
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
p
filename = paste(directory, marker,'_pertot_throughtime_Stenobrachius.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')

p <- test %>% ggplot(aes(x = time, y = depth))+ 
  geom_point(aes(color=fct_reorder(hour_char, hour), size=sum_ASVs, shape=ESP))+
  scale_y_reverse()+
  scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
p
filename = paste(directory, marker,'_sumASVs_throughtime_Stenobrachius.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')


p <- test %>% ggplot(aes(x = time, y = depth))+ 
  geom_point(aes(color=fct_reorder(hour_char, hour), size=sum_reads, shape=ESP))+
  scale_y_reverse()+
  scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
p
filename = paste(directory, marker,'_sumreads_throughtime_Stenobrachius.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')

p <- test %>% ggplot(aes(x = time, y = depth))+ 
  geom_point(aes(color=fct_reorder(hour_char, hour), size=sum_per_tot, shape=ESP))+
  scale_y_reverse()+
  scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
p
filename = paste(directory, marker,'_pertot_throughtime_Stenobrachius.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')


# take max of each plate and then proportion of that max?
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(hour=hour(time)) %>%
  mutate(hour_char=as.character(hour)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Stenobrachius') %>%
  filter(reads>1) %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(sum_ASVs = sum(count)) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  group_by(PlateID) %>%
  mutate(max_read_plate = max(sum_reads)) %>%
  mutate(prop_read_plate = sum_reads/max_read_plate*100) %>%
  ungroup() %>%
  distinct(SampleID,.keep_all=TRUE)

p <- test %>% ggplot(aes(x = time, y = depth))+ 
  geom_point(aes(color=fct_reorder(hour_char, hour), size=prop_read_plate, shape=ESP))+
  scale_y_reverse()+
  scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
p
filename = paste(directory, marker,'_pertot_throughtime_Stenobrachius_byplate.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')


#look at distribution of proportion of reads by cast?
#histogram at different time periods?
p <-test %>% 
  filter(depth<50) %>%
  ggplot(aes(x = prop_read_plate, fill = diel))+ 
  geom_histogram(position='dodge')
p

p <-test %>% 
  filter((depth>50)&(depth<=100)) %>%
  ggplot(aes(x = prop_read_plate, fill = diel))+ 
  geom_histogram(position='dodge')
p


p <-test %>% 
  filter((depth>200)&(depth<=500)) %>%
  ggplot(aes(x = prop_read_plate, fill = diel))+ 
  geom_histogram(position='dodge')
p


# Lipolagus
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(hour=hour(time)) %>%
  mutate(hour_char=as.character(hour)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Lipolagus') %>%
  filter(reads>1) %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(sum_ASVs = sum(count)) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID,.keep_all=TRUE)


p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_per_tot, shape=ESP, color=PlateID))+ 
  geom_point() +
  facet_wrap(~diel, ncol=1, nrow=3)
p

p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_reads, shape=ESP, color=PlateID))+ 
  geom_point() +
  facet_wrap(~diel, ncol=1, nrow=3)
p


p <- test %>% ggplot(aes(x = time, y = depth))+ 
  geom_point(aes(color=fct_reorder(hour_char, hour), size=sum_reads, shape=ESP))+
  scale_y_reverse()+
  scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
p

# Merluccius
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(hour=hour(time)) %>%
  mutate(hour_char=as.character(hour)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Genus == 'Merluccius') %>%
  filter(reads>1) %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(sum_ASVs = sum(count)) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID,.keep_all=TRUE)

p <- test %>% ggplot(aes(x = time, y = depth))+ 
  geom_point(aes(color=fct_reorder(hour_char, hour), size=per_tot, shape=ESP))+
  scale_y_reverse()+
  scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
p

p <- test %>% ggplot(aes(x = depth, y = sum_per_tot))+ 
  geom_point() %>%
  facet_wrap(~ diel)
  #geom_point(aes(color=fct_reorder(hour_char, hour), size=per_tot, shape=ESP))+
  facet_wrap(~diel) %>%
  #scale_y_reverse()+
  scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
#geom_point(aes(color=diel, size=sum_reads, shape=diel))
  
p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_per_tot, shape=ESP, color=PlateID))+ 
    geom_point() +
    facet_wrap(~diel, ncol=1, nrow=3)
p


# p <- test %>% ggplot(aes(x = depth, y = prop_read_plate))+ 
#   #geom_point(aes(color=fct_reorder(hour_char, hour), size=prop_read_plate, shape=ESP))+
#   geom_point(aes(color=diel, size=prop_read_plate, shape=ESP))+
#   stat_summary(aes(y = prop_read_plate,group=1), fun.y=mean, colour="red", geom="line",group=1)+
#   #geom_line(aes(color=diel)) +
#   scale_y_reverse()+
#   scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)
# p





# Myctophid Reads and ASV Diversity

test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  filter(Family == 'Myctophidae') %>%
  filter(reads>1) %>%
  mutate(count=1) %>%
  group_by(SampleID) %>%
  mutate(sum_ASVs = sum(count)) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID,.keep_all=TRUE)


p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_per_tot, shape=ESP, color=PlateID))+ 
  geom_point() +
  facet_wrap(~diel, ncol=1, nrow=4)
p

p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_reads, shape=ESP, color=PlateID))+ 
  geom_point() +
  facet_wrap(~diel, ncol=1, nrow=4)
p

#number of samples at each depth
p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, fill=diel))+ 
  geom_histogram() +
  facet_wrap(~diel, ncol=1, nrow=4)
p

p <- test %>%
  mutate(depth_bin = case_when(depth <=50 ~ "0-50m",
                               (depth >50 & depth <=100) ~ "50-100m",
                               (depth >100 & depth <=300) ~ "50-100m",
                               (depth >300 & depth <=400) ~ "50-100m",
                               (depth >400 & depth <=500) ~ "50-100m",
                               (depth >500 & depth <=750) ~ "50-100m", TRUE ~ "unknown"
                               )) %>%
  filter(depth>=0) %>%
  ggplot(aes(x = depth_bin, y = sum_reads)) +
  geom_boxplot(aes(fill=diel)) +
  geom_point()
p


p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=diel))
p

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=diel, size=depth, shape=ESP))
p
filename = paste(directory, marker,'_NumASVs_NumReads_byDiel_Myctophidae.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=PlateID, size=depth, shape=ESP))
p
filename = paste(directory, marker,'_NumASVs_NumReads_byPlateID_Myctophidae.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_per_tot))+ 
  geom_point(aes(color=PlateID, size=depth, shape=ESP))
p
filename = paste(directory, marker,'_NumASVs_PerReads_byPlateID_Myctophidae.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')


# Anchovy Reads and ASV Diversity -------
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

p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_per_tot, shape=ESP, color=PlateID))+ 
  geom_point() +
  facet_wrap(~diel, ncol=1, nrow=4)
p

p <- test %>%
  filter(depth >=0 ) %>%
  ggplot(aes(x = depth, y = sum_reads, shape=ESP, color=PlateID))+ 
  geom_point() +
  facet_wrap(~diel, ncol=1, nrow=4)
p


p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=diel))
p

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_reads))+ 
  geom_point(aes(color=depth))
p

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
filename = paste(directory, marker,'_NumASVs_NumReads_byPlateID_Engraulis.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')

p <- test %>% ggplot(aes(x = sum_ASVs, y = sum_per_tot))+ 
  geom_point(aes(color=PlateID, size=depth, shape=ESP))
p
filename = paste(directory, marker,'_NumASVs_PerReads_byPlateID_Engraulis.png', sep='')
print(var)
filename
ggsave(filename,height = 8, width =10, units = 'in')


