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
    ungroup() %>%
    filter(reads>0)
  return(new_potu)
}

#can't return multiple objects, have to do it as a list
limit_by_meta <- function(meta, taxa, potu) {
  new_potu <- left_join(meta %>% select(SampleID), potu) %>%
    filter(reads>0)
  new_taxa <- left_join(new_potu %>% select(ASV), taxa)
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

out <- limit_by_meta(meta_filt, taxa_new, potu_filt) 
potu_filt <- out[[1]]
taxa_filt <- out[[2]]


# Plot --------------------------------------------------------------------

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


