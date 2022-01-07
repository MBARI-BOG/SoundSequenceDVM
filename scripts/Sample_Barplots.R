#040221
#kpitz

# Set directory to save plots
directory <- 'figures/Limited_barplots/'


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




# Create Seasonal Variables -----------------------------------------------

meta <- samp.c %>% 
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  mutate(month =  month(time)) %>%
  mutate(day =  day(time)) %>%
  mutate(year =  year(time)) %>%
  mutate(jday = yday(time)) %>%
  mutate(month_char = as.character(month)) %>%
  mutate(year_char = as.character(year)) %>%
  mutate(depth_bin = case_when(depth <=50 ~ "0-50m",
                               depth >50 & depth <=100 ~ "50-100m",
                               depth >100 & depth <=200 ~ "100-200m",
                               depth >200 & depth <=300 ~ "200-300m",
                               depth >300 & depth <=400 ~ "300-400m",
                               depth >400 & depth <=500 ~ "400-500m",
                               depth >400 & depth <=600 ~ "500-600m",
                               depth >600 & depth <=750 ~ "600-750m", TRUE ~ "unknown"
  )) 

# Bar by Taxa_level ------------------------------------------------------------
taxas <- c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
#### Percent Reads
for (val in taxas) {
  taxa_level = sym(val)
  top_taxa <- potu.c %>%
    full_join(species_label) %>%
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
  
  filename = paste(directory, marker,'_top10',taxa_level,'_bar.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 6, width =20, units = 'in')
}
####  Raw Reads
for (val in taxas) {
  taxa_level = sym(val)
  top_taxa <- potu.c %>%
    full_join(species_label) %>%
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
  
  filename = paste(directory, marker,'_top10',taxa_level,'_bar_rawreads.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 6, width =10, units = 'in')
}
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
  ggsave(filename,height = 8, width =10, units = 'in')
  filename = paste(directory, marker,'_top20sp_bar_',var,'_comp.svg', sep='')
  print(var)
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  
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

#stat summary and geom smooth showing distribution of organisms:
#percent of total reads

library(stringr)
#how does # of ASVs correspond with # of reads? Per_tot? w/depth and time
#merge ASVS with same taxonomy
test <- inner_join(potu.c, samp.c %>% select(SampleID, depth, local_time, time_label,diel, PlateID),  by = c("SampleID")) %>%
  mutate(time = mdy_hm(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         str_detect(SampleID, 'Bongo')==TRUE ~'Bongo',
                         TRUE~'CTD')) %>%
  inner_join(species_label,  by = c("ASV")) %>%
  #filter(Genus == 'Stenobrachius') %>%
  filter(Genus == 'Diaphus') %>%
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


