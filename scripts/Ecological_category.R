#081122
#kpitz

# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- 'figures/Ecological_category/'
# Set directory to retrieve data
data_directory = "Data/filtered_seq_data/"

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
library(viridis)

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


# Make tibble of lowest taxonomic annotation level for ASVs
# Group together different unassigned categories as 'Unknown'
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


# Modify metadata -----------------------------------------------

# create time variables, manual depth bins, sample type variable
meta <- samp.c %>% 
  mutate(time = ymd_hms(local_time)) %>%
  mutate(time_since = as.numeric(time)) %>%
  mutate(ESP = case_when(str_detect(SampleID, 'SC')==TRUE ~'ESP',
                         TRUE~'CTD')) %>%
  mutate(month =  month(time)) %>%
  mutate(hour = hour(time)) %>%
  mutate(day =  day(time)) %>%
  mutate(year =  year(time)) %>%
  mutate(jday = yday(time)) %>%
  mutate(month_char = as.character(month)) %>%
  mutate(year_char = as.character(year)) %>%
  mutate(depth_bin = case_when(depth <=25 ~ "0-25m",
                               depth >25 & depth <=75 ~ "25-75m",
                               #depth >50 & depth <=75 ~ "50-75m",
                               depth >75 & depth <=100 ~ "75-100m",
                               depth >100 & depth <=150 ~ "100-150m",
                               depth >150 & depth <=200 ~ "150-200m",
                               depth >200 & depth <=250 ~ "200-250m",
                               depth >250 & depth <=300 ~ "250-300m",
                               depth >300 & depth <=400 ~ "300-400m",
                               depth >400 & depth <=500 ~ "400-500m",
                               depth >400 & depth <=600 ~ "500-600m",
                               depth >600 & depth <=750 ~ "600-750m", TRUE ~ "unknown"
  )) 

# Species composition within Ecological Category ------------------

#cats <- c("epipelagic","mesopelagic", "bathypelagic","cosmopolitan")
cats <- c("benthopelagic or nearshore bottom", "epipelagic or benthic coastal", "meso and bathypelagic", "epipelagic","mesopelagic", "bathypelagic","cosmopolitan", "Atlantic species",
          "benthopelagic", "errant classification", "Exclude", "nearshore bottom", 
          "nearshore surface")
tax_levels = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

#organize by overall depth distribution across data set?
#merge by taxonomy and take average %reads in depth bin
#so for all the reads in a depth bin, what % were this thing?
#then also do mean read % by depth

#Kingdom,Phylum,Class,Order,Family,Genus,Species

for (val in cats) {
  var = sym(val)
  merged_data <- tax.c %>% full_join(sp_desig, by=tax_levels) %>%
    left_join(potu.c) %>%
    group_by(Kingdom,Phylum,Class,Order,Family,Genus,Species, SampleID) %>%
    mutate(sum_reads = sum(reads)) %>%
    mutate(sum_per_tot = sum(per_tot)) %>%
    ungroup() %>%
    distinct(Kingdom,Phylum,Class,Order,Family,Genus,Species, SampleID, .keep_all=TRUE) %>%
    filter(Ecological_Category == var) %>%
    #filter(Ecol_Cat2 == var) %>%
    select(-reads, -per_tot) %>%
    left_join(meta %>% select(SampleID,depth_bin, depth)) %>%
    mutate(count =1) %>%
    group_by(Kingdom,Phylum,Class,Order,Family,Genus,Species,depth_bin) %>%
    mutate(mean_reads = mean(sum_reads)) %>%
    mutate(mean_per_tot = mean(sum_per_tot)) %>%
    mutate(dsum_reads = sum(sum_reads)) %>%
    mutate(dsum_per_tot = sum(sum_per_tot)) %>%
    mutate(sum_count = sum(count)) %>%
    ungroup() %>%
    distinct(Kingdom,Phylum,Class,Order,Family,Genus,Species, depth_bin, .keep_all=TRUE) %>%
    select(-sum_reads, -sum_per_tot, -SampleID, -count, -ASV)
  
  
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- merged_data %>%
    #filter(depth>-1) %>%
    ggplot(aes(x = fct_reorder(depth_bin, desc(depth)), y = mean_reads)) +
    geom_bar(stat = "identity", aes(fill = Common_name))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="Depth Bin",y="Mean reads per sample", title = var)+
    theme_minimal() +
    coord_flip() +
    guides(fill=guide_legend(ncol=1)) +
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
      axis.title.y = element_text(size=8),
      axis.title.x = element_text(size=8),
      plot.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = .25),
      plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))
  bp_top
  filename = paste(plot_directory, marker,'_',var,'_EcolCat_sp_bydepth_mean_reads.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 6, width =5, units = 'in')
  
  # assign text colour
  textcol <- "grey40"
  print("Begin plotting...")
  bp_top <- merged_data %>%
    #filter(depth>-1) %>%
    ggplot(aes(x = fct_reorder(depth_bin, desc(depth)), y = mean_per_tot)) +
    geom_bar(stat = "identity", aes(fill = Common_name))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="Depth Bin",y="Mean percent reads per sample", title = var)+
    theme_minimal() +
    coord_flip() +
    guides(fill=guide_legend(ncol=1)) +
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
      axis.title.y = element_text(size=8),
      axis.title.x = element_text(size=8),
      plot.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = .25),
      plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))
  bp_top
  filename = paste(plot_directory, marker,'_',var,'_EcolCat_sp_bydepth_mean_percentreads.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 6, width =5, units = 'in')
  
  
}
