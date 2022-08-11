#040221
#kpitz

# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- 'figures/Fig1/'
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


# Group data by Ecological Category --------------------------------------------

Ecol_data <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot)


# Ecological Components by Depth - Percent Reads --------------------------

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- full_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  #filter(depth<600) %>%
  filter(diel %in% c('day', 'night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
                                    'cosmopolitan')) %>%
  ggplot(aes(x = -depth, y = sum_per_tot)) +
  geom_point(aes(color = Ecological_Category, shape=ESP), alpha=0.7) +
  scale_shape_manual(values= c(16,8)) +
  geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5, span=0.6) +
  facet_grid(. ~ diel, margins=FALSE)+
  scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_x_continuous(breaks=c(-890,-800,-700,-600,-500,-400,-300,-200,-100,0), limits=c(-890,0))+
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
filename = paste(plot_directory, marker,'_EcolCat_scatter_Diel_depthy.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_EcolCat_scatter_Diel_depthy.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =7, units = 'in')


# Presence Absence --------------------------------------------------------

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_per_tot, sum_reads) %>%
  #filter(sum_reads ==0) %>%
  mutate(present = case_when(sum_per_tot<0.1  ~ 0,
                            TRUE ~ 1)) %>%
  #mutate(present = 1) %>%
  select(SampleID, Ecological_Category, present) %>%
  left_join(meta %>% select(diel, SampleID, depth))


cats <- c("epipelagic","mesopelagic", "bathypelagic","cosmopolitan")
for (val in cats) {
  cat = sym(val)
  p <- df %>% filter(Ecological_Category == cat) %>%
    ggplot(aes(y=present, x=-depth,color=diel, shape=diel))+
    geom_point()+
    geom_smooth()+
    coord_flip()+
    labs(y='Percent Present in eDNA Sample (>0.1% reads)')
  filename = paste(plot_directory, marker,'_',cat,'_PA_01.png', sep='')
  print(filename)
  ggsave(filename,height = 5, width =6, units = 'in')
  # filename = paste(plot_directory, marker,'_',cat,'_PA_01.svg', sep='')
  # print(filename)
  # ggsave(filename,height = 5, width =6, units = 'in')
}


# Ridge Plot of Detections ------------------------------------------------





# Number of Samples by Depth and Diel -------------------------------------

df <- meta %>%
  mutate(count=1) %>%
  group_by(diel, depth_bin) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(diel, depth_bin, sum_count, .keep_all=TRUE) %>%
  arrange(depth_bin, diel) %>%
  ggplot(aes(y=fct_reorder(depth_bin, -depth),x=diel, fill = sum_count))+
  geom_tile()+
  geom_text(aes(label=sum_count), colour = "white", check_overlap = TRUE)+  
  scale_fill_viridis()+
  labs(y='Depth', x='', fill='Number of Samples' )
filename = paste(plot_directory, marker,'_','sample_heatmap.png', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
filename = paste(plot_directory, marker,'_','sample_heatmap.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
df  

# depths split evenly:
bp_top <- meta %>%
  select(depth, diel, ESP, SampleID) %>%
  mutate(count=1) %>%
  mutate( ints = cut(depth ,breaks=25, include.lowest = TRUE)) %>%
  group_by(ints, diel) %>% 
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(ints,diel,.keep_all = TRUE) %>%
  ggplot(aes(y= fct_reorder(ints, -depth), x=diel, fill=sum_count)) +
  #facet_grid(. ~ diel, margins=FALSE, scales = "free")+
  geom_tile()+
  geom_text(aes(label=sum_count), colour = "white", check_overlap = TRUE)+  
  #scale_y_reverse()+
  scale_fill_viridis()+
  labs(y="Depth Bins",x="")+
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
filename = paste(plot_directory, marker,'_','sample_heatmap_depthsplit.png', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
filename = paste(plot_directory, marker,'_','sample_heatmap_depthsplit.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
df  

# Histogram plot by depth
bp_top <- meta %>%
  select(depth, diel, ESP, SampleID) %>%
  mutate(count=1) %>%
  filter(diel %in% c("day", "night")) %>%
  #mutate( ints = cut(depth ,breaks=25, include.lowest = TRUE)) %>%
  #group_by(ints, diel) %>% 
  #mutate(sum_count = sum(count)) %>%
  #ungroup() %>%
  #distinct(ints,diel,.keep_all = TRUE) %>%
  ggplot(aes(x= depth,color=diel, fill=diel)) +
  #facet_grid(. ~ diel, margins=FALSE, scales = "free")+
  geom_histogram(aes(y=..density..),alpha=0.5, bins=40)+
  geom_density(alpha=0.5)+
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  #labs(y="Depth",x="Number of Samples")+
  labs(title="Sampling density with depth")+
  theme_minimal() +
  #guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=8,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=8,colour=textcol),
    axis.title.y = element_text(size=8),
    axis.title.x = element_text(size=8),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))

bp_top
filename = paste(plot_directory, marker,'_','sample_hist.png', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
filename = paste(plot_directory, marker,'_','sample_hist.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
df  





# SCRAP



  mutate( ints = cut(depth ,breaks = 25, labels=TRUE)) %>% 
  
  group_by(ints, diel, ESP) %>% 
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(ints,diel, ESP,.keep_all = TRUE)


# merged by diel
bp_top <- meta %>%
  select(depth, diel, ESP, SampleID) %>%
  mutate(count=1) %>%
  mutate( ints = cut(depth ,breaks = 25, labels=FALSE)) %>% 
  group_by(ints, diel, ESP) %>% 
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(ints,diel, ESP,.keep_all = TRUE) %>%
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
filename = paste(plot_directory, marker,'_','sample_heatmap_depthsplit.png', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
filename = paste(plot_directory, marker,'_','sample_heatmap_depthsplit.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
df  

# Box plot -----------------------------------------------------

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE)

textcol <- "grey40"
print("Begin plotting...")

cats <- c("epipelagic","mesopelagic", "bathypelagic","cosmopolitan")
for (val in cats) {
  cat = sym(val)
  bp_top <- full_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  filter(depth>-1) %>%
  #filter(depth<600) %>%
  filter(diel %in% c('day', 'night')) %>%
  # filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
  #                                   'cosmopolitan')) %>%
  filter(Ecological_Category %in% c(cat)) %>%
  #ggplot(aes(x = -depth, y = sum_per_tot)) +
  ggplot(aes(y = fct_reorder(depth_bin, -depth), x = sum_per_tot, fill=diel)) +
  geom_boxplot() +
  facet_grid(. ~ Ecological_Category) +
  labs(y='Depth', x='Percent Reads')
  bp_top
  filename = paste(plot_directory, marker,'_',cat,'_boxplot_diel_depth_perreads.png', sep='')
  #print('Plot of top 20 Genus average by month:')
  print(filename)
  ggsave(filename,height = 6, width =5, units = 'in')
}

