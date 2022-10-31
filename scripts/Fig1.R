#040221
#kpitz

# Create Figure 1 for Manuscript (eDNA portion):
# day versus night line plots of Ecological Category Percent Total Reads
# alongside heatmap of sampling effort

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
library(wesanderson)


# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- 'figures/Fig1/'
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


# Ecological Categories to focus on  ------------------

cats <- c("epipelagic","mesopelagic", "benthopelagic","cosmopolitan")

# Group data by Ecological Category --------------------------------------------

tax_levels = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

Ecol_data <- tax.c %>% full_join(sp_desig, by=tax_levels) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot) %>%
  # average replicate filters
  left_join(meta %>% select(SampleID, FilterID)) %>%
  # get mean values of replicate sequenced filters (some ESP samples)
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE)


# Ecological Components by Depth - Percent Reads --------------------------

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- left_join(Ecol_data, meta %>% select(SampleID, depth, diel, ESP),  by = c("SampleID")) %>% #join with metadata
  filter(diel %in% c('day', 'night')) %>%
  # filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
  #                                   'cosmopolitan')) %>%
  filter(Ecological_Category %in% cats) %>%
  ggplot(aes(x = -depth, y = sum_per_tot)) +
  geom_point(aes(color = Ecological_Category, shape=ESP), alpha=0.7) +
  scale_shape_manual(values= c(16,8)) +
  geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5, span=0.6) +
  facet_grid(. ~ diel, margins=FALSE)+
  scale_color_manual(values=paletteEcolCat )+
  scale_fill_manual(values=paletteEcolCat)+
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


# Add Acoustic Data -------------------------------------------------------

# assign text colour
textcol <- "grey40"
print("Begin plotting...")
bp_top <- left_join(Ecol_data, meta %>% select(SampleID, depth, diel, ESP),  by = c("SampleID")) %>% #join with metadata
  filter(diel %in% c('day', 'night')) %>%
  # filter(Ecological_Category %in% c('mesopelagic', 'epipelagic', 'benthopelagic',
  #                                   'cosmopolitan')) %>%
  filter(Ecological_Category %in% cats) %>%
  ggplot(aes(x = -depth, y = sum_per_tot)) +
  geom_point(aes(color = Ecological_Category, shape=ESP), alpha=0.7) +
  scale_shape_manual(values= c(16,8)) +
  geom_smooth(aes(color = Ecological_Category, fill=Ecological_Category), alpha=0.5, span=0.6) +
  facet_grid(. ~ diel, margins=FALSE)+
  # #add acoustic data
  # geom_point(data=ac_sum, aes(x=-depth, y=mean,shape=diel, size=diel), alpha=0.7)+
  # geom_errorbar(data=ac_sum,aes(ymin=mean-SD, ymax =mean+SD), width=.2)+
  # geom_line(data=ac_sum,aes(x=-depth, y=mean, color=diel))+
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

+ ggplot(data = ac_sum, aes(x = -depth, y=mean, color=diel, group=diel), alpha=0.4)+
  geom_point(aes(shape=diel, size=diel), alpha=0.7)+
  #geom_errorbar(aes(ymin=Mean-SD, ymax =Mean+SD), width=.2,
  #              position=position_dodge(0.05))+
  geom_errorbar(aes(ymin=mean-SD, ymax =mean+SD), width=.2)+
  geom_line()+theme_minimal() +
  coord_flip()

filename = paste(plot_directory, marker,'_EcolCat_scatter_Diel_depthy_Acoustic.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_EcolCat_scatter_Diel_depthy_Acoustic.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =7, units = 'in')






# Number of Samples by Depth and Diel -------------------------------------

# define new depth bins that are evenly split:
meta2 <- meta %>% 
  mutate(depth_bin2 = case_when(depth <=100 ~ "0-100m",
                               depth >100 & depth <=200 ~ "100-200m",
                               depth >200 & depth <=300 ~ "200-300m",
                               depth >300 & depth <=400 ~ "300-400m",
                               depth >400 & depth <=500 ~ "400-500m",
                               depth >400 & depth <=600 ~ "500-600m",
                               depth >600 & depth <=750 ~ "600-700m", TRUE ~ "unknown"
  )) 

df <- meta2 %>%
  #drop replicates
  distinct(FilterID, .keep_all=TRUE) %>%
  mutate(count=1) %>%
  filter(diel %in% c('day', 'night')) %>%
  group_by(diel, depth_bin2) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(diel, depth_bin2, sum_count, .keep_all=TRUE) %>%
  arrange(depth_bin2, diel) %>%
  ggplot(aes(y=fct_reorder(depth_bin2, -depth),x=diel, fill = sum_count))+
  geom_tile()+
  geom_text(aes(label=sum_count), colour = "white", check_overlap = TRUE)+  
  #scale_fill_viridis(option="cividis")+
  #scale_color_grey()+
  labs(y='Depth', x='', fill='Number of Samples' )
df
filename = paste(plot_directory, marker,'_','sample_heatmap.png', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
filename = paste(plot_directory, marker,'_','sample_heatmap.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
df  

# Compare with Presence Absence --------------------------------------------------------

# look above a threshold (0.5%) and also anything > 0% reads is present

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_per_tot, sum_reads) %>%
  # average replicate filters
  left_join(meta %>% select(SampleID, FilterID)) %>%
  # get mean values of replicate sequenced filters (some ESP samples)
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  mutate(present = case_when(sum_per_tot<0.1  ~ 0,
                             TRUE ~ 1)) %>%
  #mutate(present = 1) %>%
  select(SampleID, Ecological_Category, present) %>%
  left_join(meta %>% select(diel, SampleID, depth)) %>%
  filter(diel %in% c('day', 'night'))


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
}

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_per_tot, sum_reads) %>%
  # average replicate filters
  left_join(meta %>% select(SampleID, FilterID)) %>%
  # get mean values of replicate sequenced filters (some ESP samples)
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  mutate(present = case_when(sum_per_tot<0  ~ 0,
                             TRUE ~ 1)) %>%
  #mutate(present = 1) %>%
  select(SampleID, Ecological_Category, present) %>%
  left_join(meta %>% select(diel, SampleID, depth)) %>%
  filter(diel %in% c('day', 'night'))


for (val in cats) {
  cat = sym(val)
  p <- df %>% filter(Ecological_Category == cat) %>%
    ggplot(aes(y=present, x=-depth,color=diel, shape=diel))+
    geom_point()+
    geom_smooth()+
    coord_flip()+
    labs(y='Percent Present in eDNA Sample (>0.0% reads)')
  filename = paste(plot_directory, marker,'_',cat,'_PA_0.png', sep='')
  print(filename)
  ggsave(filename,height = 5, width =6, units = 'in')
  
}

