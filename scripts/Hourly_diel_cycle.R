#082422
#kpitz

# figure tracking reads in binned depth through time (hourly)

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

# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- 'figures/Hourly_diel_cycle/'
# Set directory to retrieve data
data_directory = "Data/filtered_seq_data/"

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

# Ecological Reads in Depth Bins through Time -----------------------------

# take means of sequenced replicate filters (duplicate FilterID)

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot) %>%
  filter(Ecological_Category == "mesopelagic") %>%
  left_join(meta %>% select(SampleID, FilterID, depth_bin, local_time,Sampling_method, diel, hour, ESP, depth_bin2)) %>%
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour == 24, 0)) %>%
  filter(depth_bin2 %in% c("0-100m","100-200m","200-300m")) %>%
  # get mean values of replicate sequenced filters (some ESP samples)
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  #filter(depth_bin2 == "200-300m") %>%
  group_by(local_time, depth_bin2) %>%
  mutate(mean_per_tot = mean(sum_per_tot)) %>%
  ungroup()
# add in 24 hours of data on either side to get smoothed
test <- df %>% select(hour, depth_bin2, sum_per_tot, SampleID, ESP) %>%
  mutate(hour = hour +24)
test2 <- df %>% select(hour, depth_bin2, sum_per_tot, SampleID, ESP) %>%
  mutate(hour = hour -24)

# assign text colour
textcol <- "grey40"
print("Begin plotting...")


df %>% full_join(test) %>%
  full_join(test2) %>%
  ggplot(aes(x=hour, y=sum_per_tot, color=depth_bin2, fill=depth_bin2, linetype=depth_bin2))+
  geom_point(aes(shape=ESP))+
  geom_smooth(span=0.3)+
  coord_cartesian(xlim = c(0, 23), expand = FALSE)+
  #scale_color_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  #scale_fill_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", title='Mesopelagic', color="depth bin", fill="depth bin", linetype="depth bin")+
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(size=8,face="bold"),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))


filename = paste(plot_directory, marker,'_Mesopelagic_byhour.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_Mesopelagic_byhour.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

# just 200-300m

df %>% full_join(test) %>%
  full_join(test2) %>%
  filter(depth_bin2=="200-300m") %>%
  #ggplot(aes(x=hour, y=sum_per_tot, color=depth_bin2, fill=depth_bin2, linetype=depth_bin2))+
  ggplot(aes(x=hour, y=sum_per_tot))+
  #geom_point(aes(shape=ESP))+
  geom_point(aes(shape=ESP, color=diel, fill=diel))+
  geom_smooth(span=0.3)+
  coord_cartesian(xlim = c(0, 23), expand = FALSE)+
  #scale_color_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  #scale_fill_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", title='Mesopelagic, >200-300m samples')+
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(size=8,face="bold"),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))


filename = paste(plot_directory, marker,'_Mesopelagic_byhour_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_Mesopelagic_byhour_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

#mesopelagic through time

df %>% full_join(test) %>%
  full_join(test2) %>%
  filter(depth_bin2=="200-300m") %>%
  #ggplot(aes(x=hour, y=sum_per_tot, color=depth_bin2, fill=depth_bin2, linetype=depth_bin2))+
  ggplot(aes(x=local_time, y=sum_per_tot))+
  #geom_point(aes(shape=ESP))+
  geom_point(aes(shape=ESP, color=diel, fill=diel))+
  geom_smooth(span=0.2)+
  #coord_cartesian(xlim = c(0, 23), expand = FALSE)+
  #scale_color_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  #scale_fill_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", title='Mesopelagic, >200-300m samples')+
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(size=8,face="bold"),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))


filename = paste(plot_directory, marker,'_Mesopelagic_throughtime_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_Mesopelagic_throughtime_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')


# Just look at narrower depth bin -----------

# take means of sequenced replicate filters (duplicate FilterID)

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot) %>%
  filter(Ecological_Category == "mesopelagic") %>%
  left_join(meta %>% select(SampleID, FilterID, depth_bin, local_time,Sampling_method, diel, hour, ESP, depth_bin2)) %>%
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour == 24, 0)) %>%
  #filter(depth_bin %in% c("00_0-25m", "01_25-75m","03_75-100m","04_100-150m", "05_150-200m", "06_200-250m", "07_250-300m", "08_300-400m")) %>%
  filter(depth_bin %in% c("00_0-25m", "01_25-75m","03_75-100m","04_100-150m", "05_150-200m", "06_200-250m", "07_250-300m", "08_300-400m")) %>%
  # get mean values of replicate sequenced filters (some ESP samples)
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  #filter(depth_bin2 == "200-300m") %>%
  group_by(local_time, depth_bin) %>%
  mutate(mean_per_tot = mean(sum_per_tot)) %>%
  ungroup()
# add in 24 hours of data on either side to get smoothed
test <- df %>% select(hour, depth_bin, sum_per_tot, SampleID, ESP) %>%
  mutate(hour = hour +24)
test2 <- df %>% select(hour, depth_bin, sum_per_tot, SampleID, ESP) %>%
  mutate(hour = hour -24)

# assign text colour
textcol <- "grey40"
print("Begin plotting...")

df %>% full_join(test) %>%
  full_join(test2) %>%
  ggplot(aes(x=hour, y=sum_per_tot, color=depth_bin, fill=depth_bin, linetype=depth_bin))+
  geom_point(aes(shape=ESP))+
  geom_smooth(span=0.3, alpha=0.4,se=FALSE)+
  coord_cartesian(xlim = c(0, 23), expand = FALSE)+
  #scale_color_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  #scale_fill_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", title='Mesopelagic', color="depth bin", fill="depth bin", linetype="depth bin")+
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(size=8,face="bold"),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))


filename = paste(plot_directory, marker,'_Mesopelagic_byhour_smallbins.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_Mesopelagic_byhour_smallbins.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

bins = c("00_0-25m", "01_25-75m","03_75-100m","04_100-150m", "05_150-200m", "06_200-250m", "07_250-300m", "08_300-400m")
bins[1]

for (bin in bins) {
  print(bin)
  plot_title = paste('Mesopelagic', bin, sep=', ')
  df %>% full_join(test) %>%
    full_join(test2) %>%
    filter(depth_bin==bin) %>%
    #ggplot(aes(x=hour, y=sum_per_tot, color=depth_bin, fill=depth_bin, linetype=depth_bin))+
    ggplot(aes(x=hour, y=sum_per_tot))+
    #geom_point(aes(shape=ESP))+
    geom_point(aes(shape=ESP, color=diel, fill=diel))+
    geom_smooth(span=0.3)+
    coord_cartesian(xlim = c(0, 23), expand = FALSE)+
    #scale_color_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
    #scale_fill_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
    theme_minimal() +
    guides(fill=guide_legend(ncol=1)) +
    labs(y="Percent Total Reads",x="Hour (24-hour)", title=plot_title)+
    theme(
      #legend
      legend.position="right",legend.direction="vertical",
      legend.text=element_text(size=8,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(size=8,face="bold"),
      plot.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = .25),
      plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))
  
  filename = paste(plot_directory, marker,'_Mesopelagic_byhour_',bin,'.png', sep='')
  #filename = paste(plot_directory, marker,'_Mesopelagic_byhour_200300.png', sep='')
  print(filename)
  ggsave(filename,height = 5, width =8, units = 'in')
  filename = paste(plot_directory, marker,'_Mesopelagic_byhour_',bin,'.svg', sep='')
  #filename = paste(plot_directory, marker,'_Mesopelagic_byhour_200300.svg', sep='')
  print(filename)
  ggsave(filename,height = 5, width =8, units = 'in')
  
  
}
df %>% full_join(test) %>%
  full_join(test2) %>%
  filter(depth_bin==bins[1]) %>%
  #ggplot(aes(x=hour, y=sum_per_tot, color=depth_bin, fill=depth_bin, linetype=depth_bin))+
  ggplot(aes(x=hour, y=sum_per_tot))+
  #geom_point(aes(shape=ESP))+
  geom_point(aes(shape=ESP, color=diel, fill=diel))+
  geom_smooth(span=0.3)+
  coord_cartesian(xlim = c(0, 23), expand = FALSE)+
  #scale_color_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  #scale_fill_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", title='Mesopelagic, >200-300m samples')+
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(size=8,face="bold"),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))


filename = paste(plot_directory, marker,'_Mesopelagic_byhour_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_Mesopelagic_byhour_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

#mesopelagic through time

df %>% full_join(test) %>%
  full_join(test2) %>%
  filter(depth_bin=="200-300m") %>%
  #ggplot(aes(x=hour, y=sum_per_tot, color=depth_bin, fill=depth_bin, linetype=depth_bin))+
  ggplot(aes(x=local_time, y=sum_per_tot))+
  #geom_point(aes(shape=ESP))+
  geom_point(aes(shape=ESP, color=diel, fill=diel))+
  geom_smooth(span=0.2)+
  #coord_cartesian(xlim = c(0, 23), expand = FALSE)+
  #scale_color_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  #scale_fill_manual(values=c("chartreuse4", "darkorange", "darkblue")) +
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", title='Mesopelagic, >200-300m samples')+
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(size=8,face="bold"),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"))


filename = paste(plot_directory, marker,'_Mesopelagic_throughtime_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_Mesopelagic_throughtime_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
