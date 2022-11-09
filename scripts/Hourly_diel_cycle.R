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
library(wesanderson)
library(cowplot)
#library(gridExtra)
#library(patchwork) #putting graphs together on same plot

# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- 'figures/Hourly_diel_cycle/'
# Set directory to retrieve data
data_directory = "Data/filtered_seq_data/"

paletteDayNight <- c(wes_palette("Chevalier1", type = "discrete")[2], wes_palette("Darjeeling2", type = "discrete")[2], 'darkgrey')
paletteEcolCat <- c('blue3', 'darkorchid', 'deepskyblue','chartreuse' )

# ticks at 0,5,10,15,20 hours
x_ticks_hr <- c(0,5,10,15,20)

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

# alternate way of plotting 24-hours (7am to 6am
# as hours since 7am)
meta %<>% mutate(hour_int = as.numeric(hour)) %>%
  mutate(hr_since_7 = NaN) %>%
  mutate(hr_since_7 = case_when(hour_int >=7 ~ hour_int-7,
                                hour_int <7 ~ 18 + hour_int, TRUE ~ NaN))


# Import Acoustic Data ----------------------------------------------------

filepath <- 'data/acoustic_data/Spring canon acoustics summaries_200_300m.csv'
ac_sum <- read_csv(filepath) 

# Plot Acoustic Data -----------------------------------------
glimpse(ac_sum)

#geom_smooth

x_ticks_min <- x_ticks_hr *60*60

p1 <- ac_sum %>%
  mutate(local_time2 = hm(local_time))%>%
  mutate(duration = as.duration(local_time2)) %>%
  ggplot(aes(x = duration, y=mean)) +
  geom_smooth(aes(ymin=mean-SD, ymax =mean+SD),  stat="identity")+
  geom_errorbar(aes(ymin=mean-SD, ymax =mean+SD), width=.2)+
  coord_cartesian(xlim = c(0, 86400), expand = FALSE)+
  scale_x_continuous(breaks=x_ticks_min, labels= x_ticks_hr) +
  theme_minimal() +
  labs(y=expression(paste("Mean Scattering 200-300m "," (dB re 1", m^-1, ')')))
  #ylab("Mean Scattering 200-300 m (dB re 1m^-1)")

p1
filename = paste(plot_directory, 'Acoustic_DayNight_duration_smooth.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')



# Plot Ecological Groups in Depth Bins through Time 200-300m --------------

# 59 unique samples in the 200-300m bin (80 including sequenced duplicates)
# should be 59*4 = 236  observations across all Ecol Cats

## Original ---------
df <- tax.c %>% 
  # create merged reads and percent reads for each ecological group
  left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot) %>%
  filter(Ecological_Category %in% c("mesopelagic", "epipelagic", "cosmopolitan", "benthopelagic")) %>%
  #filter(Ecological_Category %in% c("mesopelagic", "epipelagic", "cosmopolitan")) %>%
  left_join(meta %>% select(SampleID, FilterID, depth_bin, local_time,Sampling_method, diel, hour, ESP, depth_bin2)) %>%
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour == 24, 0)) %>%
  filter(depth_bin2 %in% c("200-300m")) %>%
  # get mean values of replicate sequenced filters (some ESP samples)
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) 
  # now have 236 observations from 59 unique filters

# add in 24 hours of data on either side to get smoothed
test <- df %>% select(Ecological_Category, hour, depth_bin2, sum_per_tot, SampleID, ESP) %>%
  mutate(hour = hour +24)
test2 <- df %>% select(Ecological_Category, hour, depth_bin2, sum_per_tot, SampleID, ESP) %>%
  mutate(hour = hour -24)

p2 <- df %>% full_join(test) %>%
  full_join(test2) %>%
  #filter(depth_bin2=="200-300m") %>%
  ggplot(aes(x=hour, y=sum_per_tot, color=Ecological_Category, fill=Ecological_Category))+
  #geom_point(aes(shape=ESP))+
  #don't repeat 0hr and 24hr points:
  geom_point(data = df %>% full_join(test) %>%
               full_join(test2) %>% filter(hour!=24), aes(color=Ecological_Category, fill=Ecological_Category))+
  geom_smooth(span=0.3)+
  #geom_rug(sides='b', alpha=0.1, size=4) +
  coord_cartesian(xlim = c(0, 24), expand = FALSE)+
  scale_x_continuous(breaks=x_ticks_hr) +
  scale_color_manual(values=paletteEcolCat )+
  scale_fill_manual(values=paletteEcolCat)+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", fill='Ecological Category', color='Ecological Category')+
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

p2
filename = paste(plot_directory, marker,'_AllEco_byhour_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_AllEco_byhour_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

#Just mesopelagic
p <- df %>% full_join(test) %>%
  full_join(test2) %>%
  #filter(depth_bin2=="200-300m") %>%
  filter(Ecological_Category == 'mesopelagic') %>%
  ggplot(aes(x=hour, y=sum_per_tot, color=Ecological_Category, fill=Ecological_Category))+
  #geom_point(aes(shape=ESP))+
  #don't repeat 0hr and 24hr points:
  geom_point(data = df %>% full_join(test) %>%
               full_join(test2) %>% filter(hour!=24), aes(color=Ecological_Category, fill=Ecological_Category))+
  geom_smooth(span=0.3)+
  #geom_rug(sides='b', alpha=0.1, size=4) +
  coord_cartesian(xlim = c(0, 24), expand = FALSE)+
  scale_x_continuous(breaks=x_ticks_hr) +
  scale_color_manual(values=paletteEcolCat )+
  scale_fill_manual(values=paletteEcolCat)+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", fill='Ecological Category', color='Ecological Category')+
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

p
filename = paste(plot_directory, marker,'_Meso_byhour_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_Meso_byhour_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

## Split out Ecol Cats ---------------

p2 <- df %>% full_join(test) %>%
  full_join(test2) %>%
  #filter(depth_bin2=="200-300m") %>%
  ggplot(aes(x=hour, y=sum_per_tot, color=Ecological_Category, fill=Ecological_Category))+
  #geom_point(aes(shape=ESP))+
  #don't repeat 0hr and 24hr points:
  geom_point(data = df %>% full_join(test) %>%
               full_join(test2) %>% filter(hour!=24), aes(color=Ecological_Category, fill=Ecological_Category))+
  geom_smooth(span=0.3)+
  #geom_rug(sides='b', alpha=0.1, size=4) +
  coord_cartesian(xlim = c(0, 24), expand = FALSE)+
  scale_x_continuous(breaks=x_ticks_hr) +
  scale_color_manual(values=paletteEcolCat )+
  scale_fill_manual(values=paletteEcolCat)+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", fill='Ecological Category', color='Ecological Category')+
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

#Horizonal split
p2h <- p2 + facet_grid( ~ Ecological_Category, margins = FALSE)
p2h
filename = paste(plot_directory, marker,'_SplitEco_horiz_byhour_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')
filename = paste(plot_directory, marker,'_SplitEco_horiz_byhour_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =6, units = 'in')

#Box split
p2b <- p2 + facet_wrap( Ecological_Category ~., scales="free")
p2b
filename = paste(plot_directory, marker,'_SplitEco_box_byhour_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')
filename = paste(plot_directory, marker,'_SplitEco_box_byhour_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

#Vertical split
p2v <- p2 + facet_wrap( Ecological_Category ~.,  scales="free", ncol = 1)
p2v
filename = paste(plot_directory, marker,'_SplitEco_vert_byhour_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')
filename = paste(plot_directory, marker,'_SplitEco_vert_byhour_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

# 2 Ecol cats each:
p2me <- df %>% full_join(test) %>%
  full_join(test2) %>%
  #filter(depth_bin2=="200-300m") %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  ggplot(aes(x=hour, y=sum_per_tot, color=Ecological_Category, fill=Ecological_Category))+
  #geom_point(aes(shape=ESP))+
  #don't repeat 0hr and 24hr points:
  geom_point(data = df %>% full_join(test) %>%
               filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
               full_join(test2) %>% filter(hour!=24), aes(color=Ecological_Category, fill=Ecological_Category))+
  geom_smooth(span=0.3)+
  #geom_rug(sides='b', alpha=0.1, size=4) +
  coord_cartesian(xlim = c(0, 24), expand = FALSE)+
  scale_x_continuous(breaks=x_ticks_hr) +
  scale_color_manual(values=paletteEcolCat )+
  scale_fill_manual(values=paletteEcolCat)+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", fill='Ecological Category', color='Ecological Category')+
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
p2me

p2bc <- df %>% full_join(test) %>%
  full_join(test2) %>%
  #filter(depth_bin2=="200-300m") %>%
  filter(Ecological_Category %in% c('benthopelagic', 'cosmopolitan')) %>%
  ggplot(aes(x=hour, y=sum_per_tot, color=Ecological_Category, fill=Ecological_Category))+
  #geom_point(aes(shape=ESP))+
  #don't repeat 0hr and 24hr points:
  geom_point(data = df %>% full_join(test) %>%
               filter(Ecological_Category %in% c('benthopelagic', 'cosmopolitan')) %>%
               full_join(test2) %>% filter(hour!=24), aes(color=Ecological_Category, fill=Ecological_Category))+
  geom_smooth(span=0.3)+
  #geom_rug(sides='b', alpha=0.1, size=4) +
  coord_cartesian(xlim = c(0, 24), expand = FALSE)+
  scale_x_continuous(breaks=x_ticks_hr) +
  scale_color_manual(values=paletteEcolCat )+
  scale_fill_manual(values=paletteEcolCat)+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Percent Total Reads",x="Hour (24-hour)", fill='Ecological Category', color='Ecological Category')+
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
p2bc


# Sampling Effort --------

# 59 unique samples in the 200-300m bin (80 including sequenced duplicates)

p3 <- meta %>% 
  filter(depth_bin2 == "200-300m") %>%
  #remove ESP replicates
  distinct(local_time, depth, SAMPLING_bottle, .keep_all = TRUE) %>%
  mutate(count=1) %>%
  ggplot(aes(x=hour,y=count, color=diel, fill=diel))+
  geom_bar(stat='identity')+
  coord_cartesian(xlim = c(0, 24), expand = FALSE)+
  scale_y_continuous(breaks=c(2,4,6,8,10)) +
  scale_x_continuous(breaks=x_ticks_hr) +
  scale_color_manual(values=paletteDayNight )+
  scale_fill_manual(values=paletteDayNight)+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Observations",x="Hour (24-hour)")+
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
p3
filename = paste(plot_directory, marker,'_samplingeffort_byhour_200-300m.png', sep='')
print(filename)
ggsave(filename,height = 2, width =8, units = 'in')
filename = paste(plot_directory, marker,'_samplingeffort_byhour_200-300m.svg', sep='')
print(filename)
ggsave(filename,height = 2, width =8, units = 'in')

# Join plots --------------

## Vertical combined -------
#Put all together
pf1 <- p1 + theme_minimal() + theme(text = element_text(size=10), 
                                    axis.title = element_text(size = 8),
                                    legend.position = 'none',
                                    axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
                                    axis.line = element_line(color = "black", size=0.2),
                                    axis.ticks = element_line(color = "black", size=0.2),
)
# two vertical plots with 2 ecol cats each
# first mesopelagic + epipelagic ; then cosmopolitan and benthopelagic

pf2 <- p2me + theme_minimal() + theme(text = element_text(size=10), 
                                     legend.position = 'none',
                                     axis.title.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     #axis.ticks.x=element_blank(),
                                     panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank(),
                                     strip.background = element_blank(), 
                                     strip.text.x = element_blank(),
                                     #axis.title.y=element_blank(),
                                     #axis.text.y=element_blank(),
                                     #axis.ticks.y=element_blank(),
                                     #axis.line.x  = element_line(color = "black"), 
                                     #axis.line.y  = element_line(color = "black")
                                     axis.line = element_line(color = "black", size=0.2),
                                     axis.ticks = element_line(color = "black", size=0.2),
)

pf4 <- p2bc + theme_minimal() + theme(text = element_text(size=10), 
                                      legend.position = 'none',
                                      axis.title.x=element_blank(),
                                      axis.text.x=element_blank(),
                                      #axis.ticks.x=element_blank(),
                                      panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(),
                                      strip.background = element_blank(), 
                                      strip.text.x = element_blank(),
                                      #axis.title.y=element_blank(),
                                      #axis.text.y=element_blank(),
                                      #axis.ticks.y=element_blank(),
                                      #axis.line.x  = element_line(color = "black"), 
                                      #axis.line.y  = element_line(color = "black")
                                      axis.line = element_line(color = "black", size=0.2),
                                      axis.ticks = element_line(color = "black", size=0.2),
)

pf3 <- p3 + theme_minimal() + theme(text = element_text(size=10), 
                                    legend.position = 'none',
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
                                    axis.line = element_line(color = "black", size=0.2),
                                    axis.ticks = element_line(color = "black", size=0.2),
)

library(cowplot)

# legend
legend <- get_legend(p2)

#aligned <- align_plots(pf1, pf2,pf3, align = "v")
aligned <- align_plots(pf1, pf2, pf3, align = "v")

#ggdraw(aligned[[1]])
#plot_grid(pf1, pf2, ncol = 1, align = "v")

plot_grid(pf1, pf2,pf4, pf3,ncol = 1, align = "v", rel_heights = c(1,1,1, 0.5))

filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp_splitv_comb.png', sep='')
filename
ggsave(filename, width=4, height=8, units = 'in')
filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp_splitv_comb.svg', sep='')
filename
ggsave(filename, width=4, height=8, units = 'in')

## Vertical --------------
#Put all together
pf1 <- p1 + theme_minimal() + theme(text = element_text(size=13), 
                                    legend.position = 'none',
                                    axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
)
pf2 <- p2v + theme_minimal() + theme(text = element_text(size=13), 
                                     legend.position = 'none',
                                     axis.title.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank(),
                                     panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank(),
                                     strip.background = element_blank(), 
                                     strip.text.x = element_blank(),
                                     #axis.title.y=element_blank(),
                                     #axis.text.y=element_blank(),
                                     #axis.ticks.y=element_blank(),
                                     #axis.line.x  = element_line(color = "black"), 
                                     #axis.line.y  = element_line(color = "black")
)
pf3 <- p3 + theme_minimal() + theme(text = element_text(size=13), 
                                    legend.position = 'none',
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
)

library(cowplot)

# legend
legend <- get_legend(p2)

#aligned <- align_plots(pf1, pf2,pf3, align = "v")
aligned <- align_plots(pf1, pf2, pf3, align = "v")

#ggdraw(aligned[[1]])
#plot_grid(pf1, pf2, ncol = 1, align = "v")

plot_grid(pf1, pf2, pf3,ncol = 1, align = "v", rel_heights = c(1.2,4.5, 0.6))

filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp_splitv.png', sep='')
filename
ggsave(filename, width=4, height=15, units = 'in')
filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp_splitv.svg', sep='')
filename
ggsave(filename, width=4, height=15, units = 'in')
## Box --------------

#Put all together
pf1 <- p1 + theme_minimal() + theme(text = element_text(size=13), 
                                    legend.position = 'none',
                                    axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
)
pf2 <- p2b + theme_minimal() + theme(text = element_text(size=13), 
                                     legend.position = 'none',
                                     #axis.title.x=element_blank(),
                                     #axis.text.x=element_blank(),
                                     #axis.ticks.x=element_blank(),
                                     panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank(),
                                     #axis.title.y=element_blank(),
                                     #axis.text.y=element_blank(),
                                     #axis.ticks.y=element_blank(),
                                     #axis.line.x  = element_line(color = "black"), 
                                     #axis.line.y  = element_line(color = "black")
)
pf3 <- p3 + theme_minimal() + theme(text = element_text(size=13), 
                                    legend.position = 'none',
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
)

library(cowplot)

# legend
legend <- get_legend(p2)

#aligned <- align_plots(pf1, pf2,pf3, align = "v")
aligned <- align_plots(pf1, pf2, pf3, align = "v")


plot_grid(pf1, pf3, pf2,ncol = 1, align = "v", rel_heights = c(1.2,0.5,2.5))

filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp_splitb.png', sep='')
filename
ggsave(filename, width=4, height=12, units = 'in')
filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp_splitb.svg', sep='')
filename
ggsave(filename, width=4, height=12, units = 'in')

## Horizontal --------------
#Put all together
pf1 <- p1 + theme_minimal() + theme(text = element_text(size=13), 
                                    legend.position = 'none',
                                    axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
                                    )
pf2 <- p2h + theme_minimal() + theme(text = element_text(size=13), 
                                    legend.position = 'none',
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
)
pf3 <- p3 + theme_minimal() + theme(text = element_text(size=13), 
                                    legend.position = 'none',
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
)

library(cowplot)

# legend
legend <- get_legend(p2)

#aligned <- align_plots(pf1, pf2,pf3, align = "v")
aligned <- align_plots(pf1, pf3, pf2, align = "v")

#ggdraw(aligned[[1]])
#plot_grid(pf1, pf2, ncol = 1, align = "v")

plot_grid(pf1, pf3, pf2,ncol = 1, align = "v", rel_heights = c(1,0.3, 1))

filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp_splith.png', sep='')
filename
ggsave(filename, width=8, height=10, units = 'in')
filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp_splith.svg', sep='')
filename
ggsave(filename, width=8, height=10, units = 'in')

## Original --------------
#Put all together
pf1 <- p1 + theme_minimal() + theme(text = element_text(size=13), 
                                    legend.position = 'none',
                                    axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
)
pf2 <- p2 + theme_minimal() + theme(text = element_text(size=13), 
                                     legend.position = 'none',
                                     axis.title.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                     #axis.title.y=element_blank(),
                                     #axis.text.y=element_blank(),
                                     #axis.ticks.y=element_blank(),
                                     #axis.line.x  = element_line(color = "black"), 
                                     #axis.line.y  = element_line(color = "black")
)
pf3 <- p3 + theme_minimal() + theme(text = element_text(size=13), 
                                    legend.position = 'none',
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    #axis.line.x  = element_line(color = "black"), 
                                    #axis.line.y  = element_line(color = "black")
)

library(cowplot)

# legend
legend <- get_legend(p2)

aligned <- align_plots(pf1, pf2,pf3, align = "v")
#ggdraw(aligned[[1]])
#plot_grid(pf1, pf2, ncol = 1, align = "v")

plot_grid(pf1, pf2,pf3 ,ncol = 1, align = "v", rel_heights = c(1, 1, 0.3))

filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp.png', sep='')
filename
ggsave(filename, width=8, height=10, units = 'in')
filename = paste(plot_directory, 'Acoustic_Ecolcat_Samp.svg', sep='')
filename
ggsave(filename, width=8, height=10, units = 'in')

## Exploration below ----------

# for each unique time point, plot boxplot
df %>% #full_join(test) %>%
  #full_join(test2) %>%
  #filter(depth_bin2=="200-300m") %>%
  mutate(hour = as.factor(hour)) %>%
  ggplot(aes(x=hour, y=sum_per_tot, color=Ecological_Category, fill=Ecological_Category))+
  geom_boxplot(alpha=0.3) 

# number of observations

p <- meta %>% 
  filter(depth_bin2 == "200-300m") %>%
  #remove ESP replicates
  distinct(local_time, depth, .keep_all = TRUE) %>%
  mutate(count=1) %>%
  ggplot(aes(x=hour,y=count, color=ESP, fill=ESP))+
  geom_bar(stat='identity')+
  scale_y_continuous(breaks=c(2,4,6,8,10))
p

p <- meta %>% 
  filter(depth_bin2 == "200-300m") %>%
  #remove ESP replicates
  distinct(local_time, depth, .keep_all = TRUE) %>%
  mutate(count=1) %>%
  ggplot(aes(x=hour,y=count, color=diel, fill=diel))+
  geom_bar(stat='identity')+
  scale_y_continuous(breaks=c(2,4,6,8,10)) +
  scale_color_manual(values=paletteDayNight )+
  scale_fill_manual(values=paletteDayNight)+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  labs(y="Number of Observations",x="Hour (24-hour)", title='Sampling Effort')+
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
p
filename = paste(plot_directory, marker,'_samplingeffort_byhour_200-300m.png', sep='')
print(filename)
ggsave(filename,height = 2, width =8, units = 'in')
filename = paste(plot_directory, marker,'_samplingeffort_byhour_200-300m.svg', sep='')
print(filename)
ggsave(filename,height = 2, width =8, units = 'in')


df <-  tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot) %>%
  filter(Ecological_Category %in% c("mesopelagic", "epipelagic", "cosmopolitan", "benthopelagic")) %>%
  #filter(Ecological_Category %in% c("mesopelagic", "epipelagic", "cosmopolitan")) %>%
  left_join(meta %>% select(SampleID, FilterID, depth_bin, local_time,Sampling_method, diel, hour, ESP, depth_bin2)) %>%


df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot) %>%
  filter(Ecological_Category %in% c("mesopelagic", "epipelagic", "cosmopolitan", "benthopelagic")) %>%
  #filter(Ecological_Category %in% c("mesopelagic", "epipelagic", "cosmopolitan")) %>%
  left_join(meta %>% select(SampleID, FilterID, depth_bin, local_time,Sampling_method, diel, hour, ESP, depth_bin2)) %>%
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour == 24, 0)) %>%
  filter(depth_bin2 %in% c("200-300m")) %>%
  # get mean values of replicate sequenced filters (some ESP samples)
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  group_by(Ecological_Category,local_time, depth_bin2) %>%
  mutate(mean_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct( Ecological_Category,local_time, depth_bin2, .keep_all=TRUE)



+
  #geom_point(aes(shape=ESP))+
  geom_point(aes(shape=ESP, color=Ecological_Category, fill=Ecological_Category))+
  geom_smooth(span=0.3)+
  geom_rug(sides='b', alpha=0.1, size=4) +
  coord_cartesian(xlim = c(0, 23), expand = FALSE)+
  scale_color_manual(values=paletteEcolCat )+
  scale_fill_manual(values=paletteEcolCat)+
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

# Plot Acoustic Data -----------------------------------------
glimpse(ac_sum)

# p <- ac_sum %>%
#   mutate(local_time = hm(local_time))%>%
#   mutate(duration = as.duration(local_time)) %>%
#   ggplot(aes(x = duration, y=mean))+
#   geom_point(alpha=0.7)+
#   geom_errorbar(aes(ymin=mean-SD, ymax =mean+SD), width=.2)+
#   geom_line()+theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# p
# filename = paste(plot_directory, 'Acoustic_DayNight_duration.png', sep='')
# #print('Plot of top 20 Genus average by month:')
# print(filename)
# ggsave(filename,height = 5, width =8, units = 'in')
# 
# p <- ac_sum %>%
#   mutate(local_time2 = hm(local_time))%>%
#   mutate(duration = as.duration(local_time2)) %>%
#   ggplot(aes(x = fct_reorder(local_time, duration), y=mean))+
#   geom_point(alpha=0.7)+
#   geom_errorbar(aes(ymin=mean-SD, ymax =mean+SD), width=.2)+
#   geom_line()+theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# p
# filename = paste(plot_directory, 'Acoustic_DayNight_point.png', sep='')
# #print('Plot of top 20 Genus average by month:')
# print(filename)
# ggsave(filename,height = 5, width =8, units = 'in')

#geom_smooth


# Mesopelagic in Depth Bins through Time 200-300m-----------------------------

# take means of sequenced replicate filters (duplicate FilterID)


df <- tax.c %>% 
  # create merged reads and percent reads for each ecological group
  left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot) %>%
  filter(Ecological_Category %in% c("mesopelagic", "epipelagic", "cosmopolitan", "benthopelagic")) %>%
  #filter(Ecological_Category %in% c("mesopelagic", "epipelagic", "cosmopolitan")) %>%
  left_join(meta %>% select(SampleID, FilterID, depth_bin, local_time,Sampling_method, diel, hour, ESP, depth_bin2)) %>%
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour == 24, 0)) %>%
  filter(depth_bin2 %in% c("200-300m")) %>%
  # get mean values of replicate sequenced filters (some ESP samples)
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  group_by(Ecological_Category,local_time, depth_bin2) %>%
  mutate(mean_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct( Ecological_Category,local_time, depth_bin2, .keep_all=TRUE)



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


# mesopelagic eDNA and acoustic data 200-300m  ----------------------------

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


filename = paste(plot_directory, marker,'_Mesopelagic_throughtime_200300_ac.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_Mesopelagic_throughtime_200300_ac.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

# Cosmopolitan in Depth Bins through Time 200-300m-----------------------------

# take means of sequenced replicate filters (duplicate FilterID)

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot) %>%
  filter(Ecological_Category == "cosmopolitan") %>%
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


filename = paste(plot_directory, marker,'_cosmopolitan_byhour.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_cosmopolitan_byhour.svg', sep='')
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
  labs(y="Percent Total Reads",x="Hour (24-hour)", title='cosmopolitan, >200-300m samples')+
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


filename = paste(plot_directory, marker,'_cosmopolitan_byhour_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_cosmopolitan_byhour_200300.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

# Epipelagic in Depth Bins through Time 200-300m-----------------------------

# take means of sequenced replicate filters (duplicate FilterID)

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot) %>%
  filter(Ecological_Category == "epipelagic") %>%
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
  labs(y="Percent Total Reads",x="Hour (24-hour)", title='epipelagic', color="depth bin", fill="depth bin", linetype="depth bin")+
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


filename = paste(plot_directory, marker,'_epipelagic_byhour.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_epipelagic_byhour.svg', sep='')
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
  labs(y="Percent Total Reads",x="Hour (24-hour)", title='epipelagic, >200-300m samples')+
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


filename = paste(plot_directory, marker,'_epipelagic_byhour_200300.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_epipelagic_byhour_200300.svg', sep='')
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
