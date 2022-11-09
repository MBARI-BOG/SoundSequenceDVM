#040221
#kpitz

# Create boxplot figure of mesopelagic and epipelagic percent reads across depth bins
# Add acoustic variability across night/day

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
library(wesanderson)
library(cowplot)
#library(viridis)

# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- 'figures/Boxplot/'
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


# Import Acoustic data -------------------------------------------

filepath <- 'data/acoustic_data/Spring canon acoustics summaries_bydepth_overall.csv'
ac_sum <- read_csv(filepath) 

# Plot Acoustic Data -----------------------------------------
glimpse(ac_sum)
p1 <- ac_sum %>%
  ggplot(aes(x = -depth, y=mean, color=diel, group=diel), alpha=0.4)+
  geom_line(alpha=0.8)+
  geom_errorbar(aes(ymin=mean-SD, ymax =mean+SD), width=.4,alpha=0.8)+
  geom_point(aes(shape=diel), alpha=0.7, size=6)+
  #coord_cartesian(ylim = c(100, -700), expand = FALSE)+
  scale_color_manual(values=paletteDayNight )+
  scale_fill_manual(values=paletteDayNight)+
  theme_minimal() +
  coord_flip(xlim = c(-660, 60), expand = FALSE)+
  scale_x_continuous(n.breaks = 6)+
  #scale_x_continuous(breaks=x_ticks_hr) +
  xlab('Depth (m)')+
  ylab('Mean Scattering (dB)')+
  labs(fill='Diel', color='Diel', shape='Diel')

p1
filename = paste(plot_directory, 'Acoustic_DayNight_point.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 6, width =5, units = 'in')
filename = paste(plot_directory, 'Acoustic_DayNight_point.svg', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 6, width =5, units = 'in')

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
  #more detailed depth bins
  mutate(depth_bin = case_when(depth <=25 ~ "00_0-25m",
                               depth >25 & depth <=50 ~ "01_25-50m",
                               depth >50 & depth <=75 ~ "02_50-75m",
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

# Also include different depth bin

meta %<>% mutate(depth_bin2 = case_when(depth <=100 ~ "0-100m",
                                        depth >100 & depth <=200 ~ "100-200m",
                                        depth >200 & depth <=300 ~ "200-300m",
                                        depth >300 & depth <=400 ~ "300-400m",
                                        depth >400 & depth <=500 ~ "400-500m",
                                        depth >400 & depth <=600 ~ "500-600m",
                                        depth >600 & depth <=750 ~ "600-700m", TRUE ~ "unknown"
)) 


# Create Diel Boxplot of Mesopelagic and Epipelagic EcolCats --------------

# For the following analysis, take average values of replicate sequenced filters
# 245 unique samples
# 207 from day or night periods
meta %>% distinct(FilterID, .keep_all = TRUE) %>% filter(diel %in% c('day', 'night')) %>% glimpse()

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  left_join(meta %>% select(SampleID, FilterID)) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot, FilterID) %>%
  # average by unique filter
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic'))

stats <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
  mutate(count=1) %>%
  group_by(diel, Ecological_Category, depth_bin2) %>%
  mutate(median = median(sum_per_tot)) %>%
  mutate(max = max(sum_per_tot)) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(diel, Ecological_Category, depth_bin2, median, sum_count, max) %>%
  # make values for segments, end is current median, start is median from previous depth
  mutate(xend = NaN) %>% # will hold depth bin values
  #mutate(yend = NaN) %>% # will hold per tot values
  mutate(xend = case_when(depth_bin2 == "0-100m" ~ "0-100m",
                          depth_bin2 == "100-200m" ~ "0-100m",
                          depth_bin2 == "200-300m" ~ "100-200m",
                          depth_bin2 == "300-400m" ~ "200-300m",
                          depth_bin2 == "400-500m" ~ "300-400m",
                          depth_bin2 == "500-600m" ~ "400-500m",
                          depth_bin2 == "600-700m" ~ "500-600m", TRUE ~ "unknown"))
#make yend df
# xend as index, merge back in to get median value
stats2 <- stats %>%
  left_join(stats %>% select(depth_bin2,median, Ecological_Category, diel) %>% rename(yend=median) %>% rename(xend=depth_bin2))

# assign text colour
textcol <- "grey40"
print("Begin plotting...")

p2 <- left_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  #filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
  #add in stats df
  full_join(stats2) %>%
  #now plot
  ggplot(aes(x=Ecological_Category, y=sum_per_tot, fill=diel)) +
  geom_boxplot(aes(x=fct_rev(depth_bin2)), alpha=0.5)+
  geom_point(aes(y=median, x=depth_bin2, color=diel)) +
  # add in count of number of samples?
  geom_text(data=.%>%filter(Ecological_Category=='mesopelagic') %>% distinct(depth_bin2,diel, .keep_all = TRUE), aes(label=sum_count, y=106, x=depth_bin2, color=diel), position = position_dodge(width = 0.9))+
  #geom_text(data=. %>% filter(diel=='night') %>%filter(Ecological_Category=='mesopelagic'),
  #          aes(label=sum_count, y=102, x=depth_bin2, color=diel), position = position_dodge(width = 0.9)) +
  #geom_text(data=. %>% filter(diel=='day') %>%filter(Ecological_Category=='mesopelagic'),
  #          aes(label=sum_count, y=112, x=depth_bin2, color=diel), position = position_dodge(width = 0.9)) +
  geom_segment(aes(y=median,yend=yend, x=depth_bin2,xend=xend, color=diel)) +
  facet_grid(. ~ Ecological_Category, scales = "free") +
  #formatting...
  scale_color_manual(values=paletteDayNight )+
  scale_fill_manual(values=paletteDayNight)+
  #scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = -1)+
  #scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = -1)+
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
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8,colour=textcol),
    axis.text.y=element_text(size=8,colour=textcol),
    #axis.title.y = element_text(size=6),
    plot.background=element_blank(),
    #panel.border=element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())

p2

filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')



# Join plots --------------

#Put all together
pf1 <- p1  + theme_minimal() + theme(text = element_text(size=10), 
                                    legend.position = 'none',
                                    axis.line = element_line(color = "black"),
                                    axis.ticks = element_line(color = "black"),
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    axis.title.y=element_blank(),
                                    axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()
)
pf2 <- p2  + theme_minimal() + theme(text = element_text(size=10), 
                                    legend.position = 'none',
                                    axis.line = element_line(color = "black"),
                                    axis.ticks = element_line(color = "black"),
                                    #axis.title.x=element_blank(),
                                    #axis.text.x=element_blank(),
                                    #axis.ticks.x=element_blank(),
                                    #axis.title.y=element_blank(),
                                    #axis.text.y=element_blank(),
                                    #axis.ticks.y=element_blank(),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank()
)


plot_grid(pf2, pf1,ncol = 2, align = "h", axis = "bt", rel_widths = c(1, 0.5))

filename = paste(plot_directory, 'Acoustic_EpiMeso_Boxplot.png', sep='')
filename
ggsave(filename, width=8, height=6, units = 'in')
filename = paste(plot_directory, 'Acoustic_EpiMeso_Boxplot.svg', sep='')
filename
ggsave(filename, width=8, height=6, units = 'in')





# Exploration Below --------




# # swarm/jitter plot
# 
# bp_top <- full_join(df, meta,  by = c("SampleID")) %>% #join with metadata
#   mutate(hour = as.integer(hour)) %>%
#   mutate(hour = replace(hour, hour==24,0)) %>%
#   filter(diel %in% c('day', 'night')) %>%
#   filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
#   select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
#   #add in stats df
#   full_join(stats2) %>%
#   #now plot
#   ggplot(aes(x=Ecological_Category, y=sum_per_tot, fill=diel)) +
#   geom_boxplot(aes(x=fct_rev(depth_bin2)), alpha=0.3)+
#   geom_point(aes(x=fct_rev(depth_bin2), color=diel),size=0.4, alpha=0.7, position=position_jitterdodge()) +
#   #geom_jitter(aes(x=fct_rev(depth_bin2), color=diel),size=0.4, alpha=0.9) +
#   geom_point(aes(y=median, x=depth_bin2, color=diel)) +
#   # add in count of number of samples?
#   geom_text(data=. %>% filter(diel=='night') %>%filter(Ecological_Category=='mesopelagic'),aes(label=sum_count, y=102, x=depth_bin2, color=diel)) +
#   geom_text(data=. %>% filter(diel=='day') %>%filter(Ecological_Category=='mesopelagic'),aes(label=sum_count, y=112, x=depth_bin2, color=diel)) +
#   geom_segment(aes(y=median,yend=yend, x=depth_bin2,xend=xend, color=diel)) +
#   facet_grid(. ~ Ecological_Category, scales = "free") +
#   #formatting...
#   scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = -1)+
#   scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = -1)+
#   coord_flip() +
#   theme_minimal() +
#   guides(fill=guide_legend(ncol=2)) +
#   labs(y="Percent Total Reads",x="Depth Bin")+
#   theme(
#     #legend
#     legend.position="right",legend.direction="vertical",
#     legend.text=element_text(colour=textcol,size=8,face="bold"),
#     legend.key.height=grid::unit(0.3,"cm"),
#     legend.key.width=grid::unit(0.3,"cm"),
#     legend.title=element_text(colour=textcol,size=8,face="bold"),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8,colour=textcol),
#     axis.text.y=element_text(size=8,colour=textcol),
#     #axis.title.y = element_text(size=6),
#     plot.background=element_blank(),
#     panel.border=element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(size = .25),
#     plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
#     plot.title=element_blank())
# 
# bp_top
# 
# filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2_jitter.png', sep='')
# print(filename)
# ggsave(filename,height = 5, width =8, units = 'in')
# filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2_jitter.svg', sep='')
# print(filename)
# ggsave(filename,height = 5, width =8, units = 'in')








# Look at replicate sequenced filters ------

dup_samps <- meta %>% group_by(FilterID) %>% filter(n()>1) %>% ungroup()

#number of samples by diel group:
test <- meta %>% distinct(FilterID, .keep_all = TRUE) %>% group_by(diel) %>%
  mutate(count = n())%>%
  ungroup() %>%
  distinct(diel, count)

# #deep samples by diel group:
# test <- meta %>% filter(depth_bin =='09_400-500m')
# test <- meta %>% filter(depth_bin =='08_300-400m')
# test <- meta %>% filter(depth_bin =='11_600-750m')

# plot
df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot)

dup_samps %>% select(SampleID, FilterID, local_time, depth) %>%
  left_join(df) %>%
  # unite(label,depth,local_time,FilterID,SampleID,sep='_' ) %>%
  ggplot(aes(x=SampleID, y=sum_per_tot))+
  #geom_bar(aes( y=0.5),stat='identity', fill = "grey",alpha=0.8, width=20)+
  geom_bar(stat='identity', aes(fill = Ecological_Category))+
  coord_flip()+
  #scale_x_reverse()+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  labs(x='',y="Percent Total Reads", title='Replicate Filters')+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="bottom",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=8,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=6,colour=textcol),
    axis.title.y = element_text(size=6),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())
  
filename = paste(plot_directory, marker,'_EcolCat_replicate_ESP_comp.png', sep='')
print(filename)
ggsave(filename,height = 10, width =10, units = 'in')
filename = paste(plot_directory, marker,'_EcolCat_replicate_ESP_comp.svg', sep='')
print(filename)
ggsave(filename,height = 10, width =10, units = 'in')

# values as points on top of each other

dup_samps %>% select(SampleID, FilterID, local_time, depth) %>%
  left_join(df) %>%
  # unite(label,depth,local_time,FilterID,SampleID,sep='_' ) %>%
  ggplot(aes(x=FilterID, y=sum_per_tot,color= Ecological_Category))+
  #geom_bar(aes( y=0.5),stat='identity', fill = "grey",alpha=0.8, width=20)+
  #geom_bar(stat='identity', aes(fill = Ecological_Category))+
  geom_point()+
  geom_line()+
  coord_flip()+
  #scale_x_reverse()+
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
  labs(x='',y="Percent Total Reads", title='Replicate Filters')+
  theme_minimal() +
  guides(fill=guide_legend(ncol=2)) +
  theme(
    #legend
    legend.position="bottom",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=8,face="bold"),
    legend.key.height=grid::unit(0.3,"cm"),
    legend.key.width=grid::unit(0.3,"cm"),
    legend.title=element_text(colour=textcol,size=8,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=6,colour=textcol),
    axis.title.y = element_text(size=6),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())

filename = paste(plot_directory, marker,'_EcolCat_replicate_ESP_points.png', sep='')
print(filename)
ggsave(filename,height = 10, width =10, units = 'in')
filename = paste(plot_directory, marker,'_EcolCat_replicate_ESP_points.svg', sep='')
print(filename)
ggsave(filename,height = 10, width =10, units = 'in')

# Boxplot of Ecological Categories through depth ------------

# For the following analysis, take average values of replicate sequenced filters

df <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  left_join(meta %>% select(SampleID, FilterID)) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot, FilterID) %>%
  # average by unique filter
  group_by(Ecological_Category, FilterID) %>%
  mutate(sum_reads = mean(sum_reads)) %>%
  mutate(sum_per_tot = mean(sum_per_tot)) %>%
  ungroup() %>%
  distinct(FilterID, Ecological_Category, .keep_all=TRUE)

#want to plot median value on top of bar plot, show number of samples in bin

stats <- full_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin, diel) %>%
  # drop duplicates
  distinct(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin, diel, .keep_all = TRUE) %>%
  mutate(count=1) %>%
  group_by(diel, Ecological_Category, depth_bin) %>%
  mutate(median = median(sum_per_tot)) %>%
  mutate(max = max(sum_per_tot)) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(diel, Ecological_Category, depth_bin, median, sum_count, max) %>%
  # make values for segments, end is current median, start is median from previous depth
  mutate(xend = NaN) %>% # will hold depth bin values
  mutate(xend = case_when(depth_bin == "00_0-25m" ~ "00_0-25m",
                          depth_bin == "01_25-50m" ~ "00_0-25m",
                          depth_bin == "02_50-75m" ~ "01_25-50m",
                          depth_bin == "03_75-100m" ~ "02_50-75m",
                          depth_bin == "04_100-150m" ~ "03_75-100m",
                          depth_bin == "05_150-200m" ~ "04_100-150m",
                          depth_bin == "06_200-250m" ~ "05_150-200m",
                          depth_bin == "07_250-300m" ~ "06_200-250m",
                          depth_bin == "08_300-400m" ~ "07_250-300m",
                          depth_bin == "09_400-500m" ~ "08_300-400m",
                          depth_bin == "10_500-600m" ~ "09_400-500m",
                          depth_bin == "11_600-750m" ~ "10_500-600m", TRUE ~ "unknown"))
# make yend df
# xend as index, merge back in to get median value
stats2 <- stats %>%
  left_join(stats %>% select(depth_bin,median, Ecological_Category, diel) %>% rename(yend=median) %>% rename(xend=depth_bin))

# assign text colour
textcol <- "grey40"
print("Begin plotting...")

bp_top <- full_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  # drop duplicates
  distinct(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin, diel, .keep_all = TRUE) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin, diel) %>%
  #add in stats df
  full_join(stats2) %>%
  #now plot
  ggplot(aes(x=Ecological_Category, y=sum_per_tot, fill=diel)) +
  geom_boxplot(aes(x=fct_rev(depth_bin)), alpha=0.5)+
  geom_point(aes(y=median, x=depth_bin, color=diel)) +
  geom_segment(aes(y=median,yend=yend, x=depth_bin,xend=xend, color=diel)) +
  # add in count of number of samples?
  geom_text(data=. %>% filter(diel=='night') %>%filter(Ecological_Category=='mesopelagic'),aes(label=sum_count, y=102, x=depth_bin, color=diel)) +
  geom_text(data=. %>% filter(diel=='day') %>%filter(Ecological_Category=='mesopelagic'),aes(label=sum_count, y=112, x=depth_bin, color=diel)) +
  facet_grid(. ~ Ecological_Category, scales = "free") +
  #formatting...
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

filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')

# Alternative depth bin boxplot -------------

#want to plot median value on top of bar plot, show number of samples in bin

stats <- full_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
  mutate(count=1) %>%
  group_by(diel, Ecological_Category, depth_bin2) %>%
  mutate(median = median(sum_per_tot)) %>%
  mutate(max = max(sum_per_tot)) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(diel, Ecological_Category, depth_bin2, median, sum_count, max) %>%
  # make values for segments, end is current median, start is median from previous depth
  mutate(xend = NaN) %>% # will hold depth bin values
  #mutate(yend = NaN) %>% # will hold per tot values
  mutate(xend = case_when(depth_bin2 == "0-100m" ~ "0-100m",
                          depth_bin2 == "100-200m" ~ "0-100m",
                          depth_bin2 == "200-300m" ~ "100-200m",
                          depth_bin2 == "300-400m" ~ "200-300m",
                          depth_bin2 == "400-500m" ~ "300-400m",
                          depth_bin2 == "500-600m" ~ "400-500m",
                          depth_bin2 == "600-700m" ~ "500-600m", TRUE ~ "unknown"))
#make yend df
# xend as index, merge back in to get median value
stats2 <- stats %>%
  left_join(stats %>% select(depth_bin2,median, Ecological_Category, diel) %>% rename(yend=median) %>% rename(xend=depth_bin2))

# assign text colour
textcol <- "grey40"
print("Begin plotting...")

bp_top <- full_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
  #add in stats df
  full_join(stats2) %>%
  #now plot
  ggplot(aes(x=Ecological_Category, y=sum_per_tot, fill=diel)) +
  geom_boxplot(aes(x=fct_rev(depth_bin2)), alpha=0.3)+
  geom_point(aes(y=median, x=depth_bin2, color=diel)) +
  # add in count of number of samples?
  geom_text(data=. %>% filter(diel=='night') %>%filter(Ecological_Category=='mesopelagic'),aes(label=sum_count, y=102, x=depth_bin2, color=diel)) +
  geom_text(data=. %>% filter(diel=='day') %>%filter(Ecological_Category=='mesopelagic'),aes(label=sum_count, y=112, x=depth_bin2, color=diel)) +
  geom_segment(aes(y=median,yend=yend, x=depth_bin2,xend=xend, color=diel)) +
  facet_grid(. ~ Ecological_Category, scales = "free") +
  #formatting...
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
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8,colour=textcol),
    axis.text.y=element_text(size=8,colour=textcol),
    #axis.title.y = element_text(size=6),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())

bp_top

filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')


# swarm/jitter plot

bp_top <- full_join(df, meta,  by = c("SampleID")) %>% #join with metadata
  mutate(hour = as.integer(hour)) %>%
  mutate(hour = replace(hour, hour==24,0)) %>%
  filter(diel %in% c('day', 'night')) %>%
  filter(Ecological_Category %in% c('mesopelagic', 'epipelagic')) %>%
  select(time, depth, sum_per_tot, Ecological_Category, hour, depth_bin2, diel) %>%
  #add in stats df
  full_join(stats2) %>%
  #now plot
  ggplot(aes(x=Ecological_Category, y=sum_per_tot, fill=diel)) +
  geom_boxplot(aes(x=fct_rev(depth_bin2)), alpha=0.3)+
  geom_point(aes(x=fct_rev(depth_bin2), color=diel),size=0.4, alpha=0.7, position=position_jitterdodge()) +
  #geom_jitter(aes(x=fct_rev(depth_bin2), color=diel),size=0.4, alpha=0.9) +
  geom_point(aes(y=median, x=depth_bin2, color=diel)) +
  # add in count of number of samples?
  geom_text(data=. %>% filter(diel=='night') %>%filter(Ecological_Category=='mesopelagic'),aes(label=sum_count, y=102, x=depth_bin2, color=diel)) +
  geom_text(data=. %>% filter(diel=='day') %>%filter(Ecological_Category=='mesopelagic'),aes(label=sum_count, y=112, x=depth_bin2, color=diel)) +
  geom_segment(aes(y=median,yend=yend, x=depth_bin2,xend=xend, color=diel)) +
  facet_grid(. ~ Ecological_Category, scales = "free") +
  #formatting...
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
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8,colour=textcol),
    axis.text.y=element_text(size=8,colour=textcol),
    #axis.title.y = element_text(size=6),
    plot.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())

bp_top

filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2_jitter.png', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')
filename = paste(plot_directory, marker,'_EcolCat_boxplot_through_depth2_jitter.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =8, units = 'in')





