#Aitchenson_distance_plot
#Katie Pitz
#040721

# Set directory to save plots
plot_dir <- 'figures/Aitchenson_distance/'
results_directory <- 'Qiime_Results/'


# Load Libraries -----------------------------------------------------------------
library(qiime2R)
library(tidyverse)
library(lubridate) #for date modifications
library(vegan)
library(ggthemes)
library(magrittr)
#for heatmap section
library(gridExtra)
library(RColorBrewer) #colors for plotting
library(forcats) #working with factor data (reordering)
library(patchwork) #putting graphs together on same plot

library(ggdendro)

marker = sym("12S")

# Import Data -------------------------------------------------------------
data_directory = "data/filtered_seq_data/"

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


# Modify Metadata ---------------------------------------------------------
library(magrittr)
library(lubridate)
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


# Import DEICODE Results --------------------------------------------------

library(qiime2R)

#Import Qiime2 Results
file = paste(results_directory,"ordination.qza",sep="")
print(file)
pco<-read_qza(file)
pco$uuid
#look at data
head(pco$data$ProportionExplained)
pco$data$Vectors[1:5, 1:4]

#create proportion explained labels
label.PC1 <- paste("PC1: ", round(pco$data$ProportionExplained$PC1, 3)*100,"%")
label.PC1
label.PC2 <- paste("PC2: ", round(pco$data$ProportionExplained$PC2, 3)*100,"%")
label.PC2
label.PC3 <- paste("PC3: ", round(pco$data$ProportionExplained$PC3, 3)*100,"%")
label.PC3

#Join with sample data
pcscores <- left_join(pco$data$Vectors, meta, by="SampleID")

#format loading scores
loadings <- as.data.frame(pco$data$Species)
loadings$ASV <- loadings$FeatureID

#join on OTU, adding taxa info
loadings <- left_join(loadings, tax.c, by="ASV")

#export pcscores
file = paste(results_directory, "pcscores_",marker,"_Dada2_Qiime2.csv",sep="")
print(file)
write.csv(pcscores, file)
file = paste(results_directory, "loadings_",marker,"_Dada2_Qiime2.csv",sep="")
print(file)
write.csv(loadings, file)
pcscores[1:5, 1:9]  #long because of sample data
head(loadings)

# Ait distance:

file = paste(results_directory,"distance.qza",sep="")
print(file)
ait <-read_qza(file)
ait
aitmat <- as.matrix(ait$data)
aitdist <- as.dist(aitmat)



# Distance between depths -------------------------------------------------

#transform aitchenson distance matrix into format for 
#plotting distance between depths
# get into long format to get distance from sample to sample

ait_tib <-as_tibble(aitmat,rownames = "SampleID") %>%
  pivot_longer(-SampleID,names_to = "SampleID2", values_to = "distance")

# add in metadata for both 1 and 2:
lim_meta <- select(meta,c('SampleID','depth_bin', 'depth', 'diel', 'day', 'time'))
ait_tib %<>% left_join(lim_meta, by="SampleID")
ait_tib %<>% left_join(lim_meta %>% mutate(SampleID2 = SampleID), by="SampleID2", suffix= (c("1","2"))) %>%
  select(-SampleID22) %>%
  mutate(interaction_type = paste(diel1,'-',diel2),sep='')


#plot

#Number of samples in depth bins:
p <- lim_meta %>%
  mutate(count=1) %>%
  filter(diel !='transition') %>%
  group_by(depth_bin) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  distinct(depth_bin,.keep_all=TRUE) %>%
  ggplot(aes(x= fct_reorder(depth_bin, depth), y=sum_count)) +
  geom_bar(stat='identity') +
  labs(y='Number of Samples', x='Depth bin')
filename = paste(plot_dir, 'Aitdist_',marker,'_numberofsamples.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')
p

# distribution of distance in interaction types between depths
p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='0-50m') %>%
  ggplot(aes(x=fct_reorder(depth_bin2, depth2), y=distance)) +
  geom_boxplot(aes(color=interaction_type))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

filename = paste(plot_dir, 'Aitdist_',marker,'_surface_to_alldepths.png', sep='')
filename
ggsave(filename,height = 6, width =12, units = 'in')
p

p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='200-300m') %>%
  ggplot(aes(x=fct_reorder(depth_bin2, depth2), y=distance)) +
  geom_boxplot(aes(color=interaction_type))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
filename = paste(plot_dir, 'Aitdist_',marker,'_200300_to_alldepths.png', sep='')
filename
ggsave(filename,height = 6, width =12, units = 'in')
p

p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='100-200m') %>%
  ggplot(aes(x=fct_reorder(depth_bin2, depth2), y=distance)) +
  geom_boxplot(aes(color=interaction_type))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
filename = paste(plot_dir, 'Aitdist_',marker,'_100200_to_alldepths.png', sep='')
filename
ggsave(filename,height = 6, width =12, units = 'in')
p

# test if distribution is the same of distance values
data1 <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='0-50m') %>%
  filter(depth_bin2=='300-400m') %>%
  filter(interaction_type == 'day - day') %>%
  pull(distance)
data2 <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='200-300m') %>%
  filter(interaction_type == 'day - night') %>%
  pull(distance)
ks.test(data1, data2)

p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='0-50m') %>%
  filter(depth_bin2=='0-50m') %>%
  ggplot(aes(x=diel1, y=distance)) +
  geom_boxplot(aes(color=diel2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

filename = paste(plot_dir, 'Aitdist_',marker,'_allsurface_bydiel.png', sep='')
filename
ggsave(filename,height = 8, width =8, units = 'in')

p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='0-50m') %>%
  filter(depth_bin2=='0-50m') %>%
  mutate(day1_char = as.character(day1)) %>%
  mutate(day2_char = as.character(day2)) %>%
  ggplot(aes(x=fct_reorder(day1_char, time1), y=distance)) +
  geom_boxplot(aes(color=fct_reorder(day2_char, time2)))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

filename = paste(plot_dir, 'Aitdist_',marker,'_allsurface_bytime.png', sep='')
filename
ggsave(filename,height = 8, width =8, units = 'in')
p



p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='0-50m') %>%
  #filter(diel1=='day' ) %>%
  #ggplot(aes(x=fct_reorder(depth_bin, depth), y=distance)) +
  ggplot(aes(x=fct_reorder(depth_bin2, depth2), y=distance)) +
  #geom_boxplot(aes(color=diel2))+
  geom_boxplot(aes(color=diel2))+
  #geom_point(alpha=0.1)+
  facet_wrap(~diel1, ncol=2, nrow=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
filename = paste(plot_dir, 'Aitdist_',marker,'_surface1_facetdiel1.png', sep='')
filename
ggsave(filename,height = 8, width =12, units = 'in')
p


p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='0-50m') %>%
  #filter(diel1=='day' ) %>%
  #ggplot(aes(x=fct_reorder(depth_bin, depth), y=distance)) +
  ggplot(aes(x=fct_reorder(depth_bin2, depth2), y=distance)) +
  #geom_boxplot(aes(color=diel2))+
  geom_boxplot(aes(color=diel2))+
  #geom_point(alpha=0.1)+
  #facet_wrap(~diel1, ncol=2, nrow=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
filename = paste(plot_dir, 'Aitdist_',marker,'_surface1.png', sep='')
filename
ggsave(filename,height = 8, width =8, units = 'in')
p




p <- ait_tib %>%
  filter(depth1>=0) %>%
  filter(depth2>=0) %>%
  filter(depth_bin1=='0-50m') %>%
  filter(diel1=='day' ) %>%
  #ggplot(aes(x=fct_reorder(depth_bin, depth), y=distance)) +
  ggplot(aes(x=fct_reorder(depth_bin2, depth2), y=distance)) +
  #geom_boxplot(aes(color=diel2))+
  geom_boxplot(aes(color=diel2))+
  #geom_point(alpha=0.1)+
  #facet_wrap(~diel1, ncol=2, nrow=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='300-400m') %>%
  #filter(diel1=='day' ) %>%
  #ggplot(aes(x=fct_reorder(depth_bin, depth), y=distance)) +
  ggplot(aes(x=fct_reorder(depth_bin2, depth2), y=distance)) +
  #geom_boxplot(aes(color=diel2))+
  geom_boxplot(aes(color=diel2))+
  #geom_point(alpha=0.1)+
  facet_wrap(~diel1, ncol=2, nrow=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
filename = paste(plot_dir, 'Aitdist_',marker,'_300400day_1_facetdiel1.png', sep='')
filename
ggsave(filename,height = 8, width =12, units = 'in')
p

p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='300-400m') %>%
  #filter(diel1=='day' ) %>%
  #ggplot(aes(x=fct_reorder(depth_bin, depth), y=distance)) +
  ggplot(aes(x=fct_reorder(depth_bin2, depth2), y=distance)) +
  #geom_boxplot(aes(color=diel2))+
  geom_boxplot(aes(color=diel1))+
  #geom_point(alpha=0.1)+
  #facet_wrap(~diel1, ncol=2, nrow=3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
filename = paste(plot_dir, 'Aitdist_',marker,'_300400day_1_colordiel1.png', sep='')
filename
ggsave(filename,height = 8, width =12, units = 'in')
p

#calculate mean distances
mean1 <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='300-400m') %>%
  filter(depth_bin2=='0-50m') %>%
  #filter(diel1=='day') %>%
  group_by(diel1) %>%
  mutate(mean_dist = mean(distance)) %>%
  ungroup() %>%
  distinct(diel1,mean_dist)

p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='300-400m') %>%
  filter(depth_bin2=='0-50m') %>%
  #filter(diel1=='day' ) %>%
  #ggplot(aes(x=fct_reorder(depth_bin, depth), y=distance)) +
  ggplot(aes(x=distance, color=diel1, fill=diel1)) +
  geom_freqpoly(size=2, alpha=0.8)+
  geom_vline(data=mean1, aes(xintercept=mean_dist, color=diel1), linetype = "longdash")+
  labs(color='300-400m', x='Distance to surface 0-50m samples (all)')
  #geom_vline(xintercept=2.06)

filename = paste(plot_dir, 'Aitdist_',marker,'_300400day_tosurface_all.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')
p

p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='200-300m') %>%
  filter(depth_bin2=='0-50m') %>%
  #filter(diel1=='day' ) %>%
  #ggplot(aes(x=fct_reorder(depth_bin, depth), y=distance)) +
  ggplot(aes(x=distance, color=diel1, fill=diel1)) +
  geom_freqpoly(size=2, alpha=0.8)+
  geom_vline(data=mean1, aes(xintercept=mean_dist, color=diel1), linetype = "longdash")+
  labs(color='200-300m', x='Distance to surface 0-50m samples (all)')
#geom_vline(xintercept=2.06)

filename = paste(plot_dir, 'Aitdist_',marker,'_200300day_tosurface_all.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')
p

p <- ait_tib %>%
  drop_na() %>%
  filter(diel1 !='transition') %>%
  filter(diel2 !='transition') %>%
  filter(depth_bin1=='400-500m') %>%
  filter(depth_bin2=='0-50m') %>%
  #filter(diel1=='day' ) %>%
  #ggplot(aes(x=fct_reorder(depth_bin, depth), y=distance)) +
  ggplot(aes(x=distance, color=diel1, fill=diel1)) +
  geom_freqpoly(size=2, alpha=0.8)+
  geom_vline(data=mean1, aes(xintercept=mean_dist, color=diel1), linetype = "longdash")+
  labs(color='400-500m', x='Distance to surface 0-50m samples (all)')
#geom_vline(xintercept=2.06)

filename = paste(plot_dir, 'Aitdist_',marker,'_400500day_tosurface_all.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')
p

# Plot Organisms in layers at day or night ----------------------------------------
taxa_level = sym('Family')

top_taxa <- potu.c %>%
  inner_join(meta %>% select(SampleID, depth_bin),  by = c("SampleID")) %>%
  filter(depth_bin == '50-100m') %>%
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
bp_top <- inner_join(potu.c, meta,  by = c("SampleID")) %>% #join with metadata
  inner_join(species_label,  by = c("ASV")) %>%  #join with taxonomy
  right_join(top_taxa) %>% #limit to top taxa
  filter(depth_bin == '50-100m') %>%
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
bp_top


# Cluster -----------------------------------------------------------------

aitclust <- hclust(aitdist, method= "ward.D2")
plot(aitclust)
aitdend <- as.dendrogram(aitclust, hang = -1, lwd = 3, lty = 3, sub = "")
#dend_cols <- select(samp.c,c('SampleID','date'))
plot(aitdend)


#Can plot with ggplot2

dendr <- dendro_data(aitclust, type="rectangle") 

#Define cluster based on similarity percent - good for plotting
clust <- cutree(aitclust, h = 6)               # find 'cut' clusters (k) or choose level (h) numeric scalar or vector with heights where the tree should be cut.
clust2 <- cutree(aitclust, h = 7)
clust3 <- cutree(aitclust, h = 10)  
clust4 <- cutree(aitclust, h = 16)  

clust.df <- data.frame(label = names(clust), cluster_6 = clust, cluster_7 = clust2, cluster_10 = clust3, cluster_16 = clust4)
sapply(clust.df, mode)
clust.df$label <- as.character(clust.df$label)
clust.df <- clust.df[order(clust.df$label),]

#join sample data with cluster df and simprofCLUSTERS df:
tree_scores <- as_tibble(clust.df) %>% rename(SampleID = label) %>% left_join(pcscores, by='SampleID')
tree_data <- as_tibble(dendr$labels) %>% 
  rename(SampleID = label) %>% 
  left_join(tree_scores, by='SampleID') %>% 
  arrange(depth_bin)

# In PCA space:
tree_data %>%ggplot(aes(x= PC1, y=PC2, color=cluster_16))+
  geom_point()
tree_data %>%ggplot(aes(x= PC1, y=PC2, color=cluster_10))+
  geom_point()
tree_data %>%ggplot(aes(x= PC1, y=PC2, color=cluster_7))+
  geom_point()

#By depth_bin
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  #geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=year_char),alpha=1, size=1)+
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=fct_reorder(depth_bin,depth)),alpha=1, size=1)+
  geom_text(data=tree_data, aes(x=x, y=y, label=diel, hjust=-.4), size=1) +
  geom_text(data=tree_data, aes(x=x, y=y, label=PlateID, hjust=-6), size=2) +
  coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  scale_color_brewer(palette='Set1')+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

filename = paste(plot_dir, 'Aitdist_',marker,'_hclust_depth.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')




#Nice plots for paper
#month
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=fct_reorder(depth_bin,depth)),alpha=1, size=1)+
  coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  scale_color_brewer(palette='Paired')+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

filename = paste(plot_dir, 'Aitdist_',marker,'_hclust_depth_lim.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')

#type (ESP)
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=ESP),alpha=1, size=1)+
  coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  scale_color_brewer(palette='Paired')+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

filename = paste(plot_dir, 'Aitdist_',marker,'_hclust_ESP_lim.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')

#PlateID
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=PlateID),alpha=1, size=1)+
  coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  scale_color_brewer(palette='Paired')+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

filename = paste(plot_dir, 'Aitdist_',marker,'_hclust_PlateID_lim.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')

#Day?
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=day),alpha=1, size=1)+
  geom_text(data=tree_data, aes(x=x, y=y, label=SampleID, hjust=-.4), size=1) +
  coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  #scale_fill_viridis() +
  #scale_color_brewer(palette='Paired')+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

filename = paste(plot_dir, 'Aitdist_',marker,'_hclust_time_lim.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')

#Diel
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=diel),alpha=1, size=1)+
  geom_text(data=tree_data, aes(x=x, y=y, label=depth, hjust=-.4), size=1) +
  coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  #scale_fill_viridis() +
  scale_color_brewer(palette='Paired')+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

filename = paste(plot_dir, 'Aitdist_',marker,'_hclust_diel_lim.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')

#By Year
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=fct_reorder(year_char,year)),alpha=1, size=1)+
  coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  scale_colour_manual(values= c( "#6baed6", "#3182bd", "#08519c", 
                                 "#a1d99b", "#74c476", "#31a354", "#006d2c", 
                                 "#fdae6b", "#fd8d3c", "#e6550d", "#a63603","#9ecae1"))+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

filename = paste(directory, 'Aitdist_',marker,'_hclust_year_lim.svg', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')
filename = paste(directory, 'Aitdist_',marker,'_hclust_year_lim.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')
#By Year/BW
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=year),alpha=1, size=1)+
  coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

filename = paste(directory, 'Aitdist_',marker,'_hclust_year_lim_scale.svg', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')
filename = paste(directory, 'Aitdist_',marker,'_hclust_year_lim_scale.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')
#By PC score
ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=PC1),alpha=1, size=1)+
  geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-2.5,ymax=y-3.5, color=PC2),alpha=1, size=1)+
  #geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-4,ymax=y-5, color=PC3),alpha=1, size=1)+
  coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  scale_colour_gradient2()+
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

filename = paste(directory, 'Aitdist_',marker,'_hclust_PCs_lim_scale.svg', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')
filename = paste(directory, 'Aitdist_',marker,'_hclust_PCs_lim_scale.png', sep='')
filename
ggsave(filename,height = 4, width =6, units = 'in')



# SCRAP -----

# Start Loop over Markers -------------------------------------------------


# Loop over script for each maker
pcs <- c("PC1", "PC2", "PC3")
markers <- c("12S","COI","18S","12SnoCL", "COInoH")
#marker = sym("18S")
for (val in markers) {
  marker = sym(val)
  
  # Import Data -------------------------------------------------------------
  data_directory = "Data/Filtered_Seq_Data/"
  #ASV table
  print('ASV table')
  file = paste("A2W_",marker,"_otu_Filtered.csv", sep='')
  filepath = paste(data_directory, file, sep='')
  print(filepath)
  otu.c <- read_csv(filepath) %>% rename('ASV' = 'X1')
  
  #taxa table
  print('taxa table')
  file = paste("A2W_",marker,"_taxa_Filtered.csv", sep='')
  filepath = paste(data_directory, file, sep='')
  print(filepath)
  tax.c <- read_csv(filepath) %>% rename('ASV' = 'X1')
  
  #metadata table
  print('metadata table')
  #file = "A2W_12S_meta_Filtered.csv"
  file = paste("A2W_",marker,"_meta_Filtered.csv", sep='')
  filepath = paste(data_directory, file, sep='')
  print(filepath)
  samp.c <- read_csv(filepath) %>% rename('SampleID' = 'sample_name')
  
  #OTU table long format with percent total reads
  potu.c <- otu.c %>%
    tidyr::pivot_longer( -ASV, names_to ='SampleID',values_to = 'reads' ) %>%
    group_by(SampleID) %>%
    mutate(per_tot = reads / sum(reads) *100) %>%
    arrange(-reads)
  head(potu.c)
  
  # Create Seasonal Variables -----------------------------------------------
  
  #edit sample metadata to get date columns
  samp.c %<>% mutate(DATE_TIME = as_datetime(SAMPLING_date_time)) %>%
    # extract date components
    mutate(month =  month(DATE_TIME)) %>%
    mutate(day =  day(DATE_TIME)) %>%
    mutate(year =  year(DATE_TIME)) %>%
    mutate(jday = yday(DATE_TIME)) %>%
    mutate(month_2 = format.Date(DATE_TIME, "%m")) %>%  #two digit month
    #mutate(ym = ym(DATE_TIME)) %>%   #lubridate will be adding this function soon; have to work around it for now.
    #change column data types
    mutate(month_char = as.character(month)) %>%
    mutate(year_char = as.character(year)) %>%
    mutate(month_2 = as.character(month_2)) %>%
    #year/month for labeling
    unite(ym, year_char, month_2, sep='/', remove=FALSE) %>%
    #create seasonal categories
    # 12,1,2; 3,4,5; 6,7,8; 9,10,11
    mutate(season = "Dec-Feb") %>% #create new column
    mutate(season = case_when(month<=11 & month >=9  ~"Sep-Nov",
                              TRUE ~ season)) %>%
    mutate(season = case_when(month<=8 & month >=6  ~"June-Aug",
                              TRUE ~ season)) %>%
    mutate(season = case_when(month<=5 & month >=3  ~"March-May",
                              TRUE ~ season)) %>%
    #create "winters" so can group late fall and early spring of same winter
    mutate(winter = "Not winter") %>% #create new column
    mutate(year_before = year-1) %>%
    mutate(year_before = as.character(year_before)) %>%
    mutate(winter = case_when(month<=3  ~year_before,
                              TRUE ~ winter)) %>%
    mutate(winter = case_when(month>=11  ~year_char,
                              TRUE ~ winter)) %>%
    mutate(winter = as.character(winter))
  
  # Import Qiime2 RPCA Results ----------------------------------------------
  
  #Import Qiime2 Results
  file = paste("Qiime_Results/",marker,"/ordination.qza",sep="")
  print(file)
  pco<-read_qza(file)
  pco$uuid
  #look at data
  head(pco$data$ProportionExplained)
  pco$data$Vectors[1:5, 1:4]
  
  #create proportion explained labels
  label.PC1 <- paste("PC1: ", round(pco$data$ProportionExplained$PC1, 3)*100,"%")
  label.PC1
  label.PC2 <- paste("PC2: ", round(pco$data$ProportionExplained$PC2, 3)*100,"%")
  label.PC2
  label.PC3 <- paste("PC3: ", round(pco$data$ProportionExplained$PC3, 3)*100,"%")
  label.PC3
  
  #Join with sample data
  pcscores <- left_join(pco$data$Vectors, samp.c, by="SampleID")
  
  #format loading scores
  loadings <- as.data.frame(pco$data$Species)
  loadings$ASV <- loadings$FeatureID
  
  #join on OTU, adding taxa info
  loadings <- left_join(loadings, tax.c, by="ASV")
  
  pcscores[1:5, 1:9]  #long because of sample data
  head(loadings)
  

# Import Aitchenson Distance Matrix ---------------------------------------
  file = paste("Qiime_Results/",marker,"/distance.qza",sep="")
  print(file)
  ait <-read_qza(file)
  ait
  aitmat <- as.matrix(ait$data)
  aitdist <- as.dist(aitmat)
  

# Cluster -----------------------------------------------------------------
  
  aitclust <- hclust(aitdist, method= "ward.D2")
  plot(aitclust)
  aitdend <- as.dendrogram(aitclust, hang = -1, lwd = 3, lty = 3, sub = "")
  #dend_cols <- select(samp.c,c('SampleID','date'))
  plot(aitdend)
  
  
  #Can plot with ggplot2
  
  dendr <- dendro_data(aitclust, type="rectangle") 
  
  #Define cluster based on similarity percent - good for plotting
  clust <- cutree(aitclust, h = 6)               # find 'cut' clusters (k) or choose level (h) numeric scalar or vector with heights where the tree should be cut.
  clust2 <- cutree(aitclust, h = 7)
  clust3 <- cutree(aitclust, h = 10)  
  clust4 <- cutree(aitclust, h = 16)  
  
  clust.df <- data.frame(label = names(clust), cluster_6 = clust, cluster_7 = clust2, cluster_10 = clust3, cluster_16 = clust4)
  sapply(clust.df, mode)
  clust.df$label <- as.character(clust.df$label)
  clust.df <- clust.df[order(clust.df$label),]
  
  #join sample data with cluster df and simprofCLUSTERS df:
  tree_scores <- as_tibble(clust.df) %>% rename(SampleID = label) %>% left_join(pcscores, by='SampleID')
  tree_data <- as_tibble(dendr$labels) %>% 
    rename(SampleID = label) %>% 
    left_join(tree_scores, by='SampleID') %>% 
    arrange(month)

  
  #By MONTH
  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
    #geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=year_char),alpha=1, size=1)+
    geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=fct_reorder(month_char,month)),alpha=1, size=1)+
    geom_text(data=tree_data, aes(x=x, y=y, label=year, hjust=-.4), size=1) +
    geom_text(data=tree_data, aes(x=x, y=y, label=libraryID, hjust=-6), size=2) +
    coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
    scale_colour_manual(values= c( "#6baed6", "#3182bd", "#08519c", 
                                   "#a1d99b", "#74c476", "#31a354", "#006d2c", 
                                   "#fdae6b", "#fd8d3c", "#e6550d", "#a63603","#9ecae1"))+
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  
  filename = paste(directory, 'Aitdist_',marker,'_hclust_month.png', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  #By Year
  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
    #geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=year_char),alpha=1, size=1)+
    geom_text(data=tree_data, aes(x=x, y=y, label=cluster_6, hjust=40), size=2) +
    geom_hline(yintercept=6, color="blue", linetype = "longdash", alpha=0.8) +
    geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=fct_reorder(year_char,year)),alpha=1, size=1)+
    geom_text(data=tree_data, aes(x=x, y=y, label=year, hjust=-.4), size=1) +
    coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
    scale_colour_manual(values= c( "#6baed6", "#3182bd", "#08519c", 
                                   "#a1d99b", "#74c476", "#31a354", "#006d2c", 
                                   "#fdae6b", "#fd8d3c", "#e6550d", "#a63603","#9ecae1"))+
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  
  filename = paste(directory, 'Aitdist_',marker,'_hclust_year.png', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  
  #Nice plots for paper
  #month
  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
    geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=fct_reorder(month_char,month)),alpha=1, size=1)+
    coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
    scale_colour_manual(values= c( "#6baed6", "#3182bd", "#08519c", 
                                   "#a1d99b", "#74c476", "#31a354", "#006d2c", 
                                   "#fdae6b", "#fd8d3c", "#e6550d", "#a63603","#9ecae1"))+
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  
  filename = paste(directory, 'Aitdist_',marker,'_hclust_month_lim.svg', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  filename = paste(directory, 'Aitdist_',marker,'_hclust_month_lim.png', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  #By Year
  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
    geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=fct_reorder(year_char,year)),alpha=1, size=1)+
    coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
    scale_colour_manual(values= c( "#6baed6", "#3182bd", "#08519c", 
                                   "#a1d99b", "#74c476", "#31a354", "#006d2c", 
                                   "#fdae6b", "#fd8d3c", "#e6550d", "#a63603","#9ecae1"))+
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  
  filename = paste(directory, 'Aitdist_',marker,'_hclust_year_lim.svg', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  filename = paste(directory, 'Aitdist_',marker,'_hclust_year_lim.png', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  #By Year/BW
  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
    geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=year),alpha=1, size=1)+
    coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  
  filename = paste(directory, 'Aitdist_',marker,'_hclust_year_lim_scale.svg', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  filename = paste(directory, 'Aitdist_',marker,'_hclust_year_lim_scale.png', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  #By PC score
  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
    geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=PC1),alpha=1, size=1)+
    geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-2.5,ymax=y-3.5, color=PC2),alpha=1, size=1)+
    #geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-4,ymax=y-5, color=PC3),alpha=1, size=1)+
    coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
    scale_colour_gradient2()+
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  
  filename = paste(directory, 'Aitdist_',marker,'_hclust_PCs_lim_scale.svg', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  filename = paste(directory, 'Aitdist_',marker,'_hclust_PCs_lim_scale.png', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  
  # p0 <- ggplot() + 
  #   geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  #   geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=PC1),alpha=1, size=1)+
  #   geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-2.5,ymax=y-3.5, color=PC2),alpha=1, size=1)+
  #   geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-4,ymax=y-5, color=PC3),alpha=1, size=1)+
  #   coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  #   scale_colour_gradient2()+
  #   theme(axis.line.y=element_blank(),
  #         axis.ticks.y=element_blank(),
  #         axis.text.y=element_blank(),
  #         axis.title.y=element_blank(),
  #         panel.background=element_rect(fill="white"),
  #         panel.grid=element_blank())
  # p1 <- ggplot() + 
  #   #geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
  #   geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=PC1),alpha=1, size=1)+
  #   geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-2.5,ymax=y-3.5, color=PC2),alpha=1, size=1)+
  #   geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-4,ymax=y-5, color=PC3),alpha=1, size=1)+
  #   coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  #   scale_colour_gradient2()+
  #   theme(axis.line.y=element_blank(),
  #         axis.ticks.y=element_blank(),
  #         axis.text.y=element_blank(),
  #         axis.title.y=element_blank(),
  #         panel.background=element_rect(fill="white"),
  #         panel.grid=element_blank())
  # 
  # p2 <- ggplot()+
  #   geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=fct_reorder(month_char,month)),alpha=1, size=1)+
  #   coord_flip() + scale_y_reverse(expand=c(.5, 0)) + 
  #   scale_colour_manual(values= c( "#6baed6", "#3182bd", "#08519c", 
  #                                  "#a1d99b", "#74c476", "#31a354", "#006d2c", 
  #                                  "#fdae6b", "#fd8d3c", "#e6550d", "#a63603","#9ecae1"))+
  #   theme(axis.line.y=element_blank(),
  #         axis.ticks.y=element_blank(),
  #         axis.text.y=element_blank(),
  #         axis.title.y=element_blank(),
  #         panel.background=element_rect(fill="white"),
  #         panel.grid=element_blank())
  #   p2
  #   grid.arrange(p0,p1, p2, nrow = 1)
  
  # #Plot distribution of samples in clustered groups
  # ggplot(econdatalong, aes(x=Country, y=value))+
  #   geom_bar(stat='identity', fill="forest green")+
  #   facet_wrap(~measure,  ncol=1)
  # 
  # p <- tree_data %>% ggplot(aes(x=month, color=fct_reorder(month_char,month)))+ 
  #   geom_histogram(binwidth=1,stat="count") +
  #   facet_wrap('cluster_6') +
  #   scale_x_continuous(breaks = seq(1, 12, by = 1)) +
  #   scale_colour_manual(values= c( "#6baed6", "#3182bd", "#08519c", 
  #                                  "#a1d99b", "#74c476", "#31a354", "#006d2c", 
  #                                  "#fdae6b", "#fd8d3c", "#e6550d", "#a63603","#9ecae1"))
  # p
  # 
  # p <- tree_data %>% ggplot(aes(x=year, color=fct_reorder(year_char,year)))+ 
  #   geom_histogram(binwidth=1,stat="count") +
  #   facet_wrap('cluster_6') +
  #   scale_x_continuous(breaks = seq(2008, 2018, by = 1)) +
  #   scale_colour_manual(values= c( "#6baed6", "#3182bd", "#08519c", 
  #                                  "#a1d99b", "#74c476", "#31a354", "#006d2c", 
  #                                  "#fdae6b", "#fd8d3c", "#e6550d", "#a63603","#9ecae1")) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # p
  # p <- tree_data %>% ggplot(aes(x=libraryID, color=libraryID,))+ 
  #   geom_histogram(stat="count") +
  #   facet_wrap('cluster_6') +
  #   #scale_x_continuous(breaks = seq(2008, 2018, by = 1)) +
  #   scale_colour_manual(values= c( "#6baed6", #"#3182bd", "#08519c", 
  #                                  "#a1d99b", #"#74c476", "#31a354", "#006d2c", 
  #                                  "#fdae6b"#, "#fd8d3c", "#e6550d", "#a63603","#9ecae1"
  #                                  )) +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # p
  

# Seasonality -------------------------------------------------------------


# Distance between each timepoint -----------------------------------------

  #Plot distance between adjacent points
  #transform aitchenson distance matrix into format for plotting change over time
  aitmat <- as.matrix(ait$data)
  aitmat <-as_tibble(aitmat,rownames = "SampleID")
  #just a couple columns
  test2 <- select(samp.c,c('SampleID','DATE_TIME'))
  aitmat <- left_join(aitmat, test2, by="SampleID")
  aitmat <- arrange(aitmat, DATE_TIME)
  aitmat <- select(aitmat, -SampleID)
  #for plotting dist between points, want date to next date?
  df <- aitmat
  df %<>%
    pivot_longer(-DATE_TIME, names_to = "SampleID", values_to = "distance")
  df <-  rename(df,'DATE_TIME1' = 'DATE_TIME')
  #want date of recipient sample
  df <- left_join(df, test2, by="SampleID")
  
  #Plot date1 to date, with value distance
  #sort by both dates
  df <- arrange(df, DATE_TIME1,DATE_TIME)
  df
  #just look at first date to everything else:
  #df <- filter(df, date1 == '2008-01-03 15:44:11')
  
  df <- filter(df, DATE_TIME1 != DATE_TIME)
  #no dates that are less than the sample data
  df <- filter(df, DATE_TIME > DATE_TIME1)
  df
  
  #Just keep the first unique value (closest date to current date)
  dates <- distinct(df, DATE_TIME1, .keep_all = TRUE)
  dates
  # Plot
  p<- ggplot(data=dates, aes(x=DATE_TIME, y=distance)) +
    geom_path()+
    geom_point()
  p
  filename = paste(directory, 'Aitdist_',marker,'_distance_between_samps.png', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  

# Distance from Reference Time Point --------------------------------------

  #Distance from first time point (reference time point)
  
  #transform aitchenson distance matrix into format for plotting change over time
  aitmat <- as.matrix(ait$data)
  aitmat <-as_tibble(aitmat,rownames = "SampleID")
  #just a couple columns
  test2 <- select(samp.c,c('SampleID','DATE_TIME'))
  aitmat <- left_join(aitmat, test2, by="SampleID")
  aitmat <- arrange(aitmat, DATE_TIME)
  aitmat <- select(aitmat, -SampleID)
  
  df <- aitmat
  #df %>% select(id, num1, num2, cat) %>%
  #  pivot_longer(., cols = c(num1,num2), names_to = "Var", values_to = "Val")
  
  df %<>%
    pivot_longer(-DATE_TIME, names_to = "SampleID", values_to = "distance")
  
  df <-  rename(df,'DATE_TIME1' = 'DATE_TIME')
  #want date of recipient sample
  df <- left_join(df, test2, by="SampleID")
  
  #Plot date1 to date, with value distance
  #sort by both dates
  df <- arrange(df, DATE_TIME1,DATE_TIME)
  df
  #just look at first date to everything else:
  df <- filter(df, DATE_TIME1 == first(df$DATE_TIME1))
  #df <- filter(df, DATE_TIME1 == '2008-01-03 15:44:00')
  
  
  p<- ggplot(data=df, aes(x=DATE_TIME, y=distance)) +
    geom_path()+
    geom_point()
  p
  
  filename = paste(directory, 'Aitdist_',marker,'_distance_from_firstsamp.png', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  

# Plot distance as a function of time in between samples ------------------

  aitmat <- as.matrix(ait$data)
  aitmat <-as_tibble(aitmat,rownames = "SampleID")
  #just a couple columns
  test2 <- select(samp.c,c('SampleID','DATE_TIME'))
  aitmat <- left_join(aitmat, test2, by="SampleID")
  aitmat <- arrange(aitmat, DATE_TIME)
  aitmat <- select(aitmat, -SampleID)
  #for plotting dist between points, want date to next date?
  #aitmat %>%
  #  gather(var, val, 2:ncol(aitmat)) %>%
  #  spread(Series.Description, val)
  df <- aitmat
  #df %>% select(id, num1, num2, cat) %>%
  #  pivot_longer(., cols = c(num1,num2), names_to = "Var", values_to = "Val")
  
  df %<>%
    pivot_longer(-DATE_TIME, names_to = "SampleID", values_to = "distance")
  
  df <-  rename(df,'DATE_TIME1' = 'DATE_TIME')
  #want date of recipient sample
  df <- left_join(df, test2, by="SampleID")
  
  #Plot date1 to date, with value distance
  #sort by both dates
  df <- arrange(df, DATE_TIME1,DATE_TIME)
  #just look at first date to everything else:
  #df <- filter(df, date1 == '2008-01-03 15:44:11')
  #remove distance from sample sample to each other (same date)
  df <- filter(df, DATE_TIME1 != DATE_TIME)
  
  df$diff_in_days<- difftime(df$DATE_TIME ,df$DATE_TIME1 , units = c("days"))
  #df$diff_in_days<- difftime(df$date ,df$date1 , units = c("weeks"))
  #make absolute value of diff in days, right now it's negative if before
  df$diff_in_days<-abs(df$diff_in_days)
  #round to nearest integer
  df$diff_in_days<-round(df$diff_in_days,digits=1 )
  df
  #round to nearest 10 days/weeks
  #df$diff_in_days<-round(df$diff_in_days,digits=-1 )
  #round to nearest 100 days
  #df$diff_in_days<-round(df$diff_in_days,digits=-2 )
  
  p<- ggplot(data=df, aes(x=diff_in_days, y=distance)) +
    #geom_path()+
    #geom_point(alpha=0.08, color='cornflowerblue') +
    geom_point(alpha=0.15, color='cornflowerblue') +
    #stat_summary(aes(y = distance,group=1), fun.y=mean, colour="red", geom="line",group=1)+
    #stat_summary(aes(y = distance,group=1), fun.y=mean, colour="red", geom="point",group=1, alpha=0.8, size=1)
    stat_summary_bin(aes(y = distance),fun.y=mean, colour="red", geom="line", binwidth=30) +
    #stat_summary_bin(aes(y = distance),fun.y=mean, colour="red", geom="line", binwidth=7)
    stat_summary_bin(aes(y = distance),colour='red',
                     fun.data = mean_se,
                     binwidth=30, size= 0.4)
  #p + geom_vline(xintercept = 30)
  #p + geom_vline(xintercept = 60)
  p <- p + geom_vline(xintercept = 365, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 730, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 1095, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 1460, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 1825, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 2190, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 2555, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 2920, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 3285, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 3650, linetype='dashed',
                      color = "black", size=0.8)
  p <- p + geom_vline(xintercept = 4015, linetype='dashed',
                      color = "black", size=0.8)
  
  
  
  
  #stat_summary(aes(y = distance,group=1), fun=mean, colour="red", geom="line",group=1)
  p+labs(x='Time Difference in Days' , y='RPCA Aitchenson Distance' )+ ggtitle(marker)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.text = element_text(size = 12, colour = "black"),
          legend.title = element_text(face = "bold", size=12),
          legend.key = element_rect(fill = FALSE, colour = FALSE),
          legend.key.size = unit(0.1,"line")
    )+ guides(color = guide_legend(override.aes = list(size=5, shape=15)))
  
  filename = paste(directory, 'Aitdist_',marker,'_distance_over_time_30daybin_error.png', sep='')
  filename
  ggsave(filename,height = 4, width =6, units = 'in')
  

}
  
  