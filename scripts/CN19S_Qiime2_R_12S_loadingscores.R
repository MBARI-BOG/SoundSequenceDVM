## written Katie Pitz
## 11/23/20
## Plot DEICODE loading score results alongside read# heatmap

marker = '12S'
prefix = 'CN19S'
data_directory = "data/filtered_seq_data/"
# Set directory to save plots
plot_dir <- 'figures/Loadingscores/'
results_directory <- 'Qiime_Results/'

# Load Libraries -----------------------------------------------------------------
library(qiime2R)
library(tidyverse)
library(lubridate) #for date modifications
library(ggthemes)
library(magrittr)
#for heatmap section
library(gridExtra)
library(RColorBrewer) #colors for plotting
library(forcats) #working with factor data (reordering)
library(patchwork) #putting graphs together on same plot

# Import Data -------------------------------------------------------------

#ASV table
print('ASV table')
file = paste(prefix,"_",marker,"_otu_filtered.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
otu.c <- read_csv(filepath) %>% rename_with(.cols = 1, ~"ASV")

#taxa table
print('taxa table')
file = paste(prefix,"_",marker,"_taxa_filtered.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
tax.c <- read_csv(filepath) %>% rename_with(.cols = 1, ~"ASV")

#metadata table
print('metadata table')
file = paste(prefix,"_",marker,"_meta_filtered.csv", sep='')
filepath = paste(data_directory, file, sep='')
print(filepath)
samp.c <- read_csv(filepath) %>% rename('SampleID' = 'sample_name')

#OTU table long format with percent total reads
potu.c <- otu.c %>%
  tidyr::pivot_longer( -ASV, names_to ='SampleID',values_to = 'reads' ) %>%
  group_by(SampleID) %>%
  mutate(per_tot = reads / sum(reads) *100) %>%
  ungroup() %>%
  arrange(-reads)
head(potu.c)

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

# Import Qiime2 Results ---------------------------------------------------

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


# Plot high loading scores -----------------

#for 12S color by Family
pcs <- c("PC1", "PC2", "PC3")
for (val in pcs) {
  pc = sym(val)
  #get top 25 positive and negative loading scores
  df_l <-as_tibble(loadings)
  df_l %<>% arrange(!!pc)
  top_l <- top_n(df_l,25, !!pc)
  bot_l <- top_n(df_l, -25, !!pc)
  tot_l <- full_join(top_l, bot_l)
  tot_l  #contains top and bottom 25 ASVs
  
  #Plot these loading scores as barplot:
  
  #barplot of top 50 loadings scores (by species)
  bp_species <- tot_l %>%
    mutate(round_PC= round(!!pc,5) ) %>%
    unite( ASV_label,round_PC, FeatureID,Order,Family,Genus,Species, sep = "_", remove=FALSE) %>%
    mutate(ASV_label = fct_reorder(ASV_label, !!pc)) %>%
    ggplot(aes(x = !!pc, y = ASV_label)) +
    geom_bar(stat = "identity", aes(fill = Family))+
    geom_text(aes(label=Species),position = position_stack(vjust = 0.5), size=2) +
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x = paste(pc," loading score", sep=""))+
    theme_minimal() +
    theme(legend.position="right",axis.title.y = element_blank(), axis.text.y = element_blank(),
          legend.text=element_text(colour='grey40',size=7,face="bold"),
          legend.key.height=grid::unit(0.3,"cm"),
          legend.key.width=grid::unit(0.3,"cm"),
          legend.title=element_text(colour="grey40",size=7,face="bold"),
          plot.background=element_blank(),
          panel.border=element_blank(),
          plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
          plot.title=element_blank())
  bp_species
  filename = paste(plot_dir, '12S_',pc,'_toploadingsp_bar.png', sep='')
  print('Plot of top 50 loadings scores by species:')
  filename
  ggsave(filename,height = 6, width =6, units = 'in')
  
  #barplot of top 50 loading scores (by genus)
  bp_final <- tot_l %>%
    mutate(round_PC= round(!!pc,5) ) %>%
    unite( ASV_label,round_PC, FeatureID,Order,Family,Genus,Species, sep = "_", remove=FALSE) %>%
    mutate(ASV_label = fct_reorder(ASV_label, !!pc)) %>%
    ggplot(aes(x = !!pc, y = ASV_label)) +
    geom_bar(stat = "identity", aes(fill = Family))+
    geom_text(aes(label=Genus),position = position_stack(vjust = 0.5), size=2) +
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x = paste(pc," loading score", sep=""))+
    theme_minimal() +
    theme(legend.position="right",axis.title.y = element_blank(), axis.text.y = element_blank(),
          legend.text=element_text(colour='grey40',size=7,face="bold"),
          legend.key.height=grid::unit(0.3,"cm"),
          legend.key.width=grid::unit(0.3,"cm"),
          legend.title=element_text(colour="grey40",size=7,face="bold"),
          plot.background=element_blank(),
          panel.border=element_blank(),
          plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
          plot.title=element_blank())
  bp_final
  filename = paste(plot_dir, '12S_',pc,'_toploadinggenus_bar.png', sep='')
  print('Plot of top 50 loadings scores by genus:')
  filename
  ggsave(filename,height = 6, width =6, units = 'in')
  
  #plot these species across dataset by depth
  #ylabels
  meta %>%
    arrange(depth) %>%
    #distinct(depth_bin, .keep_all = TRUE) %>%
    pull(depth_bin) -> depth_labels
  
  meta %>%
    arrange(depth) %>%
    #distinct(depth_bin, .keep_all = TRUE) %>%
    pull(SampleID) -> depth_ticks 
  
  #depth_labels
  #depth_ticks
  
  bp_top3 <- tot_l %>%
    mutate(PC_sign  = case_when(!!pc>0 ~'Pos', !!pc<0 ~ 'Neg')) %>%
    mutate(Family = case_when(Family=='unassigned' | Family =='unknown'| Family =='s_'| Family =='no_hit' ~as.character('Unknown'),
                              TRUE ~ as.character(Family))) %>%
    left_join(potu.c) %>%
    left_join(meta %>% select(SampleID, depth, depth_bin)) %>%
    ggplot(aes(x = fct_reorder(SampleID, depth), y = per_tot)) +
    geom_bar(stat = "identity", aes(fill = Family))+
    scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(x="",y="Percent Total Reads")+
    scale_x_discrete(breaks = depth_ticks, labels = depth_labels, name = "",drop = FALSE)+
    #scale_y_discrete(breaks = ASV_ticks, labels = ASV_labels, name = "",expand=c(0,0))+
    theme_minimal() +
    facet_grid(rows = vars(PC_sign)) +
    theme(
      #legend
      legend.position="right",legend.direction="vertical",
      legend.text=element_text(colour="grey40",size=7,face="bold"),
      legend.key.height=grid::unit(0.3,"cm"),
      legend.key.width=grid::unit(0.3,"cm"),
      legend.title=element_text(colour="grey40",size=7,face="bold"),
      axis.text.x=element_text(size=4,colour="grey40", angle = 90),
      axis.text.y=element_text(size=7,colour="grey40"),
      plot.background=element_blank(),
      panel.border=element_blank(),
      plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
      plot.title=element_blank())
  bp_top3
  filename = paste(plot_dir, '12S_',pc,'_toploadingsplit_bar.png', sep='')
  filename
  ggsave(filename,height = 5, width =10, units = 'in')
  
  # mean percent abundance of taxa over depth:
  # smoothed line:
  p <- tot_l %>%
    left_join(potu.c) %>%
    left_join(meta %>% select(SampleID, depth, depth_bin)) %>%
    group_by(SampleID, Family) %>%
    mutate(per_tot = sum(per_tot)) %>%
    ungroup() %>%
    distinct(SampleID, Family, .keep_all = TRUE) %>%
    ggplot(aes(x=depth,y=per_tot, color=Family)) +
    geom_point()+
    geom_smooth()+
    scale_color_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    labs(y='Percent Total Reads', x='Depth (m)')
  filename = paste(plot_dir, '12S_',pc,'_toploading_smoothlines.png', sep='')
  filename
  ggsave(filename,height = 8, width =10, units = 'in')
  p
  
  #By depth bin, boxplot?:
  #to generate list of depths to bin over:
  #seq(0, 750, by=10)
  
  p <- tot_l %>%
    #mutate(PC_sign  = case_when(!!pc>0 ~'Pos', !!pc<0 ~ 'Neg')) %>%
    mutate(Family = case_when(Family=='unassigned' | Family =='unknown'| Family =='s_'| Family =='no_hit' ~as.character('Unknown'),
                              TRUE ~ as.character(Family))) %>%
    #filter(Family %in% c('Engraulidae','Merlucciidae', 'Myctophidae', 'Sebastidae')) %>%
    left_join(potu.c) %>%
    left_join(meta %>% select(SampleID, depth, depth_bin)) %>%
    mutate(depth_cut = cut(depth,seq(0, 750, by=50),include.lowest = TRUE )) %>%
    group_by(SampleID, Family) %>%
    mutate(per_tot = sum(per_tot)) %>%
    ungroup() %>%
    distinct(SampleID, Family, .keep_all = TRUE) %>%
    group_by(depth_cut, Family) %>%
    mutate(mean_per_tot = mean(per_tot)) %>%
    ungroup() %>%
    distinct(SampleID, Family, .keep_all = TRUE) %>%
    ggplot(aes(x = depth_cut, y = per_tot, fill=Family)) +
    geom_boxplot()+
    scale_color_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
    scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1)+
    labs(x='Depth Bin (m)', y='Percent Total Reads')
    #geom_point(aes(x=depth,y=mean_per_tot, color=Family))+
    #geom_smooth(aes(x=depth,y=mean_per_tot, color=Family))
  filename = paste(plot_dir, '12S_',pc,'_toploading_boxplot_50m.png', sep='')
  filename
  ggsave(filename,height = 8, width =15, units = 'in')
  p
  
}

