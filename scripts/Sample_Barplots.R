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

markers <- c("12S")

for (val in markers) {
  marker = sym(val)
  
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
  
  # Create Seasonal Variables -----------------------------------------------
  
  # samp.c %<>% mutate(DATE_TIME = as_datetime(SAMPLING_date_time)) %>%
  #   # extract date components
  #   mutate(month =  month(DATE_TIME)) %>%
  #   mutate(day =  day(DATE_TIME)) %>%
  #   mutate(year =  year(DATE_TIME)) %>%
  #   mutate(jday = yday(DATE_TIME)) %>%
  #   mutate(DATE = as_date(ymd_hms(SAMPLING_date_time))) %>%  # in order to ignore time, change to lubridate 'date' item
  #   mutate(ymd = ymd(DATE)) %>%
  #   mutate(month_2 = format.Date(DATE_TIME, "%m")) %>%  #two digit month
  #   #change column data types
  #   mutate(month_char = as.character(month)) %>%
  #   mutate(year_char = as.character(year)) %>%
  #   mutate(month_2 = as.character(month_2)) %>%
  #   #year/month for labeling
  #   unite(ym, year_char, month_2, sep='/', remove=FALSE) %>%
  #   #create seasonal categories
  #   # 12,1,2; 3,4,5; 6,7,8; 9,10,11
  #   mutate(season = "Dec-Feb") %>% #create new column
  #   mutate(season = case_when(month<=11 & month >=9  ~"Sep-Nov",
  #                             TRUE ~ season)) %>%
  #   mutate(season = case_when(month<=8 & month >=6  ~"June-Aug",
  #                             TRUE ~ season)) %>%
  #   mutate(season = case_when(month<=5 & month >=3  ~"March-May",
  #                             TRUE ~ season))
  
  
  # Lowest Taxonomic Annotation ---------------------------------------------
  
  head(tax.c)
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
      mutate(Genus = case_when(Species=='s_'| Species =='unassigned' | Species =='unknown'~as.character('unknown'),
                               TRUE ~ as.character(Species))) %>%
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
      mutate(Genus = case_when(Species=='s_'| Species =='unassigned' | Species =='unknown'~as.character('unknown'),
                               TRUE ~ as.character(Species))) %>%
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
      mutate(Genus = case_when(Species=='s_'| Species =='unassigned' | Species =='unknown'~as.character('unknown'),
                               TRUE ~ as.character(Species))) %>%
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
      mutate(Genus = case_when(Species=='s_'| Species =='unassigned' | Species =='unknown'~as.character('unknown'),
                               TRUE ~ as.character(Species))) %>%
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
      mutate(Genus = case_when(Species=='s_'| Species =='unassigned' | Species =='unknown'~as.character('unknown'),
                               TRUE ~ as.character(Species))) %>%
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
      mutate(Genus = case_when(Species=='s_'| Species =='unassigned' | Species =='unknown'~as.character('unknown'),
                               TRUE ~ as.character(Species))) %>%
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
    # #Genus
    # top_taxa <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    #   inner_join(gen_label,  by = c("ASV")) %>%
    #   mutate(Genus = case_when(Genus=='g_'| Genus =='unassigned' | Genus =='unknown'~as.character('unknown'),
    #                            TRUE ~ as.character(Genus))) %>%
    #   mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
    #                             TRUE ~ as.character(Phylum))) %>%
    #   filter(Phylum == var) %>%  #limit to a certain group
    #   group_by(Genus) %>%
    #   mutate(sum_per_tot = sum(per_tot)) %>%
    #   arrange(-sum_per_tot) %>%
    #   distinct(Genus,.keep_all = TRUE ) %>%
    #   select(Kingdom, Phylum, Class, Order, Family,Genus, sum_per_tot) %>%
    #   ungroup() %>%
    #   filter(Genus != var) %>% #remove general phylum level annot
    #   filter(Genus != 'Dinophyceae') %>%
    #   filter(Genus != 'Dinophyceae sp. CCMP1878') %>% #remove unclassified dino
    #   filter(Genus != 'Dinophyceae sp. UDMS0803') %>% #remove unclassified dino
    #   filter(Genus != 'unknown') %>% #remove unclassified genera
    #   select(Genus, sum_per_tot) %>%
    #   top_n(10)
    # print(top_taxa)
    # 
    # 
    # bp_top <- inner_join(potu.c, samp.c,  by = c("SampleID")) %>%
    #   inner_join(gen_label,  by = c("ASV")) %>%
    #   mutate(Genus = case_when(Genus=='g_'| Genus =='unassigned' | Genus =='unknown'~as.character('unknown'),
    #                            TRUE ~ as.character(Genus))) %>%
    #   mutate(Phylum = case_when(Class=='Dinophyceae' ~as.character('Dinophyceae'),
    #                             TRUE ~ as.character(Phylum))) %>%
    #   right_join(top_taxa) %>% #limit to most abundant taxa
    #   unite(Label, Class, Order, Family, Genus, sep="_", remove='False') %>%
    #   #fct_reorder(name, desc(val))
    #   ggplot(aes(x = fct_reorder(SampleID, desc(depth)), y = per_tot)) +
    #   geom_bar(stat = "identity", aes(fill = Label))+
    #   scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1)+
    #   labs(x="",y=paste("Percent ",var," Reads", sep=""))+
    #   #scale_x_discrete(breaks = year_ticks, labels = year_labels, name = "",drop = FALSE)+
    #   theme_minimal() +
    #   theme(
    #     #legend
    #     legend.position="bottom",legend.direction="vertical",
    #     legend.text=element_text(colour=textcol,size=5,face="bold"),
    #     legend.key.height=grid::unit(0.3,"cm"),
    #     legend.key.width=grid::unit(0.3,"cm"),
    #     legend.title=element_text(colour=textcol,size=5,face="bold"),
    #     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=7,colour=textcol),
    #     #axis.text.x=element_text(size=7,colour=textcol),
    #     axis.text.y=element_text(size=6,colour=textcol),
    #     axis.title.y = element_text(size=6),
    #     plot.background=element_blank(),
    #     panel.border=element_blank(),
    #     panel.grid.minor = element_blank(),
    #     panel.grid.major = element_line(size = .25),
    #     plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    #     plot.title=element_blank())
    # bp_top 
    # 
    # 
    # filename = paste(directory, marker,'_top10gen_bar_',var,'_comp.png', sep='')
    # print(var)
    # filename
    # ggsave(filename,height = 8, width =10, units = 'in')
    # filename = paste(directory, marker,'_top10gen_bar_',var,'_comp.svg', sep='')
    # print(var)
    # filename
    # ggsave(filename,height = 8, width =10, units = 'in')
    
  }
}




