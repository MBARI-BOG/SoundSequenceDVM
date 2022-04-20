#040221
#kpitz

# Set directory to save plots
directory <- 'figures/Model/'


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
marker = sym("12S")

# Import Data -------------------------------------------------------------
data_directory = "Data/filtered_seq_data/"

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


# Import species designations

filepath = "/Users/kpitz/Projects/CN19all/Canon2019Spring_Taxa_Categories_KBB_KP.csv"

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


# Group data by Ecological Category --------------------------------------------

Ecol_data <- tax.c %>% left_join(sp_desig) %>%
  left_join(potu.c) %>%
  group_by(Ecological_Category, SampleID) %>%
  mutate(sum_reads = sum(reads)) %>%
  mutate(sum_per_tot = sum(per_tot)) %>%
  ungroup() %>%
  distinct(SampleID, Ecological_Category, .keep_all=TRUE) %>%
  select(SampleID, Ecological_Category, sum_reads, sum_per_tot)

# Create Loess Model ------------------------------------------------------

test <- Ecol_data %>% 
  filter(Ecological_Category == 'epipelagic') %>%
  left_join(meta %>% select(SampleID, diel, depth)) %>%
  filter(diel == 'night') %>%
  select(Ecological_Category, depth, sum_per_tot) %>%
  arrange(depth)

# In base R:

loessMod10 <- loess(sum_per_tot ~ depth, data=test, span=0.30)
smoothed10 <- predict(loessMod10) 

loessMod20 <- loess(sum_per_tot ~ depth, data=test, span=0.50)
smoothed20 <- predict(loessMod20) 

loessMod30 <- loess(sum_per_tot ~ depth, data=test, span=0.75)
smoothed30 <- predict(loessMod30) 

# Plot it
plot(test$sum_per_tot, x=test$depth, type="l", main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
plot(test$sum_per_tot, x=test$depth,main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
lines(smoothed10, x=test$depth, col="red")
lines(smoothed20, x=test$depth, col="purple")
lines(smoothed30, x=test$depth, col="green")

# In tidyverse:

test <- Ecol_data %>% 
  #filter(Ecological_Category == 'epipelagic') %>%
  filter(Ecological_Category %in% c('epipelagic', 'mesopelagic')) %>%
  left_join(meta %>% select(SampleID, diel, depth)) %>%
  filter(diel == 'night') %>%
  select(Ecological_Category, depth, sum_per_tot) %>%
  arrange(depth)


#split by Ecol_cat
sdata <-split(test, test$Ecological_Category)
#calculate the predicted values for each country
data.1.with.pred <- lapply(sdata, function(df){
  df$pred.response  <-stats::predict(stats::loess(sum_per_tot ~ depth, span = .60, data=df))
  df
})

#merge back into 1 dataframe
data.1.with.pred <-dplyr::bind_rows(data.1.with.pred )

ggplot(data.1.with.pred, aes(x=depth, y=pred.response, color=Ecological_Category)) +
  geom_point(aes(x=depth, y=sum_per_tot), size = 3, alpha=0.3) + 
  geom_line( size=0.5, alpha= 1) 

# split by time period

test <- Ecol_data %>% 
  filter(Ecological_Category == 'epipelagic') %>%
  #filter(Ecological_Category %in% c('epipelagic', 'mesopelagic')) %>%
  left_join(meta %>% select(SampleID, diel, depth)) %>%
  filter(diel %in% c('night', 'day')) %>%
  #filter(diel == 'night') %>%
  select(Ecological_Category, depth, sum_per_tot, diel) %>%
  arrange(depth)


#split by Ecol_cat
sdata <-split(test, test$diel)
#calculate the predicted values for each country
data.1.with.pred <- lapply(sdata, function(df){
  df$pred.response  <-stats::predict(stats::loess(sum_per_tot ~ depth, span = .60, data=df))
  df
})

#merge back into 1 dataframe
data.1.with.pred <-dplyr::bind_rows(data.1.with.pred )

ggplot(data.1.with.pred, aes(x=-depth, y=pred.response, color=diel)) +
  geom_point(aes(y=sum_per_tot), size = 3, alpha=0.3) + 
  geom_line( size=0.5, alpha= 1) +
  coord_flip()
