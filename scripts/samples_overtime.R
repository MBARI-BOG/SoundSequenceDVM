#032824
#kpitz

# Samples through time by depth and plateID (methodID and sampling type)

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
#library(cowplot)
library(tibble)


# Set Constants -----------------------------------------------------------------

marker = sym("12S")

# Set directory to save plots
plot_directory <- 'figures/Plate_stats/samples_overtime/'
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


# Plot --------------------------------------------------------------------

p1 <- meta %>%
  mutate(Method = case_when(PlateID %in% c('BT', 'CE') ~ 'Method 2: Plates BT, CE',
                            PlateID %in% c('JJ', 'RR') ~ 'Method 1: Plates JJ, RR')) %>%
  ggplot(aes(x=time, y=depth, shape=SAMPLING_platform, color=Method)) +
  geom_point(alpha=0.7)+
  scale_y_reverse()+
  scale_color_manual(values=c('darkorange', 'darkblue'))+
  labs(x='Day', y='Depth (meters)', title='Samples Taken Through Time by Method') +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"))
p1
filename = paste(plot_directory, 'samples_overtime_bymethod.png', sep='')
filename
ggsave(filename, width=8, height=2.5, units = 'in')
filename = paste(plot_directory, 'samples_overtime_bymethod.svg', sep='')
filename
ggsave(filename, width=8, height=2.5, units = 'in')
p1

