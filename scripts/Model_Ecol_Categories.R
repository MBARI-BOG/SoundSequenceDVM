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

# Loess Model in base R with diff span values ------------------------------------------------------

# Individual

test <- Ecol_data %>% 
  filter(Ecological_Category == 'epipelagic') %>%
  left_join(meta %>% select(SampleID, diel, depth)) %>%
  filter(diel == 'night') %>%
  select(Ecological_Category, depth, sum_per_tot) %>%
  arrange(depth)

# In base R:

loessMod10 <- loess(sum_per_tot ~ depth, data=test, span=0.30)
smoothed10 <- predict(loessMod10, se=TRUE) 

loessMod20 <- loess(sum_per_tot ~ depth, data=test, span=0.50)
smoothed20 <- predict(loessMod20, se=TRUE) 

loessMod25 <- loess(sum_per_tot ~ depth, data=test, span=0.60)
smoothed25 <- predict(loessMod25, se=TRUE) 

loessMod30 <- loess(sum_per_tot ~ depth, data=test, span=0.75)
smoothed30 <- predict(loessMod30, se=TRUE) 

# Plot it
plot(test$sum_per_tot, x=test$depth,main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
lines(smoothed10$fit, x=test$depth, col="red")
lines(x=test$depth,smoothed10$fit - qt(0.975,smoothed10$df)*smoothed10$se, lty=2, col="red")
lines(x=test$depth,smoothed10$fit + qt(0.975,smoothed10$df)*smoothed10$se, lty=2, col="red")

lines(smoothed20$fit, x=test$depth, col="purple")
lines(x=test$depth,smoothed20$fit - qt(0.975,smoothed20$df)*smoothed20$se, lty=2, col="purple")
lines(x=test$depth,smoothed20$fit + qt(0.975,smoothed20$df)*smoothed20$se, lty=2, col="purple")

lines(smoothed25$fit, x=test$depth, col="orange")
lines(x=test$depth,smoothed25$fit - qt(0.975,smoothed25$df)*smoothed25$se, lty=2, col="orange")
lines(x=test$depth,smoothed25$fit + qt(0.975,smoothed25$df)*smoothed25$se, lty=2, col="orange")

lines(smoothed30$fit, x=test$depth, col="green")
lines(x=test$depth,smoothed30$fit - qt(0.975,smoothed30$df)*smoothed30$se, lty=2, col="green")
lines(x=test$depth,smoothed30$fit + qt(0.975,smoothed30$df)*smoothed30$se, lty=2, col="green")

#lines(smoothed20, x=test$depth, col="purple")
#lines(smoothed25, x=test$depth, col="orange")
#lines(smoothed30, x=test$depth, col="green")

#add legend to plot
legend('topright',
       col = c('red', 'purple', 'orange','green'),
       lwd = 2,
       c('span = 0.3', 'span = 0.5', 'span = 0.6', 'span = 0.75'))

# Open a pdf file
filename = paste(directory, marker,'_LOESS_spans_epipelagic_night_withse.pdf', sep='')
pdf(filename) 
# 2. Create a plot
plot(test$sum_per_tot, x=test$depth,main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
lines(smoothed10$fit, x=test$depth, col="red")
lines(x=test$depth,smoothed10$fit - qt(0.975,smoothed10$df)*smoothed10$se, lty=2, col="red")
lines(x=test$depth,smoothed10$fit + qt(0.975,smoothed10$df)*smoothed10$se, lty=2, col="red")

lines(smoothed20$fit, x=test$depth, col="purple")
lines(x=test$depth,smoothed20$fit - qt(0.975,smoothed20$df)*smoothed20$se, lty=2, col="purple")
lines(x=test$depth,smoothed20$fit + qt(0.975,smoothed20$df)*smoothed20$se, lty=2, col="purple")

lines(smoothed25$fit, x=test$depth, col="orange")
lines(x=test$depth,smoothed25$fit - qt(0.975,smoothed25$df)*smoothed25$se, lty=2, col="orange")
lines(x=test$depth,smoothed25$fit + qt(0.975,smoothed25$df)*smoothed25$se, lty=2, col="orange")

lines(smoothed30$fit, x=test$depth, col="green")
lines(x=test$depth,smoothed30$fit - qt(0.975,smoothed30$df)*smoothed30$se, lty=2, col="green")
lines(x=test$depth,smoothed30$fit + qt(0.975,smoothed30$df)*smoothed30$se, lty=2, col="green")
#add legend to plot
legend('topright',
       col = c('red', 'purple', 'orange','green'),
       lwd = 2,
       c('span = 0.3', 'span = 0.5', 'span = 0.6', 'span = 0.75'))
# Close the pdf file
dev.off()


# Loop over span values across Ecol Cat and Diel --------------------------

# with standard error

ecols = c('epipelagic', 'mesopelagic', 'cosmopolitan', 'benthopelagic')
diels = c('day', 'night')

for (val in ecols) {
  ecol_level = sym(val)
  for (val2 in diels) {
    diel_level = sym(val2)
    test <- Ecol_data %>% 
      filter(Ecological_Category == ecol_level) %>%
      left_join(meta %>% select(SampleID, diel, depth)) %>%
      filter(diel == diel_level) %>%
      select(Ecological_Category, depth, sum_per_tot) %>%
      arrange(depth)
    #LOESS MODEL
    loessMod10 <- loess(sum_per_tot ~ depth, data=test, span=0.30)
    smoothed10 <- predict(loessMod10, se=TRUE) 
    loessMod20 <- loess(sum_per_tot ~ depth, data=test, span=0.50)
    smoothed20 <- predict(loessMod20, se=TRUE) 
    loessMod25 <- loess(sum_per_tot ~ depth, data=test, span=0.60)
    smoothed25 <- predict(loessMod25, se=TRUE) 
    loessMod30 <- loess(sum_per_tot ~ depth, data=test, span=0.75)
    smoothed30 <- predict(loessMod30, se=TRUE) 
    # Plot and Save
    # Open a pdf file
    filename = paste(directory, marker,'_LOESS_spans_',ecol_level, '_', diel_level, '_withse.pdf', sep='')
    pdf(filename) 
    title = paste('Loess Smoothing and Prediction : ',ecol_level, ' ', diel_level, sep='')
    # 2. Create a plot
    #plot(test$sum_per_tot, x=test$depth, type="l", main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
    plot(test$sum_per_tot, x=test$depth,main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
    lines(smoothed10$fit, x=test$depth, col="red")
    lines(x=test$depth,smoothed10$fit - qt(0.975,smoothed10$df)*smoothed10$se, lty=2, col="red")
    lines(x=test$depth,smoothed10$fit + qt(0.975,smoothed10$df)*smoothed10$se, lty=2, col="red")
    
    lines(smoothed20$fit, x=test$depth, col="purple")
    lines(x=test$depth,smoothed20$fit - qt(0.975,smoothed20$df)*smoothed20$se, lty=2, col="purple")
    lines(x=test$depth,smoothed20$fit + qt(0.975,smoothed20$df)*smoothed20$se, lty=2, col="purple")
    
    lines(smoothed25$fit, x=test$depth, col="orange")
    lines(x=test$depth,smoothed25$fit - qt(0.975,smoothed25$df)*smoothed25$se, lty=2, col="orange")
    lines(x=test$depth,smoothed25$fit + qt(0.975,smoothed25$df)*smoothed25$se, lty=2, col="orange")
    
    lines(smoothed30$fit, x=test$depth, col="green")
    lines(x=test$depth,smoothed30$fit - qt(0.975,smoothed30$df)*smoothed30$se, lty=2, col="green")
    lines(x=test$depth,smoothed30$fit + qt(0.975,smoothed30$df)*smoothed30$se, lty=2, col="green")
    #add legend to plot
    legend('topright',
           col = c('red', 'purple', 'orange','green'),
           lwd = 2,
           c('span = 0.3', 'span = 0.5', 'span = 0.6', 'span = 0.75'))
    # Close the pdf file
    dev.off()
  }
}

# Without SE:

for (val in ecols) {
  ecol_level = sym(val)
  for (val2 in diels) {
    diel_level = sym(val2)
    test <- Ecol_data %>% 
      filter(Ecological_Category == ecol_level) %>%
      left_join(meta %>% select(SampleID, diel, depth)) %>%
      filter(diel == diel_level) %>%
      select(Ecological_Category, depth, sum_per_tot) %>%
      arrange(depth)
    #LOESS MODEL
    loessMod10 <- loess(sum_per_tot ~ depth, data=test, span=0.30)
    smoothed10 <- predict(loessMod10, se=FALSE) 
    loessMod20 <- loess(sum_per_tot ~ depth, data=test, span=0.50)
    smoothed20 <- predict(loessMod20, se=FALSE) 
    loessMod25 <- loess(sum_per_tot ~ depth, data=test, span=0.60)
    smoothed25 <- predict(loessMod25, se=FALSE) 
    loessMod30 <- loess(sum_per_tot ~ depth, data=test, span=0.75)
    smoothed30 <- predict(loessMod30, se=FALSE) 
    # Plot and Save
    # Open a pdf file
    filename = paste(directory, marker,'_LOESS_spans_',ecol_level, '_', diel_level, '.pdf', sep='')
    pdf(filename) 
    title = paste('Loess Smoothing and Prediction : ',ecol_level, ' ', diel_level, sep='')
    # 2. Create a plot
    #plot(test$sum_per_tot, x=test$depth, type="l", main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
    plot(test$sum_per_tot, x=test$depth,main=title, xlab="depth", ylab="sum_per_tot")
    lines(smoothed10, x=test$depth, col="red")
    lines(smoothed20, x=test$depth, col="purple")
    lines(smoothed25, x=test$depth, col="orange")
    lines(smoothed30, x=test$depth, col="green")
    #add legend to plot
    legend('topright',
           col = c('red', 'purple', 'orange','green'),
           lwd = 2,
           c('span = 0.3', 'span = 0.5', 'span = 0.6', 'span = 0.75'))
    # Close the pdf file
    dev.off()
  }
}


# In base R without loop:

loessMod10 <- loess(sum_per_tot ~ depth, data=test, span=0.30)
smoothed10 <- predict(loessMod10) 

loessMod20 <- loess(sum_per_tot ~ depth, data=test, span=0.50)
smoothed20 <- predict(loessMod20) 

loessMod25 <- loess(sum_per_tot ~ depth, data=test, span=0.60)
smoothed25 <- predict(loessMod25) 

loessMod30 <- loess(sum_per_tot ~ depth, data=test, span=0.75)
smoothed30 <- predict(loessMod30) 

# Plot it
plot(test$sum_per_tot, x=test$depth, type="l", main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
plot(test$sum_per_tot, x=test$depth,main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
lines(smoothed10, x=test$depth, col="red")
lines(smoothed20, x=test$depth, col="purple")
lines(smoothed25, x=test$depth, col="orange")
lines(smoothed30, x=test$depth, col="green")

#add legend to plot
legend('topright',
       col = c('red', 'purple', 'orange','green'),
       lwd = 2,
       c('span = 0.3', 'span = 0.5', 'span = 0.6', 'span = 0.75'))

# Open a pdf file
filename = paste(directory, marker,'_LOESS_spans_epipelagic_night.pdf', sep='')
pdf(filename) 
# 2. Create a plot
#plot(test$sum_per_tot, x=test$depth, type="l", main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
plot(test$sum_per_tot, x=test$depth,main="Loess Smoothing and Prediction : Epipelagic Night", xlab="depth", ylab="sum_per_tot")
lines(smoothed10, x=test$depth, col="red")
lines(smoothed20, x=test$depth, col="purple")
lines(smoothed25, x=test$depth, col="orange")
lines(smoothed30, x=test$depth, col="green")
#add legend to plot
legend('topright',
       col = c('red', 'purple', 'orange','green'),
       lwd = 2,
       c('span = 0.3', 'span = 0.5', 'span = 0.6', 'span = 0.75'))
# Close the pdf file
dev.off()


# Find optimum span based on minimizing sum of squares errors (http://r-statistics.co/Loess-Regression-With-R.html)

# define function that returns the SSE
calcSSE <- function(x){
  loessMod <- try(loess(sum_per_tot ~ depth, data=test, span=x), silent=T)
  res <- try(loessMod$residuals, silent=T)
  if(class(res)!="try-error"){
    if((sum(res, na.rm=T) > 0)){
      sse <- sum(res^2)  
    }
  }else{
    sse <- 99999
  }
  return(sse)
}

# Run optim to find span that gives min SSE, starting at 0.5
optim(par=c(0.5), calcSSE, method="SANN")

#$par
#[1] 0.08663962

# try this span value:
loessMod30 <- loess(sum_per_tot ~ depth, data=test, span=0.08663962)
smoothed30 <- predict(loessMod30) 

# Plot it
plot(test$sum_per_tot, x=test$depth, type="l", main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
plot(test$sum_per_tot, x=test$depth,main="Loess Smoothing and Prediction", xlab="depth", ylab="sum_per_tot")
lines(smoothed30, x=test$depth, col="red")

# Too messy I think!!

### Model in Tidyverse and combine ---------------
#DAY
test <- Ecol_data %>% 
  #filter(Ecological_Category == 'epipelagic') %>%
  filter(Ecological_Category %in% c('epipelagic', 'mesopelagic', 'benthopelagic', 'cosmopolitan')) %>%
  left_join(meta %>% select(SampleID, diel, depth, ESP)) %>%
  filter(diel == 'day') %>%
  #filter(diel %in% c('night', 'day')) %>%
  select(Ecological_Category, depth, sum_per_tot, ESP) %>%
  arrange(depth)

#split by Ecol_cat
sdata <-split(test, test$Ecological_Category)
#calculate the predicted values for each country
data.1.with.pred <- lapply(sdata, function(df){
  df$pred.response  <-stats::predict(stats::loess(sum_per_tot ~ depth, span = .50, data=df), se=TRUE)$fit
  df$se  <-stats::predict(stats::loess(sum_per_tot ~ depth, span = .50, data=df), se=TRUE)$se
  df$df  <-stats::predict(stats::loess(sum_per_tot ~ depth, span = .50, data=df), se=TRUE)$df
  df
})

#merge back into 1 dataframe
day.with.pred <-dplyr::bind_rows(data.1.with.pred ) %>%
  mutate(diel = 'day')

#NIGHT
test <- Ecol_data %>% 
  #filter(Ecological_Category == 'epipelagic') %>%
  filter(Ecological_Category %in% c('epipelagic', 'mesopelagic', 'benthopelagic', 'cosmopolitan')) %>%
  left_join(meta %>% select(SampleID, diel, depth, ESP)) %>%
  filter(diel == 'night') %>%
  #filter(diel %in% c('night', 'day')) %>%
  select(Ecological_Category, depth, sum_per_tot, ESP) %>%
  arrange(depth)

#split by Ecol_cat
sdata <-split(test, test$Ecological_Category)
#calculate the predicted values for each country
data.1.with.pred <- lapply(sdata, function(df){
  df$pred.response  <-stats::predict(stats::loess(sum_per_tot ~ depth, span = .50, data=df), se=TRUE)$fit
  df$se  <-stats::predict(stats::loess(sum_per_tot ~ depth, span = .50, data=df), se=TRUE)$se
  df$df  <-stats::predict(stats::loess(sum_per_tot ~ depth, span = .50, data=df), se=TRUE)$df
  df
})

#merge back into 1 dataframe
night.with.pred <-dplyr::bind_rows(data.1.with.pred ) %>%
  mutate(diel = 'night')

merged_predictions <- day.with.pred %>% add_row(night.with.pred)

# TEST
plotA <- merged_predictions %>%
  filter(diel == 'night') %>%
  ggplot(aes(x=-depth, y=pred.response, color=Ecological_Category, fill=Ecological_Category)) +
  geom_point(aes( y=sum_per_tot, shape=ESP), size = 2, alpha=0.6) + 
  #scale_shape_manual(values=c(21, 24))+
  scale_shape_manual(values=c(21, 8))+
  geom_line(  alpha= 1, size=1) +
  #geom_line(aes(y=pred.response -qt(0.975,df)*se))+
  #geom_line(aes(y=pred.response +qt(0.975,df)*se))+
  #geom_ribbon(data=subset(data.1.with.pred, 200 <= depth & depth <= 300), 
  #            aes(ymin=pred.response-qt(0.975,df)*se,ymax=pred.response+qt(0.975,df)*se), alpha=0.3,linetype='dashed') +
  geom_ribbon(aes(ymin=pred.response-qt(0.975,df)*se,ymax=pred.response+qt(0.975,df)*se), alpha=0.5,linetype = 0) +
  coord_flip()+
  scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  #xlim(-890,0)+
  #scale_x_continuous(minor_breaks= c(0, -50, -100, -200, -300, -400, -500, -600, -700, -800, -890))+
  #scale_x_continuous(minor_breaks= c(0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 890))+
  scale_x_continuous(breaks=c(-890,-800,-700,-600,-500,-400,-300,-200,-100,0), limits=c(-890,0))+
  labs(x="Depth (m)",y="Percent Total Reads")+
  theme_minimal() 
plotA





plotA <- merged_predictions %>%
  filter(diel == 'night') %>%
  ggplot(aes(x=-depth, y=pred.response, color=Ecological_Category, fill=Ecological_Category)) +
  geom_point(aes( y=sum_per_tot, shape=ESP), size = 2, alpha=0.6) + 
  #scale_shape_manual(values=c(21, 24))+
  scale_shape_manual(values=c(21, 8))+
  geom_line(  alpha= 1, size=1) +
  #geom_line(aes(y=pred.response -qt(0.975,df)*se))+
  #geom_line(aes(y=pred.response +qt(0.975,df)*se))+
  #geom_ribbon(data=subset(data.1.with.pred, 200 <= depth & depth <= 300), 
  #            aes(ymin=pred.response-qt(0.975,df)*se,ymax=pred.response+qt(0.975,df)*se), alpha=0.3,linetype='dashed') +
  geom_ribbon(aes(ymin=pred.response-qt(0.975,df)*se,ymax=pred.response+qt(0.975,df)*se), alpha=0.5,linetype = 0) +
  coord_flip()+
  scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_x_continuous(breaks=c(-890,-800,-700,-600,-500,-400,-300,-200,-100,0), limits=c(-890,10))+
  scale_y_continuous(breaks=c(0,20,40,60,80,100))+
  labs(x="Depth (m)",y="Percent Total Reads")+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  #scale_x_continuous(breaks=c(-890, -700, -500, -300, -100, 0)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=10,face="bold"),
    legend.key.height=grid::unit(0.5,"cm"),
    legend.key.width=grid::unit(0.5,"cm"),
    legend.title=element_text(colour=textcol,size=10,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12,colour=textcol),
    #axis.text.x=element_text(size=7,colour=textcol),
    axis.text.y=element_text(size=12,colour=textcol),
    axis.title.y = element_text(size=12),
    #plot.background=element_blank(),
    #panel.border=element_blank(),
    #panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())

plotA
filename = paste(directory, marker,'_EcolCat_LOESS_scatter_Nights5.png', sep='')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')
filename = paste(directory, marker,'_EcolCat_LOESS_scatter_Nights5.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

#DAY
plotA <- merged_predictions %>%
  filter(diel == 'day') %>%
  ggplot(aes(x=-depth, y=pred.response, color=Ecological_Category, fill=Ecological_Category)) +
  geom_point(aes( y=sum_per_tot, shape=ESP), size = 2, alpha=0.6) + 
  scale_shape_manual(values=c(21, 8))+
  geom_line(  alpha= 1, size=1) +
  geom_ribbon(aes(ymin=pred.response-qt(0.975,df)*se,ymax=pred.response+qt(0.975,df)*se), alpha=0.5,linetype = 0) +
  coord_flip()+
  scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_x_continuous(breaks=c(-890,-800,-700,-600,-500,-400,-300,-200,-100,0), limits=c(-890,10))+
  scale_y_continuous(breaks=c(0,20,40,60,80,100))+
  labs(x="Depth (m)",y="Percent Total Reads")+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
  theme(
    #legend
    legend.position="right",legend.direction="vertical",
    legend.text=element_text(colour=textcol,size=10,face="bold"),
    legend.key.height=grid::unit(0.5,"cm"),
    legend.key.width=grid::unit(0.5,"cm"),
    legend.title=element_text(colour=textcol,size=10,face="bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12,colour=textcol),
    axis.text.y=element_text(size=12,colour=textcol),
    axis.title.y = element_text(size=12),
    panel.grid.major = element_line(size = .25),
    plot.margin=margin(0.1,0.1,0.1,0.1,"cm"),
    plot.title=element_blank())

plotA
filename = paste(directory, marker,'_EcolCat_LOESS_scatter_Days5.png', sep='')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')
filename = paste(directory, marker,'_EcolCat_LOESS_scatter_Days5.svg', sep='')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')



### Plot by Ecological category, night and day together:
#basic:
plotA <- merged_predictions %>%
  #filter(diel == 'night') %>%
  ggplot(aes(x=depth, y=pred.response, color=diel, fill=diel, linestyle=diel)) +
  geom_point(aes( y=sum_per_tot), size = 2, alpha=0.3) + 
  geom_line(  alpha= 1, size=1) +
  facet_grid(Ecological_Category ~ ., margins=FALSE)

filename = paste(directory, marker,'_EcolCat_bydiel_lines5.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')

#highlighted between lines:
plotA <- merged_predictions %>%
  #filter(diel == 'night') %>%
  ggplot(aes(x=-depth, y=pred.response, color=diel, fill=diel, linestyle=diel)) +
  geom_point(aes( y=sum_per_tot), size = 2, alpha=0.3) + 
  geom_line(  alpha= 1, size=1) +
  #geom_polygon(aes(group = diel), alpha = 0.3)+
  geom_ribbon(aes(ymin=0,ymax=pred.response), alpha=0.3, linetype=0) +
  scale_fill_manual(values = c('chocolate1', 'royalblue3', 'darkgreen')) +
  scale_color_manual(values = c('chocolate1', 'royalblue3', 'darkgreen')) +
  facet_grid(Ecological_Category ~ ., margins=FALSE)+
  #facet_grid(. ~ Ecological_Category, margins=FALSE)+
  coord_flip()

plotA
filename = paste(directory, marker,'_EcolCat_bydiel_lines5_highlighted.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 8, width =5, units = 'in')






+
  #geom_line(aes(y=pred.response -qt(0.975,df)*se))+
  #geom_line(aes(y=pred.response +qt(0.975,df)*se))+
  #geom_ribbon(data=subset(data.1.with.pred, 200 <= depth & depth <= 300), 
  #            aes(ymin=pred.response-qt(0.975,df)*se,ymax=pred.response+qt(0.975,df)*se), alpha=0.3,linetype='dashed') +
  geom_ribbon(aes(ymin=pred.response-qt(0.975,df)*se,ymax=pred.response+qt(0.975,df)*se), alpha=0.5,linetype = 0) +
  #coord_flip()+
  scale_color_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  scale_fill_manual(values = c('blue3', 'darkorchid', 'deepskyblue','chartreuse')) +
  xlim(-890,0)+
  labs(x="Depth (m)",y="Percent Total Reads")+
  theme_minimal() +
  guides(fill=guide_legend(ncol=1)) +
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

plotA
filename = paste(directory, marker,'_EcolCat_LOESS_scatter_Nights6.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 5, width =5, units = 'in')




# ggplot(merged_predictions, aes(x=-depth, y=pred.response, color=Ecological_Category, fill=Ecological_Category, linetype=diel)) +
#   geom_point(aes( y=sum_per_tot), size = 3, alpha=0.3) + 
#   geom_line(  alpha= 1, size=2) +
#   geom_line(aes(y=pred.response -qt(0.975,df)*se),linetype='dashed')+
#   geom_line(aes(y=pred.response +qt(0.975,df)*se),linetype='dashed')+
#   #geom_ribbon(data=subset(data.1.with.pred, 200 <= depth & depth <= 300), 
#   #            aes(ymin=pred.response-qt(0.975,df)*se,ymax=pred.response+qt(0.975,df)*se), alpha=0.3,linetype='dashed') +
#   geom_ribbon(aes(ymin=pred.response-qt(0.975,df)*se,ymax=pred.response+qt(0.975,df)*se), alpha=0.3,linetype='dashed') +
#   coord_flip()
# 
# 
# filename = paste(directory, marker,'Ecolcats_LOESS_day_05.png', sep='')
# #print('Plot of top 20 Genus average by month:')
# print(filename)
# ggsave(filename,height = 8, width =10, units = 'in')

### Model in Tidyverse Ecol Cats together ---------------

test <- Ecol_data %>% 
  #filter(Ecological_Category == 'epipelagic') %>%
  filter(Ecological_Category %in% c('epipelagic', 'mesopelagic')) %>%
  left_join(meta %>% select(SampleID, diel, depth)) %>%
  filter(diel == 'day') %>%
  #filter(diel %in% c('night', 'day')) %>%
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

ggplot(data.1.with.pred, aes(x=-depth, y=pred.response, color=Ecological_Category)) +
  geom_point(aes( y=sum_per_tot), size = 3, alpha=0.3) + 
  geom_line( size=0.5, alpha= 1) +
  coord_flip()

filename = paste(directory, marker,'Ecolcats_LOESS_day_06.png', sep='')
#print('Plot of top 20 Genus average by month:')
print(filename)
ggsave(filename,height = 8, width =10, units = 'in')


### Model in Tidyverse: One EcolCat over day/night ---------------

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
  df$pred.response  <-stats::predict(stats::loess(sum_per_tot ~ depth, span = .50, data=df))
  df
})

#merge back into 1 dataframe
data.1.with.pred <-dplyr::bind_rows(data.1.with.pred )

ggplot(data.1.with.pred, aes(x=-depth, y=pred.response, color=diel)) +
  geom_point(aes(y=sum_per_tot), size = 3, alpha=0.3) + 
  geom_line( size=0.5, alpha= 1) +
  coord_flip()
