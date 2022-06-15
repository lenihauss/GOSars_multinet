#### R Script CVOO size distribution plots
# Author: H.Hauss (hhauss@geomar.de) & Nele Weigt
# Read in Ecotaxa files, calculate abundance in size bins
# plot abundance vs. size

library(data.table)
library(tidyr)
library(stringr)
library(raster)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(mltools)
library(Hmisc)

## path on Nele's PC
setwd("C:/Users/Nele/Desktop/6. Semester/BACHELORARBEIT/Ecotaxa_Daten/Datensätze") 

## path on Leni's PC
## working directory
setwd("V:/Daten/Students/NeleWeigt/ecotaxa_export") 
file_list <- list.files(pattern="*.tsv") # create list of all .tsv files in that folder
print(file_list)
##read in all ecotaxa tsvs and 
data_raw <- do.call(rbind,lapply(file_list,read.csv, header=TRUE, quote = '"', sep = "\t"))[ ,c('sample_id','object_annotation_category','object_annotation_hierarchy', 'object_date', 'object_time', 'object_depth_min', 'object_depth_max', 'object_area', 'sample_volconc')]
data_raw <- dplyr::filter(data_raw, !grepl('m105_mn01', sample_id))
memory.limit()
## filter living
data <- data_raw[-grep("not-living", data_raw$object_annotation_hierarchy),] 

rm(data_raw)

## split string to assign cruise, haul and net numbers
data$cruise  <-   unlist(lapply(strsplit(as.character(data$sample_id), "_"), '[[', 1))
data$haul_id <-   unlist(lapply(strsplit(as.character(data$sample_id), "_"), '[[', 2))
data$net_id  <-   unlist(lapply(strsplit(as.character(data$sample_id), "_"), '[[', 3))

## assign variables
data$date <- data$object_date
data$time <- data$object_time
data$depth_min <- data$object_depth_min
data$depth_max <- data$object_depth_max
data$depth_mid <- (data$object_depth_max +data$object_depth_min)/2
data$spec_id <- data$object_annotation_category

## neue Variable
data$object_area_mm2 <- data$object_area* 0.00011236

## convert from ml to m3 for concentrated volume
data$sample_volconc_m3 <- (data$sample_volconc)#/1000000 

# Define the size bins
data["size_bin"] <- bin_data (data$object_area_mm, 
                              bins = c(0, 0.064, 0.128, 0.256, 0.512, 
                                       1, 2, 4, 8, 16, 32, 64, 128, Inf), 
                              binType = "explicit")

## aggregate objects into samples (layers) to yield abundance
layers <- aggregate(object_area_mm2 ~ (cruise+sample_id+haul_id+net_id+date+
                                       time+depth_min+depth_max+depth_mid+
                                       sample_volconc_m3+size_bin), 
                                       data, length)

layers$count <- layers$object_area_mm2

layers$abundance_m3 <- layers$count /layers$sample_volconc_m3
write.table(layers, file = "layers.txt")
layers_allfractions <- aggregate(abundance_m3 ~ (cruise+haul_id+net_id+date+time+
                                                depth_min+depth_max+depth_mid+
                                                size_bin), 
                                                layers, sum)

layers_allfractions$size_bin <- as.character(layers_allfractions$size_bin)

layers_allfractions$mid_size[layers_allfractions$size_bin=="[0, 0.064)"]    <- 0.032
layers_allfractions$mid_size[layers_allfractions$size_bin=="[0.064, 0.128)"]<- 0.096
layers_allfractions$mid_size[layers_allfractions$size_bin=="[0.128, 0.256)"]<- 0.192
layers_allfractions$mid_size[layers_allfractions$size_bin=="[0.256, 0.512)"]<- 0.384
layers_allfractions$mid_size[layers_allfractions$size_bin=="[0.512, 1)"]    <- 0.756
layers_allfractions$mid_size[layers_allfractions$size_bin=="[1, 2)"]        <- 1.5
layers_allfractions$mid_size[layers_allfractions$size_bin=="[2, 4)"]        <- 3
layers_allfractions$mid_size[layers_allfractions$size_bin=="[4, 8)"]        <- 6
layers_allfractions$mid_size[layers_allfractions$size_bin=="[8, 16)"]       <- 12
layers_allfractions$mid_size[layers_allfractions$size_bin=="[16, 32)"]      <- 24
layers_allfractions$mid_size[layers_allfractions$size_bin=="[32, 64)"]      <- 48
layers_allfractions$mid_size[layers_allfractions$size_bin=="[64, 128)"]     <- 96
layers_allfractions$mid_size[layers_allfractions$size_bin=="[128, Inf]"]    <- 192

str(layers_allfractions)

Size <- subset(layers_allfractions, mid_size>0.2)

## fit linear regression separately for each depth
## explore is there could be a day/night pattern or differences between cruises
## color by depth, make color scale light for surface and dark for deeper

## Plot
ggplot(Size, aes(x=mid_size, y=abundance_m3, color=depth_mid)) +
  geom_point(size=2) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2')

# 1. Regression durch alles
cor(Size$mid_size, Size$abundance_m3)

Reg1 <- lm(abundance_m3 ~ mid_size, data=Size)
summary(Reg1)    
coef(Reg1)
attributes(Reg1)
confint(Reg1)

# add regression line
abline(Reg1) ##klappt nicht
confint(Reg1, level = 0.99)
summary(Reg1)
anova(Reg1)

ggplot(Size, aes(x=mid_size, y=abundance_m3, color=depth_mid)) +
  geom_point(size=2) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  geom_smooth(method='lm')


## Subsets für Tiefe

#50m
depth_50  <- subset(Size, depth_mid== (50), select = c(cruise, haul_id, net_id, date, time, 
                                                                     depth_mid,
                                                                     abundance_m3, mid_size))
p1 <- ggplot(depth_50, aes(x=mid_size, y=abundance_m3)) +
      geom_point(size=2) +
      scale_x_continuous(trans = 'log2') +
      scale_y_continuous(trans = 'log2') +
      geom_smooth(method='lm')

#150m
depth_150 <- subset(Size, depth_mid==(150), select = c(cruise, haul_id, net_id, date, time, 
                                                                     depth_mid,
                                                                     abundance_m3, mid_size))
p2 <- ggplot(depth_150, aes(x=mid_size, y=abundance_m3)) +
      geom_point(size=2) +
      scale_x_continuous(trans = 'log2') +
      scale_y_continuous(trans = 'log2') +
      geom_smooth(method='lm')

#250m
depth_250 <- subset(Size, depth_mid==(250), select = c(cruise, haul_id, net_id, date, time, 
                                                                      depth_mid,
                                                                      abundance_m3, mid_size))
p3 <- ggplot(depth_250, aes(x=mid_size, y=abundance_m3)) +
      geom_point(size=2) +
      scale_x_continuous(trans = 'log2') +
      scale_y_continuous(trans = 'log2') +
      geom_smooth(method='lm')

#450m
depth_450 <- subset(Size, depth_mid==(450), select = c(cruise, haul_id, net_id, date, time, 
                                                                      depth_mid, 
                                                                      abundance_m3, mid_size))
p4 <- ggplot(depth_450, aes(x=mid_size, y=abundance_m3)) +
      geom_point(size=2) +
      scale_x_continuous(trans = 'log2') +
      scale_y_continuous(trans = 'log2') +
      geom_smooth(method='lm')

#700m
depth_700 <- subset(Size, depth_mid==(700), select = c(cruise, haul_id, net_id, date, time, 
                                                                      depth_mid, 
                                                                      abundance_m3, mid_size))
p5 <- ggplot(depth_700, aes(x=mid_size, y=abundance_m3)) +
      geom_point(size=2) +
      scale_x_continuous(trans = 'log2') +
      scale_y_continuous(trans = 'log2') +
      geom_smooth(method='lm')

#800m
depth_800 <- subset(Size, depth_mid==(800), select = c(cruise, haul_id, net_id, date, time, 
                                                                      depth_mid, 
                                                                      abundance_m3, mid_size))
p6 <- ggplot(depth_800, aes(x=mid_size, y=abundance_m3)) +
      geom_point(size=2) +
      scale_x_continuous(trans = 'log2') +
      scale_y_continuous(trans = 'log2') +
      geom_smooth(method='lm')

library(ggpubr)
all <- ggarrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2)
all
ggsave("V:/Daten/Students/NeleWeigt/size_distribution.pdf", width = 9, height = 9, bg = "transparent")