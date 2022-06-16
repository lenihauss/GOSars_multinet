#### R Script CVOO size distribution plots
# Author: H.Hauss (hhauss@geomar.de)
# Read in Ecotaxa files, plot size distribution of different taxa
library(ggplot2)
library(ggpubr)
## working directory
setwd("V:/Daten/Cruises/GOSars_multinet") 

## read ecotaxa files and join data-sets lrg, med und sml
file_list <- list.files(pattern="*.tsv") # create list of all .tsv files in that folder
print(file_list)
##read in all ecotaxa tsvs and 
data_raw <- do.call(rbind,lapply(file_list,read.csv, header=TRUE, quote = '"', sep = "\t"))[ ,c('sample_id','object_annotation_category','object_annotation_hierarchy', 'object_depth_min', 'object_depth_max', 'object_area', 'sample_volconc')]

## filter living
data <-data_raw[-grep("not-living", data_raw$object_annotation_hierarchy),] 
rm(data_raw)

## split the string
data$sample_id <-as.character(data$sample_id)  
str(data)   ## check format of variables (numeric, character, factor, etc.)

## split string to assign net numbers
data$net_id  <- unlist(lapply(strsplit(as.character(data $sample_id), "_"), '[[', 4))

## rename some variables
data$depth_min <- data $object_depth_min
data$depth_max <- data $object_depth_max
data$depth_mid <- (data $object_depth_max +data $object_depth_min)/2
data$spec_id <- data $object_annotation_category
data$annotation_hierarchy<- data $object_annotation_hierarchy
#convert area from pixel to square mm
data$area_mm2 <- data $object_area* 0.00011236

## check in Lehette & Hernandez Leon if these are the correct values for the taxon!
##to be changed! dummy values!!
data$biomass_ug <- with(data, 
                        ifelse(grepl("Calanoida", annotation_hierarchy),45.25*data$area_mm2^1.59,
                               ifelse(grepl("Annelida", annotation_hierarchy),43.38*data $area_mm2^1.54,
                                      ifelse(grepl("Chaetognatha", annotation_hierarchy),23.45*data $area_mm2^1.19,
                                             ifelse(grepl("Gastropoda", annotation_hierarchy),43.38*data $area_mm2^1.54,
                                                    44.78*data$area_mm2^1.56)))))
str(data)
# subset by species
calanus <- data[which(data$spec_id=="Calanus"),]
metridia <- data[which(data$spec_id=="Metridinidae"),]
oithona <- data[which(data$spec_id=="Oithona"),]
oncaea <- data[which(data$spec_id=="Oncaea"),]
krill <- data[which(data$spec_id=="Euphausiacea"),]
chaetognatha <- data[which(data$spec_id=="Chaetognatha"),]

# Basic histogram

p1<-ggplot(calanus, aes(x=biomass_ug)) + 
  geom_histogram(color="black", fill="white")+
  labs(title = "Calanus",x = "Biomass (µg)", y = "Count")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2<-ggplot(krill, aes(x=biomass_ug)) + 
  geom_histogram(color="black", fill="white")+
  labs(title = "Krill",x = "Biomass (µg)", y = "Count")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p3<-ggplot(metridia, aes(x=biomass_ug)) + 
  geom_histogram(color="black", fill="white")+
  labs(title = "Metridinidae",x = "Biomass (µg)", y = "Count")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p4<-ggplot(oithona, aes(x=biomass_ug)) + 
  geom_histogram(color="black", fill="white")+
  labs(title = "Oithonidae",x = "Biomass (µg)", y = "Count")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p5<-ggplot(oithona, aes(x=biomass_ug)) + 
  geom_histogram(color="black", fill="white")+
  labs(title = "Oncaeidae",x = "Biomass (µg)", y = "Count")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p6<-ggplot(chaetognatha, aes(x=biomass_ug)) + 
  geom_histogram(color="black", fill="white")+
  labs(title = "Chaetognatha",x = "Biomass (µg)", y = "Count")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

all <- ggarrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2)
all
ggsave("size_distribution.pdf", width = 9, height = 9, bg = "transparent")