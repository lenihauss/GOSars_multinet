#### R Script GO Sars distribution plots
# Author: H.Hauss (hhauss@geomar.de)
# Read ecotaxa output, calculate individual biomass, aggregate to depth 
#layers, plot vertical distribution

library(ggplot2)
library(dplyr)
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
 
layers <- do.call("rbind", by(data, data[1:5], with, 
                              data.frame(sample_id = sample_id[1], 
                                         net_id = net_id[1],
                                         depth_min = depth_min[1], 
                                         depth_max = depth_max[1],
                                         depth_mid = depth_mid[1],
                                         sample_volconc = sample_volconc[1],
                                         spec_id = spec_id[1], 
                                         count = length(biomass_ug), biomass_ug = sum(biomass_ug))))

#filter for main groups
layers <- filter(layers, spec_id %in%  c("Calanoida", "Calanus", "Metridinidae", "Centropages", "Oithona", "Oncaea", "Euphausiacea", "Chaetognatha", "Actinopterygii"))

layers$abundance_m3 <- layers$count/layers$sample_volconc
layers$biomass_ug_m3 <- layers$biomass_ug/layers$sample_volconc

abundance <- aggregate(abundance_m3 ~ (net_id+spec_id+depth_min+depth_max+depth_mid), layers, sum)
biomass   <- aggregate(biomass_ug_m3 ~ (net_id+spec_id+depth_min+depth_max+depth_mid), layers, sum)

##Plots:

## vertical distribution barplot
#abundance
p1 <-ggplot(data=abundance, aes(x=depth_mid, y=abundance_m3 , fill=spec_id, width=(depth_max-depth_min))) +
  geom_col() +
  coord_flip() + 
  scale_x_reverse(limits=c(250,0), breaks=seq(0,250,50)) +
  labs(title = "Abundance",x = "Depth (m)", y = "Abundance (ind/m^3")+
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##biomass
p2 <- ggplot(data=biomass, aes(x=depth_mid, y=(biomass_ug_m3)/1000 , fill=spec_id, width=(depth_max-depth_min))) +
  geom_col() +
  coord_flip() + 
  scale_x_reverse(limits=c(250,0), breaks=seq(0,250,50)) +
  labs(title = "Biomass",x = "Depth (m)", y = "Biomass (mg/m^3)", fill="Group")+
  theme(legend.position=c(.8,.75), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


all <- ggarrange(p1, p2, ncol=2, nrow=1)
all
ggsave("Abundance_Biomass.png", width = 9, height = 9, bg = "transparent")<<<<<<< HEAD
