#### R Script CVOO distribution plots
# Author: H.Hauss (hhauss@geomar.de)
# Read ecotaxa output, calculate individual biomass, aggregate to depth 
#layers, plot

library(stringr)
library(raster)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)

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
 

layers <- aggregate(biomass_ug ~ (sample_id+net_id+spec_id+depth_min+depth_max+depth_mid+
                                 sample_volconc), data, length) 
layers$count <- layers$spec_id
layers$abundance_m3 <- layers$count /layers$sample_volconc

biomass <- aggregate(biomass_ug ~ (sample_id+sample_volconc), data, sum)
biomass$biomass_m3 <- biomass$biomass /biomass$sample_volconc
str(biomass)

layers_allfractions <- aggregate(abundance_m3 ~ (cruise+haul_id+net_id+date+time+depth_min+depth_max+depth_mid), layers, sum)

biomass <- aggregate(biomass_m3 ~ (cruise+haul_id+net_id+date+time+depth_min+depth_max+depth_mid), biomass, sum)

##Plots:

## October_2012
cvoo_night_20121024 <- subset(biomass, cruise =="msm22" & haul_id == "mn01")
cvoo_day_20121024<- subset(biomass, cruise =="msm22" & haul_id == "mn02")

p1 <- ggplot(data=cvoo_day_20121024, aes(x=depth_mid, y=biomass_m3 , width=(depth_max-depth_min)))+
  coord_flip() + 
  scale_y_continuous(limits=c(-60000, 60000), breaks=seq(-60000, 60000, 20000))+    
  scale_x_reverse(limits=c(1000,0), breaks=seq(0,1000,200)) +
  geom_hline(aes(yintercept=0)) +
  geom_bar(data=cvoo_night_20121024, aes(y=biomass_m3*(-1)), fill='black', stat="identity") + 
  geom_bar(data=cvoo_day_20121024, aes(y=biomass_m3),fill='grey', stat="identity")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "October_2012",x = "Depth (m)", y = "")


## November_2012
cvoo_night_20121024_1 <- subset(biomass, cruise =="msm22" & haul_id == "mn34")
cvoo_day_20121024_1<- subset(biomass, cruise =="msm22" & haul_id == "mn35")

p2 <- ggplot(data=cvoo_day_20121024_1, aes(x=depth_mid, y=biomass_m3 , width=(depth_max-depth_min)))+
  coord_flip() + 
  scale_y_continuous(limits=c(-60000, 60000), breaks=seq(-60000, 60000, 20000))+    
  scale_x_reverse(limits=c(1000,0), breaks=seq(0,1000,200)) +
  geom_hline(aes(yintercept=0)) +
  geom_bar(data=cvoo_night_20121024_1, aes(y=biomass_m3*(-1)), fill='black', stat="identity") + 
  geom_bar(data=cvoo_day_20121024_1, aes(y=biomass_m3),fill='grey', stat="identity")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "November_2012",x = "", y = "")


## May_2013
cvoo_lrg_night_052013 <- subset(biomass, cruise =="m097" & haul_id == "mn01")
cvoo_lrg_day_052013 <- subset(biomass, cruise =="m097" & haul_id == "mn02")

p3 <- ggplot(data=cvoo_lrg_day_052013, aes(x=depth_mid, y=biomass_m3 , width=(depth_max-depth_min)))+
  coord_flip() + 
  scale_y_continuous(limits=c(-60000, 60000), breaks=seq(-60000, 60000, 20000))+    
  scale_x_reverse(limits=c(1000,0), breaks=seq(0,1000,200)) +
  geom_hline(aes(yintercept=0)) +
  geom_bar(data=cvoo_lrg_night_052013, aes(y=biomass_m3*(-1)), fill='black', stat="identity") + 
  geom_bar(data=cvoo_lrg_day_052013, aes(y=biomass_m3),fill='grey', stat="identity")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "May_2013",x = "", y = "")


## March_2014
cvoo_lrg_night_032014 <- subset(biomass, cruise =="m105" & haul_id == "mn06")
cvoo_lrg_day_032014 <- subset(biomass, cruise =="m105" & haul_id == "mn05")

p4 <- ggplot(data=cvoo_lrg_day_032014, aes(x=depth_mid, y=biomass_m3 , width=(depth_max-depth_min)))+
  coord_flip() + 
  scale_y_continuous(limits=c(-60000, 60000), breaks=seq(-60000, 60000, 20000))+    
  scale_x_reverse(limits=c(1000,0), breaks=seq(0,1000,200)) +
  geom_hline(aes(yintercept=0)) +
  geom_bar(data=cvoo_lrg_night_032014, aes(y=biomass_m3*(-1)), fill='black', stat="identity") + 
  geom_bar(data=cvoo_lrg_day_032014, aes(y=biomass_m3),fill='grey', stat="identity")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "March_2014",x = "Depth (m)", y = "")


## April_2014
cvoo_lrg_night_042014 <- subset(biomass, cruise =="m106" & haul_id == "mn01")
cvoo_lrg_day_042014 <-subset(biomass, cruise =="m106" & haul_id == "mn02")

p5 <- ggplot(data=cvoo_lrg_night_042014, aes(x=depth_mid, y=biomass_m3 , width=(depth_max-depth_min)))+
  coord_flip() + 
  scale_y_continuous(limits=c(-60000, 60000), breaks=seq(-60000, 60000, 20000))+    
  scale_x_reverse(limits=c(1000,0), breaks=seq(0,1000,200)) +
  geom_hline(aes(yintercept=0)) +
  geom_bar(data=cvoo_lrg_night_042014, aes(y=biomass_m3*(-1)), fill='black', stat="identity") + 
  geom_bar(data=cvoo_lrg_day_042014, aes(y=biomass_m3),fill='grey', stat="identity")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "April_2014",x = "", y = "")


## September_2015
cvoo_lrg_night_092015 <- subset(biomass, cruise =="m119" & haul_id == "mn01")
cvoo_lrg_day_092015 <- subset(biomass, cruise =="m119" & haul_id == "mn03")

p6 <- ggplot(data=cvoo_lrg_day_092015, aes(x=depth_mid, y=biomass_m3 , width=(depth_max-depth_min)))+
  coord_flip() + 
  scale_y_continuous(limits=c(-60000, 60000), breaks=seq(-60000, 60000, 20000))+    
  scale_x_reverse(limits=c(1000,0), breaks=seq(0,1000,200)) +
  geom_hline(aes(yintercept=0)) +
  geom_bar(data=cvoo_lrg_night_092015, aes(y=biomass_m3*(-1)), fill='black', stat="identity") + 
  geom_bar(data=cvoo_lrg_day_092015, aes(y=biomass_m3),fill='grey', stat="identity")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "September_2015",x = "", y = "")


## August_2016
cvoo_lrg_night_082016 <- subset(biomass, cruise =="m130" & haul_id == "mn01")
cvoo_lrg_day_082016 <- subset(biomass, cruise =="m130" & haul_id == "mn02")

p7 <- ggplot(data=cvoo_lrg_night_082016, aes(x=depth_mid, y=biomass_m3 , width=(depth_max-depth_min)))+
  coord_flip() + 
  scale_y_continuous(limits=c(-60000, 60000), breaks=seq(-60000, 60000, 20000))+    
  scale_x_reverse(limits=c(1000,0), breaks=seq(0,1000,200)) +
  geom_hline(aes(yintercept=0)) +
  geom_bar(data=cvoo_lrg_day_082016, aes(y=biomass_m3*(-1)), fill='black', stat="identity") + 
  geom_bar(data=cvoo_lrg_night_082016, aes(y=biomass_m3),fill='grey', stat="identity")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "August_2016",x = "Depth (m)", y = bquote('Biomass ('*µg / m^3*')'))

##February_2018
cvoo_lrg_night_022018 <- subset(biomass, cruise =="pos520" & haul_id == "mn08")
cvoo_lrg_day_022018 <- subset(biomass, cruise =="pos520" & haul_id == "mn09")

p8 <-ggplot(data=cvoo_lrg_night_022018, aes(x=depth_mid, y=biomass_m3 , width=(depth_max-depth_min)))+
  coord_flip() + 
  scale_y_continuous(limits=c(-60000, 60000), breaks=seq(-60000, 60000, 20000))+    
  scale_x_reverse(limits=c(1000,0), breaks=seq(0,1000,200)) +
  geom_hline(aes(yintercept=0)) +
  geom_bar(data=cvoo_lrg_day_022018, aes(y=biomass_m3*(-1)), fill='black', stat="identity") + 
  geom_bar(data=cvoo_lrg_night_022018, aes(y=biomass_m3),fill='grey', stat="identity")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "February_2018",x = "", y = bquote('Biomass ('*µg / m^3*')'))

##February_2019
cvoo_lrg_night_022019 <- subset(biomass, cruise =="pos532" & haul_id == "mn05")
cvoo_lrg_day_022019 <- subset(biomass, cruise =="pos532" & haul_id == "mn06")

p9 <-ggplot(data=cvoo_lrg_night_022019, aes(x=depth_mid, y=biomass_m3 , width=(depth_max-depth_min)))+
  coord_flip() + 
  scale_y_continuous(limits=c(-60000, 60000), breaks=seq(-60000, 60000, 20000))+    
  scale_x_reverse(limits=c(1000,0), breaks=seq(0,1000,200)) +
  geom_hline(aes(yintercept=0)) +
  geom_bar(data=cvoo_lrg_day_022019, aes(y=biomass_m3*(-1)), fill='black', stat="identity") + 
  geom_bar(data=cvoo_lrg_night_022019, aes(y=biomass_m3),fill='grey', stat="identity")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "February_2019",x = "", y = bquote('Biomass ('*µg / m^3*')'))

grid.arrange(p1,p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, top = "Over-all Biomass")
 