## SIMPER analysis 
## Script based on O'Brien et al. (2020) paper - the one about phylosymbiosis
##Note: The SIMPER analysis calculates the contribution of each species (%) to the dissimilarity between each two group

library(vegan)
library(phyloseq)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(rbiom)
library(ranacapa)
library(microbiome)
###################################################################################################################################################################################################################
#rarefied data 4714 - at family level
#import
asv_table_rarefied4714_f <- read.csv("rarefied4714_familylevel.csv", sep = ',', header = T, row.names = 1, strip.white = T) # Family table created using qiime2
asv_table_rarefied4714_f <- subset(asv_table_rarefied4714_f, !sample.name == "TC3" & !sample.name == "TG1" &  !sample.name == "HE3" & !sample.name == "HG1"  #removed these based on the rarefied data
                         & !species == "TH.Turbo.bruneus" & !species == "A.water.TH") #better to remove this species since it only have two samples)

levels(as.factor(asv_table_rarefied4714_f$species_allknown))
asv_table_rarefied4714_f$species_allknown<- factor(asv_table_rarefied4714_f$species_allknown)                       
levels(asv_table_rarefied4714_f$species_allknown) #Perfect

#subset metadata from the .csv file - need to do this because when we merged the table it has to be in the same number of roles and column
meta_data <- asv_table_rarefied4714_f[425:439]
meta_data$species_allknown <- as.factor(meta_data$species_allknown)
meta_data$species_allknown

asv_table_rarefied4714_f2 <- as.matrix(asv_table_rarefied4714_f[, -c(425:439)]) # Remove metadata columns from the asv_table_f

#SIMPER analysis - with family level rarefied 4714 - SIMPER uses Bray-Curtis method
simper.test_rarefied4714_f <- simper(asv_table_rarefied4714_f2 , meta_data$species_allknown)
sum.simper_rarefied4714_f  <- summary(simper.test_rarefied4714_f )

# Pull out the comparison you want to look at if you want - the three main species present in both locations
?SIMPER()
Cellana_rarefied4714_f <- sum.simper_rarefied4714_f$HC.Cellana.toreuma.HK_TB.Cellana.toreuma.TH
Echinolittorina_rarefied4714_f <- sum.simper_rarefied4714_f$HG.Echinolittorina.malaccana.HK_TD.Echinolittorina.malaccana.TH
Mytilisepta_rarefied4714_f <- sum.simper_rarefied4714_f$HB.Mytilisepta.virgata.HK_TF.Mytilisepta.virgata.TH

View(Cellana_rarefied4714_f)
View(Echinolittorina_rarefied4714_f)
View(Mytilisepta_rarefied4714_f)

write.csv(Cellana_rarefied4714_f , file = "Cellana_rarefied4714_f_HKTH_simper.csv")
write.csv(Echinolittorina_rarefied4714_f, file = "Echinolittorina_rarefied4714_f_HKTH_simper.csv")
write.csv(Mytilisepta_rarefied4714_f, file = "Mytilisepta_rarefied4714_f_simper.csv")

###################################################################################################################################################################################################################
#rarefied data 4714 - at genus level
#import
asv_table_rarefied4714_g <- read.csv("rarefied4714_genuslevel.csv", sep = ',', header = T, row.names = 1, strip.white = T) # Family table created using qiime2
asv_table_rarefied4714_g <- subset(asv_table_rarefied4714_g, !sample.name == "TC3" & !sample.name == "TG1" &  !sample.name == "HE3" & !sample.name == "HG1"  #removed these based on the rarefied data
                                   & !species == "TH.Turbo.bruneus" & !species == "A.water.TH") #better to remove this species since it only have two samples)

levels(as.factor(asv_table_rarefied4714_g$species_allknown))
asv_table_rarefied4714_g$species_allknown<- factor(asv_table_rarefied4714_g$species_allknown)                       
levels(asv_table_rarefied4714_g$species_allknown) #Perfect

#subset metadata from the .csv file - need to do this because when we merged the table it has to be in the same number of roles and column
meta_data <- asv_table_rarefied4714_g[842:856]
meta_data$species_allknown <- as.factor(meta_data$species_allknown)
meta_data$species_allknown

asv_table_rarefied4714_g2 <- as.matrix(asv_table_rarefied4714_g[, -c(842:856)]) # Remove metadata columns from the asv_table_g

#SIMPER analysis - with family level rarefied 4714 - SIMPER uses Bray-Curtis method
simper.test_rarefied4714_g <- simper(asv_table_rarefied4714_g2 , meta_data$species_allknown)
sum.simper_rarefied4714_g  <- summary(simper.test_rarefied4714_g )

# Pull out the comparison you want to look at if you want - the three main species present in both locations
#?SIMPER()
Cellana_rarefied4714_g <- sum.simper_rarefied4714_g$HC.Cellana.toreuma.HK_TB.Cellana.toreuma.TH
Echinolittorina_rarefied4714_g <- sum.simper_rarefied4714_g$HG.Echinolittorina.malaccana.HK_TD.Echinolittorina.malaccana.TH
Mytilisepta_rarefied4714_g <- sum.simper_rarefied4714_g$HB.Mytilisepta.virgata.HK_TF.Mytilisepta.virgata.TH

View(Cellana_rarefied4714_g)
View(Echinolittorina_rarefied4714_g)
View(Mytilisepta_rarefied4714_g)

write.csv(Cellana_rarefied4714_g , file = "Cellana_rarefied4714_g_HKTH_simper.csv")
write.csv(Echinolittorina_rarefied4714_g, file = "Echinolittorina_rarefied4714_g_HKTH_simper.csv")
write.csv(Mytilisepta_rarefied4714_g, file = "Mytilisepta_rarefied4714_g_simper.csv")





