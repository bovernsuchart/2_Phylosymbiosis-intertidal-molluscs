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
#raw data - at family level
#import
asv_table_f_30 <- read.csv("raw_taxa_fam.csv", sep = ',', header = T, row.names = 1, strip.white = T) # Family table created using qiime2
asv_table_f_30 <- subset(asv_table_f_30, !sample.name == "TC3" & !sample.name == "TG1" &  !sample.name == "HE3" & !sample.name == "HG1"  #removed these based on the rarefied data
                         & !species == "TH.Turbo.bruneus" & !species == "A.water.TH") #better to remove this species since it only have two samples)

levels(as.factor(asv_table_f_30$species))
asv_table_f_30$species<- factor(asv_table_f_30$species)                       
levels(asv_table_f_30$species) #Perfect

#subset metadata from the .csv file - need to do this because when we merged the table it has to be in the same number of roles and column
meta_data <- asv_table_f_30[458:465]
meta_data$species <- as.factor(meta_data$species)
meta_data$species

asv_table_f2_30 <- as.matrix(asv_table_f_30[, -c(458:465)]) # Remove metadata columns from the asv_table_f

#SIMPER analysis - with family level rawdata - SIMPER uses Bray-Curtis method
simp.testc <- simper(asv_table_f2_30, meta_data$species)
sum.simpc <- summary(simp.testc)

# Pull out the comparison you want to look at if you want - the three main species present in both locations
Cellana <- sum.simpc$HC.Cellana.toreuma.HK_TB.Cellana.toreuma.TH
Echinolittorina <- sum.simpc$HG.Echinolittorina.trochoides.HK_TD.Echinolittorina.trochoides.TH
Mytilisepta <- sum.simpc$HB.Brachidontes.variabilis_TF.Unknown_bivalves_TF


View(Cellana)

write.csv(Cellana , file = "CellanaHKTH_simp.csv")
write.csv(Echinolittorina, file = "EchinolittorinaHKTH_simp.csv")
write.csv(Mytilisepta, file = "MytiliseptaHKTH_simp.csv")


##########################################################################################################################################################################################################################################################################################################################################################################################################################
#Another method #################################################################################################################################################################################################################################
#import data (important to import with read.table because the first column can be set as a vector - you can try with read.delim() to import and you will notice that at later section you can't work with this format): 

asv <- read.table("raw_feature_table.txt", sep = '\t', row.names = 1, header = T, strip.white = T)    #"\t" mean tab-delimited; " " mean space-delimited; "," for comma-delimited -normally used in .csv file
map <- read.table("metadata_phylosymbiosis.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
tax <- read.table("taxonomy_silva.tsv", sep = '\t', row.names = 1, header = T, strip.white = T)
# Separate taxonomy into different columns #you can ignore the error. It just mean that the empty value is treated as NA
tax <- separate(tax, Taxon, c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ";" , remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")
tax$Confidence <- NULL

######################################################################################################################################d#######################################################################################################################################d#######################################################################################################################################
#Prefer this method (Not fully tried yet) Following the naming system in  O'Brien et al.(2021)) - "1. R_daynight_Rarefaction.R" 
#Start with renaming first - Replace all NA's with the taxa indicator only or like in the kingdom case - change the name
#Always check with the unique names first for each column and check again at the end

tax$Kingdom[is.na(tax$Kingdom)] <- " k__unknown" 
tax$Kingdom[tax$Kingdom == " k__"] <- " k__unknown"
tax$Phylum[is.na(tax$Phylum)] <- " p__unknown"
tax$Phylum[tax$Phylum == " p__"] <- " p__unknown"
tax$Class[is.na(tax$Class)] <- " c__unknown"
tax$Class[tax$Class == " c__"] <- " c__unknown"
tax$Order[is.na(tax$Order)] <- " o__unknown"
tax$Order[tax$Order == " o__"] <- " o__unknown"
tax$Family[is.na(tax$Family)] <- " f__unknown"
tax$Family[tax$Family == " f__"] <- " f__unknown"
tax$Genus[is.na(tax$Genus)] <- " g__unknown"
tax$Genus[tax$Genus == " g__"] <- " g__unknown"
tax$Species[is.na(tax$Species)] <- " s__unknown"
tax$Species[tax$Genus == " s__"] <- " s__unknown"
unique(tax$Kingdom)
unique(tax$Phylum)
unique(tax$Class)
unique(tax$Order)
unique(tax$Family)
unique(tax$Genus)
unique(tax$Species)

#OPTIONAL - I can even remove the taxonomic label following the method in O'Brien et al.(2021): "Relative_abundance.R" file
#Don't do this is better for this dataset since there are many with unknown names even at the family levels
#tax$Domain <- gsub("d__", "", tax$Domain)
#tax$Kingdom <- gsub("k__", "", tax$Kingdom)
#tax$Phylum <- gsub("p__", "", tax$Phylum)
#tax$Class <- gsub("c__", "", tax$Class)
#tax$Order <- gsub("o__", "", tax$Order)
#tax$Family <- gsub("f__", "", tax$Family)
#tax$Genus <- gsub("g__", "", tax$Genus)
#tax$Species <- gsub("s__", "", tax$Species)
#unique(tax$Genus)
#unique(tax$Family)

# 'Phyloseq-ize' the data - following the Savary et al. (2021) way of naming
otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)
phy.all

family_data <- aggregate_taxa(phy.all, "Family") #Stothart & Newman (2021) use an unrarefied table (ANCOM-BC has bulit-in normalization procedures)
tax_table(family_data)


#To proceed - follow the R script from - Escalas et al. (2021) - I haven't tried it.


out <- df_res %>% group_by(y_var) %>% 
  dplyr::summarize(avg_p = round(mean(p_value, na.rm = TRUE),3), 
                   sd_p = round(sd(p_value, na.rm = TRUE),3), 
                   avg_f = round(mean(statistic, na.rm = TRUE),1),
                   sd_f = round(sd(statistic, na.rm = TRUE),1),
                   num_sign_0.05 = sum(p_value < 0.05, na.rm = TRUE),
                   num_sign_0.01 = sum(p_value < 0.01, na.rm = TRUE)) %>% 
  data.frame()%>% arrange(dplyr::desc(avg_f))






