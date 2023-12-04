#Microbiome actual R thesis version - 26.08.2022
#Rarefied data - rarefied_feature_table_4714.txt" 

## Load packages
library(permute)
library(dplyr)
library(lattice)
library(vegan)
library(cluster)
library(randomcoloR)
library(pairwiseAdonis)
library(ape)
library(ggplot2)
library(gg3D) #gg3D is a package created to extend ggplot2 to produce 3D plots

##Rarefied data- 4714###################################################################################################################################################################################################################################################################################################
#Import data 
#Metadata ###############################################################################################################################################################################################################################################################################################################
metadata <- read.table("metadata_phylosymbiosis.txt", sep = '\t', row.names = 1, header = T, strip.white = T) 
metadata <- subset(metadata, !sample.name == "TC3" & !sample.name == "TG1" &  !sample.name == "HE3" & !sample.name == "HG1" #removed these based on the rarefied data
                   & !species_allknown == "TH.Turbo.bruneus"
                   & !species_allknown == "A.water.TH") #better to remove this species since it only have two samples)

#Setting some data as categorical factor & re-arrangement
metadata$species_allknown <- factor(metadata$species_allknown, #re-order the levels to make it nicer in comparison
                                    levels = c( "HB.Mytilisepta.virgata.HK", 
                                                "HC.Cellana.toreuma.HK",
                                                "HD.Siphonaria.laciniosa", 
                                                "HE.Saccostrea.cucullata.HK",
                                                "HF.Lunella.granulata" ,
                                                "HG.Echinolittorina.malaccana.HK",
                                                "HH.Reishia.clavigera",            
                                                "HJ.Barbatia.virescens.HJ",
                                                "HK.Monodonta.labio_HK",
                                                "TB.Cellana.toreuma.TH",
                                                "TC.Saccostrea.mordax.TH", 
                                                "TD.Echinolittorina.malaccana.TH", 
                                                "TF.Mytilisepta.virgata.TH",          
                                                "TG.Tenguella.musiva",              
                                                "TJ.Echinolittorina.radiata.TJ" ,   
                                                "Z.Isognomon.nucleus",
                                                "A.water.TH"))
asv_rarefied4714_metadata$feeding_strategies <- as.factor(asv_rarefied4714_metadata$feeding_strategies)

#ASV ###############################################################################################################################################################################################################################################################################################################
asv_table_rarefied4714<- read.delim("rarefied_feature_table_4714.txt", sep = '\t', header = T, row.names = 1, strip.white = T)
asv_table_rarefied4714_transpose <- t(asv_table_rarefied4714)      #t() here means transpose the matrix

#the %in% operator returns a logical vector that tells us whether or not there is a match between the first object and the second and subset() does only for those that present in the first.
#This checking is important because sometimes we want to might rarefied to sampling depth that might removed some samples. So we want to know what are those samples. 
#always check what samples that was subset
asv_table_rarefied4714_transpose_relevantsamples <- subset(asv_table_rarefied4714_transpose, (rownames(asv_table_rarefied4714_transpose) %in% rownames(metadata))) 

#Check is all the samples have the same number of sequences (or sampling depth)
apply(asv_table_rarefied4714_transpose_relevantsamples, 1, sum)

#Remove any taxa that is not present in the any of the samples first
asv_table_transpose_relevantsamples_relavanttaxa <- asv_table_rarefied4714_transpose_relevantsamples[,(colSums(asv_table_rarefied4714_transpose_relevantsamples) > 0 )]

#Check again - is all the samples have the same number of sequences (or sampling depth)
apply(asv_table_transpose_relevantsamples_relavanttaxa, 1, sum)

#make a new name to the ASV table
asv_table_rarefied4714_new <- asv_table_transpose_relevantsamples_relavanttaxa

#order the metadata  following the ASV table
metadata <- metadata[match(rownames(asv_table_rarefied4714_new ), rownames(metadata)),]

##########################################################################################################################################################################################################################################################################################################################################################################################################################
## Rarefied4714 data - Statistical analysis based on Bray-Curtis dist matrix ----
#PERMANOVA 

#Verification steps - that only taxa present in either one of the samples
#Calculate column totals
total <- colSums(asv_table_rarefied4714_new)
asv_w_total <- rbind(asv_table_rarefied4714_new, total)

#re-order columns by total - check the last column
asv_sort <- asv_w_total[,order(-asv_w_total[which(rownames(asv_w_total) == 'total'), ])] #Here essentially it is re-ordering the total by descending order
asv_sort <- asv_sort[-75,] #remove the "total" row

#Know which one to combine, here I didn't use any of asv_sort because I didn't removed any samples
asv_rarefied4714_metadata <- cbind(asv_table_rarefied4714_new , metadata)

#Check if all the data is correct after combined
#View(asv_rarefied4714_metadata)

#1. PERMANOVA#######################################################################################################################################################################################################################################################################################
#Note: 1. always check if the column is right! Make sure not to include empty count taxa 
#2. The result always changes because of the random permutation and so the significance return will always be different. 
adonis.microbiome_species <- adonis2(asv_rarefied4714_metadata[,1:6468] ~ asv_rarefied4714_metadata$species_allknown, method = "bray")
adonis.microbiome_species

#                                           Df SumOfSqs      R2      F Pr(>F)    
#asv_rarefied4714_metadata$species_allknown 15   22.159 0.64395 6.9931  0.001 ***
#Residual                                   58   12.252 0.35605                  
#Total                                      73   34.411 1.00000           


adonis.microbiome_species <- adonis2(asv_rarefied4714_metadata[,1:6468] ~ asv_rarefied4714_metadata$location, method = "bray")
adonis.microbiome_species

#                                           Df SumOfSqs      R2      F Pr(>F)    
#asv_rarefied4714_metadata$location          1    1.676 0.0487 3.6863  0.001 ***
#Residual                                   72   32.735 0.9513                  
#Total                                      73   34.411 1.0000     

adonis.microbiome_location_species <- adonis2(asv_rarefied4714_metadata[,1:6468] ~ asv_rarefied4714_metadata$location + asv_rarefied4714_metadata$species_allknown, 
                                             method = "bray")
adonis.microbiome_location_species 

#                                           Df SumOfSqs      R2      F Pr(>F)    
#asv_rarefied4714_metadata$location          1    1.676 0.04870 7.9338  0.001 ***
#asv_rarefied4714_metadata$species_allknown 14   20.483 0.59524 6.9259  0.001 ***
#Residual                                   58   12.252 0.35605                  
#Total                                      73   34.411 1.00000  

adonis.microbiome_species_location <- adonis2(asv_rarefied4714_metadata[,1:6468] ~ asv_rarefied4714_metadata$species_allknown + asv_rarefied4714_metadata$location, 
                                              method = "bray")
adonis.microbiome_species_location

#                                           Df SumOfSqs      R2      F Pr(>F)    
#asv_rarefied4714_metadata$species_allknown 15   22.159 0.64395 6.9931  0.001 ***
#Residual                                   58   12.252 0.35605                  
#Total                                      73   34.411 1.00000 

adonis.microbiome_feeding <- adonis2(asv_rarefied4714_metadata[,1:6468] ~ asv_rarefied4714_metadata$feeding_strategies, 
                                              method = "bray")
adonis.microbiome_feeding
#                                           Df SumOfSqs      R2      F Pr(>F)    
#asv_rarefied4714_metadata$feeding_strategies  2    3.124 0.09079 3.5449  0.001 ***
#Residual                                     71   31.287 0.90921                  
#Total                                        73   34.411 1.00000     

adonis.microbiome_feeding_location <- adonis2(asv_rarefied4714_metadata[,1:6468] ~ asv_rarefied4714_metadata$feeding_strategies*asv_rarefied4714_metadata$location, 
                                     method = "bray")
adonis.microbiome_feeding_location
#                                                                                Df SumOfSqs      R2      F Pr(>F)    
#asv_rarefied4714_metadata$feeding_strategies                                     2    3.124 0.09079 3.9297  0.001 ***
#asv_rarefied4714_metadata$location                                               1    1.689 0.04908 4.2486  0.001 ***
#asv_rarefied4714_metadata$feeding_strategies:asv_rarefied4714_metadata$location  2    2.567 0.07460 3.2289  0.001 ***
#Residual                                                                        68   27.031 0.78553                  
#Total   

adonis.microbiome_location_feeding <- adonis2(asv_rarefied4714_metadata[,1:6468] ~ asv_rarefied4714_metadata$location + asv_rarefied4714_metadata$feeding_strategies, 
                                              method = "bray")
adonis.microbiome_location_feeding 
#                                           Df SumOfSqs      R2      F Pr(>F)    
#asv_rarefied4714_metadata$location            1    1.676 0.04870 3.9637  0.001 ***
#asv_rarefied4714_metadata$feeding_strategies  2    3.137 0.09117 3.7097  0.001 ***
#Residual                                     70   29.598 0.86013                  
#Total                                        73   34.411 1.00000    


