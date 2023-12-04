#Microbiome actual R thesis version -  27.08.2022
#Rarefied data - rarefied_feature_table_4714.txt" - Comparison made for all samples from both location

## Load packages
library(permute)
library(dplyr)
library(lattice)
library(vegan)
library(tidyr)
library(cluster)
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
                   & !species_allknown == "A.water.TH" ) #better to remove this species since it only have two samples)

# Setting some data as categorical factor:
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
                                                "Z.Isognomon.nucleus"))

metadata$group <- factor(metadata$group)
metadata$zone <- factor(metadata$zone)

#ASV ###############################################################################################################################################################################################################################################################################################################
asv_table_rarefied4714<- read.delim("rarefied_feature_table_4714.txt", sep = '\t', header = T, row.names = 1, strip.white = T)
asv_table_rarefied4714_transpose <- t(asv_table_rarefied4714)      #t() here means transpose the matrix
#the %in% operator returns a logical vector that tells us whether or not there is a match between the first object and the second and subset() does only for those that present in the first.
#This checking is important because sometimes we want to might rarefied to sampling depth that might removed some samples. So we want to know what are those samples. 

asv_table_rarefied4714_transpose_relevantsamples <- subset(asv_table_rarefied4714_transpose, (rownames(asv_table_rarefied4714_transpose) %in% rownames(metadata))) 

#Check is all the samples have the same number of sequences (or sampling depth)
apply(asv_table_rarefied4714_transpose_relevantsamples, 1, sum)

#Remove any taxa that is not present in the any of the samples first
asv_table_transpose_relevantsamples_relavanttaxa <- asv_table_rarefied4714_transpose_relevantsamples[,(colSums(asv_table_rarefied4714_transpose_relevantsamples) > 0 )]

#Check again - is all the samples have the same number of sequences (or sampling depth)
apply(asv_table_transpose_relevantsamples_relavanttaxa, 1, sum)

#make a new name to the ASV table
asv_table_rarefied4714_new <- asv_table_transpose_relevantsamples_relavanttaxa

#order the ASV table following the metadata
asv_table_rarefied4714_new <- asv_table_rarefied4714_new[match(rownames(metadata), rownames(asv_table_rarefied4714_new)),]

#Get the intra vs. inter-specific data 
#Step 1: Calculate the distance matrix
dist.mat <- as.matrix(vegdist(asv_table_rarefied4714_new, method = "bray"))

#Step 2 - arrange the data into long list - there are two methods - O'brien et al. (2020) & Pat Schloss's youtube method

#Method 1: O'Brien et al. (2020)
dist.melt <- data.frame(X=colnames(dist.mat)[col(dist.mat)], Y=rownames(dist.mat)[row(dist.mat)], ASV_dist=c(dist.mat))

#Method 2: Pat Schloss's method
#Making into a nice column-based data matrix - See Pat Schloss's video: https://www.youtube.com/watch?v=EXNOgmUyPfY&t=3s
#sample <- dist.mat %>%
#  as.matrix() %>%
#  as_tibble(rownames = "X") %>%
#  pivot_longer(-X, names_to = "Y", values_to = "ASV_dist")
#sample

#Step 3: Create a new column and assign intra- or inter-specific variation - Following O'Brien et al. (2020) R script
#I will with work with O'Brien method dist matrix file - #under what does each of the component in scripts mean
rownames(asv_table_rarefied4714_new)

dist.melt$intra.sp <- ifelse(dist.melt$X %in% c("HB1", "HB2", "HB3", "HB4", "HB5") & dist.melt$Y %in% c("HB1", "HB2", "HB3", "HB4", "HB5") |
                             dist.melt$X %in% c("HC1", "HC2", "HC3", "HC4", "HC5") & dist.melt$Y %in% c("HC1", "HC2", "HC3", "HC4", "HC5") |
                             dist.melt$X %in% c("HD1", "HD2", "HD3", "HD4", "HD5") & dist.melt$Y %in% c("HD1", "HD2", "HD3", "HD4", "HD5") |
                             dist.melt$X %in% c("HE1", "HE2", "HE4", "HE5") & dist.melt$Y %in% c("HE1", "HE2", "HE4", "HE5") |
                             dist.melt$X %in% c("HF1", "HF2", "HF3", "HF4", "HF5") & dist.melt$Y %in% c("HF1", "HF2", "HF3", "HF4", "HF5") |                  
                             dist.melt$X %in% c("HG2", "HG3", "HG4", "HG5") & dist.melt$Y %in% c("HG2", "HG3", "HG4", "HG5") |                  
                             dist.melt$X %in% c("HH1", "HH2", "HH3", "HH4") & dist.melt$Y %in% c("HH1", "HH2", "HH3", "HH4") |          
                             dist.melt$X %in% c("HJ1", "HJ2", "HJ3", "HJ4", "HJ5") & dist.melt$Y %in% c("HJ1", "HJ2", "HJ3", "HJ4", "HJ5") |
                             dist.melt$X %in% c("HK1", "HK2", "HK3", "HK4", "HK5") & dist.melt$Y %in% c("HK1", "HK2", "HK3", "HK4", "HK5") |
                             dist.melt$X %in% c("TB1", "TB2", "TB3", "TB4", "TB5") & dist.melt$Y %in% c("TB1", "TB2", "TB3", "TB4", "TB5") |                    
                             dist.melt$X %in% c("TC1", "TC2", "TC4", "TC5") & dist.melt$Y %in% c("TC1", "TC2", "TC4", "TC5") |
                             dist.melt$X %in% c("TD1", "TD2", "TD3", "TD4", "TD5") & dist.melt$Y %in% c("TD1", "TD2", "TD3", "TD4", "TD5") |                     
                             dist.melt$X %in% c("TF1", "TF2", "TF3", "TF4", "TF5") & dist.melt$Y %in% c("TF1", "TF2", "TF3", "TF4", "TF5") |                     
                             dist.melt$X %in% c("TG2", "TG3", "TG4", "TG5") & dist.melt$Y %in% c("TG2", "TG3", "TG4", "TG5") |
                             dist.melt$X %in% c("TJ1", "TJ2", "TJ3", "TJ4", "TJ5") & dist.melt$Y %in% c("TJ1", "TJ2", "TJ3", "TJ4", "TJ5") |                     
                             dist.melt$X %in% c("BRM1", "BRM2", "BRM3", "BRM4", "BRM5") & dist.melt$Y %in% c("BRM1", "BRM2", "BRM3", "BRM4", "BRM5"),                     
                            "intraspecific", "interspecific")


# Remove 0's and arcsine tranformation to normalise data - O'Brien method 
dist.melt <- dist.melt[!(dist.melt$ASV_dist ==0),]

#dist.melt$arcsin.t <- asin(sqrt(dist.melt$ASV_dist))*180/pi - O'Brien arcsin transformation reference _ I don't plan to do it 

# T-test to see if there is a significant difference between interspecific and intraspecific variation
t_test <- t.test(ASV_dist~intra.sp, dist.melt, var.equal = FALSE)
t_test #RESULT: t = 28.581, df = 274.24, p-value < 2.2e-16

p <- ggplot(dist.melt, aes(x=intra.sp  , y=ASV_dist, fill=intra.sp)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) + 
  ylim(0, NA) +
  ylab("Bray-curtis dissimilarity (arcsine transformed)\n") + xlab("Variation") +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "") + 
  theme_bw() +
  theme(axis.text=element_text(size=15, colour = "black"),axis.title=element_text(size = 18),
        legend.text=element_text(size=15)) 
p + annotate("text", x = 1, y = 0.25, label = "t=28.581, p < 0.001") 


