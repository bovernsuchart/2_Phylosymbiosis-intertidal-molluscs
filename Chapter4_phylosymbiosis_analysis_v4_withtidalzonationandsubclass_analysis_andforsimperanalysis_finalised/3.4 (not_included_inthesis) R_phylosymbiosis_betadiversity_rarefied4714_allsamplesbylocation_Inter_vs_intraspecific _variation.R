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


#HK samples #########################################################################################################################################################################################################################################################################################################################################################################################################################
#Filter samples by HK
metadata_HK <- subset(metadata, location == "HK")#subset the metadata first
asv_table_wide_rarefied4714_HK <- subset(asv_table_rarefied4714_new,rownames(asv_table_rarefied4714_new)  %in% (rownames(metadata_HK) )) 
View(asv_table_wide_rarefied4714_HK )

#Get the intra vs. inter-specific data 
#Step 1: Calculate the distance matrix
dist.mat.HK <- as.matrix(vegdist(asv_table_wide_rarefied4714_HK, method = "bray"))

#Step 2 - arrange the data into long list - there are two methods - O'brien et al. (2020) & Pat Schloss's youtube method

#Method 1: O'Brien et al. (2020)
dist.melt.HK <- data.frame(X=colnames(dist.mat.HK)[col(dist.mat.HK)], Y=rownames(dist.mat.HK)[row(dist.mat.HK)], ASV_dist=c(dist.mat.HK))

#Method 2: Pat Schloss's method
#Making into a nice column-based data matrix - See Pat Schloss's video: https://www.youtube.com/watch?v=EXNOgmUyPfY&t=3s
#sample <- dist.mat %>%
#  as.matrix() %>%
#  as_tibble(rownames = "X") %>%
#  pivot_longer(-X, names_to = "Y", values_to = "ASV_dist")
#sample

#Step 3: Create a new column and assign intra- or inter-specific variation - Following O'Brien et al. (2020) R script
#I will with work with O'Brien method dist matrix file - #under what does each of the component in scripts mean

dist.melt.HK$intra.sp <- ifelse(dist.melt.HK$X %in% c("HB1", "HB2", "HB3", "HB4", "HB5") & dist.melt.HK$Y %in% c("HB1", "HB2", "HB3", "HB4", "HB5") |
                               dist.melt.HK$X %in% c("HC1", "HC2", "HC3", "HC4", "HC5") & dist.melt.HK$Y %in% c("HC1", "HC2", "HC3", "HC4", "HC5") |
                               dist.melt.HK$X %in% c("HD1", "HD2", "HD3", "HD4", "HD5") & dist.melt.HK$Y %in% c("HD1", "HD2", "HD3", "HD4", "HD5") |
                               dist.melt.HK$X %in% c("HE1", "HE2", "HE4", "HE5") & dist.melt.HK$Y %in% c("HE1", "HE2", "HE4", "HE5") |
                               dist.melt.HK$X %in% c("HF1", "HF2", "HF3", "HF4", "HF5") & dist.melt.HK$Y %in% c("HF1", "HF2", "HF3", "HF4", "HF5") |                  
                               dist.melt.HK$X %in% c("HG2", "HG3", "HG4", "HG5") & dist.melt.HK$Y %in% c("HG2", "HG3", "HG4", "HG5") |                  
                               dist.melt.HK$X %in% c("HH1", "HH2", "HH3", "HH4") & dist.melt.HK$Y %in% c("HH1", "HH2", "HH3", "HH4") |          
                               dist.melt.HK$X %in% c("HJ1", "HJ2", "HJ3", "HJ4", "HJ5") & dist.melt.HK$Y %in% c("HJ1", "HJ2", "HJ3", "HJ4", "HJ5") |
                               dist.melt.HK$X %in% c("HK1", "HK2", "HK3", "HK4", "HK5") & dist.melt.HK$Y %in% c("HK1", "HK2", "HK3", "HK4", "HK5"),
                               "intraspecific", "interspecific") 
                             
# Remove 0's and arcsine tranformation to normalise data - O'Brien method 
dist.melt.HK <- dist.melt.HK[!(dist.melt.HK$ASV_dist ==0),] 

#dist.melt.HK$arcsin.t <- asin(sqrt(dist.melt.HK$ASV_dist))*180/pi #O'Brien arcsin transformation reference _ I don't plan to do it 


# T-test to see if there is a significant difference between interspecific and intraspecific variation
t_test <- t.test(ASV_dist~intra.sp, dist.melt.HK, var.equal = FALSE)
t_test #Result: t = 21.547, df = 156.26, p-value < 2.2e-16

p <- ggplot(dist.melt.HK, aes(x=intra.sp  , y=ASV_dist, fill=intra.sp)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) + 
  ylim(0, NA) +
  ylab("Bray-curtis dissimilarity") + xlab("Variation") +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "") + 
  theme_bw() +
  theme(axis.text=element_text(size=15, colour = "black"),axis.title=element_text(size = 18),
        legend.text=element_text(size=15)) 
p + annotate("text", x = 1, y = 0.2, label = "t=21.547, p < 0.001") 
 
#TH samples #########################################################################################################################################################################################################################################################################################################################################################################################################################
#Filter samples by TH
metadata_TH <- subset(metadata, location == "TH")#subset the metadata first
asv_table_wide_rarefied4714_TH <- subset(asv_table_rarefied4714_new,rownames(asv_table_rarefied4714_new)  %in% (rownames(metadata_TH) )) 

#Get the intra vs. inter-specific data 
#Step 1: Calculate the distance matrix
dist.mat.TH <- as.matrix(vegdist(asv_table_wide_rarefied4714_TH, method = "bray"))

#Step 2 - arrange the data into long list - there are two methods - O'brien et al. (2020) & Pat Schloss's youtube method

#Method 1: O'Brien et al. (2020)
dist.melt.TH <- data.frame(X=colnames(dist.mat.TH)[col(dist.mat.TH)], Y=rownames(dist.mat.TH)[row(dist.mat.TH)], ASV_dist=c(dist.mat.TH))

#Method 2: Pat Schloss's method
#Making into a nice column-based data matrix - See Pat Schloss's video: https://www.youtube.com/watch?v=EXNOgmUyPfY&t=3s
#sample <- dist.mat %>%
#  as.matrix() %>%
#  as_tibble(rownames = "X") %>%
#  pivot_longer(-X, names_to = "Y", values_to = "ASV_dist")
#sample

#Step 3: Create a new column and assign intra- or inter-specific variation - Following O'Brien et al. (2020) R script
#I will with work with O'Brien method dist matrix file - #under what does each of the component in scripts mean

dist.melt.TH$intra.sp <- ifelse(dist.melt.TH$X %in% c("TB1", "TB2", "TB3", "TB4", "TB5") & dist.melt.TH$Y %in% c("TB1", "TB2", "TB3", "TB4", "TB5") |                    
                                  dist.melt.TH$X %in% c("TC1", "TC2", "TC4", "TC5") & dist.melt.TH$Y %in% c("TC1", "TC2", "TC4", "TC5") |
                                  dist.melt.TH$X %in% c("TD1", "TD2", "TD3", "TD4", "TD5") & dist.melt.TH$Y %in% c("TD1", "TD2", "TD3", "TD4", "TD5") |                     
                                  dist.melt.TH$X %in% c("TF1", "TF2", "TF3", "TF4", "TF5") & dist.melt.TH$Y %in% c("TF1", "TF2", "TF3", "TF4", "TF5") |                     
                                  dist.melt.TH$X %in% c("TG2", "TG3", "TG4", "TG5") & dist.melt.TH$Y %in% c("TG2", "TG3", "TG4", "TG5") |
                                  dist.melt.TH$X %in% c("TJ1", "TJ2", "TJ3", "TJ4", "TJ5") & dist.melt.TH$Y %in% c("TJ1", "TJ2", "TJ3", "TJ4", "TJ5") |                     
                                  dist.melt.TH$X %in% c("BRM1", "BRM2", "BRM3", "BRM4", "BRM5") & dist.melt.TH$Y %in% c("BRM1", "BRM2", "BRM3", "BRM4", "BRM5"),                     
                                "intraspecific", "interspecific")

# Remove 0's and arcsine tranformation to normalise data - O'Brien method 
dist.melt.TH <- dist.melt.TH[!(dist.melt.TH$ASV_dist ==0),]

#dist.melt.TH$arcsin.t <- asin(sqrt(dist.melt.TH$ASV_dist))*180/pi  #O'Brien arcsin transformation reference _ I don't plan to do it 

# T-test to see if there is a significant difference between interspecific and intraspecific variation
t_test.TH <- t.test(ASV_dist~intra.sp, dist.melt.TH, var.equal = FALSE)
t_test.TH #Result: t = 18.427, df = 117.69, p-value < 2.2e-16

p <- ggplot(dist.melt.TH, aes(x=intra.sp  , y=ASV_dist, fill=intra.sp)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) + 
  ylim(0, NA) +
  ylab("Bray-curtis dissimilarity") + xlab("Variation") +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "") + 
  theme_bw() +
  theme(axis.text=element_text(size=15, colour = "black"),axis.title=element_text(size = 18),
        legend.text=element_text(size=15)) 
p + annotate("text", x = 1, y = 0.2, label = "t=18.427, p < 0.001") 

#Combined figures 
#create a new column with location naming
dist.melt.HK$location <- "HK"
dist.melt.TH$location <- "TH"

#combined the data
dist.df <- rbind(dist.melt.HK,dist.melt.TH)

p <- ggplot(dist.df, aes(x=location, y=ASV_dist, fill=intra.sp)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) + 
  ylim(0, NA) +
  ylab("Bray-curtis dissimilarity") + xlab("Location") +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "") + 
  theme_bw() +
  theme(axis.text=element_text(size=15, colour = "black"),axis.title=element_text(size = 18),
        legend.text=element_text(size=15)) 
p + annotate("text", x = 1, y = 0.1, label = "t=21.547\n p < 0.001") + annotate("text", x = 2, y = 0.1, label = "t=18.427\n p < 0.001") 

