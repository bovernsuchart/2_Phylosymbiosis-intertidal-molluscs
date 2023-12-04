#Microbiome actual R thesis version - 01.11.2022 - subclass subclass effect
#Rarefied data - rarefied_feature_table_4714.txt"

## Load packages
library(ggplot2)
library(ggtext)
library(FSA)
library(phyloseq)
library(forcats)
library(reshape2)
library(dplyr) 
library(carData)
library(car)
library(multcomp)
library(permute)
library(lattice)
library(vegan)
library(performance)
library(cowplot)
library(vegan)

###Rarefied data-4714###############################################################################################################################################################################################################################################################################################################################################################################################################################################################################
#Import data 
#Metadata ###############################################################################################################################################################################################################################################################################################################
metadata <- read.table("metadata_phylosymbiosis.txt", sep = '\t', row.names = 1, header = T, strip.white = T) 
metadata <- subset(metadata, !sample.name == "TC3" & !sample.name == "TG1" &  !sample.name == "HE3" & !sample.name == "HG1" #removed these based on the rarefied data
                   & !species_allknown == "TH.Turbo.bruneus" & !species_allknown == "A.water.TH") #better to remove this species since it only have two samples)

#re-order the levels
levels(as.factor(metadata$species_allknown))
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
                                                "TC.Saccostrea.sp.TH", 
                                                "TD.Echinolittorina.malaccana.TH", 
                                                "TF.Mytilisepta.virgata.TH",          
                                                "TG.Tenguella.musiva",              
                                                "TJ.Echinolittorina.radiata.TJ" ,   
                                                "Z.Isognomon.nucleus"))
levels(as.factor(metadata$species_allknown))
# Setting some data as categorical factor:
metadata$species_allknown <- factor(metadata$species_allknown)

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

###Observed ASV ##############################################################################################################################################################################################################
## Plot and calculate statistics for species richness (in this case, ASV richness) 
ASVrichness <- apply(asv_table_rarefied4714_new  >0, 1, sum) 

#check if the order is right! in ASV richness and metadata
rich_meta_rarefied4714<- cbind(ASVrichness, metadata) #cbind is to combined side-side from column. 

#Check if all the data is correct after combined
View(rich_meta_rarefied4714)

summarystat_subclass_richness <- summarise(
  group_by(rich_meta_rarefied4714, location, subclass), # Grouping by the factor
  meanASVrichness = mean(ASVrichness, na.rm = TRUE),
  sdASVrichness  = sd(ASVrichness, na.rm = TRUE))
summarystat_subclass_richness
# Plot as ASVrichness boxplot - By day.time
ASVrichness_rarefied4714 <- ggplot(rich_meta_rarefied4714, aes(x=subclass , y=ASVrichness, fill=subclass )) +
                        geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
                        #scale_fill_discrete(labels = c("D1_0915", "D1_1230", "D1_1612", "D1_1740")) +
                        geom_jitter(color="black", size=2) +
                        #scale_x_discrete(name = "subclass", 
                                   #      labels = c("High" = "High",
                                    #                "Mid" = "Mid",
                                     #               "Low" = "Low")) +
                        scale_y_continuous(name = "ASV richness",  limits = c(0,800)) + 
                        theme_bw() +
                        theme(axis.text.y =element_text(size=11),
                              axis.title.y=element_text(size = 13),
                              legend.position = "none",
                              #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                              axis.text.x = element_blank(), 
                              axis.title.x =element_blank(),
                              strip.background = element_rect(fill = "white", colour = "black")
                              #,legend.text=element_text(size=11), 
                              #legend.title=element_text(size=13)
                              )  
ASVrichness_rarefied4714
                
ASVrichness_rarefied4714_wrap <- ASVrichness_rarefied4714 + facet_wrap(~location, drop = TRUE, scales = "free_x")
ASVrichness_rarefied4714_wrap 

#ggsave("ASVrichness_rarefied4714_finalised.pdf", width = 15, height = 6, device = "pdf", dpi = "print")

#All samples from both location ############################################################################################################################################################################################################################
#Create model - rarefied4714 data 
#Calculate statistics By subclass

asv_rarefied4714.lm <- lm(ASVrichness ~ subclass, rich_meta_rarefied4714) #create a model
anova(asv_rarefied4714.lm) 
#Result: 
#                 Df Sum Sq Mean Sq F value   Pr(>F)   
#subclass   4   80909   20227  1.3295 0.2678
#Residuals 69 1049799   15214       

#Normality test - rarefied4714 data 
par(mfrow=c(2,2))
plot(asv_rarefied4714.lm) #1. Some data point seem to have some outliers, but there Q-Q plot is quite good; 2. Homogeneity of variance seem to be met

asv_rarefied4714.lm.residuals <- residuals(object = asv_rarefied4714.lm)
shapiro.test(x = asv_rarefied4714.lm.residuals) #result: W = 0.89329, p-value = 1.307e-05  #violated

#Homogeneity of Variance - rarefied4714 data 
asv_rarefied4714.lm.levene <- leveneTest(asv_rarefied4714.lm, center =mean)
asv_rarefied4714.lm.levene   #Result: Passed: p-value = 0.005369 * # violated

##Kruskal-Wallis Method - rarefied4714 data  ############################################################################################################################################################################################################################
#Reference for this code: "Statistical_Analysis_Fig2_Strains_Nsources_Condition.R"
asv_rarefied4714.kruskal <-kruskal.test(ASVrichness ~ subclass , data = rich_meta_rarefied4714)
asv_rarefied4714.kruskal #Result: Kruskal-Wallis chi-squared = 4.982, df = 4, p-value = 0.2891
dunnTest(ASVrichness ~ subclass, data = rich_meta_rarefied4714, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#All samples from HK location  ############################################################################################################################################################################################################################
rich_meta_rarefied4714_HK <- subset(rich_meta_rarefied4714, location == "HK")#subset the metadata first

#Create model - rarefied4714 data & HK samples
#Calculate statistics By species
asv_rarefied4714_HK.lm <- lm(ASVrichness ~ subclass, rich_meta_rarefied4714_HK) #create a model
anova(asv_rarefied4714_HK.lm) 
#Result: 
#                 Df Sum Sq Mean Sq F value    Pr(>F)    
#subclass 4 192526   48132  2.8489 0.0373 *
#Residuals 37 625102   16895    

#Normality test - rarefied4714 data & HK samples
par(mfrow=c(2,2))
plot(asv_rarefied4714_HK.lm) 

asv_rarefied4714_HK.lm.residuals <- residuals(object = asv_rarefied4714_HK.lm)
shapiro.test(x = asv_rarefied4714_HK.lm.residuals) #result: W = 0.8886, p-value = 0.0006676 #Failed

#Homogeneity of Variance - rarefied4714 data & HK samples
asv_rarefied4714_HK.lm.levene <- leveneTest(asv_rarefied4714_HK.lm, center =mean)
asv_rarefied4714_HK.lm.levene   #Result: Passed: p-value = 5.356e-06  #Failed

##Kruskal-Wallis Method - rarefied4714 data & HK samples############################################################################################################################################################################################################################
#Reference for this code: "Statistical_Analysis_Fig2_Strains_Nsources_Condition.R"
asv_rarefied4714_HK.kruskal <-kruskal.test(ASVrichness ~ subclass,  data = rich_meta_rarefied4714_HK)
asv_rarefied4714_HK.kruskal #Result: Kruskal-Wallis chi-squared = 9.8015, df = 4, p-value = 0.04391
dunnTest(ASVrichness ~ subclass, data = rich_meta_rarefied4714_HK, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#All samples from TH location  ############################################################################################################################################################################################################################
rich_meta_rarefied4714_TH <- subset(rich_meta_rarefied4714, location == "TH")#subset the metadata first

#Create model - rarefied4714 data & TH samples
#Calculate statistics By species
asv_rarefied4714_TH.lm <- lm(ASVrichness ~ subclass, rich_meta_rarefied4714_TH) #create a model
anova(asv_rarefied4714_TH.lm) 
#Result: 
#                 Df Sum Sq Mean Sq F value    Pr(>F)    
#subclass  2  37862 18930.9   2.084 0.1427
#Residuals 29 263427  9083.7 

#Normality test - rarefied4714 data & TH samples
par(mfrow=c(2,2))
plot(asv_rarefied4714_TH.lm) 

asv_rarefied4714_TH.lm.residuals <- residuals(object = asv_rarefied4714_TH.lm)
shapiro.test(x = asv_rarefied4714_TH.lm.residuals) #result:W = 0.95614, p-value = 0.2147 #Passed

#Homogeneity of Variance - rarefied4714 data & TH samples
asv_rarefied4714_TH.lm.levene <- leveneTest(asv_rarefied4714_TH.lm, center =mean)
asv_rarefied4714_TH.lm.levene   #Result: Passed: p-value = 0.217 #Passed

##Kruskal-Wallis Method - rarefied4714 data & TH samples############################################################################################################################################################################################################################
#Reference for this code: "Statistical_Analysis_Fig2_Strains_Nsources_Condition.R"
asv_rarefied4714_TH.kruskal <-kruskal.test(ASVrichness ~ subclass, data = rich_meta_rarefied4714_TH)
asv_rarefied4714_TH.kruskal #Result:Kruskal-Wallis chi-squared = 3.5303, df = 2, p-value = 0.1712
dunnTest(ASVrichness ~ subclass, data = rich_meta_rarefied4714_TH, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

### Shannon Diversity ##############################################################################################################################################################################################################
?diversity() #from Vegan package
shannon <- diversity(asv_table_rarefied4714_new, index = "shannon")

shannon_rarefied4714_meta <- cbind(shannon, metadata)

summarystat_shannon <- summarise(
  group_by(shannon_rarefied4714_meta, subclass, location), # Grouping by the factor
  meanshannon = mean(shannon, na.rm = TRUE),
  sdshannon  = sd(shannon, na.rm = TRUE))
summarystat_shannon

# Plot as shannon boxplot 
shannon_rarefied4714_subclass  <- ggplot(shannon_rarefied4714_meta, aes(x=subclass, y=shannon, fill=subclass)) +
                                geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
                                #scale_fill_discrete(labels = c("D1_0915", "D1_1230", "D1_1612", "D1_1740", "D2_0636", "D2_0945")) +
                                geom_jitter(color="black", size=2) +
                               scale_x_discrete(name = "Subclass") + 
                                 #                labels = c("High" = "High",
                                   #                                  "Mid" = "Mid",
                                    #                                 "Low" = "Low")) +
                                scale_y_continuous(name = "Shannon-Wiener Index (H')\n", limits = c(0, 6), breaks = c(0,1,2,3,4,5,6)) +
                                theme_bw() +
                                theme(axis.text=element_text(size=11),
                                      axis.title=element_text(size = 13),
                                      legend.position = "none",
                                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                      strip.background = element_rect(fill = "white", colour = "black")
                                      #,legend.text=element_text(size=11), 
                                      #legend.title=element_text(size=13)
                                      )  
                            #+ annotate("text", x = 3.5, y = 5.75, label = expression("F"["(4,13)"]*" = 3.76 "*"; p-value < 0.5")  , size= 5)
shannon_rarefied4714_subclass_wrap <- shannon_rarefied4714_subclass + facet_wrap(~location, drop = TRUE, scales = "free_x")
shannon_rarefied4714_subclass_wrap
#ggsave("Alpha_diversity_shannon_rarefied4714_finalised_subclass.pdf", width = 15, height = 6, device = "pdf", dpi = "print")

#All samples from both location  ############################################################################################################################################################################################################################
#Create model - rarefied4714 data 
#Calculate statistics By species

#Create shannon model - rarefied4714 data
shannon_rarefied4714_meta.lm <- lm(shannon ~ subclass, shannon_rarefied4714_meta) #create a model
anova(shannon_rarefied4714_meta.lm) 
#Results: 
#                 Df Sum Sq Mean Sq F value    Pr(>F)    
#subclass  4  12.56 3.14003    3.24 0.01702 *
#Residuals 69  66.87 0.96914     

#Normality test - rarefied4714 data 
par(mfrow=c(2,2))
plot(shannon_rarefied4714_meta.lm) 

shannon_rarefied4714_meta.lm.residuals <- residuals(object = shannon_rarefied4714_meta.lm)
shapiro.test(x = shannon_rarefied4714_meta.lm.residuals) #result: Passed: W = 0.98431, p-value = 0.4912


#Homogeneity of Variance - rarefied4714 data
shannon_rarefied4714_meta.lm.levene <- leveneTest(shannon_rarefied4714_meta.lm, center =mean)
shannon_rarefied4714_meta.lm.levene  #Result: p-value =2.268e-05 #failed

##Kruskal-Wallis Method - rarefied4714 data (shannon) All samples from both location ############################################################################################################################################################################################################################
#Reference for this code: "Statistical_Analysis_Fig2_Strains_Nsources_Condition.R"
shannon_rarefied4714.kruskal <-kruskal.test(shannon ~ subclass, data = shannon_rarefied4714_meta)
shannon_rarefied4714.kruskal #Result:Kruskal-Wallis chi-squared = 12.475, df = 4, p-value = 0.01415
dunnTest(shannon ~ subclass, data = shannon_rarefied4714_meta, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#All samples from HK location  ############################################################################################################################################################################################################################
#Create model - rarefied4714 data & HK samples
#Calculate statistics By species
shannon_rarefied4714_meta_HK <- subset(shannon_rarefied4714_meta, location == "HK")#subset the metadata first

#Create shannon model - rarefied4714 data & HK samples
shannon_rarefied4714_meta_HK.lm <- lm(shannon ~ subclass, shannon_rarefied4714_meta_HK) #create a model
anova(shannon_rarefied4714_meta_HK.lm) 

#Results: 
#                 Df Sum Sq Mean Sq F value    Pr(>F)    
#subclass  4  3.661 0.91528  0.9291 0.4576
#Residuals 37 36.449 0.98510 

#Normality test - rarefied4714 data & HK samples
par(mfrow=c(2,2))
plot(shannon_rarefied4714_meta_HK.lm) 

shannon_rarefied4714_meta_HK.lm.residuals <- residuals(object = shannon_rarefied4714_meta_HK.lm)
shapiro.test(x = shannon_rarefied4714_meta_HK.lm.residuals) #result: Passed: W = 0.98244, p-value = 0.7551

#Homogeneity of Variance - rarefied4714 data  & HK samples
shannon_rarefied4714_meta_HK.lm.levene <- leveneTest(shannon_rarefied4714_meta_HK.lm, center =mean)
shannon_rarefied4714_meta_HK.lm.levene  #Result: Passed: p-value = 0.003869 ** #Failed

#Kruskal-Wallis Method - rarefied4714 data (shannon) & HK samples ############################################################################################################################################################################################################################
#Reference for this code: "Statistical_Analysis_Fig2_Strains_Nsources_Condition.R"
shannon_rarefied4714.kruskal <-kruskal.test(shannon ~ subclass, data = shannon_rarefied4714_meta_HK)
shannon_rarefied4714.kruskal #Result: Kruskal-Wallis chi-squared = 4.5121, df = 4, p-value = 0.3411
dunnTest(shannon ~ subclass, data = shannon_rarefied4714_meta_HK, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

#All samples from TH location  ############################################################################################################################################################################################################################
#Create model - rarefied4714 data & TH samples
#Calculate statistics By species
shannon_rarefied4714_meta_TH <- subset(shannon_rarefied4714_meta, location == "TH")#subset the metadata first

#Create shannon model - rarefied4714 data & TH samples
shannon_rarefied4714_meta_TH.lm <- lm(shannon ~ subclass, shannon_rarefied4714_meta_TH) #create a model
anova(shannon_rarefied4714_meta_TH.lm) 
#Results: 
#                    Df Sum Sq Mean Sq F value    Pr(>F)    
#subclass  2 10.611  5.3053  5.4906 0.009503 **
#Residuals 29 28.022  0.9663          

#Normality test - rarefied4714 data & TH samples
par(mfrow=c(2,2))
plot(shannon_rarefied4714_meta_TH.lm) 

shannon_rarefied4714_meta_TH.lm.residuals <- residuals(object = shannon_rarefied4714_meta_TH.lm)
shapiro.test(x = shannon_rarefied4714_meta_TH.lm.residuals) #result:W = 0.95334, p-value = 0.1792#passed


#Homogeneity of Variance - rarefied4714 data  & TH samples
shannon_rarefied4714_meta_TH.lm.levene <- leveneTest(shannon_rarefied4714_meta_TH.lm, center =mean)
shannon_rarefied4714_meta_TH.lm.levene  #Result: p-value =  0.02454 #Failed

#Kruskal-Wallis Method - rarefied4714 data (shannon) & TH samples ############################################################################################################################################################################################################################
#Reference for this code: "Statistical_Analysis_Fig2_Strains_Nsources_Condition.R"
shannon_rarefied4714.kruskal <-kruskal.test(shannon ~ subclass, data = shannon_rarefied4714_meta_TH)
shannon_rarefied4714.kruskal #Result: Kruskal-Wallis chi-squared = 9.2066, df = 2, p-value = 0.01002
dunnTest(shannon ~ subclass, data = shannon_rarefied4714_meta_TH, method="bh",two.sided = TRUE) #This test carries out the multiple comparison of all groups

## Combine plots into one figure ----
library(cowplot) #method following Walke et al. (2021) in the mantelprocrustes R script
cowplot::plot_grid(ASVrichness_rarefied4714_wrap,
                   shannon_rarefied4714_subclass_wrap,
                   nrow = 2, ncol = 1,
                   scale = c(1,1),
                   rel_heights = c(1,1.39),
                   labels = c("a","b"))

ggsave("Alpha_diversity_combined_subclass.pdf", width = 12, height = 8, device = "pdf", dpi = "print")

