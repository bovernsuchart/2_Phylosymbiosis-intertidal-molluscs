#Microbiome actual R thesis version - 11.09.2022 species versus location
#Rarefied data - rarefied_feature_table_4714.txt" - Comparison made for all samples from both location

## Load packages
library(permute)
library(dplyr)
library(lattice)
library(vegan)
library(cluster)
library(pairwiseAdonis)
library(ape)
library(ggplot2)
library(gg3D) #gg3D is a package created to extend ggplot2 to produce 3D plots

##Rarefied data- 4714###################################################################################################################################################################################################################################################################################################
#Import data 
#Metadata ###############################################################################################################################################################################################################################################################################################################
metadata <- read.table("metadata_phylosymbiosis.txt", sep = '\t', row.names = 1, header = T, strip.white = T) 
metadata <- subset(metadata, !sample.name == "TC3" & !sample.name == "TG1" &  !sample.name == "HE3" & !sample.name == "HG1") #removed these based on the rarefied data

levels(as.factor(metadata$species_allknown))
#Subset samples present in both locations
metadata <- subset(metadata, species_allknown == "HB.Mytilisepta.virgata.HK" 
                   |  species_allknown == "HC.Cellana.toreuma.HK"   
                   |  species_allknown == "HG.Echinolittorina.malaccana.HK" 
                   |  species_allknown == "TF.Mytilisepta.virgata.TH"  
                   |  species_allknown == "TB.Cellana.toreuma.TH"   
                   |  species_allknown == "TD.Echinolittorina.malaccana.TH")    
levels(as.factor(metadata$species_allknown)) 

#Run the following for PERMANOVA    
metadata$species_allknown[metadata$species_allknown == "HB.Mytilisepta.virgata.HK"] <- "M. virgata"
metadata$species_allknown[metadata$species_allknown == "HC.Cellana.toreuma.HK"] <- "C. toreuma"
metadata$species_allknown[metadata$species_allknown == "HG.Echinolittorina.malaccana.HK"] <- "E. malaccana"

metadata$species_allknown[metadata$species_allknown == "TF.Mytilisepta.virgata.TH"] <- "M. virgata"
metadata$species_allknown[metadata$species_allknown == "TB.Cellana.toreuma.TH"] <- "C. toreuma"
metadata$species_allknown[metadata$species_allknown == "TD.Echinolittorina.malaccana.TH"] <- "E. malaccana"
levels(as.factor(metadata$species_allknown))

#Run the following for PERMDISP and Ordination Plot - re-run the metadata upload first
levels(as.factor(metadata$species_allknown))
metadata$species_allknown <- factor(metadata$species_allknown, levels = c("HB.Mytilisepta.virgata.HK",
                                                                          "TF.Mytilisepta.virgata.TH",
                                                                          "HC.Cellana.toreuma.HK",
                                                                          "TB.Cellana.toreuma.TH",
                                                                          "HG.Echinolittorina.malaccana.HK",
                                                                          "TD.Echinolittorina.malaccana.TH"), ordered = TRUE )
levels(as.factor(metadata$species_allknown))
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

#order the metadata  following the ASV table
metadata <- metadata[match(rownames(asv_table_rarefied4714_new ), rownames(metadata)),]
levels(as.factor(metadata$species_allknown))
##########################################################################################################################################################################################################################################################################################################################################################################################################################
## Rarefied4714 data - Statistical analysis based on Bray-Curtis dist matrix ----

#Verification steps - that only taxa present in either one of the samples
#Calculate column totals
total <- colSums(asv_table_rarefied4714_new)
asv_w_total <- rbind(asv_table_rarefied4714_new, total)

# re-order columns by total - check the last column
asv_sort <- asv_w_total[,order(-asv_w_total[which(rownames(asv_w_total) == 'total'), ])] #Here essentially it is re-ordering the total by descending order
asv_sort <- asv_sort[-30,] #remove the "total" row

#Know which one to combine, here I didn't use any of asv_sort because I didn't removed any samples
asv_rarefied4714_metadata <- cbind(asv_sort, metadata)

levels(as.factor(asv_rarefied4714_metadata$species_allknown  ))
asv_rarefied4714_metadata$species_allknown <- as.factor(asv_rarefied4714_metadata$species_allknown)
#1. PERMANOVA#######################################################################################################################################################################################################################################################################################
#Note: 1. always check if the column is right! Make sure not to include empty count taxa 
#2. The result always changes because of the random permutation and so the significance return will always be different. 
adonis.microbiome <- adonis2(asv_rarefied4714_metadata[,1:3646] ~ asv_rarefied4714_metadata$species_allknown*asv_rarefied4714_metadata$location, method = "bray")
adonis.microbiome

#                                                                              Df SumOfSqs      R2       F Pr(>F)    
#asv_rarefied4714_metadata$species_allknown                                     2   4.2816 0.35512 12.3793  0.001 ***
#asv_rarefied4714_metadata$location                                             1   1.3691 0.11356  7.9171  0.001 ***
#asv_rarefied4714_metadata$species_allknown:asv_rarefied4714_metadata$location  2   2.4287 0.20144  7.0221  0.001 ***
#Residual                                                                       23   3.9775 0.32989                   
#Total                                                                          26   11.1914                 1.00000 

#I can skipped this part - since it can't do two factors -unless I made it into groups. e.g C.toreuma_HK, C.toreuma_TH
#Post-hoc test: I plan to use fdr as correction method because bonferroni can be a bit strict.
#But I proceed with fdr which was used in qiime2 
pairwiseAdonis::pairwise.adonis(asv_rarefied4714_metadata[,1:3646], asv_rarefied4714_metadata$species_allknown, 
                               p.adjust.m ='fdr', sim.method = 'bray') 

#Here onwards (PERMDISP + Ordination plot) - you need to changes the species name following by group.
#2. PERMDISP######################################################################################################################################################################################################################################################################################
rarefied4714_bray <- vegdist(asv_rarefied4714_metadata[,1:3646], method = "bray") 

levels(as.factor(asv_rarefied4714_metadata$species_allknown))
bdisp <- betadisper(rarefied4714_bray, asv_rarefied4714_metadata$species_allknown, type=c("centroid"))
bdisp   #This is with PCoA remember! and produce the interesting graph which seem significant and different from nmds above. 
aov.bdisp <-anova(bdisp)
aov.bdisp       
permutest(bdisp) #Not significant : How to interpret more: See this http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

#Response: Distances
#           Df Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     5 0.34860 0.069721 9.5117    999  0.001 ***
#Residuals 23 0.16859 0.007330 

#The result is different between QIIME2 and PERMDISP which could be due to the random permutation

labs <- paste("Dimension", 1:8, "(", 
              round(100*bdisp$eig / sum(bdisp$eig), 2), "%)")

plot(bdisp, cex=1, pch=15:17,
     main="Bray-Curtis Dissimilarity index", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)  #Even the confidence here chosen is 0.68, so I guess my code above for nmds uses 0.7 is alright because we mainly only want to visualised. -not with any statistics logic included for the ellipse

#PCoA plotting #########################################################################################################################################################################################################################################################################################
#use the bray-curtis distance data above: rarefied4714_bray

# calculate principal coordinates analysis (Bray-Curtis)
bray.pcoa.eigen <- cmdscale(rarefied4714_bray, k = 2, eig = T)

# extract axis positions for each sample from cmdscale object and create a dataframe for plotting
bray.pcoa.eigen.plotting <- as.data.frame(bray.pcoa.eigen$points)
colnames(bray.pcoa.eigen.plotting) <- c("axis_1", "axis_2")
bray.pcoa.eigen.plotting$sample <- rownames(bray.pcoa.eigen.plotting)

bray.pcoa.eigen.plotting <- cbind(bray.pcoa.eigen.plotting, metadata)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
#Understand what is eigenvalue and eigenvector in here: https://www.youtube.com/watch?v=FgakZw6K1QQ
(bray.pcoa.eigen$eig[1]/(sum(bray.pcoa.eigen$eig)))*100 # 22.80187%

(bray.pcoa.eigen$eig[2]/(sum(bray.pcoa.eigen$eig)))*100 #16.67869%


#You need to change to species names at beginning above to create a nice figure below
#Version 2: create a PCoA plot - with spiders and location separated #The spider learnt from here: https://stackoverflow.com/questions/23463324/r-add-centroids-to-scatter-plot
axis_1 <- bray.pcoa.eigen.plotting[1]
axis_2 <- bray.pcoa.eigen.plotting[2]
centroids <- aggregate(cbind(axis_1,axis_2)~species,bray.pcoa.eigen.plotting, mean)

gg <- merge(bray.pcoa.eigen.plotting, centroids, by = "species" )

pcoa.meio.bray.plot <- ggplot(gg, aes(x = axis_1.x, y = axis_2.x, colour = species_allknown, shape = location)) +
  geom_point(size = 3) +
  geom_point(aes(x=axis_1.y,y=axis_2.y),size=2)+
  geom_segment(aes(x=axis_1.y, y=axis_2.y, xend=axis_1.x, yend=axis_2.x)) +
  theme_bw() + 
  xlab("PCoA 1 (22.80%)") +
  ylab("PCoA 2 (16.69)") +
  labs(color = "Species", shape = "Location") +
  scale_color_hue(labels = c("HC.Cellana.toreuma.HK" = expression(italic("C. toreuma")*" (HK)"), 
                          "TB.Cellana.toreuma.TH" = expression(italic("C. toreuma")*" (TH)"),
                          "HB.Mytilisepta.virgata.HK" = expression(italic("Mytilisepta virgata")~"(HK)"),
                          "TF.Mytilisepta.virgata.TH" = expression(italic("Mytilisepta virgata")~"(TH)"), 
                          "HG.Echinolittorina.malaccana.HK" =  expression(italic("E. malaccana")*"(HK)"),
                          "TD.Echinolittorina.malaccana.TH" =  expression(italic("E. malaccana")*"(TH)"))) +
    theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), #changing the title of legend
        legend.text.align = 0,
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())
pcoa.meio.bray.plot + annotate("text", x = -0.15, y = -0.3, 
                               label = "PERMANOVA:\nspecies*location, p-value < 0.05;\nPERMDISP:\nspecies(location), p-value < 0.05")

#ggsave(filename = "Rarefied_data_pcoa_plot_speciesvs_location", width = 20, height = 13, units = "cm", device = "pdf", dpi = "print")
