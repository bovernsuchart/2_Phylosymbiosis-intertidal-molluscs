#Microbiome actual R thesis version - 01.11.2022 - subclass from both locations
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
                   & !species_allknown == "TH.Turbo.bruneus" & !species_allknown == "A.water.TH") #better to remove this species since it only have two samples & water samples

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

#Setting some data as categorical factor:
levels(metadata$species_allknown)

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
asv_rarefied4714_metadata$subclass <-as.factor(asv_rarefied4714_metadata$subclass)
levels(asv_rarefied4714_metadata$subclass)

#Check if all the data is correct after combined
#View(asv_rarefied4714_metadata)

#1. PERMANOVA#######################################################################################################################################################################################################################################################################################
#Note: 1. always check if the column is right! Make sure not to include empty count taxa 
#2. The result always changes because of the random permutation and so the significance return will always be different. 
adonis.microbiome <- adonis2(asv_rarefied4714_metadata[,1:6468] ~ asv_rarefied4714_metadata$subclass, method = "bray")
adonis.microbiome

#                               Df SumOfSqs      R2     F Pr(>F)    
#asv_rarefied4714_metadata$subclass  4    7.701 0.2238 4.9736  0.001 ***
#Residual                           69   26.710 0.7762                  
#Total                              73   34.411 1.0000
#Post-hoc test: I plan to use fdr as correction method because bonferroni can be a bit strict.
#But I proceed with fdr which was used in qiime2 
pairwiseAdonis::pairwise.adonis(asv_rarefied4714_metadata[,1:6468], asv_rarefied4714_metadata$subclass, 
                                p.adjust.m ='fdr', sim.method = 'bray') 
#The comparison is huge - but many with significance

#2. PERMDISP######################################################################################################################################################################################################################################################################################
rarefied4714_bray <- vegdist(asv_rarefied4714_metadata[,1:6468], method = "bray") 
bdisp <- betadisper(rarefied4714_bray, asv_rarefied4714_metadata$subclass , type=c("centroid"))
bdisp   #This is with PCoA remember! and produce the interesting graph which seem significant and different from nmds above. 
aov.bdisp <-anova(bdisp)
aov.bdisp       
permutest(bdisp, pairwise = TRUE) #Not significant: How to interpret more: See this http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

#Response: Distances
#           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     4 0.51724 0.129311 42.615    999  0.001 ***
#Residuals 69 0.20937 0.003034  

#The result is different between QIIME2 and PERMDISP which could be due to the random permutation

labs <- paste("Dimension", 1:8, "(", 
              round(100*bdisp$eig / sum(bdisp$eig), 2), "%)")

plot(bdisp, cex=1, 
     main="Bray-Curtis Dissimilarity index", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=0.95, lwd=2)  #Even the confidence here chosen is 0.68, so I guess my code above for nmds uses 0.7 is alright because we mainly only want to visualised. -not with any statistics logic included for the ellipse


#PCoA plotting #########################################################################################################################################################################################################################################################################################
#use the bray-curtis distance data above: rarefied4714_bray

# calculate principal coordinates analysis (Bray-Curtis)
bray.pcoa.eigen <- cmdscale(rarefied4714_bray, k = 2, eig = T)
bray.pcoa.eigen

# extract axis positions for each sample from cmdscale object and create a dataframe for plotting
bray.pcoa.eigen.plotting <- as.data.frame(bray.pcoa.eigen$points)
colnames(bray.pcoa.eigen.plotting) <- c("axis_1", "axis_2")
bray.pcoa.eigen.plotting$sample <- rownames(bray.pcoa.eigen.plotting)

bray.pcoa.eigen.plotting <- cbind(bray.pcoa.eigen.plotting, metadata)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
#Understand what is eigenvalue and eigenvector in here: https://www.youtube.com/watch?v=FgakZw6K1QQ
(bray.pcoa.eigen$eig[1]/(sum(bray.pcoa.eigen$eig)))*100 #8.914112%

(bray.pcoa.eigen$eig[2]/(sum(bray.pcoa.eigen$eig)))*100 #7.719231%

#Version 1: create a PCoA plot - no spiders
pcoa.meio.bray.plot <- ggplot(bray.pcoa.eigen.plotting, aes(x = axis_1, y = axis_2, colour = subclass)) +
  geom_point(size = 3) +
  #stat_ellipse(level = 0.95) +
  theme_bw() + 
  xlab("PCoA 1 (8.91%)") +
  ylab("PCoA 2 (7.72)") +
  labs(color = "Tidal zone") +
  #scale_color_hue(labels = c("D1_0915", "D1_1230", 'D1_1612', "D1_1740")) +
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), #changing the title of legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())
pcoa.meio.bray.plot #+ annotate("text", x = -0.6, y = 0.6, label = "PERMANOVA:\n p-value > 0.05\n PERMDISP:\n p-value < 0.05")

#Version 2: create a PCoA plot - with spiders and location separated #The spider learnt from here: https://stackoverflow.com/questions/23463324/r-add-centroids-to-scatter-plot

axis_1 <- bray.pcoa.eigen.plotting[1]
axis_2 <- bray.pcoa.eigen.plotting[2]
centroids <- aggregate(cbind(axis_1,axis_2)~subclass ,bray.pcoa.eigen.plotting, mean)
centroids

gg <- merge(bray.pcoa.eigen.plotting, centroids, by = "subclass" )
gg

#get the color 
n <- 8
palette <- distinctColorPalette(n)
palette

pcoa.meio.bray.plot <- ggplot(gg, aes(x = axis_1.x, y = axis_2.x, colour = subclass, shape = location)) +
  geom_point(size = 3) +
  geom_point(aes(x=axis_1.y,y=axis_2.y),size=2)+
  geom_segment(aes(x=axis_1.y, y=axis_2.y, xend=axis_1.x, yend=axis_2.x)) +
  theme_bw() + 
  xlab("PCoA 1 (8.91%)") +
  ylab("PCoA 2 (7.72)") +
  labs(color = "Subclass", shape = "Location") +
  scale_color_hue() +
   theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 14), 
        axis.title.x = element_text(size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black"), #changing the title of legend
        legend.text.align = 0,
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())
pcoa.meio.bray.plot

pcoa.meio.bray.plot + annotate("text", x = -0.35, y = 0.15, label = "PERMANOVA: p-value < 0.05\n PERMDISP: p-value < 0.05")

ggsave(filename = "Rarefied_data_pcoa_plot_allsamples_subclass.pdf", width = 22, height = 16, units = "cm", device = "pdf", dpi = "print")

