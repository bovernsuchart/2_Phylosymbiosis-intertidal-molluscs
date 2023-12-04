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
                   & !species_allknown == "TH.Turbo.bruneus") #better to remove this species since it only have two samples)

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
                                                "Z.Isognomon.nucleus",
                                                "A.water.TH"))

#Setting some data as categorical factor:
metadata$species_allknown <- factor(metadata$species_allknown)
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
asv_sort <- asv_sort[-78,] #remove the "total" row

#Know which one to combine, here I didn't use any of asv_sort because I didn't removed any samples
asv_rarefied4714_metadata <- cbind(asv_sort , metadata)

#Check if all the data is correct after combined
#View(asv_rarefied4714_metadata)

#1. PERMANOVA#######################################################################################################################################################################################################################################################################################
#Note: 1. always check if the column is right! Make sure not to include empty count taxa 
#2. The result always changes because of the random permutation and so the significance return will always be different. 
adonis.microbiome <- adonis2(asv_rarefied4714_metadata[,1:7065] ~ asv_rarefied4714_metadata$species_allknown, method = "bray")
adonis.microbiome

#                                           Df SumOfSqs      R2      F Pr(>F)    
#asv_rarefied4714_metadata$species_allknown 16   23.637 0.65804 7.216  0.001 ***
#Residual                                   60   12.283 0.34196                 
#Total                                      76   35.920 1.00000       

#As a trial - even with margin - you still can't do the test
adonis.microbiome2<- adonis2(asv_rarefied4714_metadata[,1:7065] ~  asv_rarefied4714_metadata$species_allknown + asv_rarefied4714_metadata$location , by="margin", method = "bray")
adonis.microbiome2

#Post-hoc test: I plan to use fdr as correction method because bonferroni can be a bit strict.
#But I proceed with fdr which was used in qiime2 
pairwiseAdonis::pairwise.adonis(asv_rarefied4714_metadata[,1:7065], asv_rarefied4714_metadata$species_allknown, 
                                p.adjust.m ='fdr', sim.method = 'bray') 
#The comparison is huge - but many with significance

#2. PERMDISP######################################################################################################################################################################################################################################################################################
rarefied4714_bray <- vegdist(asv_rarefied4714_metadata[,1:7065], method = "bray") 
bdisp <- betadisper(rarefied4714_bray, asv_rarefied4714_metadata$species_allknown, type=c("centroid"))
bdisp   #This is with PCoA remember! and produce the interesting graph which seem significant and different from nmds above. 
aov.bdisp <-anova(bdisp)
aov.bdisp       
permutest(bdisp) #Not significant: How to interpret more: See this http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

#Response: Distances
#           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups    16 1.13557 0.070973 14.744    999  0.001 ***
#Residuals 60 0.28882 0.004814 

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
(bray.pcoa.eigen$eig[1]/(sum(bray.pcoa.eigen$eig)))*100 #8.558416%

(bray.pcoa.eigen$eig[2]/(sum(bray.pcoa.eigen$eig)))*100 #7.424511%

#Version 1: create a PCoA plot - no spiders
pcoa.meio.bray.plot <- ggplot(bray.pcoa.eigen.plotting, aes(x = axis_1, y = axis_2, colour = species_allknown)) +
  geom_point(size = 3) +
  #stat_ellipse(level = 0.95) +
  theme_bw() + 
  xlab("PCoA 1 (8.56%)") +
  ylab("PCoA 2 (7.42%)") +
  labs(color = "Species") +
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
centroids <- aggregate(cbind(axis_1,axis_2)~species_allknown,bray.pcoa.eigen.plotting, mean)
centroids

gg <- merge(bray.pcoa.eigen.plotting, centroids, by = "species_allknown" )
gg

#get the color 
n <- 18
palette <- distinctColorPalette(n)
palette

C18 <- c("#7AEA4A", "#B63DE7","#74B6DD", "#D75B72", "#77D872" ,"#D7D0E0", "#DBA2D5" ,"#9265E0", "#6DDAD6" ,"#82DFAF", "#DFB4A0", "#80887F",
         "#DCDA50" ,"#D75FC6", "#C9E9D9", "#D79754", "#7C80C8", "#D2D992"  )

pcoa.meio.bray.plot <- ggplot(gg, aes(x = axis_1.x, y = axis_2.x, colour = species_allknown, shape = location)) +
  geom_point(size = 3) +
  geom_point(aes(x=axis_1.y,y=axis_2.y),size=2)+
  geom_segment(aes(x=axis_1.y, y=axis_2.y, xend=axis_1.x, yend=axis_2.x)) +
  theme_bw() + 
  xlab("PCoA 1 (8.56%)") +
  ylab("PCoA 2 (7.42%)") +
  labs(color = "Species", shape = "Location") +
  #scale_fill_hue(labels = c()) +
  scale_color_manual(values = C18, labels = c("HB.Mytilisepta.virgata.HK" = expression(italic("M. virgata")~"(HK)"),    #you can set his as scale_color_hue() but without the values argument which is needed for scale_fill_manual
                              "HC.Cellana.toreuma.HK" = expression(italic("C. toreuma")~"(HK)"),
                             "HD.Siphonaria.laciniosa" = expression(italic("S. laciniosa")), 
                             "HE.Saccostrea.cucullata.HK" = expression(italic("Saccostrea cucullata")),
                             "HF.Lunella.granulata" = expression(italic("L. granulata")),
                             "HG.Echinolittorina.malaccana.HK" = expression(italic("E. malaccana")~"(HK)"),
                             "HH.Reishia.clavigera" = expression(italic("R. clavigera")),            
                             "HJ.Barbatia.virescens.HJ" = expression(italic("B. virescens")), 
                             "HK.Monodonta.labio_HK" = expression(italic("Monodonta labio")), 
                             "TB.Cellana.toreuma.TH" = expression(italic("C. toreuma")~"(TH)"), 
                             "TC.Saccostrea.sp.TH" = expression(italic("Saccostrea sp.")), 
                             "TD.Echinolittorina.malaccana.TH" = expression(italic("E. malaccana")~"(TH)"), 
                             "TF.Mytilisepta.virgata.TH" = expression(italic("M. virgata")~"(TH)"),          
                             "TG.Tenguella.musiva" = expression(italic("T. musiva")),                 
                             "TJ.Echinolittorina.radiata.TJ" = expression(italic("E. radiata")),   
                             "Z.Isognomon.nucleus" = expression(italic("I. nucleus")),
                             "A.water.TH"= "Water sample (TH)")) +
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

pcoa.meio.bray.plot + annotate("text", x = -0.3, y = 0.45,label = "PERMANOVA:p-value < 0.05\n PERMDISP:p-value < 0.05", size = 5 )

ggsave(filename = "Rarefied_data_pcoa_plot_allsamples_bothlocations", width = 22, height = 16, units = "cm", device = "pdf", dpi = "print")

