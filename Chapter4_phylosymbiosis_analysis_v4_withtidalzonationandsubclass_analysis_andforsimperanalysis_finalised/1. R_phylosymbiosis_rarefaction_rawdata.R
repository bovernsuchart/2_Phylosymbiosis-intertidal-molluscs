#Microbiome actual R thesis version - 26.08.2022

## Load packages
library(ggplot2)
library(ggtext)
library(phyloseq)
library(dplyr)
library(vegan)
library(tidyr)
library(reshape2)
library(rbiom)
library(devtools)
library(ranacapa)

getwd()
#Import data (important to import with read.table because the first column can be set as a vector - you can try with read.delim() to import and you will notice that at later section you can't work with this format): 
#Metadata ###############################################################################################################################################################################################################################################################################################################
metadata <- read.table("metadata_phylosymbiosis.txt", sep = '\t', row.names = 1, header = T, strip.white = T) 
metadata <- subset(metadata, !sample.name == "TC3" & !sample.name == "TG1" &  !sample.name == "HE3" & !sample.name == "HG1" #removed these based on the rarefied data
                 & !species_allknown == "TH.Turbo.bruneus") #better to remove this species since it only have two samples)

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
                                               "Z.Isognomon.nucleus",
                                               "A.water.TH"))

#this part is for facet wrap later : see this file here for the actual script flow how this is made: 4. R_phylosymbiosis_relativeabundance_barplot_family_allsamples_withoutisognomon_foreasevisualization.R
levels(metadata$species_allknown) #know the order of the levels or arrangement in your dataset
metadata$facets <- factor(metadata$species_allknown, 
                            labels = c( "italic(M.)~italic(virgata)~(HK)",
                                        "italic(C.)~italic(toreuma)~(HK)",
                                        "italic(Siphonaria)~italic(laciniosa)",
                                        "italic(Sacosstrea)~italic(cucullata)",
                                        "italic(L.)~italic(granulata)",
                                        "italic(E.)~italic(malaccana)~(HK)",
                                        "italic(R.)~italic(clavigera)",
                                        "italic(B.)~italic(virescens)",
                                        "italic(Monodonta)~italic(labio)",
                                        "italic(C.)~italic(toreuma)~(TH)",
                                        "italic(Saccostrea)~sp.",
                                        "italic(E.)~italic(malaccana)~(TH)", 
                                        "italic(M.)~italic(virgata)~(TH)",
                                        "italic(T.)~italic(musiva)",    
                                        "italic(E.)~italic(radiata)",
                                        "italic(I.)~italic(nucleus)",
                                        "Water~sample~(TH)"))

levels(metadata$species_allknown) 
levels(metadata$facets) 

# Setting some data as categorical factor:
metadata$sample.name <- factor(metadata$sample.name)
metadata$location <- factor(metadata$location)
metadata$date <- factor(metadata$date)
metadata$group <- factor(metadata$group)


#ASV###############################################################################################################################################################################################################################################################################################################
asv_table <- read.table("raw_feature_table.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
asv_table_transpose <- t(asv_table)      #t() here means transpose the matrix

#the %in% operator returns a logical vector that tells us whether or not there is a match between the first object and the second and subset() does only for those that present in the first.
#This checking is important because sometimes we want to might rarefied to sampling depth that might removed some samples. So we want to know what are those samples. 

asv_table_transpose_relevantsamples <- subset(asv_table_transpose, (rownames(asv_table_transpose) %in% rownames(metadata))) 

#Remove any taxa that is not present in the any of the samples first
asv_table_transpose_relevantsamples_relavanttaxa <- asv_table_transpose_relevantsamples[,(colSums(asv_table_transpose_relevantsamples) > 0 )]
asv_table_editted <- t(asv_table_transpose_relevantsamples_relavanttaxa)

#Taxonomy ###############################################################################################################################################################################################################################################################################################################
taxonomy <- read.table("taxonomy_silva.tsv", sep = '\t', row.names = 1, header = T, strip.white = T)
taxonomy$Confidence <- NULL


#Separate taxonomy into different columns #you can ignore the error. It just mean that the empty value is treated as NA
tax_sep <- separate(taxonomy, Taxon, c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "species_allknown"),
                    sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")

#There will be warning message but do not worry about it

# Replace all NA's with the taxa indicator only
tax_sep$Phylum[is.na(tax_sep$Phylum)] <- " p__unknown"
tax_sep$Class[is.na(tax_sep$Class)] <- " c__unknown"
tax_sep$Order[is.na(tax_sep$Order)] <- " o__unknown"
tax_sep$Family[is.na(tax_sep$Family)] <- " f__unknown"
tax_sep$Genus[is.na(tax_sep$Genus)] <- " g__unknown"
tax_sep$species[is.na(tax_sep$species)] <- " s__unknown"

# 'Phyloseq-ize' the data
asv.tab <- otu_table(as.matrix(asv_table_editted), taxa_are_rows = T)   #this is ASV sequence, but since Phyloseq was initially designed for Otu, hence why they named the data as OTU. But, it does not matter.

tax.tab <- tax_table(as.matrix(tax_sep))

sample.data <- sample_data(metadata)

physeq <- phyloseq(asv.tab, tax.tab, sample.data)
physeq

#Make some data as factor
sample_data(physeq)$species_allknown = factor(sample_data(physeq)$species_allknown)

sample_data(physeq)$location = factor(sample_data(physeq)$location)

###Just as a test##########################################################################################################################################################################################################
HK_samples <- sample_data(physeq)[sample_data(physeq)$location == "HK"] #just a test

#Subset sample from a location
subset_samples(physeq, sample_data(physeq)$location == "HK")
subset_samples(physeq, sample_data(physeq)$location == "TH")
####################################################################################################################################################################################################################################################################################################################################################################################################################


#Make Rarefaction curves - All samples
rarefaction_plot <- ggrare(physeq, step = 500, color = "species_allknown", label = NULL, se = FALSE) # The number of rarefaction depths to include between mindepth and max-depth

rarefaction_plot_full <- rarefaction_plot + guides(colour = guide_legend(title = "Species")) +
              guides(colour = guide_legend(title = "species")) + 
              scale_x_continuous(limits = c(0, 50000), breaks = c(5000, 25000, 50000)     ) +
             #scale_y_continuous(expand = c(0,0), breaks = seq(0, 1000, 100)) +
              xlab("Sampling depth") +  
              ylab("ASV richness") +
              geom_vline(xintercept =3500) +
              theme_bw()
rarefaction_plot_full


#ggsave("Rarefaction_ASV_Richness_full.pdf", width = 15, height = 10, dpi = "print")

rarefaction_plot_sample <- rarefaction_plot + facet_wrap(~species_allknown) + 
  guides(colour = guide_legend(title = "Species")) + 
  #scale_x_continuous(breaks = c(5000, 25000, 50000, 75000, 100000)) +
  xlab("Sampling depth") +  
  ylab("ASV richness") +
  theme_bw()
rarefaction_plot_sample

#ggsave("Rarefaction_ASV_Richness_sample.pdf", width = 15, height = 10, dpi = "print")
#how to italicize in facet_wrap. Follwoing this link here: https://stackoverflow.com/questions/34979931/adding-expressions-to-facet-labeller-in-ggplot-2-2-0


#I shifted the following datasets to the above: before physeq. 
levels(metadata$species_allknown) #know the order of the levels or arrangement in your dataset
#metadata$facets <- factor(metadata$species_allknown, 
#                    labels = c( "italic(M.)~italic(virgata)~(HK)",
#                            "italic(C.)~italic(toreuma)~(HK)",
#                            "italic(Siphonaria)~italic(laciniosa)",
#                            "italic(Saccostrea)~italic(cucullata)",
#                            "italic(L.)~italic(granulata)",
#                            "italic(E.)~italic(malaccana)~(HK)",
#                            "italic(R.)~italic(clavigera)",
#                            "italic(Barbatia)~italic(virescens)",
#                            "italic(Monodonta)~italic(labio)",
#                            "italic(C.)~italic(toreuma)~(TH)",
#                            "italic(Saccostrea)~sp.",
#                            "italic(E.)~italic(malaccana)~(TH)", 
#                            "italic(M.)~italic(virgata)~(TH)",
#                            "italic(T.)~italic(musiva)",    
#                            "italic(E.)~italic(radiata)",
#                            "italic(I.)~italic(nucleus)",
#                            "Water~sample~(TH)"))

rarefaction_plot_sample_location <- rarefaction_plot + facet_wrap(~location + facets, ncol = 9, nrow =2,labeller = label_parsed) +
  #guides(colour = guide_legend(title = "Sample")) + 
  xlab("Sampling depth") +  
  ylab("ASV richness") +
  theme_bw() +
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

rarefaction_plot_sample_location

ggsave("Rarefaction_ASV_Richness_species_location.pdf", width = 14, height = 5, dpi = "print")

