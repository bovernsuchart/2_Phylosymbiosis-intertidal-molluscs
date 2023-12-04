#Microbiome actual R thesis version - 05.09.2022

#The following pipeline follows - O'Brien et al.(2021): "Relative_abundance.R" file

#Loading library
library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(forcats)
library(gridExtra)
library(RColorBrewer)
library(randomcoloR)
#Import data 

#Paul O'Brien's email: So, during your QIIME2 analysis you would have to produce a file called ’taxa-bar-plots.qzv’ (or similar) 
#at some stage. When you open the file in QIIME2 viewer it gives you the relative abundance of each microbial taxa in each sample. 
#Here you can select the taxonomic level you are interested in, e.g., level 2 for phylum, and download a .csv file with the raw values 
#then use this csv as input into R (note: the metadata was already added to this file)

##############################################################################################################################
#####Taxonomic profiles of the 30 most abundant families######################################################################
##############################################################################################################################

#############################################################
#####Taxonomic profiles of the 30 most abundant families##### 
#############################################################

## Plot family level (Top 30) ----
asv_table_f_30 <- read.csv("raw_taxa_fam.csv", sep = ',', header = T, row.names = 1, strip.white = T) # Family table created using qiime2
asv_table_f_30 <- subset(asv_table_f_30, !sample.name == "TC3" & !sample.name == "TG1" &  !sample.name == "HE3" & !sample.name == "HG1"  #removed these based on the rarefied data
                         & !species == "TH.Turbo.bruneus") #better to remove this species since it only have two samples)

levels(as.factor(asv_table_f_30$species))
asv_table_f_30$species<- factor(asv_table_f_30$species, #re-order the levels to make it nicer in comparison
                                levels = c( "HB.Brachidontes.variabilis", 
                                            "HC.Cellana.toreuma.HK",
                                            "HD.Siphonaria.sirus", 
                                            "HE.Saccostrea.cucullata.HK",
                                            "HF.Lunella.granulata" ,
                                            "HG.Echinolittorina.trochoides.HK",
                                            "HH.Reishia.clavigera",            
                                            "HJ.Unknown_bivalves_HJ",
                                            "HK.Unknown_snail_HK",
                                            "TB.Cellana.toreuma.TH",
                                            "TC.Saccostrea.mordax.TH", 
                                            "TD.Echinolittorina.trochoides.TH", 
                                            "TF.Unknown_bivalves_TF",          
                                            "TG.Tenguella.musiva",              
                                            "TJ.Unknown_snail_TJ" ,   
                                            "Z.Isognomon.sp",
                                            "A.water.TH"))
asv_table_f_30$species<- factor(asv_table_f_30$species)                       
levels(asv_table_f_30$species) #Perfect

#subset metadata from the .csv file - need to do this because when we merged the table it has to be in the same number of roles and column
meta_data <- asv_table_f_30[458:465]
meta_data$species <- as.factor(meta_data$species)
meta_data$species

asv_table_f2_30 <- as.matrix(asv_table_f_30[, -c(458:465)]) # Remove metadata columns from the asv_table_f

# Calculate Top 30 abundant family #################################################################

# Calculate column totals
total_30 <- colSums(asv_table_f2_30)
asv_w_total_30 <- rbind(asv_table_f2_30, total_30)

# re-order columns by total
asv_sort_30 <- asv_w_total_30[,order(-asv_w_total_30[which(rownames(asv_w_total_30) == 'total_30'), ])] #Here essentially it is re-ordering the total by descending order

## Convert counts to percentages on data matrix & remove row total
asv_percent_30 <- asv_sort_30[-which(rownames(asv_sort_30) == 'total_30'),] / rowSums(asv_sort_30[-which(rownames(asv_sort_30) == 'total_30'),]) *100

# create new column "other" which is the sum of all the taxa not in top 30, remove other taxa
Other <- rowSums(asv_percent_30[,-c(1:30)])

# Combine data
asv_top30 <- cbind(asv_percent_30[,1:30], Other)
asv_rel_abunF_30 <- cbind(asv_top30, meta_data) #just to get the sample name; I'm doing an extra step more than Stevick et al.(2019)
asv_table_longF_30 <- melt(asv_rel_abunF_30, variable.name = "Family",id=c("sample.name", "location", "date", "time", "species", "group", "zone", "reference" ))

# Remove long names
unique(asv_table_longF_30$Family)

#First example: asv_table_longF$Family <- gsub("d__Bacteria.__.__.__.__", "Bacteria unknown", asv_table_longF$Family) #Always start with renaming first
#Second example: asv_table_longF$Family <- gsub("Unassigned.__.__.__.__", "Unassigned", asv_table_longF$Family)
asv_table_longF_30$Family <- gsub("Firmicutes.__.__.__", "Firmicutes unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Bacteroidia.__.__", "Bacteroidia unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Bacilli.__.__", "Bacilli unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("SAR86_clade.f__", "SAR86 clade unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("SAR11_clade.f__Clade_I", "SAR11 Clade unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Planctomycetales.__", "Planctomycetales unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Alphaproteobacteria.__.__", "Alphaproteobacteria unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Rhizobiales.__", "Rhizobiales unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("SAR324_clade.Marine_group_B..c__.o__.f__", "SAR324 clade (marine group B) unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Gammaproteobacteria.__.__", "Gammaproteobacteria unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Gracilibacteria.o__.f__", "Gracilibacteria unknwon", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Planctomycetales.f__uncultured", "Planctomycetales uncultured", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Planctomycetales.__", "Planctomycetales unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("Bacteria.k__.__.__.__.__", "Bacteria unknown", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("d__Bac.*f__", "\\1", asv_table_longF_30$Family)  #The order here is super important - do one by one to remove those unwanted parts from the names
asv_table_longF_30$Family <- gsub("d__Bac.*o__", "\\1", asv_table_longF_30$Family)  # The * means "to" : From d__Bac. "to" "p__" | or you can think of it as "times" : "a name"*"a name" as a vector. You can add a fullstop if you dont want to write the full name. e.g. Bac. instead of bacteria before the "*". Only with "*" you can do this
asv_table_longF_30$Family <- gsub("d__Bac.*c__", "\\1", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("d__Bac.*p__", "\\1", asv_table_longF_30$Family)
asv_table_longF_30$Family <- gsub("d__", "\\1", asv_table_longF_30$Family)

unique(asv_table_longF_30$Family) #Check again if anything is weird

#Factorizing your group first : this part is super important! Make sure you have factorise all the factor
asv_table_longF_30$sample.name <- factor(asv_table_longF_30$sample.name)
asv_table_longF_30$species <- factor(asv_table_longF_30$species)
asv_table_longF_30$date <- factor(asv_table_longF_30$date)
asv_table_longF_30$time <- factor(asv_table_longF_30$time)
asv_table_longF_30$location <- factor(asv_table_longF_30$location)

#(this part important to do first - Define levels of the taxa so that it plots in the correct order - learning from Stevick_etal_2019 and here: https://stackoverflow.com/questions/25098107/how-to-order-the-levels-of-factors-according-to-the-ordering-of-a-data-frame-an
asv_table_longF_30$Family<-factor(asv_table_longF_30$Family, levels = rev(as.character(unique(asv_table_longF_30$Family)))) #rev is to reverse the order

#Group by Species (To get the mean relative abundance) 
taxon_summaryF <- asv_table_longF_30 %>% 
  group_by(species, location, Family) %>%
  summarise(Percentage = mean(value, na.rm = T))

#(Not needed here - the above code (i.e., above taxon_summaryF) will make the Others till the end) taxon_summaryF <- taxon_summaryF %>% mutate(Family = fct_relevel(Family, "Other", after = Inf)) #after = Inf is to relevel to the end, when the number of levels is unknown

#Better to change the 0 value to NA- better for visualization in bubble plot - Epstein et al.(2019) did the same
taxon_summaryF$Percentage[taxon_summaryF$Percentage == "0"] <- NA

#Version 1: Plot family relative abundance
bubblef <- ggplot(taxon_summaryF , aes(x = species, y = Family)) +
  geom_point(aes(size = Percentage), colour = "#000000", fill = "#56B4E9", shape = 21) + 
  scale_size(range = range(taxon_summaryF$Percentage)) +
  #scale_size_continuous(range= c(1,10)) +
   scale_size_area(max_size = 12) +
  labs(x = "\nHost", y = "Relative Abundance (Family)\n") +
  guides(colour = guide_legend(title = "species")) +
  theme_bw() +
  theme(axis.text=element_text(size=11, colour = "black"),axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13))
bubblef

#Version 2: Plot family relative abundance
bubblef2 <- ggplot(taxon_summaryF , aes(x = species, y = Family)) +
  geom_point(aes(size = Percentage), colour = "#000000", fill = "#56B4E9", shape = 21, stroke = 0) + 
  scale_size(range = range(taxon_summaryF$Percentage)) +
  #scale_size_continuous(range= c(1,10)) +
  scale_size_area(max_size = 10) +
  scale_x_discrete(name = "Species", 
                   labels = c("HB.Brachidontes.variabilis" = expression(italic("M. virgata")~("HK")),
                              "HC.Cellana.toreuma.HK" = expression(italic("C. toreuma")~"(HK)"),
                              "HD.Siphonaria.sirus" = expression(italic("S. laciniosa")), 
                              "HE.Saccostrea.cucullata.HK" = expression(italic("S. cucullata")),
                              "HF.Lunella.granulata" = expression(italic("L. granulata")),
                              "HG.Echinolittorina.trochoides.HK" = expression(italic("E. malaccana")~"(HK)"),
                              "HH.Reishia.clavigera" = expression(italic("R. clavigera")),            
                              "HJ.Unknown_bivalves_HJ" = expression(italic("Barbatia virescens")), 
                              "HK.Unknown_snail_HK" = expression(italic("Monodonta labio")), 
                              "TB.Cellana.toreuma.TH" = expression(italic("C. toreuma")~"(TH)"), 
                              "TC.Saccostrea.mordax.TH" = expression(italic("S. mordax")), 
                              "TD.Echinolittorina.trochoides.TH" = expression(italic("E. malaccana")~"(TH)"), 
                              "TF.Unknown_bivalves_TF" = expression(italic("M. virgata")~"(TH)"),          
                              "TG.Tenguella.musiva" = expression(italic("T. musiva")),           
                              "TJ.Unknown_snail_TJ" = expression(italic("E. radiata")),   
                              "Z.Isognomon.sp" = expression(italic("I. nucleus")),
                              "A.water.TH"= "Water sample (TH)")) +
  labs(y = "Taxonomy (Family)\n") +
  theme_bw() +
  guides(size= guide_legend("Relative\nabundance (%)")) + #the "size" depends on the items(or variable) you put the color on
  theme(axis.text=element_text(size=11, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title=element_text(size = 13),
        legend.text=element_text(size=11), 
        legend.title=element_text(size=13),
        strip.background = element_rect(fill = "white", colour = "black"))
bubblef2 + facet_wrap(~location, drop = TRUE, scales = "free_x")


#Version 3: Plot family relative abundance - The best
bubblef3 <- ggplot(taxon_summaryF , aes(x = species, y = Family, color = species, fill = species)) +
  geom_point(aes(size = Percentage), stroke = 0.5) + 
  scale_size(range = range(taxon_summaryF$Percentage)) +
  #scale_size_continuous(range= c(1,10)) +
  scale_size_area(max_size = 9) +
  scale_x_discrete(name = "Species", 
                   labels = c("HB.Brachidontes.variabilis" = expression(italic("M. virgata")~("HK")),
                              "HC.Cellana.toreuma.HK" = expression(italic("C. toreuma")~"(HK)"),
                              "HD.Siphonaria.sirus" = expression(italic("S. laciniosa")), 
                              "HE.Saccostrea.cucullata.HK" = expression(italic("Saccostrea cucullata")),
                              "HF.Lunella.granulata" = expression(italic("L. granulata")),
                              "HG.Echinolittorina.trochoides.HK" = expression(italic("E. malaccana")~"(HK)"),
                              "HH.Reishia.clavigera" = expression(italic("R. clavigera")),            
                              "HJ.Unknown_bivalves_HJ" = expression(italic("B. virescens")), 
                              "HK.Unknown_snail_HK" = expression(italic("Monodonta labio")), 
                              "TB.Cellana.toreuma.TH" = expression(italic("C. toreuma")~"(TH)"), 
                              "TC.Saccostrea.mordax.TH" = expression(italic("Saccostrea")~"sp."), 
                              "TD.Echinolittorina.trochoides.TH" = expression(italic("E. malaccana")~"(TH)"), 
                              "TF.Unknown_bivalves_TF" = expression(italic("M. virgata")~"(TH)"),          
                              "TG.Tenguella.musiva" = expression(italic("T. musiva")),           
                              "TJ.Unknown_snail_TJ" = expression(italic("E. radiata")),   
                              "Z.Isognomon.sp" = expression(italic("I. nucleus")),
                              "A.water.TH"= "Water sample (TH)")) +
  labs(y = "Taxonomy (Family)\n") +
  theme_bw() +
  guides(size= guide_legend("Relative\nabundance (%)"),color = "none", fill = "none") + #the "size" or "fill" depends on the items(or variable) you label in the legend. Setting to "none" to remove the particular legend, but you still want the colour/fill/shape in the figure
   theme(axis.text=element_text(size=11, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title=element_text(size = 13),
        legend.text=element_text(size=11), 
        legend.title=element_text(size=13),
        strip.background = element_rect(fill = "white", colour = "black"),
        #panel.border = element_blank()
        )
bubblef3

bubblef3 + facet_wrap(~location, drop = TRUE, scales = "free_x")

ggsave(filename = "Raw_top30_abundfamily_phylosymbiosis_bubbleplot.pdf", device = "pdf", width = 28, height = 19, units = "cm", dpi = "print")


