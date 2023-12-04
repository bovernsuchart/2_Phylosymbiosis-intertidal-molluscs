#Microbiome actual R thesis version 26.08.2022: Relative abundance bar plots 

getwd()

#Loading library
library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(forcats)
library(gridExtra)
library(RColorBrewer)
library(randomcoloR)

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
asv_table_longF_30 <- melt(asv_rel_abunF_30, variable.name = "Family",id=c("sample.name"))

#remove all those at the end; make sure you what is the last row to remove!
asv_table_longF_30  <- asv_table_longF_30[1:2387,] 

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

#Define levels of the taxa so that it plots in the correct order
asv_table_longF_30$Family<-factor(asv_table_longF_30$Family, levels = unique(asv_table_longF_30$Family))

#Make a new data fram called meta_taxa and merge the metadata and long taxa data 
meta_taxa2 <-as.data.frame(c(meta_data, asv_table_longF_30))
meta_taxa2$value <- as.numeric(meta_taxa2$value)

#31 colors vector (The colors was chosen from distinctColorPalette() function in the package called randomcoloR) Reference:  https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
#Color choosing (the process is random and I selected the most attractive and distinctive ones) - You can ignore these and choose the one below if you want colors fewer than 26/31
n <- 31
palette <- distinctColorPalette(n)
palette

P31 = c("#7F3FE2","#A5E067", "#70E3E7" ,"#EB5C3B" ,"#EB6698" ,"#75ED4D" ,"#77E692", "#E9DFD3", "#A6638D", "#A75ACE",
"#9BA799", "#D539DC" ,"#DFAB41" ,"#71DFBA", "#6D6ED8" ,"#D97DE2" ,"#C3EDDE" ,"#6DACC9","#DC9FE3" ,"#CEE8AF",
"#E347B2" ,"#E5CB87" ,"#749AE1", "#DADD72" ,"#7FA661" ,"#D5725F", "#E8B2B9" ,"#C2B0DA" ,"#D3A280" ,"#D8EB3E",
"#C4D4E9")


levels(meta_taxa2$species) #know the order of the levels or arrangement in your dataset

#how to italicize in facet_wrap. Following this link here: https://stackoverflow.com/questions/34979931/adding-expressions-to-facet-labeller-in-ggplot-2-2-0

meta_taxa2$facets <- factor(meta_taxa2$species, 
                            labels = c("italic(M.)~italic(virgata)~(HK)",
                                       "italic(C.)~italic(toreuma)~(HK)",
                                       "italic(S.)~italic(laciniosa)",
                                       "italic(Saccostrea)~italic(cucullata)",
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
                                       "Water"~"Sample"))

levels(meta_taxa2$species) 
levels(meta_taxa2$facets) 

ggplot() +
  geom_bar(aes(y = value, x = sample.name, fill = Family), data = meta_taxa2, stat="identity", position = "fill") +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y= "Relative abundance (Family)\n", x="Sample", fill = "Taxa") +
  scale_fill_manual(values=P31) + 
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title=element_text(size = 13),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"), 
        legend.position = 'right') +
  guides(fill=guide_legend(ncol=1)) + 
  facet_wrap(~location + facets , drop = TRUE, scales = "free_x"
             ,ncol = 5, nrow = 4
             ,labeller = label_parsed)

ggsave(filename = "Raw_top30_abundfamily_phylosymbiosis.pdf", device = "pdf", width = 28, height = 20, units = "cm", dpi = "print")

#mean relative abundance figure##############################################################################################################################################################################################################
ggplot() +
  geom_bar(aes(y = value, x = species, fill = Family), data = meta_taxa2, stat="identity", position = "fill") +  
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + 
  labs(y= "Mean relative abundance (Family)\n", x="Sample", fill = "Taxa") +
  scale_fill_manual(values=P31) + 
  theme_bw() +
  theme(axis.text=element_text(size=11),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title=element_text(size = 13),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.2,"cm"), 
        legend.position = 'right') 
      + guides(fill=guide_legend(ncol=1))
      + facet_wrap(location~species, drop = TRUE, scales = "free_x", nrow = 1, ncol = 2) 
