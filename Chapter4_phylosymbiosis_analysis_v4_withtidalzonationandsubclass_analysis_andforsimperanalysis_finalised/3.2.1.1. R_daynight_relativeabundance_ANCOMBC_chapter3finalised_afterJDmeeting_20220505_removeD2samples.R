#ANCOM-Bias correction method - Microbiome actual R thesis version - 21.04.2022)
#Resources following: 
#1. Lin and Peddada (2020) Analysis of compositions of microbiomes with bias correction
#2. Follow this:https://github.com/FrederickHuangLin/ANCOMBC
#2.1 More information: https://github.com/FrederickHuangLin/ANCOMBC/issues/25
#3.  http://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html

getwd()

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

#Loading library
library(phyloseq)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(rbiom)
library(ranacapa)
library(microbiome)
library(DT)
library(ANCOMBC)
?ancombc

#Import data and build phyloseq-class object: Following this file pipeline: "1. R_daynight_rarefaction_O'Brien.R" & R_daynight_relativeabundance_ANCOMBC_Lin_Peddada_github_raw.R”

#import data (important to import with read.table because the first column can be set as a vector - you can try with read.delim() to import and you will notice that at later section you can't work with this format): 
asv <- read.table("raw_feature_table.txt", sep = '\t', row.names = 1, header = T, strip.white = T)    #"\t" mean tab-delimited; " " mean space-delimited; "," for comma-delimited -normally used in .csv file
asv <- asv[,1:20]
asv <- asv[-c(8,18)]
map <- read.table("metadata_daynight_v2.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
map <- map[1:20,]
map <- subset(map, sample.name != "BRM3" &  sample.name != "DRM3")
tax <- read.table("taxonomy_gg.tsv", sep = '\t', row.names = 1, header = T, strip.white = T)
# Separate taxonomy into different columns #you can ignore the error. It just mean that the empty value is treated as NA
tax <- separate(tax, Taxon, c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                    sep = ";" , remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")
tax$Confidence <- NULL

######################################################################################################################################d#######################################################################################################################################d#######################################################################################################################################
#Prefer this method Following the naming system in  O'Brien et al.(2021)) - "1. R_daynight_Rarefaction.R" 
#Start with renaming first - Replace all NA's with the taxa indicator only or like in the kingdom case - change the name
#Always check with the unique names first for each column and check again at the end

tax$Kingdom[is.na(tax$Kingdom)] <- " k__unknown" 
tax$Kingdom[tax$Kingdom == " k__"] <- " k__unknown"
tax$Phylum[is.na(tax$Phylum)] <- " p__unknown"
tax$Phylum[tax$Phylum == " p__"] <- " p__unknown"
tax$Class[is.na(tax$Class)] <- " c__unknown"
tax$Class[tax$Class == " c__"] <- " c__unknown"
tax$Order[is.na(tax$Order)] <- " o__unknown"
tax$Order[tax$Order == " o__"] <- " o__unknown"
tax$Family[is.na(tax$Family)] <- " f__unknown"
tax$Family[tax$Family == " f__"] <- " f__unknown"
tax$Genus[is.na(tax$Genus)] <- " g__unknown"
tax$Genus[tax$Genus == " g__"] <- " g__unknown"
tax$Species[is.na(tax$Species)] <- " s__unknown"
tax$Species[tax$Genus == " s__"] <- " s__unknown"
unique(tax$Kingdom)
unique(tax$Phylum)
unique(tax$Class)
unique(tax$Order)
unique(tax$Family)
unique(tax$Genus)
unique(tax$Species)

#OPTIONAL - I can even remove the taxonomic label following the method in O'Brien et al.(2021): "Relative_abundance.R" file
#Don't do this is better for this dataset since there are many with unknown names even at the family levels
#tax$Domain <- gsub("d__", "", tax$Domain)
#tax$Kingdom <- gsub("k__", "", tax$Kingdom)
#tax$Phylum <- gsub("p__", "", tax$Phylum)
#tax$Class <- gsub("c__", "", tax$Class)
#tax$Order <- gsub("o__", "", tax$Order)
#tax$Family <- gsub("f__", "", tax$Family)
#tax$Genus <- gsub("g__", "", tax$Genus)
#tax$Species <- gsub("s__", "", tax$Species)
#unique(tax$Genus)
#unique(tax$Family)

# 'Phyloseq-ize' the data - following the Savary et al. (2021) way of naming

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)
phy.all

################################################################################################################################################################################################################################################################################################################################################################################################
##ANCOM-BC Analysis ############################################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################################################################
family_data <- aggregate_taxa(phy.all, "Family") #Stothart & Newman (2021) use an unrarefied table (ANCOM-BC has bulit-in normalization procedures)
tax_table(family_data)

#Comparison 1: 1_0915 vs 1_1230  (C1)############################################################################################################################################################################################################################################################################################################################################################
?subset_samples
C1=subset_samples(family_data, day.time %in% c("1_0915", "1_1230"))
C1

res1 = ancombc(phyloseq = C1,formula = "day.time", 
             p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  #Zerocut is important it tells that samples with proportion of zeros greater than 90 % will be excluded in the analysis
             group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, #Even Savary et al.(2021) select the neg_lb as FALSE indicating his sample size is 5 or less than that.
             max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)


res1_df = data.frame(Family= row.names(res1$res$beta), Beta=res1[["res"]][["beta"]], Beta_se=res1[["res"]][["se"]], W=res1[["res"]][["W"]],pval=res1[["res"]][["p_val"]], 
                     qval=res1[["res"]][["q_val"]],DA=res1[["res"]][["diff_abn"]])  #[[]][[]] means selecting the item(or column) within list - we are creating a new data frame here by extracting each parameters and combine them all rather than doing one by one 
                                                                                    #like how Frederick did in Bioconductor guide (you can check my pipeline in the file here: 4.2 R_daynight_relativeabundance_ANCOMBC_Lin_Peddada_github_raw.R )
colnames(res1_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") #give the correct names to each of the column 
res1_df$Comparison="Comparison 1: D1_0915 vs D1_1230" 
res1_sig = subset(res1_df, res1_df$Diff_abundant == "TRUE") # Here selecting all that shows significant difference
res1_sig$Diff_more_abundant=ifelse(res1_sig$W < 0 , "1_0915", "1_1230") #Here it tells that if less than W < 0 means it comes from the first group = 1_0915 , 
                                                                        #otherwise it is the 2nd group = 1_1230. You can refer to the Lin & Peddada (2020) paper on how the statistics was calculated. But based on a query to Lin in Github: https://github.com/FrederickHuangLin/ANCOMBC/issues/25
                                                                        #If <0, it is the reference that is high, if > 0 then the other group is high in abundance.
#(Can be ignore for this scritp) Savary et al. (2021) R script - We don't need this in here since we are only focusing at the family level - But I kept it as a reference
#res1_sig$ASV=rownames(res1_sig) #check out: res1_sig
#res1_sig$Phylum = tax$Phylum[match(rownames(res1_sig),rownames(tax))] #Here it search through the data and try to match our sequences with the taxonomic level of our interest in our tax file. 
#res1_sig$Class = tax$Class[match(rownames(res1_sig),rownames(tax))] #Here it search through the data and try to match our sequences with the taxonomic level of our interest in our tax file. 
#res1_sig$Order = tax$Order[match(rownames(res1_sig),rownames(tax))] #Here it search through the data and try to match our sequences with the taxonomic level of our interest in our tax file. 
#res1_sig$Family = tax$Family[match(rownames(res1_sig),rownames(tax))] #Here it search through the data and try to match our sequences with the taxonomic level of our interest in our tax file. 
#res1_sig$Genus = tax$Genus[match(rownames(res1_sig),rownames(tax))] #Here it search through the data and try to match our sequences with the taxonomic level of our interest in our tax file. 
#res1_sig$Species = tax$Species[match(rownames(res1_sig),rownames(tax))] #Here it search through the data and try to match our sequences with the taxonomic level of our interest in our tax file. 


message("Number of DA ASVs: ", nrow(res1_sig), "\nNumber of DA ASVs enriched in D1_0915: ", nrow(subset(res1_sig, Diff_more_abundant == "1_0915" )), "\nNumber of DA ASVs enriched in  D1_1230: ", nrow(subset(res1_sig, Diff_more_abundant == "1_1230")))

#Comparison 2: 1_0915 vs 1_1612  (C2)############################################################################################################################################################################################################################################################################################################################################################
?subset_samples
C2=subset_samples(family_data, day.time %in% c("1_0915", "1_1612"))
C2
res2 = ancombc(phyloseq = C2,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res2_df = data.frame(Family= row.names(res2$res$beta),Beta=res2[["res"]][["beta"]], Beta_se=res2[["res"]][["se"]], W=res2[["res"]][["W"]],pval=res2[["res"]][["p_val"]], 
                     qval=res2[["res"]][["q_val"]],DA=res2[["res"]][["diff_abn"]])  

colnames(res2_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res2_df$Comparison="Comparison 1: D1_0915 vs D1_1612" 
res2_sig = subset(res2_df, res2_df$Diff_abundant == "TRUE") 
res2_sig$Diff_more_abundant=ifelse(res2_sig$W < 0 , "1_0915", "1_1612") 
message("Number of DA ASVs: ", nrow(res2_sig), "\nNumber of DA ASVs enriched in D1_0915: ", nrow(subset(res2_sig, Diff_more_abundant == "1_0915" )), "\nNumber of DA ASVs enriched in  D1_1612: ", nrow(subset(res2_sig, Diff_more_abundant == "1_1612")))

#Comparison 3: 1_0915 vs 1_1740  (C3)############################################################################################################################################################################################################################################################################################################################################################
C3=subset_samples(family_data, day.time %in% c("1_0915", "1_1740"))
C3
res3 = ancombc(phyloseq = C3,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res3_df = data.frame(Family= row.names(res3$res$beta),Beta=res3[["res"]][["beta"]], Beta_se=res3[["res"]][["se"]], W=res3[["res"]][["W"]],pval=res3[["res"]][["p_val"]], 
                     qval=res3[["res"]][["q_val"]],DA=res3[["res"]][["diff_abn"]])  

colnames(res3_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res3_df$Comparison="Comparison 1: D1_0915 vs D1_1740" 
res3_sig = subset(res3_df, res3_df$Diff_abundant == "TRUE") 
res3_sig$Diff_more_abundant=ifelse(res3_sig$W < 0 , "1_0915", "1_1740") 
message("Number of DA ASVs: ", nrow(res3_sig), "\nNumber of DA ASVs enriched in D1_0915: ", nrow(subset(res3_sig, Diff_more_abundant == "1_0915" )), "\nNumber of DA ASVs enriched in  D1_1740: ", nrow(subset(res3_sig, Diff_more_abundant == "1_1740")))

#Comparison 6: 1_1230 vs 1_1612  (C6)############################################################################################################################################################################################################################################################################################################################################################
C6=subset_samples(family_data, day.time %in% c("1_1230", "1_1612"))
C6
res6 = ancombc(phyloseq = C6,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res6_df = data.frame(Family= row.names(res6$res$beta), Beta=res6[["res"]][["beta"]], Beta_se=res6[["res"]][["se"]], W=res6[["res"]][["W"]],pval=res6[["res"]][["p_val"]], 
                     qval=res6[["res"]][["q_val"]],DA=res6[["res"]][["diff_abn"]])  

colnames(res6_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res6_df$Comparison="Comparison 1: D1_1230 vs D1_1612" 
res6_sig = subset(res6_df, res6_df$Diff_abundant == "TRUE") 
res6_sig$Diff_more_abundant=ifelse(res6_sig$W < 0 , "1_1230", "1_1612") 
message("Number of DA ASVs: ", nrow(res6_sig), "\nNumber of DA ASVs enriched in D1_1230: ", nrow(subset(res6_sig, Diff_more_abundant == "1_1230" )), "\nNumber of DA ASVs enriched in  D1_1612: ", nrow(subset(res6_sig, Diff_more_abundant == "1_1612")))

#Comparison 7: 1_1230 vs 1_1740  (C7)############################################################################################################################################################################################################################################################################################################################################################
C7=subset_samples(family_data, day.time %in% c("1_1230", "1_1740"))
C7
res7 = ancombc(phyloseq = C7,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res7_df = data.frame(Family= row.names(res7$res$beta),Beta=res7[["res"]][["beta"]], Beta_se=res7[["res"]][["se"]], W=res7[["res"]][["W"]],pval=res7[["res"]][["p_val"]], 
                     qval=res7[["res"]][["q_val"]],DA=res7[["res"]][["diff_abn"]])  

colnames(res7_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res7_df$Comparison="Comparison 1: D1_1230 vs D1_1740" 
res7_sig = subset(res7_df, res7_df$Diff_abundant == "TRUE") 
res7_sig$Diff_more_abundant=ifelse(res7_sig$W < 0 , "1_1230", "1_1740") 
message("Number of DA ASVs: ", nrow(res7_sig), "\nNumber of DA ASVs enriched in D1_1230: ", nrow(subset(res7_sig, Diff_more_abundant == "1_1230" )), "\nNumber of DA ASVs enriched in  D1_1740: ", nrow(subset(res7_sig, Diff_more_abundant == "1_1740")))

#Comparison 10: 1_1612 vs 1_1740  (C10)############################################################################################################################################################################################################################################################################################################################################################
C10=subset_samples(family_data, day.time %in% c("1_1612", "1_1740"))
C10
res10 = ancombc(phyloseq = C10,formula = "day.time", 
               p_adj_method = "holm",zero_cut = 0.9,lib_cut=1000,                  
               group = "day.time",struc_zero = TRUE ,neg_lb = FALSE, tol = 1e-05, 
               max_iter = 100,conserve = TRUE,alpha = 0.05, global = TRUE)

res10_df = data.frame(Family= row.names(res10$res$beta), Beta=res10[["res"]][["beta"]], Beta_se=res10[["res"]][["se"]], W=res10[["res"]][["W"]],pval=res10[["res"]][["p_val"]], 
                     qval=res10[["res"]][["q_val"]],DA=res10[["res"]][["diff_abn"]])  

colnames(res10_df)=c("Family","Beta",	"se", "W",	"pval",	"qval", "Diff_abundant") 
res10_df$Comparison="Comparison 1: D1_1612 vs D1_1740" 
res10_sig = subset(res10_df, res10_df$Diff_abundant == "TRUE") 
res10_sig$Diff_more_abundant=ifelse(res10_sig$W < 0 , "1_1612", "1_1740") 
message("Number of DA ASVs: ", nrow(res10_sig), "\nNumber of DA ASVs enriched in D1_1612: ", nrow(subset(res10_sig, Diff_more_abundant == "1_1612" )), "\nNumber of DA ASVs enriched in  D1_1740: ", nrow(subset(res10_sig, Diff_more_abundant == "1_1740")))

#######################################################################################################################################
#Combining all the result and produce a tabulated file
ANCOMresults_sig <- rbind(res1_sig,res2_sig, res3_sig,res6_sig, res7_sig,res10_sig)
ANCOMresults_sig
getwd()
write.table(ANCOMresults_sig,  "ANCOMBC_Familylevel_results_significant_taxa.txt", 
            sep = "\t", quote = F, row.names = T ) #Create a new file

ANCOMresults_df <- rbind(res1_df,res2_df, res3_df,res6_df, res7_df,res10_df)
ANCOMresults_df

write.table(ANCOMresults_df,  "ANCOMBC_Familylevel_results_alltaxa.txt", 
            sep = "\t", quote = F, row.names = T ) #Create a new file

#######################################################################################################################################
#Plot a bar graph based on the group abundance in a comparison category
ANCOMresults_plot <- ANCOMresults_sig %>% group_by(Diff_more_abundant, Comparison) %>% tally()
ANCOMresults_plot

ggplot(data=ANCOMresults_plot, aes(x=Comparison, y=n, fill = Diff_more_abundant, label = n )) + 
  geom_bar(stat="identity", position = "stack")  + 
  geom_text(position = position_stack(vjust = 0.5), color="black", size=3) + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle=90)) 
dev.off()

ggsave("Diff_abundance_graph_finalised.pdf", width = 10, height = 8, device = "pdf", dpi = "print")


#extras: code I got from : https://www.nicholas-ollberding.com/post/identifying-differentially-abundant-features-in-microbiome-data/

#filter for only qvalue
res15_qvalue <- res15_df %>%
  dplyr::filter(qval <0.05)
res15_qvalue

res14_qvalue <- res14_df %>%
  dplyr::filter(qval <0.05)

res13_qvalue <- res14_df %>%
  dplyr::filter(qval <0.05)
