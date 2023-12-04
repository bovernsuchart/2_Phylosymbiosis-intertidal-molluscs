#This protocol was developed using O'Brien et al.(2020) R Script (Some note may be from his script - Phylosymbiosis_analysis.R)

## Test for phylosymbiosis by; 
## a) comparing topologies between host phylogenetic tree and microbial dendrogram
## b) Correlate host phylogenetic distance with microbial dissimilarity

# Note this script requires the function RFmeasures (Mazel et al. 2018) which can be found at "https://github.com/FloMazel/Phylosymbiosis-Ecological-model"

#######################################################################################################################################################################################################################################################################################
######## Special ways of installing ggtree ############################################################################################################################################################################################################################################
#######################################################################################################################################################################################################################################################################################

#Bioconductor version
#library(BiocManager)
#BiocManager::install("ggtree")

#github version
#devtools::install_github("YuLab-SMU/ggtree")

#######################################################################################################################################################################################################################################################################################

library(ape)
library(vegan)
library(phangorn)
library(ggtree)
library(ggplot2)
library(GUniFrac)
library(rbiom) #better than phyloseq in calculation of unifrac - see here: https://github.com/joey711/phyloseq/issues/956
source("PhyloSymbiosis_Functions.R") #Download the script from here: https://github.com/FloMazel/Phylosymbiosis-Ecological-model 

## a) comparing topologies between host phylogenetic tree and microbial dendrogram

#host phylogeny - timetree created phylogenetic tree
host_tree <- read.nexus("Final_Nexus_Tree_JD.nxs") #the nexus tree JD made - a species from two different location were regarded as 1. *There seem to have problem in uploading the "Final_Tree_JD.tre" file with read.tree()
host_tree[["tip.label"]] #if here works then all is fine
host_tree  
plot.phylo(host_tree, direction = "rightward", align.tip.label = TRUE)
plot.phylo(host_tree, direction = "rightward", use.edge.length = FALSE, 
           edge.width = 1.5,  
           node.pos = 2,
           align.tip.label = TRUE)
ggtree(host_tree, branch.length = "none") #ink: https://yulab-smu.top/treedata-book/chapter4.html
plot.phylo(host_tree, type="fan", edge.color = "deeppink", tip.color = "black", cex=1)

#rooting a tree (this is ignore for now and only served as a reference) - #I can't do this because JD didn't include this out group - But it is important to note this should be done first.######################################################################################
#host_rooted_tree <- root(host_tree, outgroup = "Sipunculus nudus", resolve.root = TRUE) #you need to root first to have a rooted tree and then remove the outgroup.
#plot.phylo(host_rooted_tree)
#host_rooted_tree$edge.length
#is.rooted(host_rooted_tree)

#remove the outgroup
#drop_sipunculus_rooted_tree <- drop.tip(host_tree, "Sipunculus nudus")
#drop_sipunculus_rooted_tree
#plot.phylo(drop_sipunculus_rooted_tree) #see if sipunculus is removed and the tree is well. 
#drop_sipunculus_rooted_tree$edge.length
#is.rooted(drop_sipunculus_rooted_tree)
#host_rooted_tree5sp <- drop_sipunculus_rooted_tree #renaming
#host_rooted_tree5sp
#host_rooted_tree5sp$tip.label

##########################################################################################################################################################################################################################################################################################################################################################################################################################################
#As a confirmation again
plot.phylo(host_tree, direction = "rightward", align.tip.label = TRUE)
host_tree[["tip.label"]] 

#Changing the tiplabels
host_tree$tip.label <- gsub("'Isognomonnucleus'","Isognomon.nucleus", host_tree$tip.label)
host_tree$tip.label <- gsub("'Saccostreacucullata'","Saccostrea.cucullata", host_tree$tip.label)
host_tree$tip.label <- gsub("'Saccostreamordax'","Saccostrea.mordax", host_tree$tip.label)
host_tree$tip.label <- gsub("'Mytiliseptavirgata'","Mytilisepta.virgata", host_tree$tip.label)
host_tree$tip.label <- gsub("'Barbatiavirescens'","Barbatia.virescens", host_tree$tip.label)
host_tree$tip.label <- gsub("'Cellanatoreuma'","Cellana.toreuma", host_tree$tip.label)
host_tree$tip.label <- gsub("'Siphonarialaciniosa'","Siphonaria.laciniosa", host_tree$tip.label)
host_tree$tip.label <- gsub("'Tenguellamusiva'","Tenguella.musiva", host_tree$tip.label)
host_tree$tip.label <- gsub("'Reishiaclavigera'","Reishia.clavigera", host_tree$tip.label)
host_tree$tip.label <- gsub("'Echinolittorinaradiata'","Echinolittorina.radiata", host_tree$tip.label)
host_tree$tip.label <- gsub("'Echinolittorinamalaccana'","Echinolittorina.malaccana", host_tree$tip.label)
host_tree$tip.label <- gsub("'Monodontalabio'","Monodonta.labio", host_tree$tip.label)
host_tree$tip.label <- gsub("'Lunellagranulata'","Lunella.granulata", host_tree$tip.label)
host_tree[["tip.label"]]
plot.phylo(host_tree)

#Remove non-TH samples
host_tree <- drop.tip(host_tree, c("Siphonaria.laciniosa","Saccostrea.cucullata","Lunella.granulata","Reishia.clavigera","Barbatia.virescens","Monodonta.labio"))
host_tree[["tip.label"]] #perfect
plot.phylo(host_tree, direction = "rightward", align.tip.label = TRUE)

#Microbiome dendrogram###################################################################################################################################################################################################################################################
#Full Species
microbiome_tree <- read.tree("sample-clustering-upgma_TH.tre") #this is an original tree with the root word in it #you can ignore the warming or follow this link here for the solution: https://github.com/YuLab-SMU/tidytree/issues/10
microbiome_tree[["tip.label"]] #if here works then all is fine
plot.phylo(microbiome_tree) 

#check if the tree is binary and if names matches between host and microbiome tree 
#it still works even if the tree is not a binary tree
is.binary(host_tree) #the tree needs to be in binary - here with bayesian tree it is not binary 
is.binary(microbiome_tree) #the tree needs to be in binary

is.element(host_tree$tip.label, microbiome_tree$tip.label)

# Test if tree and dendrogram are congruent
RFmeasures(host_tree, microbiome_tree, nRandom =9999) #make sure the R script from Mazel et al.(2018) code is in the folder

#RF      RFrd
#stat 0.0 0.0
#pval 0.0 0.0

## b) Correlate host phylogenetic distance with microbial dissimilarity ----

# create distance matrix of host phylogenetic distances
host_tree_dist <- cophenetic.phylo(host_tree)
host_tree_dist

library(adephylo)
?distTips()
distTips(host_tree, method = "patristic" ) #it is still the same as above! 

# Create BC and UniFrac distance matrices of ASV table (collapsed by host species in qiime2)
# Import and format table

asv <- read.delim("phylosym_table_group_TH.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
asv_t <- t(asv)

# Edit row names (Skipped this because all the species names are correct - which i corrected them previously)
# Check matching row names
is.element(row.names(asv_t), row.names(host_tree_dist)) #Perfect

# align rows between dataframes
asv_t <- asv_t[match(row.names(host_tree_dist), row.names(asv_t)),]
asv_t

# standardise data - by rarefying as how it was done with the dendrogram rarefaction in qiime2: 19925 (I can actually rarefy it from QIIME2 - would better)
asv_r <- rrarefy(asv_t, min(apply(asv_t, 1, sum)))

#verify if all samples have equal sampling depth
apply(asv_r, 1, sum) #yes


#Remove any taxa that is not present in the any of the species first
#asv_r <- asv_r[,colSums(asv_r) > 0]

#Verification steps - that only taxa present in either one of the samples
#Calculate column totals
#total <- colSums(asv_r)
#asv_r_total <- rbind(asv_r, total)

# re-order columns by total - check the last column
#asv_sort <- asv_r_total[,order(-asv_r_total[which(rownames(asv_r_total) == 'total'), ])] #Here essentially it is re-ordering the total by descending order

#Bray-curtis distance calculation
#Method 1(my way):
bc_dist <- vegdist(asv_r[,1:1861], method = "bray")
bc_mat <- as.matrix(bc_dist)
bc_mat

# Test for Bray-Curtis correlation with mantel test # the mantel test our will keep changing depending on the rarefied data
mantel(host_tree_dist, bc_dist, method = "pearson", permutations = 9999, na.rm = TRUE) #changing to bc_dist to bc_mat would still be the same
#Mantel statistic r: 0.8185 
#Significance: 0.00079365 

mantel(host_tree_dist, bc_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
#Mantel statistic r: 0.3755 
#Significance: 0.049206 

###Just as an extra###########################################################################################################################################################################################################################################################################################
#Unifrac distance calculation (Weighted and unweighted) - slightly different than how I did with my Chapter 3 data - 3.2 R_daynight_betadiversity_ordination_statistics_rarefied5612_chapter3finalised_afterJDmeeting_20220505_removeD2samples_Unifrac.R

# For UniFrac
asv_tree <- read.tree("tree.nwk") #this is rooted-tree from qiime's exported-tree-rooted folder

Unifrac <- GUniFrac(asv_r, asv_tree, alpha = c(0,0.5,1)) # the tree has more OTU than OTU table, For alpha = 0 means unweighted, alpha =1 means weighted
Unifrac

#this is how we extract the unweighted/weighted uniFrac information
Unweighted_unifrac_dist <- Unifrac$unifracs[, , "d_UW"] 
Weighted_unifrac_dist <-  Unifrac$unifracs[, , "d_1"]

#Test for Weighted/unweighted correlation with mantel test # the mantel test our will keep changing depending on the rarefied data

mantel(host_tree_dist, Unweighted_unifrac_dist, method = "pearson", permutations = 9999, na.rm = TRUE)
mantel(host_tree_dist, Weighted_unifrac_dist, method = "pearson", permutations = 9999, na.rm = TRUE)

###Extra on how to use loop##############################################
#Learn how to use loop to access and open the many files in your folder - without the need to open one by one #But the folder must only contains the file where the function can open it. e.g. read.nexus() can only open .nex or .tre. files.
#reference: O'Brien et al. (2020)R's script: "Phylosymbiosis_analysis"
getwd()
host_filepath <- "/Users/bovern/Desktop/Phylosymbiosis trial/" #need to add / for concatenation with this directory - "Phylosymbiosis trial" is the folder contains the nexus files. 
host_files <- list.files(host_filepath) # list.files specifies the list of files in the folder. if I don't this, it will then just be a string of vectors. i.e. [1]file1,file2 instead of [1]file1 [2]file2
host_trees <- list() #need to specify the files that we have will be in a list
for (f in seq_along(host_files)) {                                                    #seq_along() specification to specify indicies and more convenient when for loop is used. reference:   https://r-lang.com/seq_along-function-in-r/   
  host_trees[[f]] <- read.nexus(paste0(host_filepath, file = host_files[[f]]))        #paste0() function is to concatenate the first argument and the second argument without any separator indicator. 
}
host_trees[[3]]$tip.label

# Root host phylogenies
for (i in seq(2)){
  host_trees[[i]] <- root(host_trees[[i]], "Eferg" , resolve.root = T) # set root
}

