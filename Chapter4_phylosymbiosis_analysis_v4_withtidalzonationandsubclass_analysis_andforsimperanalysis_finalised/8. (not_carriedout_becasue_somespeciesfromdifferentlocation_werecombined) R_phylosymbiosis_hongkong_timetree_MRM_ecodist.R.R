#This protocol was developed following the description from Escalas et al. (2021) & Bletz et al. (2017) & Weinstein et al. (2021)
#The MRM code was from Weinstein et al. (2021)

#######################################################################################################################################################################################################################################################################################

library(ape)
library(vegan)
library(phangorn)
library(ggtree)
library(ggplot2)
library(GUniFrac)
library(cluster) #for Gower distance metic calculation
library(rbiom) #better than phyloseq in calculation of unifrac - see here: https://github.com/joey711/phyloseq/issues/956

#Edit the names in host phylogenetic tree 
#host phylogeny - timetree created phylogenetic tree
host_tree <- read.tree("timetree_hk_v4.nwk") 
host_tree[["tip.label"]] #if here works then all is fine
host_tree  
plot.phylo(host_tree, direction = "rightward", align.tip.label = TRUE)
ggtree(host_tree, branch.length = "none") #ink: https://yulab-smu.top/treedata-book/chapter4.html
plot.phylo(host_tree, type="fan", edge.color = "deeppink", tip.color = "black", cex=1)

#rooting a tree (this is ignore for now and only served as a reference)
host_rooted_tree <- root(host_tree, outgroup = "Sipunculus nudus", resolve.root = TRUE) #you need to root first to have a rooted tree and then remove the outgroup.
plot.phylo(host_rooted_tree)
host_rooted_tree$edge.length
is.rooted(host_rooted_tree)

#remove the outgroup
drop_sipunculus_rooted_tree <- drop.tip(host_tree, "Sipunculus nudus")
drop_sipunculus_rooted_tree
plot.phylo(drop_sipunculus_rooted_tree) #see if sipunculus is removed and the tree is well. 
drop_sipunculus_rooted_tree$edge.length
is.rooted(drop_sipunculus_rooted_tree)
host_rooted_tree5sp <- drop_sipunculus_rooted_tree #renaming
host_rooted_tree5sp
host_rooted_tree5sp$tip.label

#Changing the tiplabels
host_rooted_tree5sp$tip.label <- gsub("Echinolittorina trochoides A","HG.Echinolittorina.trochoides.HK", host_rooted_tree5sp$tip.label)
host_rooted_tree5sp$tip.label <- gsub("Cellana toreuma","HC.Cellana.toreuma.HK", host_rooted_tree5sp$tip.label)
host_rooted_tree5sp$tip.label <- gsub("Saccostrea mordax","HE.Saccostrea.cucullata.HK", host_rooted_tree5sp$tip.label)
host_rooted_tree5sp$tip.label <- gsub("Siphonaria pectinata","HD.Siphonaria.sirus", host_rooted_tree5sp$tip.label)
host_rooted_tree5sp$tip.label <- gsub("Lunella granulata","HF.Lunella.granulata", host_rooted_tree5sp$tip.label)
host_rooted_tree5sp[["tip.label"]]
plot.phylo(host_rooted_tree5sp)

# create distance matrix of host phylogenetic distances -  co.phenetic.phylo() can be used - but you have to do extre step of converting the data matrix output into distance metric - using as.dist()
library(adephylo)
host_rooted_tree5sp_dist <- distTips(host_rooted_tree5sp, method = "patristic" ) #it is still the same as above! 
as.matrix(host_rooted_tree5sp_dist)

# Create BC and UniFrac distance matrices of ASV table (collapsed by host species in qiime2)
# Import and format table

asv <- read.delim("species_table_group_HK.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
asv_t <- t(asv)

# Edit row names (Skipped this because all the species names are correct - which i corrected them previously)
# Check matching row names
is.element(row.names(asv_t), row.names(as.matrix(host_rooted_tree5sp_dist))) #Perfect

# align rows between dataframes
asv_t <- asv_t[match(row.names(as.matrix(host_rooted_tree5sp_dist)), row.names(asv_t)),]

# standardise data - by rarefying as how it was done with the dendrogram rarefaction in qiime2: 19925 (I can actually rarefy it from QIIME2 - would better)
asv_r <- rrarefy(asv_t, min(apply(asv_t, 1, sum)))

#verify if all samples have equal sampling depth
apply(asv_r, 1, sum) #yes

sample_dist <- vegdist(asv_r[,1:711], method = "bray")
sample_dist 
host_rooted_tree5sp_dist #perfectly aligned for both

#load the metadata 
library(ecodist)
metadata <- read.table("metadata_phylosymbiosis_group_HK_modified.txt", sep = '\t', row.names = 1, header = T, strip.white = T)

?read.table


metadata
is.element(row.names(asv_r), row.names(metadata )) #Check if the names are present for all and matched first

# align rows between dataframes
metadata <- metadata[match(row.names(asv_r), row.names(metadata)),]
rownames(metadata)
rownames(asv_r)

#We need to use Gower distance metric for categorical variable - like how Huang et al. (2022) & Youngblut et al. (2019)
#The Gower distance take on the scale between 0 and 1 only. So, you need to factor your categorical variable! Don't rank it if isn't ordinal - See my note here: Multifactor_effect_on_microbial_composition.docx

#Calculate the Gower distance matrix for group (body structure)

metadata$group<- as.factor(metadata$group)
group_gower <- daisy(metadata[6] , metric = c("gower"))
View(as.matrix(group_gower))
View(as.matrix(sample_dist))
View(as.matrix(host_rooted_tree5sp_dist))

#Calculate the distance matrix for zone
metadata$zone<- as.factor(metadata$zone)
zone_gower <- daisy(metadata[8] , metric = c("gower"))
View(as.matrix(zone_gower))

#Multiple regression on matrices (MRM) - using Ecodist package - code from Weinstein et al. (2021)  
MRM <- MRM(sample_dist ~ group_gower + zone_gower , mrank = TRUE)
MRM

#Multiple regression on matrices (MRM) - using Ecodist package - code from Weinstein et al. (2021)  
MRM_full <- MRM(sample_dist ~ group_gower + zone_gower + host_rooted_tree5sp_dist , mrank = TRUE)
MRM_full


