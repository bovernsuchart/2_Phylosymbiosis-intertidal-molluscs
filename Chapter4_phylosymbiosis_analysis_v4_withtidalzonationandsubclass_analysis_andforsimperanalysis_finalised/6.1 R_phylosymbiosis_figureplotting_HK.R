#This protocol was developed using O'Brien et al.(2020) R Script (Some note may be from his script - Plot_phylosymbiosis_figs.R)
## Plot dendrogram against host phylogeny 

library(phytools)
library(ggtree)
library(ggplot2)
library(cowplot)
library(ape)
library(treeio)
library(dendextend) #for making tanglegram
## Mollusc ----

## Host tree

#Host phylogeny - Timetree created phylogenetic tree - JD made tree
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

#rooting a tree (this is ignore for now and only served as a reference) - #I can't do this because JD didn't include this out group - But it is important to note this should be done first.
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

#Remove non-HK samples
host_tree <- drop.tip(host_tree, c("Saccostrea.mordax","Tenguella.musiva","Echinolittorina.radiata","Isognomon.nucleus"))
host_tree[["tip.label"]] #perfect
plot.phylo(host_tree)

#Edit tips - for better visualization with tree _ see the dendrogram section for more details (sorry I worked backward on how to edit the label. i.e., get italicize)
tiplabels_host <- host_tree$tip.label  
tiplabels_host

a <- data.frame(label = tiplabels_host,
                label2 = c("S.", "M.","B.", "C.", "S.", "R.",  "E.","M.", "L." ),
                label3 = c( "cucullata","virgata","virescens","toreuma","laciniosa", "clavigera", "malaccana", "labio", "granulata"))
a

host_tree_plot <- ggtree(host_tree,branch.length = "none") %<+%  a + #the values is adjustable and depends on the how combine figure would like - so you may need to adjust accordingly 
  geom_tiplab(offset = 0.3, aes(label= paste0('italic(', label2, ')~italic(', label3, ')')), parse = TRUE) + 
  geom_nodepoint() + 
  geom_tippoint(aes(color = label3), size=3, show.legend = FALSE) + 
  #geom_nodelab(aes(label = label), nudge_x = 0.2) +
  geom_treescale(x=1, y=1, offset = 0.25) + 
  xlim(0, 6.75) + #becareful of this (you would want to test few time by changing the values so that all the information are covered)
  theme_tree() 
host_tree_plot

## Microbiome dendrogram - BrayCurtis-upgma
a_bc <- read.tree("sample-clustering-upgma_HK.tre") #this is an original tree with the root word in it #you can ignore the warming or follow this link here for the solution: https://github.com/YuLab-SMU/tidytree/issues/10
#a_bc <- root(a_bc, outgroup = "HE.Saccostrea.cucullata.HK", resolve.root = TRUE)
plot.phylo(a_bc)

#Beautify the tree try to match with the host phylogenetic tree as much as possible
#(Do this only after you see the host phylogeny) Rotate the leaves. This depends on the node. Link : http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html
#need to reinstall and load the phytools package again - each time you do this
plot.phylo(a_bc)
nodelabels() #to know that is the node names - it's named based on numbers
a_bc_rotated <- rotate(a_bc,12)
plot.phylo(a_bc_rotated) 
nodelabels()
a_bc_rotated_2 <- rotate(a_bc_rotated,11)
plot.phylo(a_bc_rotated_2) 
nodelabels() 
a_bc_rotated_3 <- rotate(a_bc_rotated_2,17)
plot.phylo(a_bc_rotated_3)
nodelabels()
#perfect

#rename/put into a new object
a_bc <- a_bc_rotated_3
tiplabels_microbiome <- a_bc$tip.label  
tiplabels_microbiome


#Label editing (only for ggtree): following the link here: https://guangchuangyu.github.io/software/ggtree/faq/
#Method 1: Following the website
#a_b <- data.frame(label = tiplabels_microbiome, 
#               label2 = c("S. ","L. ","C. ", "Siphonaria ","E. "),
#                label3 = c( "cucullata","granulata", "toreuma","sirius","trocoides"))

tiplabels_microbiome
a_b <- data.frame(label = tiplabels_microbiome, 
                  label2 = c("L.", "M.","C.", "R.", "S.", "E.", "S.", "M.", "B."),
                  label3 = c( "granulata", "labio", "toreuma", "clavigera", "laciniosa", "malaccana", "cucullata", "virgata", "virescens"))
a_b

den_a_bc <- ggtree(a_bc,branch.length = "none") %<+% a_b + 
  geom_tiplab(offset = -2.5, aes(label= paste0('italic(', label2, ')~italic(', label3, ')')), parse = TRUE) + #the values is adjustable and depends on the how combine figure would like - so you may need to adjust accordingly. Parse = T: is important for the paste0 to work
  geom_nodepoint() + 
  geom_tippoint(aes(colour = label3), size = 3, show.legend = FALSE) + 
  #geom_nodelab(nudge_x = -0.035) + 
  scale_x_reverse(limits=c(7.75,0)) + 
  theme_tree() 
den_a_bc

# Combine dendrogram and host phylogeny -you adjsut based on the saved file produced, always edits and compare with the output file
asc1 <- ggdraw() + 
  draw_plot(host_tree_plot, x = 0, y = 0.35, width = 0.5, height = 0.6) + draw_label("Host Phylogeny", size=13, x=0.42, y=0.975) +
  draw_plot(den_a_bc, x = 0.49, y = 0.35, width = 0.45, height = 0.6) + draw_label("Microbial dendrogram", size=13, x=0.62, y=0.975) +
  draw_label("nRF = 0.71, p < 0.05 \n\nMantel r = 0.36, p < 0.05", size = 11, x = 0.5, y = 0.3)
asc1  

ggsave(filename = "Phylosymbiosis_timetreehosttree_HK.pdf", device = "pdf", width = 24, height = 15, units = "cm", dpi = "print")

#(additional - I didnt do this) Making tanglegram following Gregor et al. (2021) paper R script######################################################################################################################################################################################################################################
#the tree has to be binary in order for this to work
#Host dendrogram setting
dend_host <- as.dendrogram(host_rooted_tree5sp)

#(incase your phylogenetic tree is not binary) How to make an ultrametric tree: https://brettkvo.com/how-to-create-a-ultrametric-dichotomous-phylogeny-in-r-using-vertlife-org-or-birdtree-org/
host_rooted_tree5sp.ultra <- chronos(host_rooted_tree5sp, lambda = 0)
as.dendrogram(host_rooted_tree5sp.ultra) #failed because the tree is not binary - can't make tanglegram

#microbiome dendrogram setting
dend_ps <-as.dendrogram(a_bc)

dendlist <- dendlist(dend_host, dend_ps)

tanglegram(dendlist, sort=T, main_left = "Host phylogeny", main_right = "16S", common_subtrees_color_branches = FALSE, lwd=1, lab.cex= 1, margin_inner = 7, main="", columns_width= c(5,1,5), highlight_branches_lwd=FALSE, highlight_distinct_edges=FALSE)


#Statistical test - See the "Fig1-Gregor-et-al.nb.html" file for more information
#How is the entanglement (quality of alignment, lower is better)?
dendlist %>% entanglement


#Correlation

cor.dendlist(dendlist, method_coef = "pearson")

set.seed(42)
the_cor <- cor_bakers_gamma(dend_ps, dend_ps)
the_cor2 <- cor_bakers_gamma(dend_ps, dend_host)
the_cor
the_cor2

#create the null model
dend_mixed <- dend_ps

R <- 1000
cor_bakers_gamma_results <- numeric(R)


for(i in 1:R) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
  cor_bakers_gamma_results[i] <- cor_bakers_gamma(dend_host, dend_mixed)
}

##Create plot to show results
plot(density(cor_bakers_gamma_results),
     main = "16S to host phylogeny, Baker's gamma distribution under H0",
     xlim = c(-1,1))
abline(v = 0, lty = 2)

abline(v = the_cor, lty = 2, col = 2)
abline(v = the_cor2, lty = 2, col = 4)


legend("topleft", legend = c("16S to itself", "16S to host phylogeny"), fill = c(2,4))

round(sum(the_cor2 < cor_bakers_gamma_results)/ R, 4)

title(sub = paste("One sided p-value:",
                  "to itself =",  round(sum(the_cor < cor_bakers_gamma_results)/ R, 4),
                  " ; to host phylogeny tree =",  round(sum(the_cor2 < cor_bakers_gamma_results)/ R, 4)
))

