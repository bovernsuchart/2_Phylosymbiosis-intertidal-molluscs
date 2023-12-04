#Microbiome actual R thesis version 08.09.2022: New relative abundance figure

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


#Top 30 for each species################################################################################################################################################################################################

#1: Isognomon.nucleus
asv_top30.IS <- asv_top30[1:5,] #BRM1,2,3,4,5
View(asv_top30.IS)
# Calculate column totals
total_30.IS <- colSums(asv_top30.IS )
asv_top30_total.IS <- rbind(asv_top30.IS , total_30.IS)
# re-order columns by total - Top 3:
asv_sort_30_IS <- asv_top30_total.IS[,order(-asv_top30_total.IS[which(rownames(asv_top30_total.IS) == 'total_30.IS'), ])] #Here essentially it is re-ordering the total by descending order
colnames(asv_sort_30_IS)
#write.table(asv_sort_30_IS ,  "Family_asv_sort_30_IS_isognomon.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#2: HB.Brachidontes.variabilis
asv_top30.HB <- asv_top30[6:10,]
View(asv_top30.HB)
# Calculate column totals
total_30.HB <- colSums(asv_top30.HB )
asv_top30_total.HB <- rbind(asv_top30.HB , total_30.HB)
# re-order columns by total - Top 3: 
asv_sort_30_HB <- asv_top30_total.HB[,order(-asv_top30_total.HB[which(rownames(asv_top30_total.HB) == 'total_30.HB'), ])] #Here essentially it is re-ordering the total by descending order
colnames(asv_sort_30_HB) 
#write.table(asv_sort_30_HB ,  "Family_asv_sort_30_HB_Mytilisepta_HK.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#3: HC.Cellana.toreuma.HK
asv_top30.HC <- asv_top30[11:15,]
View(asv_top30.HC)
# Calculate column totals
total_30.HC <- colSums(asv_top30.HC )
asv_top30_total.HC <- rbind(asv_top30.HC , total_30.HC)
# re-order columns by total - Top 3: 
asv_sort_30_HC <- asv_top30_total.HC[,order(-asv_top30_total.HC[which(rownames(asv_top30_total.HC) == 'total_30.HC'), ])] #HCre essentially it is re-ordering tHC total by descending order
colnames(asv_sort_30_HC)
#write.table(asv_sort_30_HC ,  "Family_asv_sort_30_HC_Cellana_HK.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#4: HD.Siphonaria sirius
asv_top30.HD <- asv_top30[16:20,]
View(asv_top30.HD)
# Calculate column totals
total_30.HD <- colSums(asv_top30.HD )
asv_top30_total.HD <- rbind(asv_top30.HD , total_30.HD)
# re-order columns by total - Top 3: 
asv_sort_30_HD <- asv_top30_total.HD[,order(-asv_top30_total.HD[which(rownames(asv_top30_total.HD) == 'total_30.HD'), ])] #HDre essentially it is re-ordering tHD total by descending order
colnames(asv_sort_30_HD)
#write.table(asv_sort_30_HD ,  "Family_asv_sort_30_HD_Siphonaria.txt",  sep = "\t", quote = F, row.names = T ) #Create a new file

#5: HE.Saccostrea.cucullata.HK
asv_top30.HE <- asv_top30[21:24,]
View(asv_top30.HE)
# Calculate column totals
total_30.HE <- colSums(asv_top30.HE )
asv_top30_total.HE <- rbind(asv_top30.HE , total_30.HE)
# re-order columns by total - Top 3: 
asv_sort_30_HE <- asv_top30_total.HE[,order(-asv_top30_total.HE[which(rownames(asv_top30_total.HE) == 'total_30.HE'), ])] #HEre essentially it is re-ordering tHE total by descending order
colnames(asv_sort_30_HE)
#write.table(asv_sort_30_HE ,  "Family_asv_sort_30_HE_cucullata.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#6: HF.Lunella.granulata
asv_top30.HF <- asv_top30[25:29,]
View(asv_top30.HF)
# Calculate column totals
total_30.HF <- colSums(asv_top30.HF )
asv_top30_total.HF <- rbind(asv_top30.HF , total_30.HF)
# re-order columns by total - Top 3:
asv_sort_30_HF <- asv_top30_total.HF[,order(-asv_top30_total.HF[which(rownames(asv_top30_total.HF) == 'total_30.HF'), ])] #HFre essentially it is re-ordering tHF total by descending order
colnames(asv_sort_30_HF)
#write.table(asv_sort_30_HF ,  "Family_asv_sort_30_HF_lunella.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#7: HG.Echinolittorina.trochoides.HK
asv_top30.HG <- asv_top30[30:33,]
View(asv_top30.HG)
# Calculate column totals
total_30.HG <- colSums(asv_top30.HG )
asv_top30_total.HG <- rbind(asv_top30.HG , total_30.HG)
# re-order columns by total - Top 3: 
asv_sort_30_HG <- asv_top30_total.HG[,order(-asv_top30_total.HG[which(rownames(asv_top30_total.HG) == 'total_30.HG'), ])] #HGre essentially it is re-ordering tHG total by descending order
colnames(asv_sort_30_HG)
#write.table(asv_sort_30_HG ,  "Family_asv_sort_30_HG_E_malaccana.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#8: HH.Reishia.clavigera
asv_top30.HH <- asv_top30[34:37,]
View(asv_top30.HH)
# Calculate column totals
total_30.HH <- colSums(asv_top30.HH )
asv_top30_total.HH <- rbind(asv_top30.HH , total_30.HH)
# re-order columns by total - Top 3: 
asv_sort_30_HH <- asv_top30_total.HH[,order(-asv_top30_total.HH[which(rownames(asv_top30_total.HH) == 'total_30.HH'), ])] #HHre essentially it is re-ordering tHH total by descending order
colnames(asv_sort_30_HH)
#write.table(asv_sort_30_HH ,  "Family_asv_sort_30_HH_Reishia.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#9: HJ_Barbatia.virescens.HJ
asv_top30.HJ <- asv_top30[38:42,]
View(asv_top30.HJ)
# Calculate column totals
total_30.HJ <- colSums(asv_top30.HJ )
asv_top30_total.HJ <- rbind(asv_top30.HJ , total_30.HJ)
# re-order columns by total - Top 3: 
asv_sort_30_HJ <- asv_top30_total.HJ[,order(-asv_top30_total.HJ[which(rownames(asv_top30_total.HJ) == 'total_30.HJ'), ])] #HJre essentially it is re-ordering tHJ total by descending order
colnames(asv_sort_30_HJ)
#write.table(asv_sort_30_HJ ,  "Family_asv_sort_30_HJ_Barbatia.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#10: HK.Monodonta.labio_HK
asv_top30.HK <- asv_top30[43:47,]
View(asv_top30.HK)
# Calculate column totals
total_30.HK <- colSums(asv_top30.HK )
asv_top30_total.HK <- rbind(asv_top30.HK , total_30.HK)
# re-order columns by total - Top 3: 
asv_sort_30_HK <- asv_top30_total.HK[,order(-asv_top30_total.HK[which(rownames(asv_top30_total.HK) == 'total_30.HK'), ])] #HKre essentially it is re-ordering tHK total by descending order
colnames(asv_sort_30_HK)
#write.table(asv_sort_30_HK ,  "Family_asv_sort_30_HK.txt_Monodonta", sep = "\t", quote = F, row.names = T ) #Create a new file

#11: TB.Cellana.toreuma.TH
asv_top30.TB <- asv_top30[48:52,]
View(asv_top30.TB)
# Calculate column totals
total_30.TB <- colSums(asv_top30.TB )
asv_top30_total.TB <- rbind(asv_top30.TB , total_30.TB)
# re-order columns by total - Top 3: 
asv_sort_30_TB <- asv_top30_total.TB[,order(-asv_top30_total.TB[which(rownames(asv_top30_total.TB) == 'total_30.TB'), ])] #TBre essentially it is re-ordering tTB total by descending order
colnames(asv_sort_30_TB)
#write.table(asv_sort_30_TB ,  "Family_asv_sort_30_TB_Cellana_TH.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#12: TC.Saccostrea.mordax.TH
asv_top30.TC <- asv_top30[53:56,]
View(asv_top30.TC)
# Calculate column totals
total_30.TC <- colSums(asv_top30.TC )
asv_top30_total.TC <- rbind(asv_top30.TC , total_30.TC)
# re-order columns by total - Top 3: 
asv_sort_30_TC <- asv_top30_total.TC[,order(-asv_top30_total.TC[which(rownames(asv_top30_total.TC) == 'total_30.TC'), ])] #TCre essentially it is re-ordering tTC total by descending order
colnames(asv_sort_30_TC)
#write.table(asv_sort_30_TC ,  "Family_asv_sort_30_TC_mordax.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#13: TD.Echinolittorina.trochoides.TH
asv_top30.TD <- asv_top30[57:61,]
View(asv_top30.TD)
# Calculate column totals
total_30.TD <- colSums(asv_top30.TD )
asv_top30_total.TD <- rbind(asv_top30.TD , total_30.TD)
# re-order columns by total - Top 3: 
asv_sort_30_TD <- asv_top30_total.TD[,order(-asv_top30_total.TD[which(rownames(asv_top30_total.TD) == 'total_30.TD'), ])] #TDre essentially it is re-ordering tTD total by descending order
colnames(asv_sort_30_TD)
#write.table(asv_sort_30_TD ,  "Family_asv_sort_30_TD_malaccana_TH.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#14: TF_Mytilisepta.virgata.TF
asv_top30.TF <- asv_top30[62:66,]
View(asv_top30.TF)
# Calculate column totals
total_30.TF <- colSums(asv_top30.TF )
asv_top30_total.TF <- rbind(asv_top30.TF , total_30.TF)
# re-order columns by total - Top 3: 
asv_sort_30_TF <- asv_top30_total.TF[,order(-asv_top30_total.TF[which(rownames(asv_top30_total.TF) == 'total_30.TF'), ])] #TFre essentially it is re-ordering tTF total by descending order
colnames(asv_sort_30_TF)
#write.table(asv_sort_30_TF ,  "Family_asv_sort_30_TF_Mytilisepta.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#15: TG.Tenguella.musiva
asv_top30.TG <- asv_top30[67:69,]
View(asv_top30.TG)
# Calculate column totals
total_30.TG <- colSums(asv_top30.TG )
asv_top30_total.TG <- rbind(asv_top30.TG , total_30.TG)
# re-order columns by total - Top 3: 
asv_sort_30_TG <- asv_top30_total.TG[,order(-asv_top30_total.TG[which(rownames(asv_top30_total.TG) == 'total_30.TG'), ])] #TGre essentially it is re-ordering tTG total by descending order
colnames(asv_sort_30_TG)
#write.table(asv_sort_30_TG ,  "Family_asv_sort_30_TG_Tenguella.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#16: TJ.Echinolittorina.radiata.TJ
asv_top30.TJ <- asv_top30[70:74,]
View(asv_top30.TJ)
# Calculate column totals
total_30.TJ <- colSums(asv_top30.TJ )
asv_top30_total.TJ <- rbind(asv_top30.TJ , total_30.TJ)
# re-order columns by total - Top 3: 
asv_sort_30_TJ <- asv_top30_total.TJ[,order(-asv_top30_total.TJ[which(rownames(asv_top30_total.TJ) == 'total_30.TJ'), ])] #TJre essentially it is re-ordering tTJ total by descending order
colnames(asv_sort_30_TJ)
#write.table(asv_sort_30_TJ ,  "Family_asv_sort_30_TJ_radiata.txt", sep = "\t", quote = F, row.names = T ) #Create a new file

#17: A.water.TH
asv_top30.water <- asv_top30[75:77,]
View(asv_top30.water)
# Calculate column totals
total_30.water <- colSums(asv_top30.water )
asv_top30_total.water <- rbind(asv_top30.water , total_30.water)
# re-order columns by total - Top 3: 
asv_sort_30_water <- asv_top30_total.water[,order(-asv_top30_total.water[which(rownames(asv_top30_total.water) == 'total_30.water'), ])] #waterre essentially it is re-ordering twater total by descending order
colnames(asv_sort_30_water)
#write.table(asv_sort_30_water ,  "Family_asv_sort_30_water.txt", sep = "\t", quote = F, row.names = T ) #Create a new file









