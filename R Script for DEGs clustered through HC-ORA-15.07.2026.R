#Import data with gene contained in 3 clusters from hierarchical (DEGs-Airway)
#My working directory-setwd("~/Training/Enrichment-DEG-ORA")


#Cluster 01
library(openxlsx)
cl1<-read.xlsx("./cl1-Dendrogram.xlsx",rowNames = TRUE)
#Inspect rows 1:6 and columns 1:6
cl1[1:6,1:6] 
#Extract gene names present in cluster 1 using rownames()
genes_cl1<-rownames(cl1)
paste0("Number of genes in cluster 1 = ",length(genes_cl1))

#Convert SYMBOL into Gene Entrez ID using mapIds() function implemented in the “org.Hs.eg.db” annotation package 
#Columns that can be used as keys (for conversion)
library(org.Hs.eg.db) 
keytypes(org.Hs.eg.db) 

#Conversion of the from ENSEMBL IDs included in the “genes_cl1” vector to Gene Entrez IDs 
genes_cl1_Entrez<-mapIds(org.Hs.eg.db, 
                          keys=genes_cl1, 
                          keytype="ENSEMBL", # input ID
                          column="ENTREZID") # output ID 

head(genes_cl1_Entrez)


#Remove duplicates
duplicated<-genes_cl1_Entrez[duplicated(genes_cl1_Entrez)] 
paste0("Number ofduplicated values (NA and duplicated Gene Entrez ID) = ",length(duplicated))

NAs<-duplicated[is.na(duplicated)] # check for NA values
paste0("Number of NAs = ",length(NAs))

duplicated_Gene_Entrez_ID<-duplicated[!is.na(duplicated)] # check for duplicated IDs
paste0("Number of duplicatedGene Entrez ID = ",length(duplicated_Gene_Entrez_ID))

#Remove the elements with duplicated names
genes_cl1_Entrez<-genes_cl1_Entrez[!duplicated(genes_cl1_Entrez) &                        
                                       !duplicated(genes_cl1_Entrez,fromLast=TRUE)]
paste0("Number of genes in genes_cl1_Entrez vectorafterremoving elements with duplicated names= ", 
       length(genes_cl1_Entrez))


#ORA using ReactomePA
library(ReactomePA)
library(ggplot2)


#Run enrichPathway function by keeping pvalueCutoff=0.05
ora<-enrichPathway(gene=genes_cl1_Entrez,pvalueCutoff=0.05, readable=T)  
head(ora@result[,c(2,3,5,6,7)])


#Re-run Reactome enrichment with a relaxed threshold
ora <- enrichPathway(gene = genes_cl1_Entrez, 
                     pvalueCutoff = 0.2, 
                     qvalueCutoff = 0.2, # Added to prevent q-value filtering
                     readable = TRUE)
head(ora@result[,c(2,3,5,6,7)])


# View the description alongside the raw and adjusted p-values
head(ora@result[, c("Description", "pvalue", "p.adjust")])


# Check if there are any rows in your enrichment result
summary(as.data.frame(ora))


#dotplot and barpolt were drawan results with pvalueCutoff=0.2
#Visualization ORA results: dotplot(enrichplot)
dotplot(ora,x="GeneRatio")+ 
  ggtitle("Dotplot ORA withgenes in cluster 1")+    
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))



#Visualization ORA results: barplot(enrichplot)
barplot(ora,x="GeneRatio")+ 
  ggtitle("Barplot ORA withgenes in cluster 1")+ 
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))

#Cluster 02
library(openxlsx)
cl2<-read.xlsx("./cl2-Dendrogram.xlsx",rowNames = TRUE)
#Inspect rows 1:6 and columns 1:6
cl2[1:6,1:6] 
#Extract gene names present in cluster 2 using rownames()
genes_cl2<-rownames(cl2)
paste0("Number of genes in cluster 2 = ",length(genes_cl2))

#Convert SYMBOL into Gene Entrez ID using mapIds() function implemented in the “org.Hs.eg.db” annotation package 
#Columns that can be used as keys (for conversion)
library(org.Hs.eg.db) 
keytypes(org.Hs.eg.db) 

#Conversion of the from ENSEMBL IDs included in the “genes_cl2” vector to Gene Entrez IDs 
genes_cl2_Entrez<-mapIds(org.Hs.eg.db, 
                         keys=genes_cl2, 
                         keytype="ENSEMBL", # input ID
                         column="ENTREZID") # output ID 

head(genes_cl2_Entrez)


#Remove duplicates
duplicated<-genes_cl2_Entrez[duplicated(genes_cl2_Entrez)] 
paste0("Number ofduplicated values (NA and duplicated Gene Entrez ID) = ",length(duplicated))

NAs<-duplicated[is.na(duplicated)] # check for NA values
paste0("Number of NAs = ",length(NAs))

duplicated_Gene_Entrez_ID<-duplicated[!is.na(duplicated)] # check for duplicated IDs
paste0("Number of duplicatedGene Entrez ID = ",length(duplicated_Gene_Entrez_ID))

#Remove the elements with duplicated names
genes_cl2_Entrez<-genes_cl2_Entrez[!duplicated(genes_cl2_Entrez) &                        
                                     !duplicated(genes_cl2_Entrez,fromLast=TRUE)]
paste0("Number of genes in genes_cl2_Entrez vectorafterremoving elements with duplicated names= ", 
       length(genes_cl2_Entrez))


#ORA using ReactomePA
library(ReactomePA)
library(ggplot2)


#Run enrichPathway function by keeping pvalueCutoff=0.05
ora2<-enrichPathway(gene=genes_cl2_Entrez,pvalueCutoff=0.05, readable=T)  
head(ora2@result[,c(2,3,5,6,7)])

#View the description alongside the raw and adjusted p-values
head(ora2@result[, c("Description", "pvalue", "p.adjust")])

# Check if there are any rows in enrichment result
summary(as.data.frame(ora2))

#Re-run Reactome enrichment with a relaxed threshold
ora2 <- enrichPathway(gene = genes_cl2_Entrez, 
                     pvalueCutoff = 0.2, 
                     qvalueCutoff = 0.2, # Added to prevent q-value filtering
                     readable = TRUE)
head(ora2@result[,c(2,3,5,6,7)])

#View the description alongside the raw and adjusted p-values
head(ora2@result[, c("Description", "pvalue", "p.adjust")])


# Check if there are any rows in enrichment result
summary(as.data.frame(ora2))


#Re-run Reactome enrichment for cluster 2 with a slightly relaxed cutoff (0.25)
ora2 <- enrichPathway(gene = genes_cl2_Entrez, 
                      pvalueCutoff = 0.25,      # Increased to catch the 0.202 pathways
                      qvalueCutoff = 0.25, 
                      readable = TRUE)

head(ora2@result[,c(2,3,5,6,7)])

#View the description alongside the raw and adjusted p-values
head(ora2@result[, c("Description", "pvalue", "p.adjust")])
# Verify that it now contains rows
nrow(as.data.frame(ora2))



#dotplot and barpolt were drawan results with pvalueCutoff=0.25
#Visualization ORA results: dotplot(enrichplot)
dotplot(ora2,x="GeneRatio")+ 
  ggtitle("Dotplot ORA withgenes in cluster 2")+    
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))

#Cluster 03
library(openxlsx)
cl3<-read.xlsx("./cl3-Dendrogram.xlsx",rowNames = TRUE)
#Inspect rows 1:6 and columns 1:6
cl3[1:6,1:6] 
#Extract gene names present in cluster 2 using rownames()
genes_cl3<-rownames(cl3)
paste0("Number of genes in cluster 3 = ",length(genes_cl3))

#Convert SYMBOL into Gene Entrez ID using mapIds() function implemented in the “org.Hs.eg.db” annotation package 
#Columns that can be used as keys (for conversion)
library(org.Hs.eg.db) 
keytypes(org.Hs.eg.db) 

#Conversion of the from ENSEMBL IDs included in the “genes_cl3” vector to Gene Entrez IDs 
genes_cl3_Entrez<-mapIds(org.Hs.eg.db, 
                         keys=genes_cl3, 
                         keytype="ENSEMBL", # input ID
                         column="ENTREZID") # output ID 

head(genes_cl3_Entrez)


#Remove duplicates
duplicated<-genes_cl3_Entrez[duplicated(genes_cl3_Entrez)] 
paste0("Number ofduplicated values (NA and duplicated Gene Entrez ID) = ",length(duplicated))

NAs<-duplicated[is.na(duplicated)] # check for NA values
paste0("Number of NAs = ",length(NAs))

duplicated_Gene_Entrez_ID<-duplicated[!is.na(duplicated)] # check for duplicated IDs
paste0("Number of duplicatedGene Entrez ID = ",length(duplicated_Gene_Entrez_ID))

#Remove the elements with duplicated names
genes_cl3_Entrez<-genes_cl3_Entrez[!duplicated(genes_cl3_Entrez) &                        
                                     !duplicated(genes_cl3_Entrez,fromLast=TRUE)]
paste0("Number of genes in genes_cl3_Entrez vectorafterremoving elements with duplicated names= ", 
       length(genes_cl3_Entrez))


#ORA using ReactomePA
library(ReactomePA)
library(ggplot2)


#Run enrichPathway function by keeping pvalueCutoff=0.05
ora3<-enrichPathway(gene=genes_cl3_Entrez,pvalueCutoff=0.05, readable=T)  
head(ora3@result[,c(2,3,5,6,7)])

#View the description alongside the raw and adjusted p-values
head(ora3@result[, c("Description", "pvalue", "p.adjust")])

# Check if there are any rows in enrichment result
summary(as.data.frame(ora3))

#Re-run Reactome enrichment with a relaxed threshold
ora3 <- enrichPathway(gene = genes_cl3_Entrez, 
                      pvalueCutoff = 0.2, 
                      qvalueCutoff = 0.2, # Added to prevent q-value filtering
                      readable = TRUE)
head(ora3@result[,c(2,3,5,6,7)])

#View the description alongside the raw and adjusted p-values
head(ora3@result[, c("Description", "pvalue", "p.adjust")])


# Check if there are any rows in enrichment result
summary(as.data.frame(ora3))



# Verify that it now contains rows
nrow(as.data.frame(ora3))



#dotplot and barpolt were drawan results with pvalueCutoff=0.2
#Visualization ORA results: dotplot(enrichplot)
dotplot(ora3,x="GeneRatio")+ 
  ggtitle("Dotplot ORA withgenes in cluster 3")+    
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))


#Visualization ORA results: barplot(enrichplot)
barplot(ora3,x="GeneRatio")+ 
  ggtitle("Barplot ORA withgenes in cluster3")+ 
  theme(plot.title = element_text(color="black", size=14, face="bold.italic"))




