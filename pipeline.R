getwd()
setwd("F:/PhD/Other than PhD/DR data")
data <- read.delim('GSE221521_193_samples_gene_expression.xls', header = T)
library(tidyverse)
data <- data %>%
  select(-c(2:198))

# Print the first column header name
#print(colnames(data)[1])

data<-data %>%
  gather(key = 'samples', value = 'counts',-gene_id)

data<-data %>%mutate(samples = gsub('\\_count', '', samples))
data<-data %>%mutate(samples = gsub('\\_', '', samples))

#remove the environmental variable
#rm(data1)

library(GEOquery)

geo_id <- "GSE221521"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
phenoData<-phenoData %>%
  select(-c(18:44))

phenoData<-phenoData %>%
  select(-c(3:16))

phenoData<-phenoData %>% 
  mutate(
    title = gsub('\\_', '', 
                 gsub('\\leukocytes, Control group ', '', 
                      gsub('\\leukocytes, DR group ', '', 
                           gsub('\\leukocytes, DM group ', '', title)))),
    description = gsub('\\leukocytes, ', '', description))


write.csv(data, "phrnodata_of_193_data_after_processing.csv", row.names = FALSE)

data<-data%>%
  inner_join(phenoData, by = c('samples' = 'title'))

filtered_data <- data %>%
  filter(description != "DM group")

#write.csv(data, "filtered_data_excluding_filtered_data.csv", row.names = FALSE)


filtered_data <- filtered_data %>%
  mutate(samples_description = paste(samples, description))
head(filtered_data)

filtered_data <- filtered_data %>%
  select(-samples, -geo_accession, -description)

#write.csv(data, "filtered_data_excluding_3_columns.csv", row.names = FALSE)


filtered_data <- filtered_data %>%
  spread(key = 'samples_description', value = 'counts')%>%
  column_to_rownames(var = 'gene_id')

#write.csv(filtered_data, "wide_format_data_for_analysis.csv", row.names = FALSE)



#write.csv(data, "193_merged_with_phenodata.csv", row.names = FALSE)
  
  
#write.csv(data, "193_data_long.csv", row.names = FALSE)
#data1 <- read.delim('193_data_long.csv', sep = ",", header = T)
#write.csv(data, "193_data_long_after_removing_different_symbols.csv", row.names = FALSE)


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#  BiocManager::install("WGCNA")

#version check
#packageVersion("WGCNA")

# detect outlier genes
# goodSamplesGenes is a function from WGCNA where it needs data rows should be samples and columns should be
# genes
gsg <- goodSamplesGenes(t(filtered_data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

# To keep only the good genes based on the gsg object returned by goodSamplesGenes. The gsg$goodGenes vector is a logical vector indicating which genes are good (TRUE) and which are not (FALSE).
filtered_data <- filtered_data[gsg$goodGenes == TRUE,]
write.csv(filtered_data, "filtered_data_removing_bad_genes.csv", row.names = FALSE)

# detect outlier samples - hierarchical clustering - method 1
#htree <- hclust(dist(t(filtered_data)), method = "average")
#plot(htree)


# pca - method 2

#pca <- prcomp(t(filtered_data))
#pca.dat <- pca$x


#pca.var <- pca$sdev^2
#pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

#pca.dat <- as.data.frame(pca.dat)

#ggplot(pca.dat, aes(PC1, PC2)) +
  #geom_point() +
  #geom_text(label = rownames(pca.dat)) +
  #labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       #y = paste0('PC2: ', pca.var.percent[2], ' %'))



#final_data<-read.delim('filtered_data_removing_bad_genes.csv', sep = ",", header = T, )
#head(final_data)


#For read first column as row names
final_data<-read.delim('filtered_data_removing_bad_genes.csv', sep = ",", header = TRUE, row.names = 1)
final_phenoData<-read.delim('phenoData.csv', sep = ",", header = TRUE, row.names = 1)


# making the rownames and column names identical
all(rownames(final_phenoData) %in% colnames(final_data))
all(rownames(final_phenoData) == colnames(final_data))


#If the phenodata rows not matched with columns of dataset then below should follow
# Get the common names between row names of final_phenoData and column names of final_data

matched_names <- intersect(rownames(final_phenoData), colnames(final_data))



# subsets final_phenoData to keep only the rows whose names are in matched_names, preserving all columns and ensuring the result remains a data frame even if it has only one row or one column. 
final_phenoData <- final_phenoData[matched_names, , drop = FALSE]

#The empty space before the comma indicates that we are not subsetting the rows, so we keep all rows.
#matched_names: This is a character vector containing the names of the columns that we want to keep.

final_data <- final_data[, matched_names]

#write.csv(final_data, "final_data_matched_with_phenodata.csv")
#write.csv(final_phenoData, "final_phenoData_matched_with_finaldata.csv")

samples_to_exclude<-c('RNA335','RJS040','RNA224','RNA441','RNA51','RNA50','RNA49','RNA124','RNA160','RNA310','RJS014')

final_data_after_exclude <- final_data[,!(colnames(final_data) %in% samples_to_exclude)]
head(final_data_after_exclude)

colData <- final_phenoData %>% 
  filter(!row.names(.) %in% samples_to_exclude)

# Convert 'description' to factor
colData$description <- as.factor(colData$description)


# create dds
dds <- DESeqDataSetFromMatrix(countData = final_data_after_exclude,
                              colData = colData,
                              design = ~ description)

head(colData)

# Run the DESeq2 analysis
real_dds <- DESeq(dds)
dds_results<-results(real_dds)
#write.csv(dds_results, file = "dds_results.csv")

## remove all genes with counts < 15 in more than 75% of samples (108*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 81,]
nrow(dds75) # 13949 genes

normalized_dds_norm<-DESeq(dds75)
dds75_results<-results(normalized_dds_norm)
head(dds75_results)

head(normalized_dds_norm)


dds75_results<-results(normalized_dds_norm)
head(dds75_results)


dds75_results_df <- as.data.frame(dds75_results)
# Save the results data frame to a CSV file
write.csv(dds75_results_df, file = "dds75_results.csv")

normalized_counts <- counts(normalized_dds_norm, normalized = TRUE)
head(normalized_counts)
write.csv(normalized_counts, "normalized_expression_data.csv")

# Save the normalized expression data for CIBERSORT input
write.table(normalized_counts, file = "normalized_expression_data_for_cibersort.txt", sep = "\t", quote = FALSE, row.names = TRUE)

log2_normalized_counts_only <- log2(normalized_counts + 0.1)
head(log2_normalized_counts_only)

#log2_normalized_counts_only_df <- as.data.frame(log2_normalized_counts_only)
#write.csv(log2_normalized_counts_only_df, file = "log2_normalized_counts_only_df.csv")


# Save metadata (colData) to CSV
write.csv(as.data.frame(colData(dds75)), "metadata.csv")

write.csv(as.data.frame(colData), "metadata.csv")



# Filter for significant DEGs based on adjusted p-value (padj <= 0.05)
significant_genes <- subset(dds75_results, padj <= 0.05)
nrow(significant_genes) #3614
write.csv(rownames(significant_genes), file = "significant_genes", row.names = FALSE)

# Filter for significant upregulated genes (log2FoldChange > 0 and padj <= 0.05)
#upregulated_genes <- subset(dds75_results, log2FoldChange > 0 & padj <= 0.05)
# Count upregulated genes (log2FoldChange > 0 and padj <= 0.05)
#num_upregulated <- sum(dds75_results$log2FoldChange > 0 & dds75_results$padj <= 0.05)


# Identify upregulated genes (log2FoldChange > 0)
upregulated_genes <- subset(significant_genes, log2FoldChange > 0.5)
nrow(upregulated_genes) #389 Upregulated genes
write.csv(rownames(upregulated_genes), file = "upregulated_genes.csv", row.names = FALSE)

# Identify downregulated genes (log2FoldChange < 0)
downregulated_genes <- subset(significant_genes, log2FoldChange < -0.5)
nrow(downregulated_genes) #109 Downregulated genes
write.csv(rownames(downregulated_genes), file = "downregulated_genes.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# Convert DESeqResults object to data frame
significant_genes <- as.data.frame(subset(dds75_results, padj <= 0.05))
nrow(significant_genes) # 3614
write.csv(rownames(significant_genes), file = "significant_genes.csv", row.names = FALSE)

# Identify upregulated genes (log2FoldChange > 0.5)
upregulated_genes <- subset(significant_genes, log2FoldChange > 0.5)
nrow(upregulated_genes) # 389 Upregulated genes
write.csv(rownames(upregulated_genes), file = "upregulated_genes.csv", row.names = FALSE)

# Identify downregulated genes (log2FoldChange < -0.5)
downregulated_genes <- subset(significant_genes, log2FoldChange < -0.5)
nrow(downregulated_genes) # 109 Downregulated genes
write.csv(rownames(downregulated_genes), file = "downregulated_genes.csv", row.names = FALSE)

# Extract top 25 upregulated genes based on log2FoldChange
top_25_upregulated <- upregulated_genes %>%
  arrange(desc(log2FoldChange)) %>%
  head(25)

# Extract top 25 downregulated genes based on log2FoldChange
top_25_downregulated <- downregulated_genes %>%
  arrange(log2FoldChange) %>%
  head(25)

# Save top 25 upregulated and downregulated genes to CSV files
write.csv(rownames(top_25_upregulated), file = "top_25_upregulated_genes.csv", row.names = FALSE)
write.csv(rownames(top_25_downregulated), file = "top_25_downregulated_genes.csv", row.names = FALSE)

# Combine top 25 upregulated and downregulated genes into one dataframe
top_50_genes <- rbind(top_25_upregulated, top_25_downregulated)

# Save combined top 50 genes to a CSV file
write.csv(rownames(top_50_genes), file = "top_50_genes.csv", row.names = FALSE)
#-----------------------------------------------------------------------

# Initialize a new plot
plot.new()

# Define colors for upregulated and downregulated genes based on log2FC thresholds
colors <- ifelse(dds75_results$padj <= 0.05 & dds75_results$log2FoldChange > 0.5, "red",    # Upregulated genes in red
                 ifelse(dds75_results$padj <= 0.05 & dds75_results$log2FoldChange < -0.5, "blue", "black"))  # Downregulated genes in blue, others in black

# Create volcano plot
plot(dds75_results$log2FoldChange, -log10(dds75_results$padj), pch=20, col=colors, 
     main="Volcano Plot of DEGs",
     xlab="log2(Fold Change)", ylab="-log10(Adjusted p-value)")

# Add significance threshold lines
abline(h = -log10(0.05), col = "gray", lty = 2)  # Adjusted p-value threshold

# Add vertical lines for log2FC thresholds
abline(v = c(0.5, -0.5), col = c("red", "blue"), lty = 2)  # Vertical lines for log2FC thresholds

# Add legend
legend("topright", legend = c("Upregulated", "Downregulated"),
       col = c("red", "blue"), pch = 20, cex = 0.8)



#---------------------------------------------------------------------------------------------


# perform variance stabilization
dds_norm <- vst(dds75)




# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices


# visualization to pick power


a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


install.packages("gridExtra")
library(gridExtra)

grid.arrange(a1, a2, nrow = 2)
grid.arrange(a1, nrow = 1)
grid.arrange(a2, nrow = 1)

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 12
temp_cor <- cor
cor <- WGCNA::cor


install.packages("doParallel")
install.packages("foreach")
library(doParallel)
library(foreach)
library(WGCNA)

# Detect the number of available cores
num_cores <- detectCores()

# Register the parallel backend
cl <- makeCluster(num_cores - 1)  # Reserve one core for the OS
registerDoParallel(cl)

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor



#Module Eigengenes 
module_eigengenes <- bwnet$MEs

head(module_eigengenes)
# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

#------------------------------------------------------------------------------------------------
# Create traits file - binarize categorical variables for DR and Control
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(description == "DR", 1, 0)) %>% 
  select(disease_state_bin)
head(traits)

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


#---------------------------------------------------------------------------------------
# Calculate module-trait correlations
module_trait_corr <- cor(module_eigengenes, traits, use = 'p')
head(module_trait_corr)
module_trait_corr_pvals <- corPvalueStudent(module_trait_corr, nSamples)


#---------------------------------------------------------------------------------------


moduleColors <- bwnet$colors
# Find the genes belonging to the top associated module
top_module_genes <- names(moduleColors[moduleColors == "purple"])

head(top_module_genes)
write.csv(top_module_genes, file = "top_module_genes.csv", row.names = FALSE, quote = FALSE)

#---------------------------------------------------------------------------------------

# Visualize module-trait association as a heatmap
heatmap_data <- merge(module_eigengenes, traits, by = 'row.names')
heatmap_data <- heatmap_data %>% column_to_rownames(var = 'Row.names')
head(heatmap_data)

#--------------------------------------------------------------------------------------------
# Perform t-tests for each module
results <- heatmap_data %>%
#  select(-sample) %>%
  group_by(disease_state_bin) %>%
  summarise(across(starts_with("ME"), list(mean = ~mean(.), sd = ~sd(.)))) %>%
  pivot_longer(-disease_state_bin, names_to = c("Module", "Metric"), names_sep = "_") %>%
  pivot_wider(names_from = Metric, values_from = value)
#---------------------------------------------------------------------------------------------






heatmap_data_long <- melt(heatmap_data, id.vars = c("sample", "disease_state_bin"),
                          variable.name = "Module", value.name = "Value")


ggplot(heatmap_data_long, aes(x = factor(disease_state_bin), y = Module, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  labs(title = "Heatmap of Module Eigengenes by Disease State",
       x = "Disease State Bin", y = "Module", fill = "Value")




# Convert 'disease_state_bin' to a factor with meaningful labels
heatmap_data_long1 <- heatmap_data_long %>%
  mutate(disease_state_bin = factor(disease_state_bin, levels = c(0, 1), labels = c("Control", "DR")))

# Plotting
 
ggplot(heatmap_data_long1, aes(x = disease_state_bin, y = Module, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-0.5, 0.5)) +
  theme_minimal() +
  labs(title = "Heatmap of Module Eigengenes by Disease State",
       x = "Disease State", y = "Module", fill = "Value")


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Assuming `heatmap_data` is already in the correct format with numeric values
# and rownames(heatmap_data) are the sample names.

# If you have sample annotations like `disease_state_bin`
annotation <- data.frame(disease_state_bin = factor(heatmap_data$disease_state_bin, labels = c("Control", "DR")))
rownames(annotation) <- rownames(heatmap_data)


# Remove non-numeric columns if necessary
heatmap_data_numeric <- heatmap_data[ , !colnames(heatmap_data) %in% c("sample", "disease_state_bin")]

pheatmap(
  heatmap_data_numeric,
  annotation_row = annotation,
  color = colorRampPalette(c("blue", "white", "red"))(50), # Custom color scale
  main = "Heatmap of Module Eigengenes by Disease State"
)



#----------------------------------------------------------------------------------------------------


# Aggregate the data by disease_state_bin
agg_heatmap_data <- heatmap_data %>%
  group_by(disease_state_bin) %>%
  summarize(across(starts_with("ME"), mean)) %>%
  as.data.frame() # Convert to data frame

# Convert disease_state_bin to factor for proper labeling
agg_heatmap_data$disease_state_bin <- factor(agg_heatmap_data$disease_state_bin, labels = c("Control", "DR"))

# Set row names as the disease_state_bin for easier plotting
rownames(agg_heatmap_data) <- agg_heatmap_data$disease_state_bin

# Remove the disease_state_bin column for heatmap plotting
agg_heatmap_data <- agg_heatmap_data[, -1]



# Plot heatmap with aggregated data
pheatmap(
  as.matrix(agg_heatmap_data),
  cluster_cols = FALSE, # Don't cluster columns
  cluster_rows = FALSE, # Don't cluster rows
  color = colorRampPalette(c("blue", "white", "red"))(50), # Custom color scale
  main = "Heatmap of Module Eigengenes by Disease State",
  labels_row = rownames(agg_heatmap_data), # Use group labels (Control, DR)
  show_rownames = TRUE, # Show group names
  show_colnames = TRUE  # Show module names
)


agg_heatmap_data_t <- t(agg_heatmap_data)
pheatmap(
  agg_heatmap_data_t,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  fontsize = 10,
  main = "Module-Trait Relationships"
)

#----------------------------------------------------------------------------------------------------



# Calculate module membership and associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.
module_membership_measure <- cor(module_eigengenes, norm.counts, use = 'p')
module_membership_measure_pvals <- corPvalueStudent(module_membership_measure, nSamples)
module_membership_measure_pvals[1:10,1:10]

#----------------------------------------------------------------------------------------


# Calculate gene significance and associated p-values
gene_signf_corr <- cor(norm.counts, traits$disease_state_bin, use = 'p')
gene_signf_corr_pvals <- corPvalueStudent(gene_signf_corr, nSamples)




# Identify top genes
top_genes <- gene_signf_corr_pvals %>% as.data.frame() %>% arrange(V1) %>% head(25)
print(top_genes)

#-------------------------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("clusterProfiler")
# library(clusterProfiler)


# Load necessary libraries
library(pheatmap)

# Assuming dds75_results is your DESeq2 results dataframe
# and log2_normalized_counts_only_df is your expression data matrix

# Filter significant genes with adjusted p-value <= 0.05
#significant_genes <- subset(dds75_results, padj <= 0.05)

# Sort significant genes by adjusted p-value in increasing order and select the top 50
top_significant_genes <- head(significant_genes[order(significant_genes$padj), ], 50)

# Extract the gene symbols (or IDs) of the top 50 significant genes
top_genes <- rownames(top_significant_genes)

# Subset expression data for these top 50 significant genes
top_exp_data <- final_data_after_exclude[top_genes, ]
top_exp_data1<-log2_normalized_counts_only[top_genes, ]
head(final_data_after_exclude)
head(top_exp_data)
annotation_col <- colData

rownames(annotation_col1) <- colnames(top_exp_data)



# Plot the heatmap
pheatmap(top_exp_data1,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Heatmap of Top 50 Significant DEGs",
         annotation_col = annotation_col)


all_degs <- rbind(upregulated_genes, downregulated_genes)

log2_normalized_counts_only_five <- log2(normalized_counts - 0.5)
deg_genes <- rownames(all_degs)
deg_exp_data3 <- log2_normalized_counts_only[deg_genes, ]


# Plot the heatmap
pheatmap(deg_exp_data3,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Heatmap of Top 50 Significant DEGs",
         annotation_col = annotation_col)

#-------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

head(all_degs)
all_degs_df <- data.frame(all_degs)
ensembl_ids <- rownames(all_degs_df)  

entrez_ids <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(entrez_ids)
# # Convert gene symbols to Entrez IDs
# gene_list <- bitr(all_degs$GeneSymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)




# GO Biological Process
ego_bp <- enrichGO(gene = entrez_ids$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)



#write.csv(ego_bp, "ego_bp.csv")

# Biological Process
dotplot(ego_bp, showCategory = 10) + ggtitle("GO Biological Processes")

#-------------------------------------------------------------------------------

# GO Cellular Component
ego_cc <- enrichGO(gene = entrez_ids$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.04,
                   qvalueCutoff = 0.04)

# Cellular Component
dotplot(ego_cc, showCategory = 10) + ggtitle("GO Cellular Components")

if (is.null(ego_bp) || nrow(ego_bp) == 0) {
  cat("No enriched GO terms were found.")
} else {
  print(head(ego_bp))
}



#-------------------------------------------------------------------------------

pvals <- all_degs_df$padj
pvals <- pvals[!is.na(pvals)]
ggplot(data.frame(padj = pvals), aes(x = padj)) + 
  geom_histogram(bins = 50, fill = "blue", color = "black") + 
  theme_minimal() + 
  xlab("Adjusted p-values") + 
  ggtitle("Distribution of Adjusted p-values")
#-------------------------------------------------------------------------------
# GO Molecular Function
ego_mf <- enrichGO(gene = entrez_ids$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# Molecular Function
dotplot(ego_mf, showCategory = 10) + ggtitle("GO Molecular Functions")

#-------------------------------------------------------------------------------
ekk <- enrichKEGG(gene = entrez_ids$ENTREZID,
                  organism = "hsa",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

#KEGG Pathways
dotplot(ekk, showCategory = 10) + ggtitle("KEGG Pathways")
#-------------------------------------------------------------------------------
# Extract the Ensembl IDs from row names
ensembl_ids <- rownames(all_degs_df)

# Map Ensembl IDs to gene symbols
gene_mapping <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# Merge gene symbols with the original data frame
all_degs_with_symbols <- merge(all_degs_df, gene_mapping, by.x = "row.names", by.y = "ENSEMBL")

# Set the row names to gene symbols
rownames(all_degs_with_symbols) <- all_degs_with_symbols$SYMBOL

# Remove the redundant columns
all_degs_with_symbols <- all_degs_with_symbols[, -which(names(all_degs_with_symbols) %in% c("Row.names", "SYMBOL"))]
head(all_degs_with_symbols)
#--------------------------------(For Box Plot)-----------------------------------------------
# Read the CSV file containing the selected genes
selected_genes <- read.csv("selected_features_elastic_net.csv")

# Read the metadata file containing type information
#metadata <- read.csv("metadata.csv")

# Assuming the column with gene IDs is named "Feature" and type is in "Type"
selected_gene_ids <- selected_genes$Feature

# Convert all_degs_df to a data frame if it's not already
all_degs_df <- data.frame(all_degs)

# Extract log2 fold changes and merge with metadata
plot_data <- data.frame(
  Gene = rep(selected_gene_ids, each = nrow(colData)),
  Log2FoldChange = unlist(lapply(selected_gene_ids, function(gene) all_degs_df[gene, "log2FoldChange"])),
  Type = rep(colData$description, times = length(selected_gene_ids))
)


# Plot the boxplot using ggplot2
ggplot(plot_data, aes(x = Gene, y = Log2FoldChange, fill = Type)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("Gene") +
  ylab("Log2 Fold Change") +
  ggtitle("Boxplot of Log2 Fold Changes for Selected Genes by Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("blue", "red"))  # Adjust colors if needed
#-------------------------------------------------------------------------------
# # Install necessary packages if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# if (!require("clusterProfiler"))
#   BiocManager::install("clusterProfiler")
# if (!require("biomaRt"))
#   BiocManager::install("biomaRt")
#   BiocManager::install("pathview")
# 
# # Load your DE
# 
# # Load necessary libraries
# library(clusterProfiler)
# library(biomaRt)
# library(ggplot2)
# library(pathview)
# 
# # Load your DE
# 
# # Load your DEGs
# degs_df <- read.csv("all_degs.csv", header = TRUE, stringsAsFactors = FALSE)
# 
# # Extract the Ensembl gene IDs (assuming they are in a column named 'ensembl_id')
# gene_list_df <- degs_df$ensembl_id
# 
# # Use biomaRt to retrieve GO annotations from Ensembl
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# # Perform GO enrichment analysis for Biological Process (BP), Cellular Component (CC), and Molecular Function (MF)
# go_bp <- enrichGO(gene = gene_list_df, 
#                   OrgDb = org.Hs.eg.db, 
#                   keyType = "ENSEMBL", 
#                   ont = "BP", 
#                   pAdjustMethod = "BH", 
#                   pvalueCutoff = 0.05, 
#                   qvalueCutoff = 0.05)
# 
# dotplot(go_bp, showCategory = 10) + ggtitle("GO Biological Process") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# go_cc <- enrichGO(gene = gene_list, 
#                   OrgDb = org.Hs.eg.db, 
#                   keyType = "ENSEMBL", 
#                   ont = "CC", 
#                   pAdjustMethod = "BH", 
#                   pvalueCutoff = 0.05, 
#                   qvalueCutoff = 0.05)
# 
# go_mf <- enrichGO(gene = gene_list, 
#                   OrgDb = org.Hs.eg.db, 
#                   keyType = "ENSEMBL", 
#                   ont = "MF", 
#                   pAdjustMethod = "BH", 
#                   pvalueCutoff = 0.05, 
#                   qvalueCutoff = 0.05)
# 
# # Plotting the dot plot for each GO category
# 
# dotplot(go_cc, showCategory = 10) + ggtitle("GO Cellular Component") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dotplot(go_mf, showCategory = 10) + ggtitle("GO Molecular Function") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# # Perform KEGG pathway enrichment analysis
# kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID,
#                           organism = 'hsa',
#                           pAdjustMethod = "BH",
#                           pvalueCutoff = 0.05,
#                           qvalueCutoff = 0.05)
# 
# # Create a dot plot for KEGG pathway enrichment
# dotplot(kegg_enrich, showCategory = 10) + 
#   ggtitle("KEGG Pathway Enrichment") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#------------------------------CIBERSORT----------------------------------------
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(tibble)
# install.packages("devtools")
# devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
# Source the CIBERSORT function
source("F:/PhD/Other than PhD/DR data/CIBERSORT.R")
source('CIBERSORT.R')
BiocManager::install("preprocessCore")
BiocManager::install("ComplexHeatmap")
library(preprocessCore)

#-------------------------------------------------------------------------------
library(tidyverse)
library(biomaRt)

normalized_counts_df <- as.data.frame(normalized_counts)

normalized_counts_df <- normalized_counts_df %>%
  rownames_to_column("ENSEMBL_ID")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                  filters = "ensembl_gene_id",
                  values = normalized_counts_df$ENSEMBL_ID,
                  mart = ensembl)

normalized_counts_df <- normalized_counts_df %>%
  left_join(gene_map, by = c("ENSEMBL_ID" = "ensembl_gene_id"))
write.csv(normalized_counts_df, "normalized_counts_df.csv", row.names = FALSE)
normalized_counts_df_symbol <- read.table('normalized_counts_dff.csv', sep = ",", header = TRUE)
normalized_counts_df_symbol <- normalized_counts_df_symbol %>%
  filter(!is.na(ENSEMBL_ID))

normalized_counts_df_symbol <- normalized_counts_df_symbol %>%
  column_to_rownames("ENSEMBL_ID")
write.table(normalized_counts_df_symbol, file = "normalized_expression_data_for_cibersort.txt", sep = "\t", quote = FALSE, row.names = TRUE)

write.csv(normalized_counts_df_symbol, "normalized_counts_df_symbol.csv")
normalized_counts_df_g <- normalized_counts_df %>%
   full_join(gene_map, by = c("ENSEMBL_ID" = "ensembl_gene_id"))

rownames(normalized_counts_df) <- normalized_counts_df$external_gene_name
normalized_counts_df <- normalized_counts_df %>%
  select(-ENSEMBL_ID, -external_gene_name)


#-------------------------------------------------------------------------------
# Run CIBERSORT
results_cibersort <- cibersort("F:/PhD/Other than PhD/DR data/LM22.txt", "F:/PhD/Other than PhD/DR data/normalized_expression_data_for_cibersort.txt", perm = 1000, QN = TRUE)

#Prepare data for plotting
cibersort_data <- as.data.frame(results_cibersort)
normalized_data <- t(normalized_counts_df_symbol)

nSamples_cib <- ncol(cibersort_data)


#------------------------------------------------------------------------------
write.csv(results_cibersort, "results_cibersort.csv")
results_cibersort_corrected<-read.table('results_cibersort.csv', sep = ",", header = TRUE, row.names = 1)

# Compute correlation matrix for immune cells
cor_matrix <- cor(results_cibersort_corrected, use = "pairwise.complete.obs")

# Convert to long format for ggplot
cor_data_long <- melt(cor_matrix)

# Create the correlation heatmap
p_A <- ggplot(cor_data_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() +
  labs(title = "Correlation between Immune Cell Infiltration Landscapes",
       x = "Immune Cell Type", y = "Immune Cell Type", fill = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the figure
ggsave("immune_cell_correlation_heatmap.png", plot = p_A, width = 10, height = 8)
#---------------------------------------------------------------------------------
# Prepare data for plotting
heatmap_data_long <- melt(cor_matrix)
heatmap_data_long$p_value <- melt(p_values)$value

# Plot the heatmap with text annotations
p1 <- ggplot(heatmap_data_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), size = 3) +  # Add text annotations
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() +
   labs(title = "",
        x = "", y = "", fill = "Correlation") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
         axis.text.y = element_text(color = "black"))

# Save the plot
ggsave("immune_cell_correlation_heatmap_with_values.png", plot = p1)

#-------------------------------------------------------------------------------

# Assuming `metadata` contains a 'Group' column with 'DR' and 'Normal'
metadata<-read.table('metadata.csv', sep = ",", header = TRUE, row.names = 1)

cibersort_data_corrected<-results_cibersort_corrected
cibersort_data_corrected$Group <- metadata$description

# Convert to long format for ggplot
cibersort_data_long <- melt(cibersort_data_corrected, id.vars = "Group")

# Create the boxplot
p_B <- ggplot(cibersort_data_long, aes(x = variable, y = value, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "",
       x = "Immune Cell Type", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black"))
  scale_fill_manual(values = c("blue", "red"))

# Save the figure
ggsave("immune_cell_infiltration_boxplot.png", plot = p_B, width = 10, height = 8)
#-------------------------------------------------------------------------------
# Select the hub genes for correlation analysis
hub_genes <- c("DRAXIN", "TBCEL", "ZNF775", "TMSB4XP4", "SLIT1", "FSD1L", "NRADDP")
hub_gene_expression <- normalized_data[, hub_genes]

#cor_matrix_hub_vs_cells <- cor(hub_gene_expression, cibersort_data, use = "pairwise.complete.obs")

cor_matrix_hub_vs_cells <- cor(hub_gene_expression, results_cibersort_corrected, use = "pairwise.complete.obs")
# Check for these columns and remove them if they exist
#columns_to_remove <- c("P-value", "Correlation", "RMSE")

# Ensure the columns exist before attempting to remove them
#cor_matrix_hub_vs_cells_cleaned <- cor_matrix_hub_vs_cells[, !colnames(cor_matrix_hub_vs_cells) %in% columns_to_remove]





# Melt the correlation matrix for easy plotting
cor_matrix_hub_vs_cells_cleaned_melted <- melt(cor_matrix_hub_vs_cells)



# Plot the heatmap with correlation values displayed
p <- ggplot(cor_matrix_hub_vs_cells_cleaned_melted, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), size = 3) +  # Display correlation values
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() +
  labs(title = "",
       x = "Immune Cell Type", y = "Genes", fill = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"))

# Save the plot
ggsave("hub_gene_vs_immune_cells_correlation_heatmap.png", plot = p)
#--------------------------------------------------------------------------------



