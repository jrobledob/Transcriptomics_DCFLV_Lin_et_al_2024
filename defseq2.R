library("DESeq2")
library(ggplot2)
library("apeglm")
library(viridis)
library(reshape2)
library(ggdendro)
library(gridExtra)
library(gtable)
library(grid)
library(dplyr)
library(stringr)
library(AnnotationHub)
#library(clusterProfiler)
require(DOSE)
#library(enrichplot)
library(ggridges)
library(pathview)
library(ggupset)
library(clusterProfiler)
library(ggrepel)
library(enrichplot)
cts<- read.csv("C:/Users/jrobl/OneDrive - University of Florida/Documents/MASTER Thesis/Transcriptomics_Tony/Adults/gene_count_matrix_D_CITRI.csv",row.names = 1)
coldata<- read.csv("C:/Users/jrobl/OneDrive - University of Florida/Documents/MASTER Thesis/Transcriptomics_Tony/Adults/Names_of_samples.csv")
coldata$Status<- factor(coldata$Status, levels = c("Healthy", "FLAVI"))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~Status)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Status_FLAVI_vs_Healthy")
res <- res[order(res$padj),]
head(res)
#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))
plotCounts(dds, gene="MSTRG.25837|LOC103516011", intgroup="Status")
plotCounts(dds, gene="MSTRG.1026", intgroup="Status")
plotCounts(dds, gene="MSTRG.18639|LOC103512116", intgroup="Status")
plotCounts(dds, gene="MSTRG.41124", intgroup="Status")
plotCounts(dds, gene="MSTRG.13023|LOC103508934", intgroup="Status")
plotCounts(dds, gene="MSTRG.32196", intgroup="Status")
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
threshold_p_value<- 0.05
fold_change_threshold<- 2
with(subset(res, padj<=threshold_p_value ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<=threshold_p_value & abs(log2FoldChange)>=fold_change_threshold), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
#count DEGs according to parameters
nrow((subset(res, padj<= threshold_p_value & abs(log2FoldChange)>=fold_change_threshold)))
#count DEGs according to parameters (UPREGULATED)
nrow((subset(res, padj<= threshold_p_value & log2FoldChange>=fold_change_threshold)))
#count DEGs according to parameters (DOWNREGULATED)
nrow((subset(res, padj<= threshold_p_value & log2FoldChange<=(-fold_change_threshold))))
#PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Status")
# Transform count data using the variance stablilizing transform
deseq2VST <- vst(dds)
# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)
# Coerce to a data frame
deseq2ResDF <- as.data.frame(res)
# Examine this data frame
head(deseq2ResDF)
# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < threshold_p_value, "Significant", NA)
# Keep only the significantly differentiated genes with threshold_p_value and fold_change_threshold as parameters
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= fold_change_threshold & abs(deseq2ResDF$log2FoldChange) >= fold_change_threshold,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]
# Convert the VST counts to long format for ggplot2
# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))
head(deseq2VST_wide)
head(deseq2VST_long)
# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap
# Convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL
# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))
# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()
# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])
# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap
# Combine the dendrogram and the heatmap
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))
# Convert both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)
# Check the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths
# Add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)
# Make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
# Arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))
# Draw the plot
grid.draw(finalGrob)
coldata$ids <- factor(coldata$ids, levels=clusterSample$labels[clusterSample$order])
colours <- c("orange", "blue")
sampleClinical <- ggplot(coldata, aes(x=ids, y=1, fill=Status)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Status", values=colours) + theme_void()
sampleClinicalGrob <- ggplotGrob(sampleClinical)
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)
# Arrange and output the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleClinicalGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
grid.draw(finalGrob)
################################################################################
################# Step 1: create dendrogram for genes ##########################
# we already did the clustering for genes in the tutorial, get the data to make a dendrogram with ggplot
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))
# construct the dendrogram in ggplot
geneDendrogram <- ggplot(geneDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) + theme_dendro()
################################################################################
################# Step 2: Re-arrange the heatmap cells #########################
# re-factor genes for ggplot2
deseq2VST$Gene <- factor(deseq2VST$Gene, levels=clusterGene$labels[clusterGene$order])
# recreate the heatmap with this new factoring
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
################################################################################
################# Step 3: convert to everything to grobs #######################
# note! before this step as mentioned you might need to alter the expand parameters in the plot scales for all the plots we do that here
# convert the heatmap to a grob
heatmapGrob <- ggplotGrob(heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)))
# convert the dendrogram to a grob
# note! we flipped the axis above so the x-axis is now displayed as what we would think of as the y-axis
geneDendrogramGrob <- ggplotGrob(geneDendrogram + scale_x_discrete(expand=c(0, 0)))
# we already have a sample Dendrogram, but here it is again
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0)))
# we already have our sample clinical plot but here it is again
sampleClinicalGrob <- sampleClinicalGrob
################################################################################
######### Step 4: align the gene dendrograms to match the heatmap ##############
# check that the both the heatmap and gene dendrogram have the same number of vertical elements
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)
# make sure every height between the two grobs is the same
maxHeight <- unit.pmax(geneDendrogramGrob$heights, heatmapGrob$heights)
geneDendrogramGrob$heights <- as.list(maxHeight)
heatmapGrob$heights <- as.list(maxHeight)
################################################################################
# Step 4b: we have a new heatmap so we need to re-align the horizontal elements #
# repeat the steps in the tutorial
# check the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths
sampleClinicalGrob$widths
# add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)
# make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)
################################################################################
############### Step 5: create a blank panel ###################################
# we can use grid graphics for this
blankPanel <- grid.rect(gp=gpar(col="white"))
################################################################################
############### Step 6: Arrange the final result ###############################
# arrange all the plots together
finalGrob_v2 <- arrangeGrob(blankPanel, sampleDendrogramGrob, blankPanel, sampleClinicalGrob, geneDendrogramGrob, heatmapGrob, ncol=2, nrow=3, widths=c(1,5), heights=c(2,.8,6))
# draw the final result
grid.draw(finalGrob_v2)


#Setting Data base
deGenes_DB <- subset(deseq2ResDF, padj<= threshold_p_value & abs(log2FoldChange)>=fold_change_threshold)
deGenes_DB$Gene_ID_BLAST_search<- str_extract(rownames(deGenes_DB), "(?<=\\|)LOC\\d+")
deGenes_DB$Gene_ID_BLAST_search<-gsub("LOC(\\d+)", "\\1", deGenes_DB$Gene_ID_BLAST_search)
deGenes_DB$Trasncrip_ID<- deGenes_DB$Gene_ID_BLAST_search
deGenes_DB$Trasncrip_ID[which(is.na(deGenes_DB$Trasncrip_ID))]<- rownames(deGenes_DB)[which(is.na(deGenes_DB$Gene_ID_BLAST_search))]
list_of_annotated_genes<- deGenes_DB$Gene_ID_BLAST_search[which(!is.na(deGenes_DB$Gene_ID_BLAST_search))]
#write.table(deGenes_DB$Gene_ID_BLAST_search[which(!is.na(deGenes_DB$Gene_ID_BLAST_search))], "C:/Users/jrobl/OneDrive - University of Florida/Documents/MASTER Thesis/Transcriptomics_Tony/Adults/annotatated_DEGs_genes.txt", row.names = F, col.names = F)
query_NCBI<- read.delim("C:/Users/jrobl/OneDrive - University of Florida/Documents/MASTER Thesis/Transcriptomics_Tony/Adults/query_NCBI.txt")
colnames(query_NCBI)[2]<- "Gene_ID_BLAST_search"
query_NCBI$Gene_ID_BLAST_search<- as.character(query_NCBI$Gene_ID_BLAST_search)
deGenes_DB_and_query_NCBI<- left_join(deGenes_DB, query_NCBI)
write.csv(deGenes_DB_and_query_NCBI, "C:/Users/jrobl/OneDrive - University of Florida/Documents/MASTER Thesis/Transcriptomics_Tony/Adults/DEGs_padj_less0_05_and_fold_less_2_with_NCBI_information.csv")












#enrichment analysis
gene_universe_df<- data.frame(deseq2ResDF)
gene_universe_df$Gene_ID_BLAST_search<- str_extract(rownames(gene_universe_df), "(?<=\\|)LOC\\d+")
gene_universe_df$Gene_ID_BLAST_search<-gsub("LOC(\\d+)", "\\1", gene_universe_df$Gene_ID_BLAST_search)
hub <- AnnotationHub()
q <- query(hub, "Diaphorina")
id <- q$ah_id[length(q)]
D_citri <- hub[[id]]
# we want the log2 fold change 
original_gene_list <- gene_universe_df$log2FoldChange[which(!is.na(gene_universe_df$Gene_ID_BLAST_search))]
# name the vector
names(original_gene_list) <- gene_universe_df$Gene_ID_BLAST_search[which(!is.na(gene_universe_df$Gene_ID_BLAST_search))]
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = D_citri,
             pAdjustMethod="none")

dotplot(gse,  showCategory=5000, split=".sign") + facet_grid(.~.sign)
dotplot(gse,  showCategory=10, split=".sign") + facet_grid(.~.sign)+theme(axis.text.y = element_text(size = 5))
ridgeplot(gse) + labs(x = "enrichment distribution")+theme(axis.text.y = element_text(size = 5))


#KEGG
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "GID", toType = "ENTREZID", OrgDb=D_citri)
# remove duplicate IDS (here I use "ENTREZID", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENTREZID")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
gene_universe_df_NA_rm<- gene_universe_df[which(!is.na(gene_universe_df$Gene_ID_BLAST_search)),]
gene_universe_df_NA_rm_dup_rm<- gene_universe_df_NA_rm[!duplicated(gene_universe_df_NA_rm[c("Gene_ID_BLAST_search")]),] 
df2 = gene_universe_df_NA_rm_dup_rm[gene_universe_df_NA_rm_dup_rm$Gene_ID_BLAST_search %in% dedup_ids$ENTREZID,]
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID
# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange
# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y
# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "dci"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
dotplot(kk2, showCategory = 100, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)+theme(axis.text.y = element_text(size = 5))
cnetplot(gse, categorySize="pvalue", foldChange=gene_list)
ridgeplot(kk2) + labs(x = "enrichment distribution")
# Produce the native KEGG plot (PNG)
setwd("C:/Users/jrobl/OneDrive - University of Florida/Documents/MASTER Thesis/Transcriptomics_Tony/Adults/KEGG_Pathways_2/")
for (i in kk2@result$ID[!(kk2@result$ID %in% c("dci00513","dci00601", "dci00511", "dci00603", "dci00531"))]) {
  dme <- pathview(gene.data=kegg_gene_list, pathway.id=i, species = kegg_organism)
}



setwd("C:/Users/jrobl/OneDrive - University of Florida/Documents/MASTER Thesis/Transcriptomics_Tony/Adults/")
#Over-Representation Analysis with ClusterProfiler
# we want the log2 fold change 
original_gene_list <- gene_universe_df$log2FoldChange[which(!is.na(gene_universe_df$Gene_ID_BLAST_search))]
# name the vector
names(original_gene_list) <- gene_universe_df$Gene_ID_BLAST_search[which(!is.na(gene_universe_df$Gene_ID_BLAST_search))]
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
# Exctract significant results (threshold_p_value)
sig_genes_df = subset(gene_universe_df[which(!is.na(gene_universe_df$Gene_ID_BLAST_search)),], padj < threshold_p_value)
# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange
# Name the vector
names(genes) <- sig_genes_df$Gene_ID_BLAST_search
# omit NA values
genes <- na.omit(genes)
# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > fold_change_threshold]
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = D_citri, 
                      keyType = "GID",
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.25, 
                      qvalueCutoff = 0.50)




#Emphasizes the genes overlapping among different gene sets.
upsetplot(go_enrich)
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)
cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)
write.csv(deGenes_DB_and_query_NCBI, "./DEGs_padj_less0_05_and_fold_less_2_with_NCBI_information.csv", row.names = F)
write.csv(go_enrich, "./_GO_enrich_DEGs_padj_less0_05_and_fold_less_2.csv")
write.csv(gse, "./GO_set_enrich_DEGs_padj_less0_05_and_fold_less_2.csv")
write.csv(kk2, "./keeg_set_enrich_DEGs_padj_less0_05_and_fold_less_2.csv")

a<- data.frame(unlist(gse@geneSets))
a$ID<- rownames(a)
gse_plus_genes<- gse@result
gse_plus_genes$all<- NA
for (i in gse_plus_genes$ID) {
  index<- which(grepl(i, a$ID))
  genes<- a[which(grepl(i, a$ID)),]$unlist.gse.geneSets.
  gse_plus_genes$all[which(gse_plus_genes$ID==i)]<- paste0(genes, sep = "/", collapse = " ")
}
write.csv(gse_plus_genes, "./GO_set_enrich_DEGs_padj_less0_05_and_fold_less_2.csv")









cts_2<- read.csv("C:/Users/jrobl/OneDrive - University of Florida/Documents/MASTER Thesis/Transcriptomics_Tony/gene_count_matrix_FLAVI.csv",row.names = 1)
coldata<- read.csv("C:/Users/jrobl/OneDrive - University of Florida/Documents/MASTER Thesis/Transcriptomics_Tony/Names_of_samples.csv")
coldata$Status<- as.factor(coldata$Status)
dds <- DESeqDataSetFromMatrix(countData = cts_2,
                              colData = coldata,
                              design= ~Status)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Status_Healthy_vs_FLAVI")
res_shrink <- lfcShrink(dds, coef="Status_Healthy_vs_FLAVI", type="apeglm")
head(results(dds, tidy=TRUE))
summary(res)
res <- res[order(res$padj),]
head(res)




#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))

plotCounts(dds, gene="MSTRG.1", intgroup="Status")
plotCounts(dds, gene="MSTRG.2", intgroup="Status")
plotCounts(dds, gene="MSTRG.3", intgroup="Status")

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))



vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Status")

