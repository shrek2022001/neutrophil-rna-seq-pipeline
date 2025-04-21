setwd("/N/slate/spekande/RNA_seq")
# 1. Load library
library(DESeq2)

# 2. Read in count matrix
counts <- read.delim("counts/gene_counts_all.txt", comment.char="#", row.names=1)
counts <- counts[,6:ncol(counts)]  # remove first 5 metadata columns
colnames(counts) <- c("SRR507859", "SRR507861", "SRR507860", "SRR507862", "SRR507863", "SRR507864")

# 3. Define conditions for each sample
sample_info <- data.frame(
  row.names = colnames(counts),
  condition = factor(c("control", "control", "TNF", "TNF", "GMCSF", "GMCSF"))
)

# 4. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)

# ✅ 5. Now run DESeq() with the object!
dds <- DESeq(dds)

# TNF vs Control
res_TNF <- results(dds, contrast = c("condition", "TNF", "control"))
write.csv(as.data.frame(res_TNF), "DE_TNF_vs_Control.csv")

# GMCSF vs Control
res_GMCSF <- results(dds, contrast = c("condition", "GMCSF", "control"))
write.csv(as.data.frame(res_GMCSF), "DE_GMCSF_vs_Control.csv")

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

library(ggplot2)

# Filter to remove NA padj values
res_TNF_clean <- res_TNF[!is.na(res_TNF$padj), ]
ggplot(res_TNF, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.4) +
  theme_minimal() +
  ggtitle("Volcano Plot: TNF vs Control")

ggsave("volcano_TNF_vs_control.png", width = 8, height = 6, dpi = 300)

summary(res_TNF$padj)
sum(!is.na(res_TNF$padj))
sum(res_TNF$padj < 0.05, na.rm = TRUE)

res_GMCSF <- results(dds, contrast = c("condition", "GMCSF", "control"))
summary(res_GMCSF)
sum(res_GMCSF$padj < 0.05, na.rm = TRUE)

library(ggplot2)

# Clean data
res_GMCSF_clean <- na.omit(res_GMCSF)
res_GMCSF_clean$significant <- ifelse(res_GMCSF_clean$padj < 0.05 & abs(res_GMCSF_clean$log2FoldChange) > 1, "Yes", "No")

# Volcano Plot
ggplot(res_GMCSF_clean, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  ggtitle("Volcano Plot: GM-CSF vs Control")
table(res_GMCSF_clean$significant)

png("volcano_GMCSF_vs_control.png", width=1000, height=800)

ggplot(res_GMCSF_clean, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  ggtitle("Volcano Plot: GM-CSF vs Control")

dev.off()
ggplot(res_GMCSF_clean, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  theme_minimal() +
  ggtitle("Volcano Plot: GM-CSF vs Control (All Genes)")
summary(res_GMCSF_clean$log2FoldChange)
summary(res_GMCSF_clean$padj)
range(-log10(res_GMCSF_clean$padj), na.rm = TRUE)






# Replace padj = 0 with the smallest non-zero padj value
min_nonzero_padj <- min(res_GMCSF_clean$padj[res_GMCSF_clean$padj > 0])
res_GMCSF_clean$padj[res_GMCSF_clean$padj == 0] <- min_nonzero_padj

# Now calculate log10 padj safely
res_GMCSF_clean$log10padj <- -log10(res_GMCSF_clean$padj)

# Plot volcano
ggplot(res_GMCSF_clean, aes(x = log2FoldChange, y = log10padj, color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  ggtitle("Volcano Plot (Fixed): GM-CSF vs Control")

head(res_GMCSF_clean[, c("log2FoldChange", "log10padj", "significant")], 10)


png("volcano_GMCSF_fixed.png", width = 1000, height = 800)

ggplot(res_GMCSF_clean, aes(x = log2FoldChange, y = log10padj, color = significant)) +
  geom_point(size = 1) +  # no alpha!
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  ggtitle("Volcano Plot: GM-CSF vs Control")

dev.off()

top_genes <- res_GMCSF_clean[res_GMCSF_clean$significant == "Yes", ]
top_genes <- top_genes[order(top_genes$padj), ]
write.csv(head(top_genes, 50), "Top_50_DEGs_GMCSF_vs_Control.csv")

vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

library(pheatmap)

top_gene_ids <- rownames(head(top_genes, 30))  # top 30
pheatmap(assay(vsd)[top_gene_ids, ], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=as.data.frame(colData(dds)))

sig_genes <- res_GMCSF_clean[res_GMCSF_clean$padj < 0.05 & abs(res_GMCSF_clean$log2FoldChange) > 1, ]
gene_list <- rownames(sig_genes)
length(gene_list)  # how many significant genes

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")


BiocManager::install("org.Hs.eg.db")


library(org.Hs.eg.db)
library(clusterProfiler)

# Remove version suffix after the dot
gene_list_clean <- gsub("\\..*", "", gene_list)

# Map Ensembl to Entrez ID
gene_df <- bitr(gene_list,
                fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

head(gene_df)
