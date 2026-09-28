library(DESeq2)
library(uwot)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(pheatmap)
library(ggsci)
library(tibble)
#library(ReactomePA)

library(limma)
library(qvalue)

library(magrittr)
library(dplyr)

library(clusterProfiler)
library(org.Dm.eg.db)

source("C:/users/User/Documents/MEGA/ggplot_theme.R")


# Load data ----
projectDir <- "D:\\RNAseq_Project\\"
metadataDir <- paste0(projectDir, "metadata\\sample_metadata.csv")
countDir <- paste0(projectDir,"rsubread")

mapping <- read.delim(paste0(projectDir, "/reference/fbgn_fbtr_fbpp_expanded_fb_2026_01.tsv.gz")) %>%
  dplyr::select(gene_ID, gene_symbol) %>% distinct(gene_ID, .keep_all = TRUE) %>% column_to_rownames(var = "gene_ID")

colData <- read.csv(metadataDir, row.names = 1)
colData$condition <- factor(colData$condition)
colData$condition <- relevel(colData$condition,ref = "WT")

count_files <- list.files(countDir, pattern = "\\.csv$", full.names = TRUE)

geneID <- read.csv(count_files[1], header = FALSE, skip = 1)[,1]
cts <- matrix(nrow=length(geneID), ncol=length(count_files), dimnames=list(geneID, NULL))
colnames(cts) <- tools::file_path_sans_ext(basename(count_files))

for (i in seq(length(count_files))){
  cts[,i] <- read.csv(count_files[i], header = FALSE, skip = 1)[,2]
}

cts <- cts[, rownames(colData)]
stopifnot(all(rownames(colData) == colnames(cts)))

# Begin differential expression analysis ----
# modeling GlurIIA WT, GlurIIA KO, and GlurIIA KO + GlurIIC_QQ OE
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ condition)
# Prefilter to include onl genes with 10+ reads in 3+ samples
dds <- dds[rowSums(counts(dds) >= 10) >= 3,]

dds <- DESeq(dds)

vsd <- vst(dds, blind = TRUE)

# Batch Correction ----
pc1 <- prcomp(t(assay(vsd)))$x[,1]
colData(dds)$PC1 <- pc1
design(dds) <- ~ PC1 + condition
dds <- DESeq(dds)

assay(vsd) <- removeBatchEffect(assay(vsd), covariates = pc1)

# DE List ----
alpha_threshold <- 0.05

KO <- results(dds, contrast = c("condition", "KO", "WT"))
KO$symbol <- mapping[rownames(KO), "gene_symbol"]
KO$padj <- qvalue(KO$pvalue)$qvalues
KO <- KO[order(KO$padj),]
KO$sig <- if_else(KO$padj<=alpha_threshold & KO$log2FoldChange >= 1, "up",
                  if_else(KO$padj<=alpha_threshold & KO$log2FoldChange <= -1, "down",
                          "n.s."))

KO_QQ<- results(dds, contrast = c("condition", "KO_QQ", "WT"))
KO_QQ$symbol <- mapping[rownames(KO_QQ), "gene_symbol"]
KO_QQ$padj <- qvalue(KO_QQ$pvalue)$qvalues
KO_QQ <- KO_QQ[order(KO_QQ$padj),]
KO_QQ$sig <- if_else(KO_QQ$padj<=alpha_threshold & KO_QQ$log2FoldChange >= 1, "up",
                  if_else(KO_QQ$padj<=alpha_threshold & KO_QQ$log2FoldChange <= -1, "down",
                          "n.s."))

KO_QQ_KO <- results(dds, contrast = c("condition", "KO_QQ", "KO"))
KO_QQ_KO$symbol <- mapping[rownames(KO_QQ_KO), "gene_symbol"]
KO_QQ_KO$padj <- qvalue(KO_QQ_KO$pvalue)$qvalues
KO_QQ_KO <- KO_QQ_KO[order(KO_QQ_KO$padj),]
KO_QQ_KO$sig <- if_else(KO_QQ_KO$padj<=alpha_threshold & KO_QQ_KO$log2FoldChange >= 1, "up",
                     if_else(KO_QQ_KO$padj<=alpha_threshold & KO_QQ_KO$log2FoldChange <= -1, "down",
                             "n.s."))

write.csv(KO, file = paste0(projectDir, "results/analysis/dge_tables/KO Over WT.csv"))
write.csv(KO_QQ, file = paste0(projectDir, "results/analysis/dge_tables/KO_QQ Over WT.csv"))
write.csv(KO_QQ_KO, file = paste0(projectDir, "results/analysis/dge_tables/KO_QQ Over KO.csv"))

DEG <- unique(c(rownames(KO[KO$padj<0.05,]),
         rownames(KO_QQ[KO_QQ$padj<0.05,])
         #rownames(KO_QQ_KO[KO_QQ_KO$padj<0.05,])
         )
)
# UMAP ----
mat <- t(assay(vsd))

umap_res <- umap(
  mat,
  n_neighbors = 2,
  n_threads = 1,
  n_sgd_threads = 1,
  seed = 42,
  min_dist = 0.1
)
umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  Condition = colData(vsd)$condition,
  Sample = rownames(colData(vsd))
)

umap_plot <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Condition)) +
  geom_point(size = 4) + theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()) + scale_color_aaas()

ggsave(plot = umap_plot, "R:/UMAP.pdf")

# GLURII Counts ----
G2A <- plotCounts(dds, gene = "FBgn0004620", returnData = TRUE)
G2C <- plotCounts(dds, gene = "FBgn0046113", returnData = TRUE)

G2A_plot <- ggplot(G2A, aes(condition, count, fill = condition)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_beeswarm(cex = 3) + scale_color_aaas() + scale_fill_aaas()

G2C_plot <- ggplot(G2C, aes(condition, count, fill = condition)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_beeswarm(cex = 3) + scale_color_aaas() + scale_fill_aaas()

ggsave(plot = G2A_plot, "R:/G2A.pdf")
ggsave(plot = G2C_plot, "R:/G2C.pdf")

# Heat Map ----
ph <- pheatmap(assay(vsd)[DEG,],
         scale = "row",
         color = colorRampPalette(c("#3B4992FF", "#FFFFFFFF", "#EE0000FF")) (50),
         show_rownames = FALSE,
         cluster_cols = TRUE,
         treeheight_row = 0,
         treeheight_col = 5,
         clustering_distance_rows = "correlation",
         cutree_rows = 4,
         cutree_cols = 3,
         angle_col = 0,
         filename = "R:/heatmap.pdf")

heatmap_table <- as.data.frame(
  t(scale(t(assay(vsd)[DEG,])))[ph$tree_row$order, ph$tree_col$order, drop = FALSE]
)
heatmap_table$symbol <- mapping[rownames(heatmap_table), "gene_symbol"]
write.csv(heatmap_table, file = paste0(projectDir, "results/analysis/dge_heatmap_fig5c.csv"))


# PCA ----
pca <- data.frame(prcomp(mat[,DEG], scale. = TRUE)$x)
var_explained <- (prcomp(mat[,DEG], scale. = TRUE)$sdev^2) / sum(prcomp(mat[,DEG], scale. = TRUE)$sdev^2)
pca$condition <- colData$condition
pca_plot <- ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point() +
  xlab(paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"))+
  ylab(paste0("PC2 (", round(var_explained[2] * 100, 1), "%)"))+
  scale_color_aaas()

ggsave(plot = pca_plot, "R:/pca.pdf")


# Volcano Plot ----
KO_label <- rbind(
  KO[KO$sig=="up",] %>% head(n=10),
  KO[KO$sig=="down",] %>% head(n = 10)
)
KO_volcano <- ggplot(KO, aes(x=log2FoldChange, y=-log10(padj), color = sig)) + geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype=2) +
  geom_vline(xintercept = -1, linetype=2) + geom_vline(xintercept = 1, linetype=2)+
  scale_color_manual(values=c(
    "up" = "#EE0000FF",
    "down" = "#3B4992FF",
    "n.s." = "#DDDDDDDD"
  )) +
  geom_text_repel(data=KO_label, aes (label = symbol), max.overlaps = Inf)


KO_QQ_label <- rbind(
  KO_QQ[KO_QQ$sig=="up",] %>% head(n=10),
  KO_QQ[KO_QQ$sig=="down",] %>% head(n = 10)
)
KO_QQ_volcano <- ggplot(KO_QQ, aes(x=log2FoldChange, y=-log10(padj), color = sig)) + geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype=2) +
  geom_vline(xintercept = -1, linetype=2) + geom_vline(xintercept = 1, linetype=2)+
  scale_color_manual(values=c(
    "up" = "#EE0000FF",
    "down" = "#3B4992FF",
    "n.s." = "#DDDDDDDD"
  )) +
  geom_text_repel(data=KO_QQ_label, aes (label = symbol), max.overlaps = Inf)


KO_QQ_KO_label <- rbind(
  KO_QQ_KO[KO_QQ_KO$sig=="up",] %>% head(n=10),
  KO_QQ_KO[KO_QQ_KO$sig=="down",] %>% head(n = 10)
)
KO_QQ_KO_volcano <- ggplot(KO_QQ_KO, aes(x=log2FoldChange, y=-log10(padj), color = sig)) + geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype=2) +
  geom_vline(xintercept = -1, linetype=2) + geom_vline(xintercept = 1, linetype=2)+
  scale_color_manual(values=c(
    "up" = "#EE0000FF",
    "down" = "#3B4992FF",
    "n.s." = "#DDDDDDDD"
  )) +
  geom_text_repel(data=KO_QQ_KO_label, aes (label = symbol), max.overlaps = Inf)

ggsave(plot = KO_volcano, "R:/KO Volcano.pdf")
ggsave(plot = KO_QQ_volcano, "R:/KO_QQ Volcano.pdf")
ggsave(plot = KO_QQ_KO_volcano, "R:/KO_QQ_KO Volcano.pdf")