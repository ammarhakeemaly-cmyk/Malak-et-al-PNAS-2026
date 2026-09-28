library(fgsea)
library(patchwork)

nmj <- read.csv("R:/Downloads/FlyBase_IDs.txt")
nmj_KO <- KO[intersect(nmj[[1]], rownames(KO)), , drop = FALSE][order(KO[intersect(nmj[[1]], rownames(KO)), ]$padj), ]
gene_list <- sort(setNames(sign(nmj_KO$log2FoldChange) * -log10(nmj_KO$pvalue), rownames(nmj_KO)), decreasing = TRUE)

fg_plot <- function(genelist, filename){
fg_KO <- fgsea(
  pathways = list(NMJ = rownames(nmj_KO)),
  stats = genelist,
  maxSize = 200
)
bar.df <- data.frame(
  x = seq_along(genelist),
  y = 1,
  stat = genelist
)
label <- paste0(
  "KQ over WT",
  "\nNES = ", sprintf("%.2f", fg_KO$NES),
  "\nP = ", formatC(fg_KO$pval, format = "e", digits = 2)
)
p1 <- plotEnrichment(rownames(nmj_KO), genelist)+ 
  annotate("text", x = Inf, y = Inf, label = label, hjust = 1.1, vjust = 1.1, size = 4)

lims <- quantile(abs(genelist), 0.90)
p2 <- ggplot(bar.df, aes(x, y, fill = stat)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-lims, lims),
    oob = scales::squish,
    guide = "none"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(t = -5, r = 5.5, b = 5.5, l = 5.5)
  )

## Combine
p1 / p2 +
  plot_layout(heights = c(20, 0.6))
combined <- p1 / p2 +
  plot_layout(heights = c(20, 0.6))

ggsave(plot=combined,filename = paste0("R:/", filename, ".pdf"))
}
fg_plot(gene_list_KO, "GSEA_KO")
fg_plot(gene_list_QQ, "GSEA_QQ")
fg_plot(gene_list_KQ, "GSEA_KQ")


fg_KO <- fgsea(
  pathways = list(NMJ = rownames(nmj_KO)),
  stats = gene_list_KO,
  maxSize = 200
)
fg_QQ <- fgsea(
  pathways = list(NMJ = rownames(nmj_KO)),
  stats = gene_list_QQ,
  maxSize = 200
)
fg_KQ <- fgsea(
  pathways = list(NMJ = rownames(nmj_KO)),
  stats = gene_list_KQ,
  maxSize = 200
)
LE_KO <- fg_KO$leadingEdge[[1]]
LE_QQ <- fg_QQ$leadingEdge[[1]]
LE_KQ <- fg_KQ$leadingEdge[[1]]

LE_combined <- unique(c(LE_KO, LE_KQ))

symbols <- mapping[rownames(assay(vsd)[LE_KO,]), "gene_symbol"]
heatmap <- pheatmap(assay(vsd)[LE_KO,],
                    scale = "row",
                    color = colorRampPalette(c("#3B4992FF", "#FFFFFFFF", "#EE0000FF")) (50),
                    show_rownames = TRUE,
                    cluster_cols = FALSE,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    cellwidth = 10,
                    cellheight = 10,
                    labels_row = symbols,
                    clustering_distance_rows = "correlation",
                    #cutree_rows = 4,
                    #cutree_cols = 3,
                    angle_col = 0,
                    filename = "R:/LE_KO.pdf")

symbols <- mapping[rownames(assay(vsd)[LE_QQ,]), "gene_symbol"]
heatmap <- pheatmap(assay(vsd)[LE_QQ,],
                    scale = "row",
                    color = colorRampPalette(c("#3B4992FF", "#FFFFFFFF", "#EE0000FF")) (50),
                    show_rownames = TRUE,
                    cluster_cols = FALSE,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    cellwidth = 10,
                    cellheight = 10,
                    labels_row = symbols,
                    clustering_distance_rows = "correlation",
                    #cutree_rows = 4,
                    #cutree_cols = 3,
                    angle_col = 0,
                    filename = "R:/LE_QQ.pdf")

symbols <- mapping[rownames(assay(vsd)[LE_KQ,]), "gene_symbol"]
heatmap <- pheatmap(assay(vsd)[LE_KQ,],
                    scale = "row",
                    color = colorRampPalette(c("#3B4992FF", "#FFFFFFFF", "#EE0000FF")) (50),
                    show_rownames = TRUE,
                    cluster_cols = FALSE,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    cellwidth = 10,
                    cellheight = 10,
                    labels_row = symbols,
                    clustering_distance_rows = "correlation",
                    #cutree_rows = 4,
                    #cutree_cols = 3,
                    angle_col = 0,
                    filename = "R:/LE_KQ.pdf")


symbols <- mapping[rownames(assay(vsd)[LE_combined,]), "gene_symbol"]
heatmap <- pheatmap(assay(vsd)[LE_combined,],
                    scale = "row",
                    color = colorRampPalette(c("#3B4992FF", "#FFFFFFFF", "#EE0000FF")) (50),
                    show_rownames = TRUE,
                    cluster_cols = FALSE,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    cellwidth = 10,
                    cellheight = 10,
                    labels_row = symbols,
                    clustering_distance_rows = "correlation",
                    #cutree_rows = 4,
                    #cutree_cols = 3,
                    angle_col = 0,
                    filename = "R:/LE_combined.pdf")
