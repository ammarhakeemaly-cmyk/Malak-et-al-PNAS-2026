library(org.Dm.eg.db)


gene2go <- AnnotationDbi::select(
  org.Dm.eg.db,
  keys = keys(org.Dm.eg.db, keytype="FLYBASE"),
  columns = "GO",
  keytype = "FLYBASE"
)
gene2go <- unique(gene2go[, c("FLYBASE", "GO", "ONTOLOGY")])

slim <- read.delim("D:/RNAseq_Project/reference/goslim_drosophila.tsv")
slim[,1] <- sub(
  "GO_",
  "GO:",
  sub(".*GO_", "GO_", slim[[1]])
)
slim[,1] <- sub(">", "", slim[,1])
colnames(slim) <- c("TERM", "NAME")
gene2slim <- subset(
  gene2go[,c(2,1,3)],
  GO %in% slim[,1]
)

slim_BP <- gene2slim[gene2slim$ONTOLOGY=="BP", c(1,2)]
slim_CC <- gene2slim[gene2slim$ONTOLOGY=="CC", c(1,2)]
slim_MF <- gene2slim[gene2slim$ONTOLOGY=="MF", c(1,2)]

gene_list_KO <- sort(setNames(sign(KO$log2FoldChange) * -log10(KO$pvalue), rownames(KO)), decreasing = TRUE)
gene_list_QQ <- sort(setNames(sign(KO_QQ$log2FoldChange) * -log10(KO_QQ$pvalue), rownames(KO_QQ)), decreasing = TRUE)
gene_list_KQ <- sort(setNames(sign(KO_QQ_KO$log2FoldChange) * -log10(KO_QQ_KO$pvalue), rownames(KO_QQ_KO)), decreasing = TRUE)

bp_KO <- GSEA(geneList = gene_list_KO,
            universe=row.names(dds),
            TERM2GENE = slim_BP,
            TERM2NAME = slim,
            pAdjustMethod = "fdr")
bp_QQ <- GSEA(geneList = gene_list_QQ,
              universe=row.names(dds),
              TERM2GENE = slim_BP,
              TERM2NAME = slim,
              pAdjustMethod = "fdr")
bp_KQ <- GSEA(geneList = gene_list_KQ,
              universe=row.names(dds),
              TERM2GENE = slim_BP,
              TERM2NAME = slim,
              pAdjustMethod = "fdr")

CC_KO <- GSEA(geneList = gene_list_KO,
              universe=row.names(dds),
              TERM2GENE = slim_CC,
              TERM2NAME = slim,
              pAdjustMethod = "fdr")
CC_QQ <- GSEA(geneList = gene_list_QQ,
              universe=row.names(dds),
              TERM2GENE = slim_CC,
              TERM2NAME = slim,
              pAdjustMethod = "fdr")
CC_KQ <- GSEA(geneList = gene_list_KQ,
              universe=row.names(dds),
              TERM2GENE = slim_CC,
              TERM2NAME = slim,
              pAdjustMethod = "fdr")

MF_KO <- GSEA(geneList = gene_list_KO,
              universe=row.names(dds),
              TERM2GENE = slim_MF,
              TERM2NAME = slim,
              pAdjustMethod = "fdr")
MF_QQ <- GSEA(geneList = gene_list_QQ,
              universe=row.names(dds),
              TERM2GENE = slim_MF,
              TERM2NAME = slim,
              pAdjustMethod = "fdr")
MF_KQ <- GSEA(geneList = gene_list_KQ,
              universe=row.names(dds),
              TERM2GENE = slim_MF,
              TERM2NAME = slim,
              pAdjustMethod = "fdr")


write_gsea_symbols <- function(gsea_res, file,
                               OrgDb = org.Dm.eg.db) {
  
  res <- gsea_res@result
  
  # Split core enrichment strings into gene vectors
  core_lists <- strsplit(res$core_enrichment, "/")
  
  # Convert all unique IDs once
  gene_map <- bitr(
    unique(unlist(core_lists)),
    fromType = "FLYBASE",
    toType = "SYMBOL",
    OrgDb = OrgDb
  )
  
  # Build lookup table
  lookup <- setNames(
    gene_map$SYMBOL,
    gene_map$FLYBASE
  )
  
  # Replace FLYBASE IDs with SYMBOLs
  res$core_enrichment <- vapply(
    core_lists,
    function(x) {
      
      syms <- lookup[x]
      
      # Keep original ID if no symbol exists
      syms[is.na(syms)] <- x[is.na(syms)]
      
      paste(syms, collapse = "/")
      
    },
    character(1)
  )
  
  write.csv(
    res,
    file = file,
    row.names = FALSE
  )
  
  invisible(res)
}

write_gsea_symbols(bp_KO, file = "R:/KO_over_WT_BP.csv")
write_gsea_symbols(CC_KO, file = "R:/KO_over_WT_CC.csv")
write_gsea_symbols(MF_KO, file = "R:/KO_over_WT_MF.csv")

write_gsea_symbols(bp_QQ, file = "R:/QQ_over_WT_BP.csv")
write_gsea_symbols(CC_QQ, file = "R:/QQ_over_WT_CC.csv")
write_gsea_symbols(MF_QQ, file = "R:/QQ_over_WT_MF.csv")

write_gsea_symbols(bp_KQ, file = "R:/QQ_over_KO_BP.csv")
write_gsea_symbols(CC_KQ, file = "R:/QQ_over_KO_CC.csv")
write_gsea_symbols(MF_KQ, file = "R:/QQ_over_KO_MF.csv")


ggsave(plot = dotplot(bp_KO,showCategory=10, split=".sign") + facet_grid(.~.sign), filename = "R:/KO_BP_Dot.pdf")
ggsave(plot = dotplot(bp_QQ,showCategory=10, split=".sign") + facet_grid(.~.sign), filename = "R:/QQ_BP_Dot.pdf")
ggsave(plot = dotplot(bp_KQ,showCategory=10, split=".sign") + facet_grid(.~.sign), filename = "R:/KQ_BP_Dot.pdf")

ggsave(plot = dotplot(CC_KO,showCategory=10, split=".sign") + facet_grid(.~.sign), filename = "R:/KO_CC_Dot.pdf")
ggsave(plot = dotplot(CC_QQ,showCategory=10, split=".sign") + facet_grid(.~.sign), filename = "R:/QQ_CC_Dot.pdf")
ggsave(plot = dotplot(CC_KQ,showCategory=10, split=".sign") + facet_grid(.~.sign), filename = "R:/KQ_CC_Dot.pdf")

ggsave(plot = dotplot(MF_KO,showCategory=10, split=".sign") + facet_grid(.~.sign), filename = "R:/KO_MF_Dot.pdf")
ggsave(plot = dotplot(MF_QQ,showCategory=10, split=".sign") + facet_grid(.~.sign), filename = "R:/QQ_MF_Dot.pdf")
ggsave(plot = dotplot(MF_KQ,showCategory=10, split=".sign") + facet_grid(.~.sign), filename = "R:/KQ_MF_Dot.pdf")


write.csv(bp_KO, file = "R:/KO_over_WT_BP.csv")
write.csv(bp_QQ, file = "R:/QQ_over_WT_BP.csv")
write.csv(bp_KQ, file = "R:/KO_over_QQ_BP.csv")

write.csv(CC_KO, file = "R:/KO_over_WT_CC.csv")
write.csv(CC_QQ, file = "R:/QQ_over_WT_CC.csv")
write.csv(CC_KQ, file = "R:/KO_over_QQ_CC.csv")

write.csv(MF_KO, file = "R:/KO_over_WT_MF.csv")
write.csv(MF_QQ, file = "R:/QQ_over_WT_MF.csv")
write.csv(MF_KQ, file = "R:/KO_over_QQ_MF.csv")
