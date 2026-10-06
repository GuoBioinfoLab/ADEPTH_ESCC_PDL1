
install.packages("pheatmap")

# library -----------------------------------------------------------------

library(ggplot2)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# load data ---------------------------------------------------------------

gene_expr <- readr::read_tsv("output_data/all_gene_expr_length.txt", show_col_types = F)
sample_info <- readr::read_tsv("input_data/sample_info_rnaseq_final.txt", show_col_types = F)

# prepare data for heatmap -----------------------------------------------------------------

#before top20
mast = c("HLA-A","HLA-B","TMSB4X","HLA-C","RPL4","RPL3","EIF4A1","APOO","HLA-E", "EEF1A1","TUBA1B","EEF2","HES4","DDX5","RPL10A","RPSA","CPA3","RPS6","AHNAK","RPL13")
macro = c("SOD3","IGFBP7","CAV1","TAGLN","LUM","MYL9","IGFBP5","COL6A2","MFAP4", "RARRES2","FBLN1","COL3A1","MMP2","CFH","SPARCL1","TM4SF1","CAVIN1","C1R","COL1A2","C1S")

mast_macro=  c("HLA-A","HLA-B","TMSB4X","HLA-C","RPL4","RPL3","EIF4A1","APOO","HLA-E", "EEF1A1","TUBA1B","EEF2","HES4","DDX5","RPL10A","RPSA","CPA3","RPS6","AHNAK","RPL13","SOD3","IGFBP7","CAV1","TAGLN","LUM","MYL9","IGFBP5","COL6A2","MFAP4", "RARRES2","FBLN1","COL3A1","MMP2","CFH","SPARCL1","TM4SF1","CAVIN1","C1R","COL1A2","C1S")

# top 30

load("D:/01-project/02-project-postdoc/05-ESCC单细胞-周建丰课题/analysis/04-rnaseq/04-gsva/marker_f2.Rdata")
marker_f2$`Macro-SPARC`$Macro %>% 
  dplyr::arrange(desc(p_val)) %>% 
  dplyr::arrange(desc(avg_log2FC)) %>% 
  head(10) %>% 
  .$gene -> macro_top10

marker_f2$`Mast-HLA-A`$Mast %>% 
  dplyr::arrange(desc(p_val)) %>% 
  dplyr::arrange(desc(avg_log2FC)) %>% 
  head(10) %>% 
  .$gene -> mast_top30



count_df <- gene_expr %>%
  # dplyr::filter(SYMBOL %in% mast_macro) %>% 
  select(-Geneid, -Length) %>%                 # 去掉 Length
  column_to_rownames("SYMBOL") %>%     # 用 SYMBOL 作为行名
  as.data.frame()

sample_ordered <- sample_info %>%
  filter(sample_name %in% colnames(count_df)) %>%
  arrange(match(sample_name, colnames(count_df)))

rownames(sample_ordered) <- sample_info$sample_name

## 确保顺序一致
stopifnot(all(colnames(count_df) == rownames(sample_ordered)))

dds <- DESeqDataSetFromMatrix(
  countData = count_df,
  colData   = sample_ordered,
  design    = ~ response
)

## 低表达过滤（画热图很重要）
dds <- dds[rowSums(counts(dds)) >= 10, ]

## VST 变换
vsd <- vst(dds, blind = TRUE)


vsd_mat <- assay(vsd)

vsd_mat_sorted <- vsd_mat[, order(colnames(vsd_mat), decreasing = F)]
  

vsd_z <- t(scale(t(vsd_mat_sorted)))


# 准备注释数据 ------------------------------------------------------------------

annotation_col <- data.frame(
  sample_name = sample_ordered$sample_name,
  time_point = sample_ordered$time_point,
  response = sample_ordered$response

)

rownames(annotation_col) <- annotation_col$sample_name

annotation_col$response <- factor(
  annotation_col$response,
  levels = c("NR", "R")
)


sample_order <- rownames(annotation_col)[
  order(annotation_col$response)
  
]

annotation_colors <- list(
  response = c(
    NR  = "#4DBBD5",
    R = "#E64B35"
  ),
  time_point = c(
    preNICB  = "#00C6FF",
    postNICB = "#FF9289"
  )
)


hm_colors <- colorRampPalette(
  c("#3B4CC0", "white", "#B40426")
)(100)

macro_top10_sparc <- c("IGFBP7", "SOD3",   "IGFBP5", "MFAP4",  "CAV1",   "DCN",    "COL6A2", "FBLN1",  "TAGLN",  "GNG11", "SPARC","SPARCL1") 

genes_use <- intersect(macro_top10_sparc, rownames(vsd_z)) # 切换基因类型
vsd_z_macro <- vsd_z[genes_use, , drop = FALSE]

vsd_z_macro_sorted <- vsd_z_macro[, sample_order]

annotation_col_sorted <- annotation_col[sample_order, , drop = FALSE]

readr::write_csv(as.data.frame(vsd_z_macro),"~/my-data/03-博后课题/01-project/01-project/xie-rstudio/output_data/heatmap_for_core_genes.csv")

pheatmap(
  vsd_z_macro,
  color               = hm_colors,
  annotation_col      = annotation_col,
  annotation_colors   = annotation_colors,
  
  show_rownames       = TRUE,   # 基因多时隐藏
  show_colnames       = TRUE,
  
  cluster_rows        = T,
  cluster_cols        = F,
  
  gaps_col            = 30,
  
  scale               = "none",
  
  fontsize_col        = 9,
  fontsize_row        = 8,
  
  treeheight_row      = 30,
  treeheight_col      = 30,
  
  border_color        = NA,
  legend              = TRUE,
)
  


# save image --------------------------------------------------------------

save.image("output_data/plot_for_degs_heatmap.RData")
load("~/my-data/03-博后课题/01-project/01-project/xie-rstudio/output_data/plot_for_degs_heatmap.RData")
