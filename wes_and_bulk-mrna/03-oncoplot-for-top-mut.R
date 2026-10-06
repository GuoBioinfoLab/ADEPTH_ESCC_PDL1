# library -----------------------------------------------------------------

library(maftools)
library(grid)
library(ggpubr)
library(ComplexHeatmap)

# load all maf file -------------------------------------------------------

.all_sample_info <- readr::read_tsv("~/wes_data/wes_sample_info.txt", show_col_types = FALSE) %>%
  dplyr::rename(Tumor_Sample_Barcode = "sample_id")

.all_48_maf_file <- read.maf("~/wes_pipline/04-modified_add_tumor_barcode_maf/all_48_filtered_file.maf", clinicalData = .all_sample_info)

# 保存为附件所需的表格

.all_48_maf_file@variant.classification.summary

readr::write_csv(.all_48_maf_file@variant.classification.summary, "~/Rproject/sczhoujianfeng/data/wes_output/all_48_sample_mut_stat.csv")


# plot waterplot ----------------------------------------------------------------


oncoplot(
  maf = .all_48_maf_file, top = 30,
  clinicalFeatures = c("response", "time_point"),
  annotationColor = annotation_colors,
  sortByAnnotation = T,
  anno_height = 1.1,
  colors = col,
  bgCol = "#EFF6FC",
  draw_titv = F,
  gene_mar = 6,
  writeMatrix = F
)

# width = 10, height=6.5

# .all_187_maf_file@gene.summary -> .aa


# 绘制美化的基因突变瀑布图 ------------------------------------------------------------

col <- c(
  "Missense_Mutation" = "#D95F02", # 主蓝（Cell 核心色）
  "Nonsense_Mutation" = "#4D7DB7", # 唯一暖色（stop，强调）
  "Splice_Site" = "#7E6FB0", # 冷紫（机制）
  "Frame_Shift_Del" = "#1F8A70", # 深青绿
  "Frame_Shift_Ins" = "#6CC3B2", # 浅青绿
  "In_Frame_Del" = "#E6AB02", # 蓝（结构改变）
  "In_Frame_Ins" = "#9CC2E6", # 浅蓝
  "Multi_Hit" = "#8E8E8E" # 中性灰
)

alter_fun <- list(
  background = alter_graphic("rect", fill = "#EFF6FC"), ## CCCCCC
  Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
  Nonsense_Mutation = alter_graphic("rect", fill = col["Nonsense_Mutation"]),
  Splice_Site = alter_graphic("rect", fill = col["Splice_Site"]),
  Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),
  Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
  In_Frame_Del = alter_graphic("rect", fill = col["In_Frame_Del"]),
  In_Frame_Ins = alter_graphic("rect", fill = col["In_Frame_Ins"]),
  Multi_Hit = alter_graphic("rect", fill = col["Multi_Hit"]),
  `0` = alter_graphic("rect", fill = col["0"])
)

# 定义图例
heatmap_legend_param <- list(
  title = "Alternations",
  at = c(
    "Missense_Mutation", "Nonsense_Mutation",
    "Splice_Site", "Frame_Shift_Del",
    "Frame_Shift_Ins", "In_Frame_Del",
    "In_Frame_Ins",
    "Multi_Hit"
  ),
  labels = c(
    "Missense_Mutation", "Nonsense_Mutation",
    "Splice_Site", "Frame_Shift_Del",
    "Frame_Shift_Ins", "In_Frame_Del",
    "In_Frame_Ins",
    "Multi_Hit"
  )
)

# 添加临床信息

annotation_colors <- list(
  response = c(
    "NR" = "#4C6A92", # 深冷蓝（non-responder）
    "R"  = "#6FA7C7" # 亮一些的蓝（responder）
  ),
  time_point = c(
    "preNICB"  = "#E1E4E8", # 很浅的冷灰（baseline）
    "postNICB" = "#8A9FB3" # 蓝灰（treatment）
  )
)


columnanno <- HeatmapAnnotation(
  response = .all_sample_info$response,
  time_point = .all_sample_info$time_point,
  show_annotation_name = T,
  col = annotation_colors
  # annotation_height = c(
  #   Group = unit(2, "mm"),
  #   Gender = unit(4, "mm"),
  #   Clinical_stage = unit(4, "mm")
  # ),
  # height = unit(4, "mm")  # 总高度强制设置
  # annotation_legend_param = list(
  #   Group = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7)),
  #   Gender = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7)),
  #   Alternations=list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5)),
  #   Clinical_stage = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
  # )
)


# HLA family gene ---------------------------------------------------------


hla_gene <- c("HLA-A", "HLA-B", "HLA-C", "ACVR2A", "B2M", "CIITA", "CTNND1", "CYLD", "ERAP1", "ERAP2", "FBXW7", "FGFR1", "IRF1", "JAK2", "KEAP1", "PSME1", "RNF43", "TAP1", "TAP2", "TAPBP")

oncoplot(
  maf = .all_48_maf_file, genes = hla_gene,
  clinicalFeatures = c("response", "time_point"),
  annotationColor = annotation_colors,
  sortByAnnotation = T,
  anno_height = 1.1,
  colors = col,
  bgCol = "#EFF6FC",
  draw_titv = F,
  gene_mar = 6,
  writeMatrix = F
)


# MMDR gene

dMMR_genes <- c(
  "POLD3", "MLH1", "MLH3", "MSH2", "MSH6",
  "PMS2", "POLD1", "POLD2", "RFC1", "RFC3",
  "RFC4", "RFC5", "RPA1", "RPA3", "EXO1", "POLE",
  "MSH3", "RPA4", "LIG1", "PCNA", "POLD4", "RFC2", "RPA2", "SSBP1", "PMS1"
)


oncoplot(
  maf = .all_48_maf_file, genes = dMMR_genes,
  clinicalFeatures = c("response", "time_point"),
  annotationColor = annotation_colors,
  sortByAnnotation = T,
  anno_height = 1.1,
  colors = col,
  bgCol = "#EFF6FC",
  draw_titv = F,
  gene_mar = 6,
  writeMatrix = F
)


# save image --------------------------------------------------------------

save.image("./RData_wes/oncoplot_for_overview_mutation.RData")
load("./RData_wes/oncoplot_for_overview_mutation.RData")
