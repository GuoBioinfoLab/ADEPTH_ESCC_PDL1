# library -----------------------------------------------------------------

library(maftools)
library(grid)
library(ggpubr)
library(magrittr)

# load all maf file -------------------------------------------------------

.all_sample_info <- readr::read_tsv("~/wes_data/wes_sample_info.txt", show_col_types = FALSE) %>%
  dplyr::rename(Tumor_Sample_Barcode = "sample_id")

.all_48_maf_file <- read.maf("~/wes_pipline/04-modified_add_tumor_barcode_maf/all_48_filtered_file.maf", clinicalData = .all_sample_info)

# 概览
pathways(maf = .all_48_maf_file, pathdb = "sigpw", plotType = "bar")
pathways(maf = .all_48_maf_file, pathdb = "smgbp", plotType = "bar")

# 治疗后 NR vs R
clin.postNR <- .all_sample_info %>%
  dplyr::filter(time_point == "postNICB") %>%
  dplyr::filter(response == "NR") %>%
  .$Tumor_Sample_Barcode
clin.postR <- .all_sample_info %>%
  dplyr::filter(time_point == "postNICB") %>%
  dplyr::filter(response == "R") %>%
  .$Tumor_Sample_Barcode

escc.postNR <- subsetMaf(maf = .all_48_maf_file, tsb = clin.postNR, isTCGA = FALSE)
escc.postR <- subsetMaf(maf = .all_48_maf_file, tsb = clin.postR, isTCGA = FALSE)

pathways(maf = escc.postNR, pathdb = "sigpw", plotType = "bar") -> .postNR
pathways(maf = escc.postR, pathdb = "sigpw", plotType = "bar") -> .postR

# 治疗前 NR vs R
clin.preNR <- .all_sample_info %>%
  dplyr::filter(time_point == "preNICB") %>%
  dplyr::filter(response == "NR") %>%
  .$Tumor_Sample_Barcode
clin.preR <- .all_sample_info %>%
  dplyr::filter(time_point == "preNICB") %>%
  dplyr::filter(response == "R") %>%
  .$Tumor_Sample_Barcode

escc.preNR <- subsetMaf(maf = .all_48_maf_file, tsb = clin.preNR, isTCGA = FALSE)
escc.preR <- subsetMaf(maf = .all_48_maf_file, tsb = clin.preR, isTCGA = FALSE)

pathways(maf = escc.preNR, pathdb = "sigpw", plotType = "bar") -> .preNR
pathways(maf = escc.preR, pathdb = "sigpw", plotType = "bar") -> .preR


# 合并绘图 （治疗前 NR vs R）--------------------------------------------------------------------

intersect(.preNR$Pathway, .preR$Pathway) -> common_pathway

.preNR %>%
  dplyr::filter(Pathway %in% common_pathway) %>%
  dplyr::mutate(type = "preNR") -> .preNR_for_plot

.preR %>%
  dplyr::filter(Pathway %in% common_pathway) %>%
  dplyr::mutate(type = "preR") -> .preR_for_plot

rbind(.preNR_for_plot, .preR_for_plot) %>%
  dplyr::select(Pathway, Fraction_mutated_samples, type) -> .for_plot

# 合并绘图 （治疗后 NR vs R）--------------------------------------------------------------------

intersect(.postNR$Pathway, .postR$Pathway) -> common_pathway_post

.postNR %>%
  dplyr::filter(Pathway %in% common_pathway_post) %>%
  dplyr::mutate(type = "postNR") -> .postNR_for_plot

.postR %>%
  dplyr::filter(Pathway %in% common_pathway_post) %>%
  dplyr::mutate(type = "postR") -> .postR_for_plot

rbind(.postNR_for_plot, .postR_for_plot) %>%
  dplyr::select(Pathway, Fraction_mutated_samples, type) -> .for_post_plot


readr::write_csv(.for_post_plot, "./data/wes_output/pathway_postNR_R.csv")

ggplot(.for_post_plot, aes(
  x = Pathway,
  y = Fraction_mutated_samples,
  fill = type
)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    size = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "postNR" = "#3C5488", # 蓝#4DBBD
      "postR"  = "#F39B7F" # 红#E64B35
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Fraction of mutated samples",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(size = 11),
    axis.line = element_line(size = 0.5, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    plot.title = element_text(
      hjust = 0,
      size = 12,
      face = "bold"
    )
  )


# 合并绘图 （治疗前 NR vs R）--------------------------------------------------------------------

intersect(.preNR$Pathway, .preR$Pathway) -> common_pathway_pre

.preNR %>%
  dplyr::filter(Pathway %in% common_pathway_pre) %>%
  dplyr::mutate(type = "preNR") -> .preNR_for_plot

.preR %>%
  dplyr::filter(Pathway %in% common_pathway_pre) %>%
  dplyr::mutate(type = "preR") -> .preR_for_plot

rbind(.preNR_for_plot, .preR_for_plot) %>%
  dplyr::select(Pathway, Fraction_mutated_samples, type) -> .for_pre_plot

readr::write_csv(.for_pre_plot, "./data/wes_output/pathway_preNR_R.csv")

ggplot(.for_pre_plot, aes(
  x = Pathway,
  y = Fraction_mutated_samples,
  fill = type
)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    size = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "preNR" = "#4A6FA5", # 蓝#4DBBD
      "preR"  = "#B24745" # 红#E64B35
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Fraction of mutated samples",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(size = 11),
    axis.line = element_line(size = 0.5, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    plot.title = element_text(
      hjust = 0,
      size = 12,
      face = "bold"
    )
  )

# width=4， height=3

# save image --------------------------------------------------------------

save.image("RData_wes/cancer_pathway.RData")
load("./RData_wes/cancer_pathway.RData")
