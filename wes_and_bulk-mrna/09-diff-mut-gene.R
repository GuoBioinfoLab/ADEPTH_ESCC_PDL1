# library -----------------------------------------------------------------

library(maftools)
library(grid)
library(ggpubr)
library(magrittr)
library(dplyr)
# install.packages("forestplot")
library(forestplot)

# load all maf file -------------------------------------------------------

.all_sample_info <- readr::read_tsv("~/wes_data/wes_sample_info.txt", show_col_types = FALSE) %>%
  dplyr::rename(Tumor_Sample_Barcode = "sample_id")

.all_48_maf_file <- read.maf("~/wes_pipline/04-modified_add_tumor_barcode_maf/all_48_filtered_file.maf", clinicalData = .all_sample_info)

# compare diff mut gene for R vs. NR ------------------------------------
# postNICB
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

postNR_R <- mafCompare(m1 = escc.postNR, m2 = escc.postR, m1Name = "postNR", m2Name = "postR", minMut = 3)

forestPlot(mafCompareRes = postNR_R, pVal = 0.59, color = c("maroon", "royalblue"), geneFontSize = 0.6)


# preNICB

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

preNR_R <- mafCompare(m1 = escc.preNR, m2 = escc.preR, m1Name = "preNR", m2Name = "preR", minMut = 3)

forestPlot(mafCompareRes = preNR_R, pVal = 0.59, color = c("maroon", "royalblue"), geneFontSize = 0.6)

readr::write_csv(preNR_R$results, "./data/wes_output/deg_oncoplot_preNR_R.csv")

readr::write_csv(postNR_R$results, "./data/wes_output/deg_oncoplot_postNR_R.csv")


# R

R_pre_post <- mafCompare(m1 = escc.preR, m2 = escc.postR, m1Name = "preR", m2Name = "postR", minMut = 3)

forestPlot(mafCompareRes = R_pre_post, pVal = 1.1, color = c("maroon", "royalblue"), geneFontSize = 0.6)


# NR

NR_pre_post <- mafCompare(m1 = escc.preNR, m2 = escc.postNR, m1Name = "preNR", m2Name = "postNR", minMut = 5)

forestPlot(mafCompareRes = NR_pre_post, pVal = 0.59, color = c("maroon", "royalblue"), geneFontSize = 0.6)


# 绘制瀑布图 -------------------------------------------------------------------

custom_colors <- c(
  "Missense_Mutation" = "#BC5659", "Nonsense_Mutation" = "#EAAA60",
  "Splice_Site" = "#B7B2D0", "Frame_Shift_Del" = "#84C3B7",
  "Frame_Shift_Ins" = "#B1DBB3", "In_Frame_Del" = "#458A74",
  "In_Frame_Ins" = "#3881BA",
  "Multi_Hit" = "#8DA0CB"
)

head(postNR_R$results, 6) -> postNR_R_deg
coOncoplot(
  m1 = escc.postR, m2 = escc.postNR,
  m1Name = "postR", m2Name = "postNR",
  genes = postNR_R_deg$Hugo_Symbol,
  sortByM1 = T,
  removeNonMutated = F,
  col = custom_colors,
  anno_height = 2,
  geneNamefont = 0.8,
  bgCol = "#EFF6FC"
)
# width = 6.5, height=1.83

head(preNR_R$results, 6) -> preNR_R_deg
coOncoplot(m1 = escc.preNR, m2 = escc.preR, m1Name = "preNR", m2Name = "preR", genes = preNR_R_deg$Hugo_Symbol, removeNonMutated = F, gene_mar = 3, titleFontSize = 1)
coOncoplot(
  m1 = escc.preR, m2 = escc.preNR,
  m1Name = "preR", m2Name = "preNR",
  genes = preNR_R_deg$Hugo_Symbol,
  sortByM1 = T,
  removeNonMutated = F,
  col = custom_colors,
  anno_height = 2,
  geneNamefont = 0.8,
  bgCol = "#EFF6FC"
)


head(R_pre_post$results, 2) -> R_post_pre_deg
coOncoplot(m1 = escc.preR, m2 = escc.postR, m1Name = "preR", m2Name = "postR", genes = R_post_pre_deg$Hugo_Symbol, removeNonMutated = F)

head(NR_pre_post$results, 3) -> NR_post_pre_deg
coOncoplot(m1 = escc.preNR, m2 = escc.postNR, m1Name = "preNR", m2Name = "postNR", genes = NR_post_pre_deg$Hugo_Symbol, removeNonMutated = F)


# save image --------------------------------------------------------------

save.image("src_wes/compare_diff_mutation.RData")
load("./src_wes/compare_diff_mutation.RData")
