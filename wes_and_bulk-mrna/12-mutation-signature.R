library(maftools)
library(Rcpp)
library(sigminer)
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyr)
library(ggplot2)


# load data ---------------------------------------------------------------

.all_sample_info <- readr::read_tsv("~/wes_data/wes_sample_info.txt", show_col_types = FALSE) %>%
  dplyr::rename(Tumor_Sample_Barcode = "sample_id")

.all_48_maf_file <- read.maf("~/wes_pipline/04-modified_add_tumor_barcode_maf/all_48_filtered_file.maf", clinicalData = .all_sample_info)




mats <- mt_tally <- sig_tally(
  .all_48_maf_file,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
  useSyn = TRUE,
  mode = "ALL"
)

# 展示突变图谱
show_catalogue(mt_tally$SBS_96 %>% t(), mode = "SBS", style = "cosmic")
show_catalogue(mt_tally$DBS_78 %>% t(), mode = "DBS", style = "cosmic")
show_catalogue(mt_tally$ID_28 %>% t(), mode = "ID", style = "cosmic")

# 提取signatures数量选择（半自动）
# SBS
est <- sig_estimate(mt_tally$SBS_96, range = 2:10, nrun = 50, verbose = TRUE)
show_sig_number_survey2(est$survey)

readr::write_csv(est$survey, "./data/wes_output/mutation_signature/nmf_figs2e.csv")

# DBS
est_dbs <- sig_estimate(mt_tally$DBS_78, range = 2:15, nrun = 20, verbose = TRUE)
show_sig_number_survey2(est_dbs$survey)

# ID
est_id <- sig_estimate(mt_tally$DBS_78, range = 2:15, nrun = 20, verbose = TRUE)
show_sig_number_survey2(est_id$survey)

# 提取上述选择数量的signatures： 半自动的结果显示选4或者3都可以
sigs <- sig_extract(mt_tally$SBS_96, n_sig = 4, nrun = 50)
# 详细碱基特征
show_sig_profile(sigs,
  mode = "SBS",
  style = "cosmic",
  bar_width = 0.9,
  free_space = "free_x",
  base_size = 13,
  font_scale = 0.9
)

readr::write_rds(sigs, "./data/wes_output/mutation_signature/figs2f.rds")

# 由于样本量太少，换种限制性方法
sigs_new <- sig_fit(
  t(mt_tally$SBS_96),
  sig_index = c("1", "6", "13", "15", "25")
)


# 绘图（每种signature在不同样本中的占比）

get_sig_exposure(sigs) # 用sigs，是手动提取的值
show_sig_exposure(sigs_new, rm_space = TRUE, style = "cosmic") # 在后面绘制饼图


# 根据已知的Signature提取与cosmic相关活动度 -----------------------------------------------------

all_cosmic_fit <- get_sig_similarity(sigs, sig_db = "SBS") # legacy 或者SBS，DBS， ID
pheatmap::pheatmap(all_cosmic_fit$similarity)

# 筛选子集绘图（余弦相似性热图）
sim_mat <- all_cosmic_fit$similarity
sim_sub <- sim_mat[, c("SBS13", "SBS6", "SBS25"), drop = FALSE]
pheatmap::pheatmap(sim_sub)

readr::write_csv(as.data.frame(sim_sub), "./data/wes_output/mutation_signature/heatmap_figs2g.csv")


all_cosmic_fit <- sig_fit(mt_tally$SBS_96 %>% t(), sig_index = c("6", "13", "25"), sig_db = "SBS")

show_sig_fit(all_cosmic_fit, palette = NULL) + ggpubr::rotate_x_text()


# 对结果绘图进行调整和显著性检验 ---------------------------------------------------------

# COSMIC在R VS NR中的差异

library(tibble)
all_cosmic_fit_clinical_info <- rownames_to_column(as.data.frame(all_cosmic_fit), var = "cosmic") %>%
  tibble::as_tibble() %>%
  tidyr::gather(key = "Tumor_Sample_Barcode", value = "signature_exposure", -cosmic) %>%
  dplyr::left_join(.all_sample_info, by = "Tumor_Sample_Barcode")


# 绘制分组箱线图 -----------------------------------------------------------------
library(ggpubr)

all_cosmic_fit_clinical_info %>%
  dplyr::filter(signature_exposure > 0) -> .filtered_for_plot

.filtered_for_plot$response <- factor(
  .filtered_for_plot$response,
  levels = c("NR", "R")
)

.filtered_for_plot %>%
  dplyr::filter(scrna_match != "no") -> scmatch_filtered_for_plot

# 绘图：preNICB vs postNICB ---------------------------------------------------------
ggplot(
  .filtered_for_plot,
  aes(x = cosmic, y = signature_exposure, fill = response)
) +
  stat_boxplot(geom = "errorbar", linewidth = 0.6) +
  geom_boxplot(outlier.shape = NA, alpha = 1) +
  geom_jitter(
    # aes(color = response),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 1,
    alpha = 0.7
  ) +
  facet_wrap(~time_point) +
  ggpubr::stat_compare_means(
    method = "t.test",
    label = "p.format",
    hide.ns = FALSE,
    # label.y = 80
  ) +
  scale_fill_manual(values = c("#1E77B4", "#E36F4D")) +
  labs(x = "", y = "Signature exposure", title = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12, colour = "black"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "Gainsboro", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13)
  )

# 绘图：R vs NR ---------------------------------------------------------
readr::write_csv(.filtered_for_plot, "./data/wes_output/mutation_signature/boxplot_figs2h.csv")

ggplot(
  .filtered_for_plot,
  aes(x = cosmic, y = signature_exposure, fill = time_point)
) +
  stat_boxplot(geom = "errorbar", linewidth = 0.6) +
  geom_boxplot(outlier.shape = NA, alpha = 1) +
  geom_jitter(
    # aes(color = response),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    size = 1,
    alpha = 0.7
  ) +
  facet_wrap(~response) +
  ggpubr::stat_compare_means(
    method = "t.test",
    label = "p.format",
    hide.ns = F
  ) +
  # facet_wrap(~SV_type, ncol = 1, scales = "free_y") + #在facet_wrap函数中添加scales = "free_y"
  # scale_fill_manual(values = c("#F18427", "#2674B2")) + #两种类型的配色
  scale_fill_manual(values = c("#2E9591", "#F5B494")) + # 四种阶段的配色
  # coord_cartesian(ylim = c(0, 200)) +  # ⭐固定y轴范围⭐
  labs(x = "", y = "Signature exposure", title = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 12, colour = "black"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "Gainsboro", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.5),

    ## === 关键：坐标轴字体大小 ===
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13)
  )


# 上述基因是否是Met vs N-Met中的差异表达基因

dMMR_genes <- c(
  "POLD3", "MLH3", "MSH6", "RPA4", "LIG1", "MLH1", "MSH2", "MSH3",
  "PCNA", "PMS2", "POLD1", "POLD2", "POLD4", "RFC1", "RFC2", "RFC3",
  "RFC4", "RFC5", "RPA1", "RPA2", "RPA3", "SSBP1", "EXO1"
)

# save image --------------------------------------------------------------

save.image("./RData_wes/mutation_signatures_analysis.RData")
load("./RData_wes/mutation_signatures_analysis.RData")
