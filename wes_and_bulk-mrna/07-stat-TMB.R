# library -----------------------------------------------------------------

library(maftools)
library(grid)
library(magrittr)
library(ggpubr)
library(ggbreak)
library(ggplot2)


# load all maf file -------------------------------------------------------

.all_sample_info <- readr::read_tsv("~/wes_data/wes_sample_info.txt", show_col_types = FALSE) %>%
  dplyr::rename(Tumor_Sample_Barcode = "sample_id")

.all_48_maf_file <- read.maf("~/wes_pipline/04-modified_add_tumor_barcode_maf/all_48_filtered_file.maf", clinicalData = .all_sample_info)

# TMB in each group-------------------------------------------

.all_tmb <- tmb(.all_48_maf_file, captureSize = 35, logScale = F)

clinical_data <- .all_48_maf_file@clinical.data
.all_tmb_clinical_data <- merge(.all_tmb, clinical_data, by = "Tumor_Sample_Barcode")

readr::write_tsv(.all_tmb_clinical_data, "~/Rproject/sczhoujianfeng/data/wes_output/all_48_tmb_clinical_info.txt")

# 绘制TMB不同条件下的散点图 ----------------------------------------------------------

annotation_colors <- list(
  time_point = c("postNICB" = "#1E77B4", "preNICB" = "#8C7ACB"),
  response = c("NR" = "#EBB261", "R" = "#E36F4D")
)

.all_tmb_clinical_data <- .all_tmb_clinical_data %>%
  arrange(response, total_perMB) %>%
  mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels = Tumor_Sample_Barcode))

# 画图

# 响应与否

library(ggplot2)
library(ggpubr)
library(rstatix)

.all_tmb_clinical_data$time_point <- factor(
  .all_tmb_clinical_data$time_point,
  levels = c("preNICB", "postNICB")
)

ggplot(
  data = .all_tmb_clinical_data,
  aes(x = response, y = total_perMB, fill = time_point) # 用 fill 区分时间点
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = time_point), alpha = 0.8, size = 1, position = position_dodge(width = 0.8)) +
  stat_compare_means(
    comparisons = list(c("NR", "R")),
    method = "wilcox.test",
    aes(group = response) # 横坐标分组比较
  ) +
  scale_fill_manual(name = "Time Point", values = annotation_colors$time_point) +
  scale_color_manual(name = "Time Point", values = annotation_colors$time_point) +
  labs(x = "", y = "TMB/MB") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(),
    panel.border = element_rect(fill = NA, size = 0.5),
    legend.position = "right", # 显示图例
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black")
  )


# 手动统计显著性差异

.all_tmb_clinical_data %>%
  filter(time_point %in% c("preNICB", "postNICB")) %>%
  group_by(response) %>%
  wilcox_test(total_perMB ~ time_point) %>%
  mutate(
    label = sprintf("p = %.3f", p),
    y.position = max(.all_tmb_clinical_data$total_perMB, na.rm = TRUE) * 1.05
  )


# width=3.5, height=3.5

save.image("RData_wes/TMB_point_plot.RData")
load("./RData_wes/TMB_point_plot.RData")
