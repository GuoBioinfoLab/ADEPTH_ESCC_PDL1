# library -----------------------------------------------------------------

library(pheatmap)
library(dplyr)
library(tidyr)

# load data ---------------------------------------------------------------

cn_postNR_matrix <- readr::read_tsv("~/wes_pipline/05-cnv-gistic2/resu_postC_NR/all_lesions.conf_90.txt", show_col_types = F) %>%
  dplyr::filter(`Amplitude Threshold` != "Actual Copy Change Given") %>%
  dplyr::select(-"...30")

cn_postR_matrix <- readr::read_tsv("~/wes_pipline/05-cnv-gistic2/resu_postC_R/all_lesions.conf_90.txt", show_col_types = F) %>%
  dplyr::filter(`Amplitude Threshold` != "Actual Copy Change Given") %>%
  dplyr::select(-"...14")

cn_preNR_matrix <- readr::read_tsv("~/wes_pipline/05-cnv-gistic2/resu_preC_NR/all_lesions.conf_90.txt", show_col_types = F) %>%
  dplyr::filter(`Amplitude Threshold` != "Actual Copy Change Given") %>%
  dplyr::select(-"...30")


# 准备数据并绘制热图 ---------------------------------------------------------------


cn_postNR_matrix %>% # change
  dplyr::filter(Descriptor %in% c("1q11", "11q13.3", "19p13.2", "19q12", "14q32.33")) %>%
  dplyr::select(
    -`Wide Peak Limits`, -`Peak Limits`, -`Region Limits`, -`Broad or Focal`, -`Amplitude Threshold`,
    -`Residual q values after removing segments shared with higher peaks`
  ) -> for_postNR_heatmap

cn_postR_matrix %>% # change
  dplyr::filter(Descriptor %in% c("1q11")) %>%
  dplyr::select(
    -`Wide Peak Limits`, -`Peak Limits`, -`Region Limits`, -`Broad or Focal`, -`Amplitude Threshold`,
    -`Residual q values after removing segments shared with higher peaks`
  ) -> for_postR_heatmap

cn_preNR_matrix %>% # change
  dplyr::filter(Descriptor %in% c("1q11", "4q13.3", "11q13.3", "2q22.1", "9p21.3")) %>%
  dplyr::select(
    -`Wide Peak Limits`, -`Peak Limits`, -`Region Limits`, -`Broad or Focal`, -`Amplitude Threshold`,
    -`Residual q values after removing segments shared with higher peaks`
  ) -> for_preNR_heatmap

# 合并Nmet & Met

for_postR_heatmap %>%
  mutate(across(
    .cols = -c(`Unique Name`, `q values`, "Descriptor"), # 除了这两列，其他都处理
    .fns = ~ ifelse(grepl("Deletion", `Unique Name`), ifelse(. == 0, 0, -.), .)
  )) %>%
  dplyr::select(-`q values`) %>%
  tidyr::unite(col = "cytoband", `Unique Name`, Descriptor, sep = "_") -> .for_postR_heatmap_changed

for_postNR_heatmap %>%
  mutate(across(
    .cols = -c(`Unique Name`, `q values`, "Descriptor"), # 除了这两列，其他都处理
    .fns = ~ ifelse(grepl("Deletion", `Unique Name`), ifelse(. == 0, 0, -.), .)
  )) %>%
  dplyr::select(-`q values`) %>%
  tidyr::unite(col = "cytoband", `Unique Name`, Descriptor, sep = "_") -> .for_postNR_heatmap_changed

for_preNR_heatmap %>%
  mutate(across(
    .cols = -c(`Unique Name`, `q values`, "Descriptor"), # 除了这两列，其他都处理
    .fns = ~ ifelse(grepl("Deletion", `Unique Name`), ifelse(. == 0, 0, -.), .)
  )) %>%
  dplyr::select(-`q values`) %>%
  tidyr::unite(col = "cytoband", `Unique Name`, Descriptor, sep = "_") -> .for_preNR_heatmap_changed


# lolipop: ---------------------------------------------------


# 准备数据

.for_preNR_heatmap_changed %>% # 改这个参数
  tidyr::pivot_longer(
    cols = -(1:1), # 保留前3列：Unique Name, Descriptor, q values
    names_to = "Sample",
    values_to = "Coverage"
  ) %>%
  dplyr::mutate(cnv_type = case_when(
    Coverage == "0" ~ "neutral",
    Coverage == "-1" ~ "loss",
    Coverage == "-2" ~ "loss",
    TRUE ~ "gain"
  )) %>%
  dplyr::group_by(cytoband, cnv_type) %>%
  dplyr::summarise(cnv_num = n(), .groups = "drop") %>%
  # dplyr::ungroup(cytoband, cnv_type) %>%
  tidyr::spread(cnv_type, cnv_num) %>%
  mutate(across((ncol(.) - 2):ncol(.), ~ replace_na(., 0))) -> .cytoband_for_lollipop

# 绘图

.cytoband_for_lollipop %>%
  arrange(desc(gain), desc(loss)) %>%
  mutate(cytoband = factor(cytoband, levels = cytoband)) -> sort_for_plot

# save source data

readr::write_csv(sort_for_plot, "./data/wes_output/cnv_lollipop_preNR.csv")

ggplot(sort_for_plot) +
  geom_linerange(aes(x = cytoband, ymin = 0, ymax = ifelse(gain > loss, gain, loss)),
    linewidth = 2, color = "grey"
  ) +
  geom_point(aes(x = cytoband, y = gain), color = "#E66D50", size = 4) +
  geom_point(aes(x = cytoband, y = loss), color = "#274752", size = 4) +
  labs(x = NULL, y = "Samples (postNR)") +
  ylim(0, 20) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9, face = "bold", color = "black")
  )

# width = 6, height=3


# 绘制瀑布图 -------------------------------------------------------------------

cnv_info_met@cytoband.summary %>%
  dplyr::filter(Cytoband %in% met_amp_del_uniq_for_60_month_sig) -> band_lst

gisticOncoPlot(
  gistic = cnv_info_met,
  top = 100,
  gene_mar = 6,
  bands = band_lst$Unique_Name,
  # clinicalData = .all_sample_info,
  # clinicalFeatures = "Group",
  # sortByAnnotation = T,
  removeNonAltered = F,
  colors = c(Amp = "#D95F02", Del = "#1B9E77")
)


# save image --------------------------------------------------------------

save.image("./RData_wes/lollipop_cnv.RData")
load("./RData_wes/lollipop_cnv.RData")
