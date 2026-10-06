# library -----------------------------------------------------------------

library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)


# load data ---------------------------------------------------------------

sam_info <- readr::read_tsv("~/wes_pipline/08-HLA-LOH-spechla/SpecHLA/wes_sample_info_abbr.txt", show_col_types = F) %>%
  dplyr::rename(Sample = "sample_name")

dpa1_info <- readr::read_tsv("~/wes_pipline/08-HLA-LOH-spechla/SpecHLA/all_dpa1_loh_merged.txt", show_col_types = F) %>%
  dplyr::mutate(Sample = str_remove(Sample, "^dedup-")) %>%
  dplyr::mutate(Sample = str_remove(Sample, "WuJun-")) %>%
  dplyr::mutate(Sample = Sample %>%
    str_replace("^glj", "gcj") %>%
    str_replace("^xsz", "wsz"))


dpb1_info <- readr::read_tsv("~/wes_pipline/08-HLA-LOH-spechla/SpecHLA/all_dpb1_loh_merged.txt", show_col_types = F) %>%
  dplyr::mutate(Sample = str_remove(Sample, "^dedup-")) %>%
  dplyr::mutate(Sample = str_remove(Sample, "WuJun-")) %>%
  dplyr::mutate(Sample = Sample %>%
    str_replace("^glj", "gcj") %>%
    str_replace("^xsz", "wsz"))


dqa1_info <- readr::read_tsv("~/wes_pipline/08-HLA-LOH-spechla/SpecHLA/all_dqa1_loh_merged.txt", show_col_types = F) %>%
  dplyr::mutate(Sample = str_remove(Sample, "^dedup-")) %>%
  dplyr::mutate(Sample = str_remove(Sample, "WuJun-")) %>%
  dplyr::mutate(Sample = Sample %>%
    str_replace("^glj", "gcj") %>%
    str_replace("^xsz", "wsz"))
dqb1_info <- readr::read_tsv("~/wes_pipline/08-HLA-LOH-spechla/SpecHLA/all_dqb1_loh_merged.txt", show_col_types = F) %>%
  dplyr::mutate(Sample = str_remove(Sample, "^dedup-")) %>%
  dplyr::mutate(Sample = str_remove(Sample, "WuJun-")) %>%
  dplyr::mutate(Sample = Sample %>%
    str_replace("^glj", "gcj") %>%
    str_replace("^xsz", "wsz"))
drb1_info <- readr::read_tsv("~/wes_pipline/08-HLA-LOH-spechla/SpecHLA/all_drb1_loh_merged.txt", show_col_types = F) %>%
  dplyr::mutate(Sample = str_remove(Sample, "^dedup-")) %>%
  dplyr::mutate(Sample = str_remove(Sample, "WuJun-")) %>%
  dplyr::mutate(Sample = Sample %>%
    str_replace("^glj", "gcj") %>%
    str_replace("^xsz", "wsz"))


# merge sam_info AND loh_info ---------------------------------------------

dpa1_info %>%
  dplyr::left_join(sam_info, by = "Sample") -> .dpa1

readr::write_csv(.dpa1, "./data/wes_output/hla_barplot/dpa1_barplot.csv")

dpb1_info %>%
  dplyr::left_join(sam_info, by = "Sample") -> .dpb1

readr::write_csv(.dpb1, "./data/wes_output/hla_barplot/dpb1_barplot.csv")

dqa1_info %>%
  dplyr::left_join(sam_info, by = "Sample") -> .dqa1

dqb1_info %>%
  dplyr::left_join(sam_info, by = "Sample") -> .dqb1

drb1_info %>%
  dplyr::left_join(sam_info, by = "Sample") -> .drb1


# width = 10， height=5

# 查看5种HLA-LOH的分布特征 --------------------------------------------------------

.dpb1 %>%
  dplyr::select(sample_id, LOH, response, time_point) -> for_stat_loh_plot

p_dpb1 <- ggplot(for_stat_loh_plot, aes(x = LOH, fill = response)) +
  geom_bar(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      R  = "#1F3A5F",
      NR = "#9E2A2B"
    ),
    labels = c(
      R  = "R",
      NR = "NR"
    )
  ) +
  labs(
    x = "LOH (DPA1)",
    y = "Number of samples",
    fill = "Response"
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )


p_dpa1 + p_dpb1 + p_dqa1 + p_dqb1 + p_drb1

# width =8， height = 4


# 堆叠条形图 （post R）-------------------------------------------------------------------

.drb1 %>%
  dplyr::filter(time_point == "preNICB") %>%
  dplyr::select(LOH, response) %>%
  dplyr::group_by(LOH, response) %>%
  dplyr::summarise(num = n(), .groups = "drop") -> for_hla_plot


# 显著性检验 -------------------------------------------------------------------
# post-hla: .dpa1
tab <- matrix(
  c(
    1, 3, # R: Y / N
    7, 13
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)
# pre-hla：.dpa1
tab <- matrix(
  c(
    1, 3, # R: Y / N
    9, 11
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)

# post-.dpb1
tab <- matrix(
  c(
    0, 4, # R: Y / N
    5, 15
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)
# pre-.dpb1
tab <- matrix(
  c(
    1, 3, # R: Y / N
    8, 12
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)

# post-dqa1
tab <- matrix(
  c(
    3, 1, # R: Y / N
    8, 12
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)
# pre-dqa1
tab <- matrix(
  c(
    3, 1, # R: Y / N
    10, 10
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)

# post-dqb1
tab <- matrix(
  c(
    3, 1, # R: Y / N
    7, 13
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)
# pre-dqb1
tab <- matrix(
  c(
    2, 2, # R: Y / N
    12, 8
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)

# post-drb1
tab <- matrix(
  c(
    2, 1, # R: Y / N
    7, 10
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)
# pre-drb1
tab <- matrix(
  c(
    3, 1, # R: Y / N
    12, 7
  ), # NR: Y / N
  nrow = 2,
  byrow = TRUE
)

rownames(tab) <- c("R", "NR")
colnames(tab) <- c("Positive", "Negative")


chisq.test(tab)

fisher.test(tab)



# 绘制堆叠条形图 -----------------------------------------------------------------

drb1_pre <- ggplot(for_hla_plot, aes(x = response, y = num, fill = LOH)) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "Y" = "#C17C54", # "#C17C54"，#6A9F58， #4C72B0
      "N" = "#DDDDDD"
    )
  ) +
  labs(
    x = NULL,
    y = "Sample number",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    plot.margin = margin(5.5, 15, 5.5, 5.5)
  )


(dpa1_pre | dpa1_post) /
  (dpb1_pre | dpb1_post) /
  (dqa1_pre | dqa1_post) /
  (dqb1_pre | dqb1_post) /
  (drb1_pre | drb1_post)


# save image --------------------------------------------------------------

save.image("RData_wes/stat_hla_resut_mhcii.RData")
load("./RData_wes/stat_hla_resut_mhcii.RData")
