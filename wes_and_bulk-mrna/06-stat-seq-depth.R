# library -----------------------------------------------------------------

library(magrittr)
library(scales)
library(ggplot2)
library(dplyr)
library(readr)
library(stringi)
library(tidyr)
library(VennDiagram)


# load data ---------------------------------------------------------------

seq_depth_normal <- readr::read_delim("data/wes_output/coverage_stat_normal.txt", delim = "\t", col_names = c("sample", "seq_depth"), show_col_types = FALSE) %>%
  dplyr::mutate(Group = "Normal")
seq_depth_pretumor <- readr::read_delim("data/wes_output/coverage_stat_pretumor.txt", delim = "\t", col_names = c("sample", "seq_depth"), show_col_types = FALSE) %>%
  dplyr::mutate(Group = "preNICB")
seq_depth_tumor <- readr::read_delim("data/wes_output/coverage_stat_tumor.txt", delim = "\t", col_names = c("sample", "seq_depth"), show_col_types = FALSE) %>%
  dplyr::mutate(Group = "postNICB")


# merge normal & tumor ----------------------------------------------------

rbind(seq_depth_normal, seq_depth_pretumor, seq_depth_tumor) -> .for_plot

readr::write_csv(.for_plot, "./data/wes_output/seq_depth.csv")

# plot for seq_depth ------------------------------------------------------

.mean_normal <- mean(seq_depth_normal$seq_depth)
.mean_pretumor <- mean(seq_depth_pretumor$seq_depth)
.mean_tumor <- mean(seq_depth_tumor$seq_depth)

nature_colors <- c(
  "Normal" = "#4E79A7",
  "preNICB" = "#59A14F",
  "postNICB" = "#E15759"
)

ggplot(
  .for_plot,
  aes(
    x = sample,
    y = seq_depth,
    color = factor(Group),
    shape = factor(Group)
  )
) +
  geom_point(size = 3, alpha = 0.9) +

  # mean lines
  geom_hline(
    yintercept = .mean_normal,
    linetype = "dashed",
    color = nature_colors["Normal"],
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = .mean_pretumor,
    linetype = "dashed",
    color = nature_colors["preNICB"],
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = .mean_tumor,
    linetype = "dashed",
    color = nature_colors["postNICB"],
    linewidth = 0.6
  ) +

  # annotations
  annotate("text",
    x = Inf, y = .mean_normal,
    label = paste("Mean normal:", round(.mean_normal, 2)),
    hjust = 3.12, vjust = -0.6,
    size = 3.5
  ) +
  annotate("text",
    x = Inf, y = .mean_pretumor,
    label = paste("Mean preNICB:", round(.mean_pretumor, 2)),
    hjust = 2.92, vjust = -0.4,
    size = 3.5
  ) +
  annotate("text",
    x = Inf, y = .mean_tumor,
    label = paste("Mean postNICB:", round(.mean_tumor, 2)),
    hjust = 2.8, vjust = -0.6,
    size = 3.5
  ) +

  # manual scales
  scale_color_manual(values = nature_colors) +
  # scale_shape_manual(values = c(16, 17, 15)) +

  labs(
    title = "Normal: 24| preNICB: 24| postNICB: 24",
    y = "Depth of coverage",
    color = "Sample type",
    shape = "Sample type"
  ) +

  # remove x axis completely
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    )
  )


# save image --------------------------------------------------------------

save.image("RData_wes/stat_seq_depth.RData")
load("./RData_wes/stat_seq_depth.RData")
