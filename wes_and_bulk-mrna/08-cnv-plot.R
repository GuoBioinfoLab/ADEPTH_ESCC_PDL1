# library -----------------------------------------------------------------

library(maftools)

our_all.lesions <- "~/wes_pipline/05-cnv-gistic2/resu_postC_NR/all_lesions.conf_90.txt"
our_amp.genes <- "~/wes_pipline/05-cnv-gistic2/resu_postC_NR/amp_genes.conf_90.txt"
our_del.genes <- "~/wes_pipline/05-cnv-gistic2/resu_postC_NR/del_genes.conf_90.txt"
our_scores.gis <- "~/wes_pipline/05-cnv-gistic2/resu_postC_NR/scores.gistic"

escc.gistic <- readGistic(
  gisticAllLesionsFile = our_all.lesions,
  gisticAmpGenesFile = our_amp.genes,
  gisticDelGenesFile = our_del.genes,
  gisticScoresFile = our_scores.gis, isTCGA = F
)

gisticChromPlot(
  gistic = escc.gistic,
  # markBands = "all",
  ref.build = "hg38",
  fdrCutOff = 0.05,
  cytobandTxtSize = 0.6,
  cytobandOffset = 0.05,
  txtSize = 0.8,
  # mutGenes = c("EPHA3"),
  color = c("#E66D50", "#2A726F")
)

# width=7.5, height=1.5


# save image --------------------------------------------------------------

save.image("RData_wes/plot_for_CNV_barplot.RData")
load("./RData_wes/plot_for_CNV_barplot.RData")
