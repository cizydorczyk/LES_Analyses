library(tidyverse)
library(ggpubr)
library(patchwork)
library(ape)
library(tidystats)
library(vegan)

geneDataHandle = "/home/conrad/les_complete/pangenome_analyses/LC8-panaroo/EDITED-gene_presence_absence.Rtab"
geneData = read.table(geneDataHandle, header=TRUE, row.names=1, stringsAsFactors=F, comment.char="", sep="\t")

# Read in isolate list and create Canada and UK lists:
isolateList = read.table(
  "/home/conrad/les_complete/les-complete-list-XNum.txt",
  header = F,
  stringsAsFactors = F
)
isolateList = isolateList$V1
isolateList = c(isolateList, "LESB58")
isolateList = isolateList[!(isolateList %in% c("PHELES19"))]

canadaList = isolateList[1:50] # all Canada isolates together
ukList = isolateList[51:205] # all UK isolates minus PHELES19
canadaList2 = isolateList[1:11] # deep branching clade Canada isolates
canadaList1 = isolateList[12:50] # non-deep branching Canada isolates

geneData = geneData[,isolateList]

geneDists = dist.gene(t(geneData), method="pairwise", pairwise.deletion=T)
geneDists = as.data.frame(as.matrix(geneDists))

# Convert matrix to long format:
outputDists = c()
namesDists = c()
colsDists = c()

start_row = 1
start_col = 1
for (i in seq(1, 206, 1)) {
  new_vec = geneDists[, start_col]
  new_vec = new_vec[-(1:start_row)]
  currentCol = colnames(geneDists)[start_col]
  currentRows = rownames(geneDists)[-(1:start_row)]
  currentCols = rep(currentCol, length(currentRows))
  outputDists = c(outputDists, new_vec)
  namesDists = c(namesDists, currentRows)
  colsDists = c(colsDists, currentCols)
  if (start_row + 1 <= nrow(geneDists) &
      start_col + 1 <= ncol(geneDists)) {
    start_row = start_row + 1
    start_col = start_col + 1
  }
}

# Generate output df:
outputDf = data.frame(Isolate1 = colsDists,
                      Isolate2 = namesDists,
                      Gene_Distance = outputDists)

# Label each comparison as intra-Canada, intra-UK, or inter-Canada/UK:
canadaUkComp = c()
for (i in seq(1, length(rownames(outputDf)), 1)) {
  isolate1 = outputDf$Isolate1[i]
  isolate2 = outputDf$Isolate2[i]
  isolate1Origin = ""
  isolate2Origin = ""
  if (isolate1 %in% canadaList) {
    isolate1Origin = "Canada"
  } else {
    isolate1Origin = "UK"
  }
  if (isolate2 %in% canadaList) {
    isolate2Origin = "Canada"
  } else {
    isolate2Origin = "UK"
  }
  canadaUkComp = c(canadaUkComp,
                   paste(isolate1Origin, "-", isolate2Origin, sep = ""))
}

outputDf$CanadaVsUK = canadaUkComp

# Label each comparison as above but split Canada into deep branching and non-deep branching:
canadaUkComp2 = c()

for (i in seq(1, length(rownames(outputDf)), 1)) {
  isolate1 = outputDf$Isolate1[i]
  isolate2 = outputDf$Isolate2[i]
  isolate1Origin = ""
  isolate2Origin = ""
  if (isolate1 %in% canadaList1) {
    isolate1Origin = "Canada1"
  } else if (isolate1 %in% canadaList2) {
    isolate1Origin = "Canada2"
  } else if (isolate1 %in% ukList) {
    isolate1Origin = "UK"
  }
  
  if (isolate2 %in% canadaList1) {
    isolate2Origin = "Canada1"
  } else if (isolate2 %in% canadaList2) {
    isolate2Origin = "Canada2"
  } else if (isolate2 %in% ukList) {
    isolate2Origin = "UK"
  }
  canadaUkComp2 = c(canadaUkComp2,
                    paste(isolate1Origin, "-", isolate2Origin, sep = ""))
}

outputDf$CanadaVsUK_2 = canadaUkComp2

################################################################################
# Calculate summary statistics per group
################################################################################

# Calculate statistics for SNP distances for CanadaVsUk comparison:
canadaVsUkGrouped = outputDf %>% group_by(CanadaVsUK) %>% summarise(
  MeanDist = mean(Gene_Distance),
  MedianDist = median(Gene_Distance),
  MinDist = min(Gene_Distance),
  MaxDist = max(Gene_Distance),
  Q1 = quantile(Gene_Distance, 0.25),
  Q3 = quantile(Gene_Distance, 0.75)
)
canadaVsUkGrouped

# Calculate statistics for SNP distances for CanadaSplitVsUK
canadaSplitVsUkGrouped = outputDf %>% group_by(CanadaVsUK_2) %>% summarise(
  MeanDist = mean(Gene_Distance),
  MedianDist = median(Gene_Distance),
  MinDist = min(Gene_Distance),
  MaxDist = max(Gene_Distance),
  Q1 = quantile(Gene_Distance, 0.25),
  Q3 = quantile(Gene_Distance, 0.75)
)
canadaSplitVsUkGrouped

################################################################################
# Plot histograms/densities
################################################################################

# # Plot histograms of CanadaVsUk SNP distance distributions:
# plot1 = ggdensity(data = outputDf,
#           x = "SNP_Distance",
#           add = "median",
#           rug = F,
#           color = "CanadaVsUK",
#           fill = "CanadaVsUK",
#           palette = c("#00AFBB", "#E7B800", "#f07dc2", "#80f07d", "#877df0", "#FC4E07")) +
#   theme(legend.title = element_blank()) +
#   xlab("SNP Distance") +
#   ylab("Density")
# plot1

# Plot histograms of CanadaVsUK_2:
outputDf$CanadaVsUK_2 = factor(
  outputDf$CanadaVsUK_2,
  levels = c(
    "Canada1-Canada1",
    "Canada1-UK",
    "UK-UK",
    "Canada2-Canada2",
    "Canada2-Canada1",
    "Canada2-UK"
  )
)
# 
# ## This first plot is just to get the legend...
# plot2a = ggdensity(
#   data = outputDf,
#   x = "Gene_Distance",
#   add = "mean",
#   rug = TRUE,
#   color = "CanadaVsUK_2",
#   fill = "CanadaVsUK_2",
#   palette = c(
#     "#00AFBB",
#     "#E7B800",
#     "#f07dc2",
#     "#80f07d",
#     "#877df0",
#     "#FC4E07"
#   ),
#   legend = "top"
# ) +
#   theme(legend.title = element_blank(),
#         legend.text = element_text(size = 13))
# plot2aLegend = get_legend(plot2a)
# 
# # And this plot is used for actual plotting:
# plot2b = ggdensity(
#   data = outputDf,
#   x = "Gene_Distance",
#   add = "median",
#   rug = F,
#   color = "CanadaVsUK_2",
#   fill = "CanadaVsUK_2",
#   palette = c(
#     "#00AFBB",
#     "#E7B800",
#     "#f07dc2",
#     "#80f07d",
#     "#877df0",
#     "#FC4E07"
#   ),
#   legend = "none"
# ) +
#   font("xy.text", size = 10) +
#   xlab("") +
#   ylab("")

# This plot includes only intra-group comparisons + UK-Canada1:
plot3 = ggdensity(
  data = outputDf,
  x = "Gene_Distance",
  add = "median",
  rug = F,
  color = "CanadaVsUK_2",
  fill = "CanadaVsUK_2",
  palette = c(
    "#00AFBB",
    "#E7B800",
    "#f07dc2",
    "#80f07d",
    "#877df0",
    "#FC4E07"
  ),
  legend = "top",
  legend.title=""
) +
  xlab("Gene Distance") +
  ylab("Density") #+
  #scale_x_continuous(breaks = get_breaks(n = 12))
plot3

# # Put plots together:
# finalPlot = (as_ggplot(plot2aLegend)) / (plot3 + inset_element(
#   plot2b,
#   left = 0.45,
#   bottom = 0.41,
#   right = 0.95,
#   top = 0.9
# )) +
#   plot_layout(heights = c(0.5, 5))
# finalPlot

################################################################################
# Make violin plots
################################################################################

# Violin plots:

# # Plot for Canada split vs UK (not much point in doing all Canada vs UK):
# violinPlot1a = ggviolin(
#   data = outputDf,
#   x = "CanadaVsUK_2",
#   y = "SNP_Distance",
#   fill = "CanadaVsUK_2",
#   palette = c(
#     "#00AFBB",
#     "#E7B800",
#     "#f07dc2",
#     "#80f07d",
#     "#877df0",
#     "#FC4E07"
#   ),
#   legend = "top"
# ) +
#   scale_x_discrete(
#     limits = c(
#       "Canada2-Canada2",
#       "Canada1-Canada1",
#       "UK-UK",
#       "Canada1-UK",
#       "Canada2-Canada1",
#       "Canada2-UK"
#     )
#   ) +
#   theme(
#     axis.text.x = element_blank(),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 12)
#   ) +
#   xlab("") +
#   ylab("SNP Distance") +
#   ylim(0, 5000)
# violinPlotLegend = get_legend(violinPlot1a)
# 
# violinPlot1b = ggviolin(
#   data = outputDf,
#   x = "CanadaVsUK_2",
#   y = "SNP_Distance",
#   fill = "CanadaVsUK_2",
#   palette = c(
#     "#00AFBB",
#     "#E7B800",
#     "#f07dc2",
#     "#80f07d",
#     "#877df0",
#     "#FC4E07"
#   ),
#   legend = "none"
# ) +
#   scale_x_discrete(
#     limits = c(
#       "Canada1-Canada1",
#       "Canada1-UK",
#       "UK-UK",
#       "Canada2-Canada2",
#       "Canada2-Canada1",
#       "Canada2-UK"
#     )
#   ) +
#   theme(
#     axis.text.x = element_blank(),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 12)
#   ) +
#   xlab("") +
#   ylab("") +
#   ylim(0, 5000)

violinPlot2 = ggviolin(
  data = outputDf,
  x = "CanadaVsUK_2",
  y = "Gene_Distance",
  fill = "CanadaVsUK_2",
  palette = c(
    "#00AFBB",
    "#E7B800",
    "#f07dc2",
    "#80f07d",
    "#877df0",
    "#FC4E07"
  ),
  legend = "top",
  legend.title="",
  add = c("boxplot")
) +
  #scale_x_discrete(limits = c("Canada1-Canada1", "Canada1-UK", "UK-UK", "Canada2-Canada2")) +
  theme(
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  xlab("") +
  ylab("Gene Distance") +
  ylim(0, 1100)

violinPlot2 = violinPlot2 + stat_compare_means(
  comparisons = list(c("Canada1-UK", "Canada1-Canada1")),
  method.args = list(alternative =
                       "two.sided"),
  label.y = c(1000),
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
    symbols = c("****", "***", "**", "*", "ns")
  )
)
violinPlot2

################################################################################
# PCoA on SNPs (fails on p-distances)
################################################################################

# ## All isolates:
# # PCoA permANOVA:
# isolateListDf = data.frame(Isolate=isolateList, Origin=c(rep("Canada", length(canadaList)), rep("UK", length(ukList))))
# permanovaRes = adonis2(geneDists ~ factor(Origin), data = isolateListDf, permutations = 10000)
# permanovaRes
# 
# # Plots:
# pcoa1 = pcoa(geneDists, correction = "cailliez")
# 
# pcoa1Df = as.data.frame(pcoa1$vectors)
# pcoa1Df$Isolate = rownames(pcoa1Df)
# pcoa1Df$Origin = c(rep("Canada", length(c(canadaList1, canadaList2))), rep("UK", length(ukList)))
# 
# pcoaPlot1 = ggscatter(
#   data = pcoa1Df,
#   x = "Axis.1",
#   y = "Axis.2",
#   color = "Origin",
#   palette = c("red", "blue"),
#   legend.title = "Country of Origin:"
# ) +
#   font("legend.text", size = 12) +
#   stat_ellipse(data=pcoa1Df, aes(x=Axis.1, y=Axis.2, color=Origin, fill=Origin), geom = "polygon", alpha=0.2)
# 
# pcoaPlot2 = ggscatter(
#   data = pcoa1Df,
#   x = "Axis.2",
#   y = "Axis.3",
#   color = "Origin",
#   palette = c("red", "blue")
# ) +
#   stat_ellipse(data=pcoa1Df, aes(x=Axis.2, y=Axis.3, color=Origin, fill=Origin), geom = "polygon", alpha=0.2)
# 
# pcoaPlot3 = ggscatter(
#   data = pcoa1Df,
#   x = "Axis.3",
#   y = "Axis.4",
#   color = "Origin",
#   palette = c("red", "blue")
# ) +
#   stat_ellipse(data=pcoa1Df, aes(x=Axis.3, y=Axis.4, color=Origin, fill=Origin), geom = "polygon", alpha=0.2)
# 
# # pcoaPlot4 = ggscatter(
# #   data = pcoa1Df,
# #   x = "Axis.4",
# #   y = "Axis.5",
# #   color = "Origin",
# #   palette = c("red", "blue")
# # )
# 
# 
# screeData = data.frame(
#   Pos = seq(1, length(pcoa1$values$Cum_corr_eig[1:50])),
#   Cumulative_Perc = (pcoa1$values$Cum_corr_eig[1:50]) * 100
# )
# screeData$Col = c(rep("lightgrey", 4), rep("darkgrey", length(screeData$Pos) -
#                                              4))
# screePlot = ggbarplot(
#   data = screeData,
#   x = "Pos",
#   y = "Cumulative_Perc",
#   color = "Col",
#   fill = "Col",
#   palette = c("lightgrey", "darkgrey"),
#   legend = "none"
# ) +
#   scale_y_continuous(breaks = get_breaks(n = 10)) +
#   ylab("Cumulative Variance (%)") +
#   xlab("Axis Number")
# 
# 
# pcoaPlots = ggarrange(
#   pcoaPlot1,
#   pcoaPlot2,
#   pcoaPlot3,
#   screePlot,
#   nrow = 2,
#   ncol = 2,
#   common.legend = T,
#   legend = "top",
#   labels = "AUTO"
# )
# 
# pcoaPlots

## No Canada2:
# permANOVA:
geneDists2 = geneDists[-(1:11), -(1:11)]
isolateListDf2 = data.frame(Isolate=c(canadaList1, ukList), Origin=c(rep("Canada", length(canadaList1)), rep("UK", length(ukList))))
permanovaRes2 = adonis2(geneDists2 ~ factor(Origin), data = isolateListDf2, permutations = 10000)
permanovaRes2

permanovaDf = data.frame(P=c("<0.0001"), X=c(-250), Y=c(300))

# Plots:
pcoa2 = pcoa(geneDists2, correction="cailliez")
pcoa2Df = as.data.frame(pcoa2$vectors)
pcoa2Df$Isolate = rownames(pcoa2Df)
pcoa2Df$Origin = c(rep("Canada", length(c(canadaList1))), rep("UK", length(ukList)))

pcoa2Plot1 = ggscatter(
  data = pcoa2Df,
  x = "Axis.1",
  y = "Axis.2",
  color = "Origin",
  palette = c("red", "blue"),
  legend.title = "Country of Origin:"
) +
  font("legend.text", size = 12) +
  stat_ellipse(data=pcoa2Df, aes(x=Axis.1, y=Axis.2, color=Origin, fill=Origin), geom = "polygon", alpha=0.2) +
  annotate("text", x=-225, y=375, label="permANOVA\np<0.0001")

pcoa2Plot2 = ggscatter(
  data = pcoa2Df,
  x = "Axis.2",
  y = "Axis.3",
  color = "Origin",
  palette = c("red", "blue")
) +
  stat_ellipse(data=pcoa2Df, aes(x=Axis.2, y=Axis.3, color=Origin, fill=Origin), geom = "polygon", alpha=0.2)

pcoa2Plot3 = ggscatter(
  data = pcoa2Df,
  x = "Axis.3",
  y = "Axis.4",
  color = "Origin",
  palette = c("red", "blue")
) +
  stat_ellipse(data=pcoa2Df, aes(x=Axis.3, y=Axis.4, color=Origin, fill=Origin), geom = "polygon", alpha=0.2)

# pcoaPlot4 = ggscatter(
#   data = pcoa1Df,
#   x = "Axis.4",
#   y = "Axis.5",
#   color = "Origin",
#   palette = c("red", "blue")
# )


screeData2 = data.frame(
  Pos = seq(1, length(pcoa2$values$Cum_corr_eig[1:50])),
  Cumulative_Perc = (pcoa2$values$Cum_corr_eig[1:50]) * 100
)
screeData2$Col = c(rep("lightgrey", 4), rep("darkgrey", length(screeData2$Pos) -
                                             4))
screePlot2 = ggbarplot(
  data = screeData2,
  x = "Pos",
  y = "Cumulative_Perc",
  color = "Col",
  fill = "Col",
  palette = c("lightgrey", "darkgrey"),
  legend = "none"
) +
  scale_y_continuous(breaks = get_breaks(n = 10)) +
  ylab("Cumulative Variance (%)") +
  xlab("Axis Number")


pcoa2Plots = ggarrange(
  pcoa2Plot1,
  pcoa2Plot2,
  pcoa2Plot3,
  screePlot2,
  nrow = 2,
  ncol = 2,
  common.legend = T,
  legend = "top",
  labels = "AUTO"
)

pcoa2Plots

################################################################################
# Write final plots/data
################################################################################
finalGeneDistPlots = ggarrange(plot3, violinPlot2, nrow=1, labels="AUTO", common.legend = T)
# svg & pdf
ggsave("/home/conrad/les_complete/gene-distances-analysis/gene-distances-distributions.svg", finalGeneDistPlots, dpi=300, width=180, height=140, units="mm", device="svg")
ggsave("/home/conrad/les_complete/gene-distances-analysis/gene-distances-distributions.pdf", finalGeneDistPlots, dpi=300, width=180, height=140, units="mm", device="pdf")

# svg & pdf
ggsave("/home/conrad/les_complete/gene-distances-analysis/gene-distances-pcoa.svg", pcoa2Plots, dpi=300, width=180, height=180, units="mm", device="svg")
ggsave("/home/conrad/les_complete/gene-distances-analysis/gene-distances-pcoa.pdf", pcoa2Plots, dpi=300, width=180, height=180, units="mm", device="pdf")


