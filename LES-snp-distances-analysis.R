library(tidyverse)
library(ggpubr)
library(patchwork)
library(ape)
library(tidystats)

################################################################################
# Process & format data
################################################################################

# Read in snp distance matrix:
snpData = read.table(
  "/home/conrad/les_complete/snp_calling/L3-snp-calling/cfml/LES.rc_masked.without_outgroups.clean.full.aln.snpdists_matrix.txt",
  sep = "\t",
  header = T,
  stringsAsFactors = F,
  comment.char = "",
  row.names = 1
)

# Set rownames equal to colnames (to keep consistent formatting since colnames are reformatted):
rownames(snpData) = colnames(snpData)
snpDataRownames = rownames(snpData)

# Read in isolate list and create Canada and UK lists:
isolateList = read.table(
  "/home/conrad/les_complete/les-complete-list-XNum.txt",
  header = F,
  stringsAsFactors = F
)
isolateList = isolateList$V1
isolateList = c(isolateList, "Reference")

canadaList = isolateList[1:50] # all Canada isolates together
ukList = isolateList[51:206] # all UK isolates
canadaList2 = isolateList[1:11] # deep branching clade Canada isolates
canadaList1 = isolateList[12:50] # non-deep branching Canada isolates

# Read in hypermutators list:
hypermutators = read.table(
  "/home/conrad/les_complete/les-phylogenetic-hypermutators-XNum.txt",
  header = F,
  stringsAsFactors = F,
  comment.char = ""
)
hypermutators = hypermutators$V1
hypermutators = c(hypermutators, "PHELES19") # remove PHELES19 b/c it's a complete outlier from rest of UK

canadaList = canadaList[!(canadaList %in% hypermutators)]
ukList = ukList[!(ukList %in% hypermutators)]
canadaList2 = canadaList2[!(canadaList2 %in% hypermutators)]
canadaList1 = canadaList1[!(canadaList1 %in% hypermutators)]

snpData = snpData[!(rownames(snpData) %in% hypermutators), ]
snpData = snpData[, !(colnames(snpData) %in% hypermutators)]

isolateList = isolateList[!(isolateList %in% hypermutators)]

# Order SNP data by isolate list (Canada on top):
snpData = snpData[isolateList, ]
snpData = snpData[isolateList]

# Convert matrix to long format:
outputDists = c()
namesDists = c()
colsDists = c()

start_row = 1
start_col = 1
for (i in seq(1, 206, 1)) {
  new_vec = snpData[, start_col]
  new_vec = new_vec[-(1:start_row)]
  currentCol = colnames(snpData)[start_col]
  currentRows = rownames(snpData)[-(1:start_row)]
  currentCols = rep(currentCol, length(currentRows))
  outputDists = c(outputDists, new_vec)
  namesDists = c(namesDists, currentRows)
  colsDists = c(colsDists, currentCols)
  if (start_row + 1 <= nrow(snpData) &
      start_col + 1 <= ncol(snpData)) {
    start_row = start_row + 1
    start_col = start_col + 1
  }
}

# Generate output df:
outputDf = data.frame(Isolate1 = colsDists,
                      Isolate2 = namesDists,
                      SNP_Distance = outputDists)

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

# write.table(outputDf, "/home/conrad/les_complete/snp_calling/L3-snp-calling/cfml/LES.melted.rc_masked.matrix.txt",
#             quote=F, sep="\t")

################################################################################
# Calculate summary statistics per group
################################################################################

# Calculate statistics for SNP distances for CanadaVsUk comparison:
canadaVsUkGrouped = outputDf %>% group_by(CanadaVsUK) %>% summarise(
  MeanDist = mean(SNP_Distance),
  MedianDist = median(SNP_Distance),
  MinDist = min(SNP_Distance),
  MaxDist = max(SNP_Distance),
  Q1 = quantile(SNP_Distance, 0.25),
  Q3 = quantile(SNP_Distance, 0.75)
)
canadaVsUkGrouped

# Calculate statistics for SNP distances for CanadaSplitVsUK
canadaSplitVsUkGrouped = outputDf %>% group_by(CanadaVsUK_2) %>% summarise(
  MeanDist = mean(SNP_Distance),
  MedianDist = median(SNP_Distance),
  MinDist = min(SNP_Distance),
  MaxDist = max(SNP_Distance),
  Q1 = quantile(SNP_Distance, 0.25),
  Q3 = quantile(SNP_Distance, 0.75)
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

## This first plot is just to get the legend...
plot2a = ggdensity(
  data = outputDf,
  x = "SNP_Distance",
  add = "mean",
  rug = TRUE,
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
  legend = "top"
) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 13))
plot2aLegend = get_legend(plot2a)

# And this plot is used for actual plotting:
plot2b = ggdensity(
  data = outputDf,
  x = "SNP_Distance",
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
  legend = "none"
) +
  font("xy.text", size = 10) +
  xlab("") +
  ylab("")

# This plot includes only intra-group comparisons + UK-Canada1:
plot3 = ggdensity(
  data = outputDf[outputDf$CanadaVsUK_2 %in% c("Canada1-Canada1", "Canada2-Canada2", "UK-UK", "Canada1-UK"), ],
  x = "SNP_Distance",
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
  legend = "none"
) +
  xlab("SNP Distance") +
  ylab("Density") +
  scale_x_continuous(breaks = get_breaks(n = 12))

# Put plots together:
finalPlot = (as_ggplot(plot2aLegend)) / (plot3 + inset_element(
  plot2b,
  left = 0.45,
  bottom = 0.41,
  right = 0.95,
  top = 0.9
)) +
  plot_layout(heights = c(0.5, 5))
finalPlot

#ggsave("/home/conrad/test-dist-plot.pdf", plot=finalPlot, device="pdf", height=180, width=180, units="mm", dpi=300)

################################################################################
# Make violin plots
################################################################################

# Violin plots:

# Plot for Canada split vs UK (not much point in doing all Canada vs UK):
violinPlot1a = ggviolin(
  data = outputDf,
  x = "CanadaVsUK_2",
  y = "SNP_Distance",
  fill = "CanadaVsUK_2",
  palette = c(
    "#00AFBB",
    "#E7B800",
    "#f07dc2",
    "#80f07d",
    "#877df0",
    "#FC4E07"
  ),
  legend = "top"
) +
  scale_x_discrete(
    limits = c(
      "Canada2-Canada2",
      "Canada1-Canada1",
      "UK-UK",
      "Canada1-UK",
      "Canada2-Canada1",
      "Canada2-UK"
    )
  ) +
  theme(
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  xlab("") +
  ylab("SNP Distance") +
  ylim(0, 5000)
violinPlotLegend = get_legend(violinPlot1a)

violinPlot1b = ggviolin(
  data = outputDf,
  x = "CanadaVsUK_2",
  y = "SNP_Distance",
  fill = "CanadaVsUK_2",
  palette = c(
    "#00AFBB",
    "#E7B800",
    "#f07dc2",
    "#80f07d",
    "#877df0",
    "#FC4E07"
  ),
  legend = "none"
) +
  scale_x_discrete(
    limits = c(
      "Canada1-Canada1",
      "Canada1-UK",
      "UK-UK",
      "Canada2-Canada2",
      "Canada2-Canada1",
      "Canada2-UK"
    )
  ) +
  theme(
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  xlab("") +
  ylab("") +
  ylim(0, 5000)

violinPlot2 = ggviolin(
  data = outputDf[outputDf$CanadaVsUK_2 %in% c("Canada1-Canada1", "Canada2-Canada2", "UK-UK", "Canada1-UK"), ],
  x = "CanadaVsUK_2",
  y = "SNP_Distance",
  fill = "CanadaVsUK_2",
  palette = c(
    "#00AFBB",
    "#E7B800",
    "#f07dc2",
    "#80f07d",
    "#877df0",
    "#FC4E07"
  ),
  legend = "none",
  add = "boxplot"
) +
  scale_x_discrete(limits = c("Canada1-Canada1", "Canada1-UK", "UK-UK", "Canada2-Canada2")) +
  theme(
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  xlab("") +
  ylab("SNP Distance") +
  ylim(0, 500)
violinPlot2 = violinPlot2 + stat_compare_means(
  comparisons = list(c("Canada1-UK", "Canada1-Canada1")),
  method.args = list(alternative =
                       "two.sided"),
  label.y = c(475),
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
    symbols = c("****", "***", "**", "*", "ns")
  )
)
finalViolinPlot = (as_ggplot(violinPlotLegend)) / (violinPlot2 + inset_element(
  violinPlot1b,
  left = 0.5,
  bottom = 0.68,
  right = 0.99,
  top = 1
)) +
  plot_layout(heights = c(0.5, 5))
finalViolinPlot

#ggsave("/home/conrad/test-violing-plot.pdf", plot=finalViolinPlot, width = 180, height = 180, units = "mm", device = "pdf", dpi = 300)

################################################################################
# Run statistical tests
################################################################################

# Wilcoxon rank sum test: are UK-Canada1 distances SMALLER than Canada1-Canada1?
# i.e., are UK isolates as distant to Canadian isolates as Canadian isolates are to one another?
wilcoxTest = wilcox.test(outputDf$SNP_Distance[outputDf$CanadaVsUK_2 == "Canada1-UK"],
                         outputDf$SNP_Distance[outputDf$CanadaVsUK_2 == "Canada1-Canada1"],
                         alternative = "two.sided")
wilcoxTest

ksTest = ks.test(outputDf$SNP_Distance[outputDf$CanadaVsUK_2 == "Canada1-UK"],
                 outputDf$SNP_Distance[outputDf$CanadaVsUK_2 == "Canada1-Canada1"],
                 alternative = "two.sided")
ksTest

kwTest = kruskal.test(list(outputDf$SNP_Distance[outputDf$CanadaVsUK_2 == "Canada1-UK"],
                           outputDf$SNP_Distance[outputDf$CanadaVsUK_2 == "Canada1-Canada1"]))
kwTest

################################################################################
################################################################################
################################################################################
# Repeat above but for p-distances
################################################################################
################################################################################
################################################################################

fasta = read.FASTA(
  "/home/conrad/les_complete/snp_calling/L3-snp-calling/cfml/LES.rc_masked.without_outgroups.clean.full.aln"
)
p_distances = dist.dna(fasta, model = "raw", pairwise.deletion = T)
pDistances = as.matrix(p_distances)
pDistances = as.data.frame(pDistances)
rownames(pDistances) = snpDataRownames
colnames(pDistances) = snpDataRownames

# remove fasta to save memory:
rm(fasta)

pDistances = pDistances[!(rownames(pDistances) %in% hypermutators), ]
pDistances = pDistances[, !(colnames(pDistances) %in% hypermutators)]

# Order SNP data by isolate list (Canada on top):
pDistances = pDistances[isolateList, ]
pDistances = pDistances[isolateList]

# Convert matrix to long format:
outputPDists = c()
namesPDists = c()
colsPDists = c()

start_row = 1
start_col = 1
for (i in seq(1, 206, 1)) {
  new_vec = pDistances[, start_col]
  new_vec = new_vec[-(1:start_row)]
  currentCol = colnames(pDistances)[start_col]
  currentRows = rownames(pDistances)[-(1:start_row)]
  currentCols = rep(currentCol, length(currentRows))
  outputPDists = c(outputPDists, new_vec)
  namesPDists = c(namesPDists, currentRows)
  colsPDists = c(colsPDists, currentCols)
  if (start_row + 1 <= nrow(pDistances) &
      start_col + 1 <= ncol(pDistances)) {
    start_row = start_row + 1
    start_col = start_col + 1
  }
}

# Generate output df:
pDistDf = data.frame(Isolate1 = colsPDists,
                     Isolate2 = namesPDists,
                     SNP_Distance = outputPDists)

# Label each comparison as intra-Canada, intra-UK, or inter-Canada/UK:
canadaUkCompP = c()
for (i in seq(1, length(rownames(pDistDf)), 1)) {
  isolate1 = pDistDf$Isolate1[i]
  isolate2 = pDistDf$Isolate2[i]
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
  canadaUkCompP = c(canadaUkCompP,
                    paste(isolate1Origin, "-", isolate2Origin, sep = ""))
}

pDistDf$CanadaVsUKP = canadaUkCompP

# Label each comparison as above but split Canada into deep branching and non-deep branching:
canadaUkComp2P = c()

for (i in seq(1, length(rownames(pDistDf)), 1)) {
  isolate1 = pDistDf$Isolate1[i]
  isolate2 = pDistDf$Isolate2[i]
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
  canadaUkComp2P = c(canadaUkComp2P,
                     paste(isolate1Origin, "-", isolate2Origin, sep = ""))
}

pDistDf$CanadaVsUK_2 = canadaUkComp2P

################################################################################
# Calculate summary statistics per group
################################################################################

# Calculate statistics for SNP distances for CanadaVsUk comparison:
canadaVsUkGroupedP = pDistDf %>% group_by(CanadaVsUKP) %>% summarise(
  MeanDist = mean(SNP_Distance),
  MedianDist = median(SNP_Distance),
  MinDist = min(SNP_Distance),
  MaxDist = max(SNP_Distance),
  Q1 = quantile(SNP_Distance, 0.25),
  Q3 = quantile(SNP_Distance, 0.75)
)
canadaVsUkGroupedP

# Calculate statistics for SNP distances for CanadaSplitVsUK
canadaSplitVsUkGroupedP = pDistDf %>% group_by(CanadaVsUK_2) %>% summarise(
  MeanDist = mean(SNP_Distance),
  MedianDist = median(SNP_Distance),
  MinDist = min(SNP_Distance),
  MaxDist = max(SNP_Distance),
  Q1 = quantile(SNP_Distance, 0.25),
  Q3 = quantile(SNP_Distance, 0.75)
)
canadaSplitVsUkGroupedP

################################################################################
# Plot histograms/densities
################################################################################

## This first plot is just to get the legend...
pDistDf$CanadaVsUK_2 = factor(
  pDistDf$CanadaVsUK_2,
  levels = c(
    "Canada1-Canada1",
    "Canada1-UK",
    "UK-UK",
    "Canada2-Canada2",
    "Canada2-Canada1",
    "Canada2-UK"
  )
)

plot2a = ggdensity(
  data = pDistDf,
  x = "SNP_Distance",
  add = "mean",
  rug = TRUE,
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
  legend = "top"
) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 13))
plot2aLegend = get_legend(plot2a)

# And this plot is used for actual plotting:
plot2b = ggdensity(
  data = pDistDf,
  x = "SNP_Distance",
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
  legend = "none"
) +
  font("xy.text", size = 10) +
  xlab("") +
  ylab("") +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    )
  )

# This plot includes only intra-group comparisons + UK-Canada1:
plot3 = ggdensity(
  data = pDistDf[pDistDf$CanadaVsUK_2 %in% c("Canada1-Canada1", "Canada2-Canada2", "UK-UK", "Canada1-UK"), ],
  x = "SNP_Distance",
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
  legend = "none"
) +
  xlab("SNP Distance") +
  ylab("Density") +
  scale_x_continuous(breaks = get_breaks(n = 12))

# Put plots together:
finalPlotP = (as_ggplot(plot2aLegend)) / (plot3 + inset_element(
  plot2b,
  left = 0.45,
  bottom = 0.41,
  right = 0.95,
  top = 0.9
)) +
  plot_layout(heights = c(0.5, 5))
finalPlotP

#ggsave()

################################################################################
# Make violin plots
################################################################################

# Violin plots:

# Plot for Canada split vs UK (not much point in doing all Canada vs UK):
violinPlot1a = ggviolin(
  data = pDistDf,
  x = "CanadaVsUK_2",
  y = "SNP_Distance",
  fill = "CanadaVsUK_2",
  palette = c(
    "#00AFBB",
    "#E7B800",
    "#f07dc2",
    "#80f07d",
    "#877df0",
    "#FC4E07"
  ),
  legend = "top"
) +
  scale_x_discrete(
    limits = c(
      "Canada2-Canada2",
      "Canada1-Canada1",
      "UK-UK",
      "Canada1-UK",
      "Canada2-Canada1",
      "Canada2-UK"
    )
  ) +
  theme(
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  xlab("") +
  ylab("SNP Distance")
violinPlotLegend = get_legend(violinPlot1a)

violinPlot1b = ggviolin(
  data = pDistDf,
  x = "CanadaVsUK_2",
  y = "SNP_Distance",
  fill = "CanadaVsUK_2",
  palette = c(
    "#00AFBB",
    "#E7B800",
    "#f07dc2",
    "#80f07d",
    "#877df0",
    "#FC4E07"
  ),
  legend = "none"
) +
  scale_x_discrete(
    limits = c(
      "Canada1-Canada1",
      "Canada1-UK",
      "UK-UK",
      "Canada2-Canada2",
      "Canada2-Canada1",
      "Canada2-UK"
    )
  ) +
  theme(
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  ) +
  xlab("") +
  ylab("")

violinPlot2 = ggviolin(
  data = pDistDf[pDistDf$CanadaVsUK_2 %in% c("Canada1-Canada1", "Canada2-Canada2", "UK-UK", "Canada1-UK"), ],
  x = "CanadaVsUK_2",
  y = "SNP_Distance",
  fill = "CanadaVsUK_2",
  palette = c(
    "#00AFBB",
    "#E7B800",
    "#f07dc2",
    "#80f07d",
    "#877df0",
    "#FC4E07"
  ),
  legend = "none",
  add = "boxplot"
) +
  scale_x_discrete(limits = c("Canada1-Canada1", "Canada1-UK", "UK-UK", "Canada2-Canada2")) +
  theme(
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  xlab("") +
  ylab("P Distance") +
  ylim(c(0, 9e-05))
violinPlot2 = violinPlot2 + stat_compare_means(
  comparisons = list(c("Canada1-UK", "Canada1-Canada1")),
  method.args = list(alternative =
                       "greater"),
  label.y = c(7.55e-5),
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
    symbols = c("****", "***", "**", "*", "ns")
  )
)
finalViolinPlotP = (as_ggplot(violinPlotLegend)) / (violinPlot2 + inset_element(
  violinPlot1b,
  left = 0.62,
  bottom = 0.65,
  right = 0.99,
  top = 1
)) +
  plot_layout(heights = c(0.5, 5))
finalViolinPlotP

#ggsave()

################################################################################
# Run statistical tests
################################################################################

# Wilcoxon rank sum test: are UK-Canada1 distances SMALLER than Canada1-Canada1?
# i.e., are UK isolates as distant to Canadian isolates as Canadian isolates are to one another?
wilcoxTestP = wilcox.test(pDistDf$SNP_Distance[pDistDf$CanadaVsUK_2 == "Canada1-UK"],
                          pDistDf$SNP_Distance[pDistDf$CanadaVsUK_2 == "Canada1-Canada1"],
                          alternative = "greater")
wilcoxTestP

ksTestP = ks.test(pDistDf$SNP_Distance[pDistDf$CanadaVsUK_2 == "Canada1-UK"],
                  pDistDf$SNP_Distance[pDistDf$CanadaVsUK_2 == "Canada1-Canada1"],
                  alternative = "less")
ksTestP

kwTestP = kruskal.test(list(pDistDf$SNP_Distance[pDistDf$CanadaVsUK_2 == "Canada1-UK"],
                            pDistDf$SNP_Distance[pDistDf$CanadaVsUK_2 == "Canada1-Canada1"]))
kwTestP

################################################################################
# PCoA on SNPs (fails on p-distances)
################################################################################

snpData2 = snpData[-(1:11), -(1:11)]

pcoa1 = pcoa(snpData2, correction = "cailliez")

pcoa1Df = as.data.frame(pcoa1$vectors)
pcoa1Df$Isolate = rownames(pcoa1Df)
pcoa1Df$Origin = c(rep("Canada", length(canadaList1)), rep("UK", length(ukList)))

pcoaPlot1 = ggscatter(
  data = pcoa1Df,
  x = "Axis.1",
  y = "Axis.2",
  color = "Origin",
  palette = c("red", "blue"),
  legend.title = "Country of Origin:"
) +
  font("legend.text", size = 12)

pcoaPlot2 = ggscatter(
  data = pcoa1Df,
  x = "Axis.2",
  y = "Axis.3",
  color = "Origin",
  palette = c("red", "blue")
)

pcoaPlot3 = ggscatter(
  data = pcoa1Df,
  x = "Axis.3",
  y = "Axis.4",
  color = "Origin",
  palette = c("red", "blue")
)

pcoaPlot4 = ggscatter(
  data = pcoa1Df,
  x = "Axis.4",
  y = "Axis.5",
  color = "Origin",
  palette = c("red", "blue")
)


screeData = data.frame(
  Pos = seq(1, length(pcoa1$values$Cum_corr_eig[1:50])),
  Cumulative_Perc = (pcoa1$values$Cum_corr_eig[1:50]) * 100
)
screeData$Col = c(rep("lightgrey", 4), rep("darkgrey", length(screeData$Pos) -
                                             4))
screePlot = ggbarplot(
  data = screeData,
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


pcoaPlots = ggarrange(
  pcoaPlot1,
  pcoaPlot2,
  pcoaPlot3,
  screePlot,
  nrow = 2,
  ncol = 2,
  common.legend = T,
  legend = "top",
  labels = "AUTO"
)
#pcoaPlots = pcoaPlots + inset_element(screePlot, left=0.77, bottom=0.51, right=0.99, top=0.78)
pcoaPlots

#ggsave("/home/conrad/test-pcoa.pdf", device="pdf", dpi=300, height=180, width=180, units="mm")



################################################################################
# Final plots
################################################################################
finalSNPsPlots = ggarrange(finalPlot, finalViolinPlot, nrow=1, labels="AUTO")
finalPdistPlots = ggarrange(finalPlotP, finalViolinPlotP, nrow=1, labels="AUTO")

ggsave("/home/conrad/les_complete/snp-distances-analysis/final-plots-SNPs.pdf", finalSNPsPlots, device="pdf", height=180, width=180, units="mm", dpi=300)
ggsave("/home/conrad/les_complete/snp-distances-analysis/final-plots-Pdists.pdf", finalPdistPlots, device="pdf", height=180, width=180, units="mm", dpi=300)

outStats = list()
outStats = add_stats(outStats, wilcoxTest)
outStats = add_stats(outStats, wilcoxTestP)
write_stats(outStats, "/home/conrad/les_complete/snp-distances-analysis/wilcox-stats.json")
