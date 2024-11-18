library(readr)
library(EnhancedVolcano)
library(ggplot2)

##############################################################################

# Define paths
data_dir <- "<DATA_DIRECTORY_PATH>"
save_path <- "<SAVE_DIRECTORY_PATH>"

##############################################################################

# Load data
AC <- read_csv(file.path(data_dir, "AC-comparison-cell_type-filtered.csv"))
CE <- read_csv(file.path(data_dir, "CE-comparison-cell_type-filtered.csv"))

##############################################################################

# AC - Inhibitory Neuron
sub_data <- AC[AC$cell_type == "Inhibitory neuron", ]
sub_data <- sub_data[sub_data$mean >= 0.1 & sub_data$qval < 0.05 & abs(sub_data$log2fc) >= 0.5, ]

filtered_genes_AC_I <- c('Meis1', 'Glra3', 'Rgs9', 'Rmst', 'Ptprd', 'Hs3st4', 'Unc13c', 'Sv2b', 'Gm48749', 'Dram1', 'Ring1')

AC_I <- EnhancedVolcano(
  sub_data,
  lab = sub_data$gene,
  x = 'log2fc',
  y = 'pval',
  pCutoffCol = 'qval',
  title = "Inhibitory neuron - AC",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 7,
  selectLab = filtered_genes_AC_I,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  legendLabSize = 10,
  max.overlaps = 20
)
ggsave(file.path(save_path, "volcano_AC_inhibitory.pdf"), plot = AC_I, device = "pdf", width = 6, height = 6, dpi = 300)

##############################################################################

# AC - Excitatory Neuron
sub_data <- AC[AC$cell_type == "Excitatory neuron", ]
sub_data <- sub_data[sub_data$mean >= 0.1 & sub_data$qval < 0.05 & abs(sub_data$log2fc) >= 0.5, ]

filtered_genes_AC_E <- c('Foxp2', 'Rorb', 'Prox1', 'Dpp10', 'Brinp3', 'Disc1', 'Sv2c', 'Synpr', 'Rfx3', 'C1ql3', 'Rmst', 'Cobl', 'Gm18870')

AC_E <- EnhancedVolcano(
  sub_data,
  lab = sub_data$gene,
  x = 'log2fc',
  y = 'pval',
  pCutoffCol = 'qval',
  title = "Excitatory neuron - AC",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 7,
  selectLab = filtered_genes_AC_E,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  legendLabSize = 10,
  max.overlaps = 20
)
ggsave(file.path(save_path, "volcano_AC_excitatory.pdf"), plot = AC_E, device = "pdf", width = 6, height = 6, dpi = 300)

##############################################################################

# CE - Inhibitory Neuron
sub_data <- CE[CE$cell_type == "Inhibitory neuron", ]
sub_data <- sub_data[sub_data$mean >= 0.1 & sub_data$qval < 0.05 & abs(sub_data$log2fc) >= 0.5, ]

filtered_genes_CE_I <- c('Rmst', 'Glra3', 'Tshz2', 'Sv2b', 'Pde7b', 'Ptprd', 'Nrg1', 'Penk', 'Cobl', 'Rgs9', 'Unc13c', 'Ttr', 'Prkcd', 'Zic4', 'Enpp2', 'Bcan')

CE_I <- EnhancedVolcano(
  sub_data,
  lab = sub_data$gene,
  x = 'log2fc',
  y = 'pval',
  pCutoffCol = 'qval',
  title = "Inhibitory neuron - CE",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 7,
  selectLab = filtered_genes_CE_I,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  legendLabSize = 10,
  max.overlaps = 20
)
ggsave(file.path(save_path, "volcano_CE_inhibitory.pdf"), plot = CE_I, device = "pdf", width = 6, height = 6, dpi = 300)

##############################################################################

# CE - Excitatory Neuron
sub_data <- CE[CE$cell_type == "Excitatory neuron", ]
sub_data <- sub_data[sub_data$mean >= 0.1 & sub_data$qval < 0.05 & abs(sub_data$log2fc) >= 0.5, ]

filtered_genes_CE_E <- c('Foxp2', 'Rorb', 'Prox1', 'Dpp10', 'Brinp3', 'Disc1', 'Sv2c', 'Synpr', 'Rfx3', 'C1ql3', 'Rmst', 'Cobl', 'Gm18870')

CE_E <- EnhancedVolcano(
  sub_data,
  lab = sub_data$gene,
  x = 'log2fc',
  y = 'pval',
  pCutoffCol = 'qval',
  title = "Excitatory neuron - CE",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 7,
  selectLab = filtered_genes_CE_E,
  drawConnectors = TRUE,
  widthConnectors = 0.75,
  legendLabSize = 10,
  max.overlaps = Inf,
  xlim = c(-2, 2)
)
ggsave(file.path(save_path, "volcano_CE_excitatory.pdf"), plot = CE_E, device = "pdf", width = 6, height = 6, dpi = 300)
