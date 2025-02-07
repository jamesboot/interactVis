# Script for implementing, testing and developing interactVis functions

#### PART 1: Load packages, data and preprocess ----

setwd('/Users/bootj/Documents/Milner_et_al/')

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(multcomp)
library(purrr)
library(tibble)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(parallel)
library(pbapply)

# Load functions
source('interactVis-main/findNeighbours.R')
source('interactVis-main/loadDB.R')
source('interactVis-main/findInteractions.R')
source('interactVis-main/annotateInteractions.R')
source('interactVis-main/interactionMatrix.R')
source('interactVis-main/interactionMeta.R')
source('interactVis-main/differentialInteraction.R')

# Load database
database <- loadDB(interaction.input = 'interactVis-main/cellphonedb-data-master/data/interaction_input.csv')

# Start analysis on all data from project GC-TM-10271

# Set base.dir
base.dir <- 'interactVis_new_outs_v2'
if (!dir.exists(base.dir)) {
  dir.create(base.dir)
}

# Create list of all sample directories
dirs <- list.dirs(path = 'spatial',
                  recursive = F)

# Create list of all sample .h5 files
files <- list.files(path = 'spatial',
                    pattern = '*.h5',
                    recursive = T,
                    full.names = T)

# Create list of sample names
samples <- gsub(
  '_spatial',
  '',
  gsub(
    'spatial/',
    '',
    dirs
  )
)

#### PART 2: Loop through all samples and perform interaction analysis ----

# For loop to go through all samples and perform all analysis
for (ITER in 13:16) {
  
  # Load spatial data
  dat <- Load10X_Spatial(
    dirs[ITER],
    filename = basename(files[ITER]),
    assay = 'Spatial',
    filter.matrix = TRUE,
    to.upper = FALSE
  )
  
  # Set the orig.ident to the sample name
  dat$orig.ident <- samples[ITER]
  
  # Filter to only spots used by Tom in Partek
  # Import the selected spots object for filtering to Tom's spots
  selectedSpots <-
    readRDS(
      'spatial/selectedSpots.RDS'
    )
  
  # Now filter each Seurat object in list for selected tissue spots
  dat <-
    subset(dat, cells = selectedSpots$`Cell name`[selectedSpots$`Sample name` == unique(dat$orig.ident)])
  
  # Normalise data
  dat <-
    NormalizeData(dat, normalization.method = "LogNormalize", verbose = FALSE)
  
  # Find neighbors
  neighbours <-
    findNeighbours(file.path(dirs[ITER], 'spatial/tissue_positions.csv'))
  
  # Find interactions
  interactions <- calculateInteractions(
    neighbours = neighbours,
    dat = dat,
    database = database,
    filter = FALSE
  )
  
  # Prepare a folder to save everything in
  dir.create(path = file.path(base.dir, samples[ITER]))
  
  # Save interaction object
  saveRDS(interactions,
          file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_interactions.rds'))
  
  # Save complex summary as csv
  write.csv(
    interactions$Complex_Summary,
    file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_complexSummary.csv')
  )
  
  # Add cluster numbers from Partek
  partekMeta <-
    read.delim('exported.txt')
  partekMetaFilt <-
    partekMeta[partekMeta$Sample.name == samples[ITER],]
  clustersMeta <- data.frame(
    row.names = partekMetaFilt$Cell.name,
    Cluster = as.factor(partekMetaFilt$Clusters..5comp.0.75res)
  )
  dat <- AddMetaData(dat, clustersMeta)
  
  # Visualise the clusters for this sample
  plt <- SpatialPlot(dat, group.by = 'Cluster', label.size = 10) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
  ggsave(
    plot = plt,
    filename = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_cluster.pdf')
  )
  
  # Annotate the Interactions results with the cluster for spot 1 and spot 2
  # Spot 1 is the receiver (where the receptor expression is taken from)
  # Spot 2 is the sender (where the ligand expression is taken from)
  diffIntDF <- annotateInteractions(SeuratObj = dat,
                                    Interactions = interactions,
                                    Attribute = 'Cluster')
  
  # How many of each interactions per cluster are there?
  # And between clusters
  # Find number of interactions between groups
  interactionSum <- diffIntDF %>%
    group_by(Sender_Anno, Receiver_Anno) %>%
    summarise(n = n())
  write.csv(interactionSum,
            file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_nInt_cluster.csv'))
  
  # Convert to matrix
  # Create empty matrix to population
  interactionSumMat <- matrix(
    data = NA,
    nrow = length(unique(interactionSum$Sender_Anno)),
    ncol = length(unique(interactionSum$Receiver_Anno)),
    dimnames = list(
      unique(interactionSum$Sender_Anno),
      unique(interactionSum$Receiver_Anno)
    )
  )
  
  # Populate matrix in for loop
  for (row in unique(interactionSum$Sender_Anno)) {
    for (col in unique(interactionSum$Receiver_Anno)) {
      # Find the n number
      n <- interactionSum %>%
        filter(Sender_Anno == row & Receiver_Anno == col) %>%
        pull(n)
      
      # If length of n == 1 add to matrix
      if (length(n) == 1) {
        # Populate with value in n
        interactionSumMat[row, col] <- n
        
      } else if (length(n) == 0) {
        # Populate with 0
        interactionSumMat[row, col] <- 0
        
      } else {
        message('Error, length of n is neither 1 or 0.')
        
      }
    }
  }
  
  # Set colours for heatmap
  colfun <- colorRampPalette(brewer.pal(8, "Blues"))(25)
  
  # Set histogram values for heatmap
  column_bp <-
    HeatmapAnnotation(Total = anno_barplot(colSums(interactionSumMat),
                                           gp = gpar(fill = '#000000')))
  row_bp <-
    rowAnnotation(Total = anno_barplot(rowSums(interactionSumMat),
                                       gp = gpar(fill = '#000000')))
  
  # Plot heatmap
  pdf(
    file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_htmap.pdf'),
    width = 5,
    height = 5
  )
  draw(Heatmap(
    interactionSumMat,
    col = colfun,
    top_annotation = column_bp,
    left_annotation = row_bp,
    name = 'n_intrtns'
  ))
  dev.off()
  
}

#### PART 3: Re-load RDS objects for each sample and summarise with Chord plots ----

# Run PART 1 before running this section
# This can be run without re-running PART 2

# Locate CSV files
nInt_files <- list.files(
  path = 'interactVis_new_outs_v2/',
  recursive = T,
  pattern = '_nInt_cluster.csv',
  full.names = T
)

# Load CSV files 
nInt_objs <- lapply(nInt_files, function(x){
  read.csv(x)[, 2:4]
})

# Names
names(nInt_objs) <- samples

# Ammend names of clusters to allow plotting
for (x in names(nInt_objs)) {
  nInt_objs[[x]][["Sender_Anno"]] <- paste0('C',nInt_objs[[x]][["Sender_Anno"]])
  nInt_objs[[x]][["Receiver_Anno"]] <- paste0('C',nInt_objs[[x]][["Receiver_Anno"]])
}

# Set grid colours
grid.col = c(
  `C1` = "#3366cc",
  `C2` = '#dc3911',
  `C3` = '#ff9900',
  `C4` = '#0d9618',
  `C5` = '#990099',
  `C6` = '#0099c5',
  `C7` = '#dd4477',
  `C8` = '#66a900',
  `C9` = '#b72e2f',
  `C10` = '#6633cc',
  `C11` = '#22a999',
  `C12` = '#306395',
  `C13` = '#aaaa11',
  `C14` = '#984499',
  `C15` = '#e67301',
  `C16` = '#8b0607',
  `C17` = '#339262',
  `C18` = '#3a3eac',
  `C19` = '#651066'
)

# Chord plot for each sample
options(scipen=1000000)
for (ITER in 1:length(nInt_objs)) {
  pdf(
    file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_chord.pdf'),
    height = 5,
    width = 5
  )
  circos.clear()
  circos.par(gap.after = 5)
  chordDiagram(
    nInt_objs[[ITER]],
    transparency = 0.5,
    grid.col = grid.col,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = mm_h(5))
  )
  for(si in get.all.sector.index()) {
    circos.axis(h = "top", labels.cex = 0.3, sector.index = si, track.index = 2, labels.facing = 'clockwise')
  }
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    xplot = get.cell.meta.data("xplot")
    circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, -1.5), cex = 0.6)
  }, bg.border = NA)
  circos.clear()
  dev.off()
}

# Plot all samples together
options(scipen=1000000)
combined_df <- do.call(rbind, nInt_objs)
pdf(
  file = paste0(base.dir, '/allSamples_chord.pdf'),
  height = 6,
  width = 6
)
circos.clear()
circos.par(gap.after = 5)
chordDiagram(combined_df, transparency = 0.5, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = mm_h(5)))
for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.3, sector.index = si, track.index = 2, labels.facing = 'clockwise')
}
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  xplot = get.cell.meta.data("xplot")
  circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, -1.5), cex = 0.6)
}, bg.border = NA)
circos.clear()
dev.off()

# Plot Exp samples together
sampleGroups <- read.csv('Tom_groups.csv')
expSamples <- sampleGroups$Sample[sampleGroups$Group == 'Exp']
combined_df <- do.call(rbind, nInt_objs[expSamples])
options(scipen=1000000)
pdf(
  file = paste0(base.dir, '/expSamples_chord.pdf'),
  height = 5,
  width = 5
)
circos.clear()
circos.par(gap.after = 5)
chordDiagram(combined_df, transparency = 0.5, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = mm_h(5)))
for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.3, sector.index = si, track.index = 2, labels.facing = 'clockwise')
}
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  xplot = get.cell.meta.data("xplot")
  circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, -1.5), cex = 0.6)
}, bg.border = NA)
circos.clear()
dev.off()

# Plot Ctrl samples together
sampleGroups <- read.csv('Tom_groups.csv')
expSamples <- sampleGroups$Sample[sampleGroups$Group == 'Cnt']
combined_df <- do.call(rbind, nInt_objs[expSamples])
options(scipen=1000000)
pdf(
  file = paste0(base.dir, '/ctrlSamples_chord.pdf'),
  height = 5,
  width = 5
)
circos.clear()
circos.par(gap.after = 5)
chordDiagram(combined_df, transparency = 0.5, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = mm_h(5)))
for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.3, sector.index = si, track.index = 2, labels.facing = 'clockwise')
}
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  xplot = get.cell.meta.data("xplot")
  circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, -1.5), cex = 0.6)
}, bg.border = NA)
circos.clear()
dev.off()
  
#### PART 4: Differential clusters analysis ----

# Run PART 1 before running this section
# This can be run without re-running PART 2 or 3

# For loop to load and process all sample Seurat objects into list
Seurat_objs <- list()
for (ITER in 1:c(length(samples))) {
  
  # Load spatial data
  Seurat_objs[[samples[ITER]]] <- Load10X_Spatial(
    dirs[ITER],
    filename = basename(files[ITER]),
    assay = 'Spatial',
    filter.matrix = TRUE,
    to.upper = FALSE
  )
  
  # Set the orig.ident to the sample name
  Seurat_objs[[samples[ITER]]]$orig.ident <- samples[ITER]
  
  # Filter to only spots used by Tom in Partek
  # Import the selected spots object for filtering to Tom's spots
  selectedSpots <-
    readRDS(
      'spatial/selectedSpots.RDS'
    )
  
  # Now filter each Seurat object in list for selected tissue spots
  Seurat_objs[[samples[ITER]]] <-
    subset(Seurat_objs[[samples[ITER]]], cells = selectedSpots$`Cell name`[selectedSpots$`Sample name` == unique(Seurat_objs[[samples[ITER]]]$orig.ident)])
  
  # Normalise data
  Seurat_objs[[samples[ITER]]] <-
    NormalizeData(Seurat_objs[[samples[ITER]]], normalization.method = "LogNormalize", verbose = FALSE)
  
  # Add cluster numbers from Partek
  partekMeta <-
    read.delim('exported.txt')
  partekMetaFilt <-
    partekMeta[partekMeta$Sample.name == samples[ITER],]
  clustersMeta <- data.frame(
    row.names = partekMetaFilt$Cell.name,
    Cluster = as.factor(partekMetaFilt$Clusters..5comp.0.75res)
  )
  Seurat_objs[[samples[ITER]]] <- AddMetaData(Seurat_objs[[samples[ITER]]], clustersMeta)

}

# Locate RDS files of interaction results
rds_files <- list.files(path = 'interactVis_new_outs_v2',
                         recursive = T,
                         pattern = '.rds',
                         full.names = T)

# Load RDS files 
rds_objs <- lapply(rds_files, function(x){
  readRDS(x)
})

# Names
names(rds_objs) <- samples

# Annotate the Interactions results with the cluster for spot 1 and spot 2
# Spot 1 is the receiver (where the receptor expression is taken from)
# Spot 2 is the sender (where the ligand expression is taken from)
AnnoInt <- list()
for (ITER in 1:c(length(samples))) {
  message(paste('Starting sample:', samples[ITER]))
  AnnoInt[[samples[ITER]]] <-
    annotateInteractions(
      SeuratObj = Seurat_objs[[samples[ITER]]],
      Interactions = rds_objs[[samples[ITER]]],
      Attribute = 'Cluster'
    )
}
# Save AnnoInt
saveRDS(AnnoInt, file = paste0(base.dir, '/AnnoInt.RDS'))

# # Loop to go through all samples and perform differential interaction
# AllDiffIntRes <- list()
# for (ITER in 1:c(length(samples))) {
#   
#   # Create Interaction Matrix
#   # Function will create for both SENDERS and RECEIVERS
#   IntMat <- interactionMatrix(AnnoInt = AnnoInt[[samples[ITER]]])
#   
#   # Create differential interaction meta data from Partek meta data
#   # Function will create for both SENDERS and RECEIVERS
#   diffIntMeta <- createMetaData(
#     SeuratObj = Seurat_objs[[samples[ITER]]],
#     InteractionMatList = IntMat,
#     Attributes = c('orig.ident', 'Cluster')
#   )
#   
#   # Perform differential interaction analysis
#   # Function will create for both SENDERS and RECEIVERS
#   # Function will also check there are enough replicates 
#   diffInt2v10 <- differentialInteraction(
#     InteractionMatList = IntMat,
#     MetaList = diffIntMeta,
#     Attribute = 'Cluster',
#     Comparison = c(2, 10)
#   )
#   
#   # Write RECEIVER results to file
#   write.csv(diffInt2v10$ReceiverResults,
#             file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_RECEIVER_wilcox.csv'))
#   
#   # Write SENDER results to file
#   write.csv(diffInt2v10$SenderResults,
#             file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_SENDER_wilcox.csv'))
#   
#   # Write result for sample to list
#   AllDiffIntRes[[samples[ITER]]] <- diffInt2v10
#   
# }

# Perform analysis on all samples merged together
# Append sample names to barcodes before merging
for (ITER in 1:c(length(AnnoInt))) {
  AnnoInt[[ITER]]$Sender_bcode <- paste0(AnnoInt[[ITER]]$Sender_bcode,
                                    '-',
                                    names(AnnoInt)[ITER])
  
  AnnoInt[[ITER]]$Receiver_bcode <- paste0(AnnoInt[[ITER]]$Receiver_bcode,
                                    '-',
                                    names(AnnoInt)[ITER])
}

# Add sample names to DF before merging
for (ITER in 1:c(length(AnnoInt))) {
  AnnoInt[[ITER]]$sample <- names(AnnoInt)[ITER]
}

# Bind all rows of all dataframes together
AnnoInt_Comb <- do.call(rbind, AnnoInt)
saveRDS(AnnoInt_Comb, file = paste0(base.dir, '/allSamplesAnnoInt.rds'))
AnnoInt_Comb <- readRDS('interactVis_new_outs_v2/allSamplesAnnoInt.rds')
AnnoInt_Comb <- as_tibble(AnnoInt_Comb)
  
# Create Interaction Matrix
# Function will create for both SENDERS and RECEIVERS
AllIntMat <- interactionMatrix(AnnoInt = AnnoInt_Comb)
saveRDS(AllIntMat, file = paste0(base.dir, '/allSamplesIntMat.rds'))

# Create meta data

# Change cell names
for (ITER in 1:length(samples)) {
  Seurat_objs[[ITER]] <- RenameCells(Seurat_objs[[ITER]],
                                     new.names = paste0(Cells(Seurat_objs[[ITER]]),
                                                        '-',
                                                        samples[ITER]))
}

# Make one large seurat object for meta data creation
AllSampSeuratObj <- merge(Seurat_objs[[1]],
                          y = c(
                            Seurat_objs[[2]],
                            Seurat_objs[[3]],
                            Seurat_objs[[4]],
                            Seurat_objs[[5]],
                            Seurat_objs[[6]],
                            Seurat_objs[[7]],
                            Seurat_objs[[8]],
                            Seurat_objs[[9]],
                            Seurat_objs[[10]],
                            Seurat_objs[[11]],
                            Seurat_objs[[12]],
                            Seurat_objs[[13]],
                            Seurat_objs[[14]],
                            Seurat_objs[[15]],
                            Seurat_objs[[16]]
                          ),
                          project = "GC-TM-10271"
)

# Create differential interaction meta data from Partek meta data
# Function will create for both SENDERS and RECEIVERS
AllDiffIntMeta <- createMetaData(
  SeuratObj = AllSampSeuratObj,
  InteractionMatList = AllIntMat,
  Attributes = c('orig.ident', 'Cluster')
)

# Add new columns to meta
# Add group
# Add glial or neuronal 
sampleGroups <- read.csv('Tom_groups.csv')
expSamples <- sampleGroups$Sample[sampleGroups$Group == 'Exp']
cntSamples <- sampleGroups$Sample[sampleGroups$Group == 'Cnt']

neuronal <- c(1,2,3,8,10,14)
glial <- c(4,5,7,12,15,16)

# Do for receiver
AllDiffIntMeta$ReceiverMeta <- AllDiffIntMeta$ReceiverMeta %>%
  mutate(group = case_when(orig.ident %in% expSamples ~ 'Exp',
                           orig.ident %in% cntSamples ~ 'Ctrl'),
         celltype = case_when(Cluster %in% neuronal ~ 'neuronal',
                              Cluster %in% glial ~ 'glial')) %>%
  mutate(celltype_group = paste(celltype, group, sep = '_'))

# Do for sender
AllDiffIntMeta$SenderMeta <- AllDiffIntMeta$SenderMeta %>%
  mutate(group = case_when(orig.ident %in% expSamples ~ 'Exp',
                           orig.ident %in% cntSamples ~ 'Ctrl'),
         celltype = case_when(Cluster %in% neuronal ~ 'neuronal',
                              Cluster %in% glial ~ 'glial')) %>%
  mutate(celltype_group = paste(celltype, group, sep = '_'))

# Now run diff interaction for neuronal comparison
neuronalComp <- differentialInteraction(
  InteractionMatList = AllIntMat,
  MetaList = AllDiffIntMeta,
  Attribute = 'celltype_group',
  Comparison = c('neuronal_Exp', 'neuronal_Ctrl')
)

# Write RECEIVER results to file
write.csv(neuronalComp$ReceiverResults,
          file = paste0(base.dir, '/', 'AllSamples_Neuronal_RECEIVER_wilcox.csv'))

# Write SENDER results to file
write.csv(neuronalComp$SenderResults,
          file = paste0(base.dir, '/', 'AllSamples_Neuronal_SENDER_wilcox.csv'))

# Save object
saveRDS(neuronalComp, file = paste0(base.dir, '/neuronalComp.rds'))

# Now run diff interaction for neuronal comparison
glialComp <- differentialInteraction(
  InteractionMatList = AllIntMat,
  MetaList = AllDiffIntMeta,
  Attribute = 'celltype_group',
  Comparison = c('glial_Exp', 'glial_Ctrl')
)

# Write RECEIVER results to file
write.csv(glialComp$ReceiverResults,
          file = paste0(base.dir, '/', 'AllSamples_Glial_RECEIVER_wilcox.csv'))

# Write SENDER results to file
write.csv(glialComp$SenderResults,
          file = paste0(base.dir, '/', 'AllSamples_Glial_SENDER_wilcox.csv'))

# Save object
saveRDS(glialComp, file = paste0(base.dir, '/glialComp.rds'))


# PART 5: CIRCOS PLOTS OF INTERACTIONS FOR SPECIFIC PATHWAYS ACROSS CLUSTERS ----
allInts <- readRDS('interactVis_new_outs_v2/allSamplesAnnoInt.rds')

# PENK interactions
selection <- grep('PENK', allInts$interaction_name)
interactionSum <- allInts[selection,] %>%
  group_by(Sender_Anno, Receiver_Anno) %>%
  summarise(n = n())

# Ammend names of clusters to allow plotting
interactionSum$Sender_Anno <- paste0('C', interactionSum$Sender_Anno)
interactionSum$Receiver_Anno <- paste0('C', interactionSum$Receiver_Anno)

# Set grid colours
grid.col = c(
  `C1` = "#3366cc",
  `C2` = '#dc3911',
  `C3` = '#ff9900',
  `C4` = '#0d9618',
  `C5` = '#990099',
  `C6` = '#0099c5',
  `C7` = '#dd4477',
  `C8` = '#66a900',
  `C9` = '#b72e2f',
  `C10` = '#6633cc',
  `C11` = '#22a999',
  `C12` = '#306395',
  `C13` = '#aaaa11',
  `C14` = '#984499',
  `C15` = '#e67301',
  `C16` = '#8b0607',
  `C17` = '#339262',
  `C18` = '#3a3eac',
  `C19` = '#651066'
)

pdf(
  file = paste0('interactVis_new_outs_v2/penk_chord.pdf'),
  height = 6,
  width = 6
)
circos.clear()
circos.par(gap.after = 5)
chordDiagram(interactionSum, transparency = 0.5, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = mm_h(5)))
for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.3, sector.index = si, track.index = 2, labels.facing = 'clockwise')
}
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  xplot = get.cell.meta.data("xplot")
  circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, -1.5), cex = 0.6)
}, bg.border = NA)
circos.clear()
dev.off()

# TAC interactions
selection <- grep('TAC1', allInts$interaction_name)
interactionSum <- allInts[selection,] %>%
  group_by(Sender_Anno, Receiver_Anno) %>%
  summarise(n = n())

# Ammend names of clusters to allow plotting
interactionSum$Sender_Anno <- paste0('C', interactionSum$Sender_Anno)
interactionSum$Receiver_Anno <- paste0('C', interactionSum$Receiver_Anno)

pdf(
  file = paste0('interactVis_new_outs_v2/tac1_chord.pdf'),
  height = 6,
  width = 6
)
circos.clear()
circos.par(gap.after = 5)
chordDiagram(interactionSum, transparency = 0.5, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = mm_h(5)))
for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.3, sector.index = si, track.index = 2, labels.facing = 'clockwise')
}
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  xplot = get.cell.meta.data("xplot")
  circos.text(mean(xlim), 1, sector.name, niceFacing = TRUE, adj = c(0.5, -1.5), cex = 0.6)
}, bg.border = NA)
circos.clear()
dev.off()






