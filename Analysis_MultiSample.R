# Script for implementing, testing and developing interactVis functions

#### PART 1: Load packages, data and preprocess ----

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

# Change working directory
setwd('/Users/jamesboot/Documents/GitHub/')

# Load functions
source('interactVis/findNeighbours.R')
source('interactVis/loadDB.R')
source('interactVis/findInteractions.R')
source('interactVis/annotateInteractions.R')
source('interactVis/interactionMatrix.R')
source('interactVis/interactionMeta.R')
source('interactVis/differentialInteraction.R')

# Load database
database <- loadDB(gene.input = 'interactVis/cellphonedb-data-4.0.0/data/gene_input_all.csv',
                   complex.input = 'interactVis/cellphonedb-data-4.0.0/data/complex_input.csv',
                   interaction.input = 'interactVis/cellphonedb-data-4.0.0/data/interaction_input.csv')

# Start analysis on all data from project GC-TM-10271

# Change working directory
setwd('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/interactVis_trial/')
base.dir <- '/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/interactVis_trial/AllSamplesTrial/'

# Create list of all sample directories
dirs <- list.dirs(path = '/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/4 Runs Combined/spatial',
                  recursive = F)

# Create list of all sample .h5 files
files <- list.files(path = '/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/4 Runs Combined/spatial',
                    pattern = '*.h5',
                    recursive = T,
                    full.names = T)

# Create list of sample names
samples <- gsub(
  '_spatial',
  '',
  gsub(
    '/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/4 Runs Combined/spatial/',
    '',
    dirs
  )
)

#### PART 2: Loop through all samples and perform interaction analysis ----
                
# For loop to go through all samples and perform all analysis
for (ITER in 1:c(length(samples))) {
  
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
      '/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/stDeconvolve_exp_v_ctrl/pt1_local/selectedSpots.RDS'
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
    neighboursList = neighbours,
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
    read.delim('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/exported.txt')
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
    filename = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_cluster.tiff'),
    height = 5,
    width = 5,
    dpi = 300,
    units = 'in'
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
    group_by(Spot1_Anno, Spot2_Anno) %>%
    summarise(n = n())
  write.csv(interactionSum,
            file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_nInt_cluster.csv'))
  
  # Convert to matrix
  # Create empty matrix to population
  interactionSumMat <- matrix(
    data = NA,
    nrow = length(unique(interactionSum$Spot1_Anno)),
    ncol = length(unique(interactionSum$Spot2_Anno)),
    dimnames = list(
      unique(interactionSum$Spot1_Anno),
      unique(interactionSum$Spot2_Anno)
    )
  )
  
  # Populate matrix in for loop
  for (row in unique(interactionSum$Spot1_Anno)) {
    for (col in unique(interactionSum$Spot2_Anno)) {
      # Find the n number
      n <- interactionSum %>%
        filter(Spot1_Anno == row & Spot2_Anno == col) %>%
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
  tiff(
    filename = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_htmap.tiff'),
    height = 5,
    width = 5,
    res = 300,
    units = 'in'
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
nInt_files <- list.files(path = '/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/interactVis_trial',
                       recursive = T,
                       pattern = '_nInt_cluster.csv',
                       full.names = T)

# Load CSV files 
nInt_objs <- lapply(nInt_files, function(x){
  read.csv(x)[, 2:4]
})

# Names
names(nInt_objs) <- samples

# Set grid colours
grid.col = c(
  `1` = "#3366cc",
  `2` = '#dc3911',
  `3` = '#ff9900',
  `4` = '#0d9618',
  `5` = '#990099',
  `6` = '#0099c5',
  `7` = '#dd4477',
  `8` = '#66a900',
  `9` = '#b72e2f',
  `10` = '#6633cc',
  `11` = '#22a999',
  `12` = '#306395',
  `13` = '#aaaa11',
  `14` = '#984499',
  `15` = '#e67301',
  `16` = '#8b0607',
  `17` = '#339262',
  `18` = '#3a3eac',
  `19` = '#651066'
)

# Chord plot for each sample
for (ITER in 1:length(nInt_objs)) {
  tiff(
    filename = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_chord.tiff'),
    height = 4,
    width = 4,
    units = 'in',
    res = 300
  )
  circos.clear()
  circos.par(gap.after = 5)
  chordDiagram(nInt_objs[[ITER]], transparency = 0.5, grid.col = grid.col)
  dev.off()
}

# Plot all samples together
combined_df <- do.call(rbind, nInt_objs)
tiff(
  filename = paste0(base.dir, '/allSamples_chord.tiff'),
  height = 4,
  width = 4,
  units = 'in',
  res = 300
)
circos.clear()
circos.par(gap.after = 5)
chordDiagram(combined_df, transparency = 0.5, grid.col = grid.col)
dev.off()

# Plot Exp samples together
sampleGroups <- read.csv('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/Tom_groups.csv')
expSamples <- sampleGroups$Sample[sampleGroups$Group == 'Exp']
combined_df <- do.call(rbind, nInt_objs[expSamples])
tiff(
  filename = paste0(base.dir, '/expSamples_chord.tiff'),
  height = 4,
  width = 4,
  units = 'in',
  res = 300
)
circos.clear()
circos.par(gap.after = 5)
chordDiagram(combined_df, transparency = 0.5, grid.col = grid.col)
dev.off()

# Plot Ctrl samples together
sampleGroups <- read.csv('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/Tom_groups.csv')
expSamples <- sampleGroups$Sample[sampleGroups$Group == 'Cnt']
combined_df <- do.call(rbind, nInt_objs[expSamples])
tiff(
  filename = paste0(base.dir, '/ctrlSamples_chord.tiff'),
  height = 4,
  width = 4,
  units = 'in',
  res = 300
)
circos.clear()
circos.par(gap.after = 5)
chordDiagram(combined_df, transparency = 0.5, grid.col = grid.col)
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
      '/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/stDeconvolve_exp_v_ctrl/pt1_local/selectedSpots.RDS'
    )
  
  # Now filter each Seurat object in list for selected tissue spots
  Seurat_objs[[samples[ITER]]] <-
    subset(Seurat_objs[[samples[ITER]]], cells = selectedSpots$`Cell name`[selectedSpots$`Sample name` == unique(Seurat_objs[[samples[ITER]]]$orig.ident)])
  
  # Normalise data
  Seurat_objs[[samples[ITER]]] <-
    NormalizeData(Seurat_objs[[samples[ITER]]], normalization.method = "LogNormalize", verbose = FALSE)
  
  # Add cluster numbers from Partek
  partekMeta <-
    read.delim('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/exported.txt')
  partekMetaFilt <-
    partekMeta[partekMeta$Sample.name == samples[ITER],]
  clustersMeta <- data.frame(
    row.names = partekMetaFilt$Cell.name,
    Cluster = as.factor(partekMetaFilt$Clusters..5comp.0.75res)
  )
  Seurat_objs[[samples[ITER]]] <- AddMetaData(Seurat_objs[[samples[ITER]]], clustersMeta)

}

# Locate RDS files of interaction results
rds_files <- list.files(path = '/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/interactVis_trial',
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
  AnnoInt[[samples[ITER]]] <-
    annotateInteractions(
      SeuratObj = Seurat_objs[[samples[ITER]]],
      Interactions = rds_objs[[samples[ITER]]],
      Attribute = 'Cluster'
    )
}

# Loop to go through all samples and perform differential interaction
AllDiffIntRes <- list()
for (ITER in 1:c(length(samples))) {
  
  # Create Interaction Matrix
  # Function will create for both SENDERS and RECEIVERS
  IntMat <- interactionMatrix(AnnoInt = AnnoInt[[samples[ITER]]])
  
  # Create differential interaction meta data from Partek meta data
  # Function will create for both SENDERS and RECEIVERS
  diffIntMeta <- createMetaData(
    SeuratObj = Seurat_objs[[samples[ITER]]],
    InteractionMatList = IntMat,
    Attributes = c('orig.ident', 'Cluster')
  )
  
  # Perform differential interaction analysis
  # Function will create for both SENDERS and RECEIVERS
  # Function will also check there are enough replicates 
  diffInt2v10 <- differentialInteraction(
    InteractionMatList = IntMat,
    MetaList = diffIntMeta,
    Attribute = 'Cluster',
    Comparison = c(2, 10)
  )
  
  # Write RECEIVER results to file
  write.csv(diffInt2v10$ReceiverResults,
            file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_RECEIVER_wilcox.csv'))
  
  # Write SENDER results to file
  write.csv(diffInt2v10$SenderResults,
            file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_SENDER_wilcox.csv'))
  
  # Write result for sample to list
  AllDiffIntRes[[samples[ITER]]] <- diffInt2v10
  
}
  


