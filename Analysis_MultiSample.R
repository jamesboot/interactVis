# Script for implementing, testing and developing interactVis functions

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(multcomp)
library(purrr)
library(tibble)
library(ComplexHeatmap)
library(RColorBrewer)

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
  diffIntDF <- annotateInteractions(SeuratObj = dat,
                                    Interactions = interactions,
                                    Attribute = 'Cluster')
  
  # Interactions will be annotated with spot 1 cluster
  # Because spot1 is where the receptor expression is taken from (i.e. the actuator of signalling)
  
  # How many of each interactions per cluster are there?
  # And between clusters
  # Plot as heatmap
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
      unique(interactionSum$Spot1_Anno)
    )
  )
  
  # Populate matrix in for loop
  for (row in unique(interactionSum$Spot1_Anno)) {
    for (col in unique(interactionSum$Spot1_Anno)) {
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
  Heatmap(
    interactionSumMat,
    col = colfun,
    top_annotation = column_bp,
    left_annotation = row_bp,
    name = 'n_intrtns'
  )
  dev.off()
  
  # DIFFERENTIAL CLUSTER ANALYSIS ----
  
  # Create Interaction Matrix
  IntMat <- interactionMatrix(AnnoInt = diffIntDF)
  
  # Create differential interaction meta data from Partek meta data
  diffIntMeta <- createMetaData(
    SeuratObj = dat,
    InteractionMat = IntMat,
    Attributes = c('orig.ident', 'Cluster')
  )
  
  # Perform differential interaction analysis
  diffInt2v10 <- differentialInteraction(
    InteractionMat = IntMat,
    Meta = diffIntMeta,
    Attribute = 'Cluster',
    Comparison = c(2, 10)
  )
  
  write.csv(diffInt2v10,
            file = paste0(base.dir, '/', samples[ITER], '/', samples[ITER], '_wilcox.csv'))
  
}

