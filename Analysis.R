# Script to try calculating interaction scores as part of the interactVis package

# Set working directory
setwd('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/interactVis_trial/')

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load spatial data
data.dir <- '/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/4 Runs Combined/spatial/NH12-297_spatial/'
dat <- Load10X_Spatial(
  data.dir,
  filename = 'NH12-297_filtered_feature_bc_matrix.h5',
  assay = 'Spatial',
  filter.matrix = TRUE,
  to.upper = FALSE
)

# Set the orig.ident to the sample name
sample <- gsub('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/4 Runs Combined/spatial/',
               '',
               data.dir)
sample <- gsub('_spatial/', '', sample)
dat$orig.ident <- sample

# Filter to only spots used by Tom in Partek
# Import the selected spots object for filtering to Tom's spots
selectedSpots <- readRDS('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/stDeconvolve_exp_v_ctrl/pt1_local/selectedSpots.RDS')

# Now filter each Seurat object in list for selected tissue spots
dat <- subset(dat, cells = selectedSpots$`Cell name`[selectedSpots$`Sample name` == unique(dat$orig.ident)])

# Normalise data
#dat <- SCTransform(dat, assay = "Spatial", verbose = FALSE)
dat <- NormalizeData(dat, normalization.method = "LogNormalize", verbose = FALSE)

# Load custom functions for InteractVis
source('findNeighbours.R')
source('loadDB.R')
source('findInteractions.R')

# Find neighbors 
neighbours <- findNeighbours(paste0(data.dir, 'spatial/tissue_positions.csv'))

# Load database
database <- loadDB(gene.input = 'cellphonedb-data-4.0.0/data/gene_input_all.csv',
                   complex.input = 'cellphonedb-data-4.0.0/data/complex_input.csv',
                   interaction.input = 'cellphonedb-data-4.0.0/data/interaction_input.csv')

# Find interactions
interactions <- calculateInteractions(neighboursList = neighbours,
                                      dat = dat,
                                      database = database,
                                      filter = FALSE)

# DATA VISUALISATION ----

# Extract all spots enriched for top occurring complex interaction
coi <- interactions$Complex_Summary$spot1_complex[1]
spots <- unique(c(interactions$Interactions$spot1[interactions$Interactions$spot1_complex == coi],
                  interactions$Interactions$spot2[interactions$Interactions$spot1_complex == coi]))

# Create dataframe to add this to meta data
select <- names(neighbours)[names(neighbours) %in% dat@assays$Spatial@counts@Dimnames[[2]]]
new_meta <- data.frame(barcodes = select,
                       interaction1 = NA)

new_meta$interaction1[new_meta$barcodes %in% spots] <- coi
new_meta$interaction1[!new_meta$barcodes %in% spots] <- 'None'

row.names(new_meta) <- new_meta$barcodes
new_meta$barcodes <- NULL

# Add new meta data ot object
dat <- AddMetaData(dat, new_meta)

# Visualise the new meta data along with genes and receptors of the complex
p1 <- SpatialPlot(dat, group.by = 'interaction1')
p1

database[[coi]]
p2 <- SpatialFeaturePlot(dat, features = database[[coi]][["Complex_Genes"]])
p2

p3 <- SpatialFeaturePlot(dat, features = database[[coi]][["Partner_Genes"]])
p3

# DIFFERENTIAL CLUSTER ANALYSIS ----

# Add cluster numbers from Partek
partekMeta <- read.delim('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/exported.txt')
partekMetaFilt <- partekMeta[partekMeta$Sample.name == 'NH12-297', ]
clustersMeta <- data.frame(row.names = partekMetaFilt$Cell.name,
                           Cluster = as.factor(partekMetaFilt$Clusters..5comp.0.75res))
dat <- AddMetaData(dat, clustersMeta)

# Visualise the clusters for this sample
p4 <- SpatialPlot(dat, group.by = 'Cluster', label.size = 10) +  
  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
p4

# Annotate the Interactions results with the cluster for spot 1 and spot 2
interactions$Interactions$Spot1_Cluster <- NA
interactions$Interactions$Spot2_Cluster <- NA
for (x in 1:nrow(test)) {
  interactions$Interactions$Spot1_Cluster[x] <- 
    partekMetaFilt$Clusters..5comp.0.75res[partekMetaFilt$Cell.name == interactions$Interactions$spot1[x]]
  interactions$Interactions$Spot2_Cluster[x] <- 
    partekMetaFilt$Clusters..5comp.0.75res[partekMetaFilt$Cell.name == interactions$Interactions$spot2[x]]
}
View(interactions$Interactions)

# Lets compare the interaction scores in two different clusters 
#
