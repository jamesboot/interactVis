# Script for implementing, testing and developing interactVis functions

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(multcomp)
library(purrr)
library(tibble)
library(edgeR)

# Change working directory
setwd('/Users/jamesboot/Documents/GitHub/')

# Load functions
source('interactVis/findNeighbours.R')
source('interactVis/loadDB.R')
source('interactVis/findInteractions.R')
source('interactVis/annotateInteractions.R')
source('interactVis/interactionMatrix.R')
source('interactVis/interactionMeta.R')

# Load database
database <- loadDB(gene.input = 'interactVis/cellphonedb-data-4.0.0/data/gene_input_all.csv',
                   complex.input = 'interactVis/cellphonedb-data-4.0.0/data/complex_input.csv',
                   interaction.input = 'interactVis/cellphonedb-data-4.0.0/data/interaction_input.csv')

# Start analysis on test data set from project GC-TM-10271

# Change working directory
setwd('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/interactVis_trial/')

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

# Find neighbors 
neighbours <- findNeighbours(paste0(data.dir, 'spatial/tissue_positions.csv'))

# Find interactions
interactions <- calculateInteractions(neighboursList = neighbours,
                                      dat = dat,
                                      database = database,
                                      filter = FALSE)

# DATA VISUALISATION ----

# Extract all spots enriched for top occurring complex interaction
coi <- interactions$Complex_Summary$spot1_complex[2]
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

# Add new meta data to object
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
diffIntDF <- annotateInteractions(SeuratObj = dat,
                                  Interactions = interactions,
                                  Attribute = 'Cluster')

# Interactions will be annotated with spot 1 cluster
# Because spot1 is where the receptor expression is taken from (i.e. the actuator of signalling) 
# How many of each interactions per cluster are there?
interactionSum <- diffIntDF %>%
  group_by(Spot1_Cluster) %>%
  summarise(n = n())
View(interactionSum)

# Which receptor ligands have the most interactions across the tissue?
recepLigSum <- diffIntDF  %>%
  group_by(spot1_complex) %>%
  summarise(n = n()) %>%
  arrange(n)

# Plot - How many total interactions does each cluster have, regardless of receptor-ligand?
p5a <-
  ggplot(interactionSum,
         aes(x = Spot1_Cluster, y = n)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ggtitle('Total number of receptor-ligand interactions per cluster')
p5a
ggsave(
  'n_interactions_cluster.tiff',
  plot = p5a,
  width = 20,
  height = 10,
  units = 'in',
  dpi = 300
)

# Plot - What is the average interaction score per cluster, regardless of receptor-ligand?
p5b <-
  ggplot(diffIntDF,
         aes(x = Spot1_Cluster, y = interaction_score)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ggtitle('Average interaction score per cluster')
p5b
ggsave(
  'av_interaction_score.tiff',
  plot = p5b,
  width = 20,
  height = 10,
  units = 'in',
  dpi = 300
)

# Create Interaction Matrix
IntMat <- interactionMatrix(AnnoInt = diffIntDF)

# Create differential interaction meta data from Partek meta data
diffIntMeta <- createMetaData(SeuratObj = dat,
                              InteractionMat = IntMat,
                              Attributes = c('orig.ident', 'Cluster'))

# Function to perform multiple t-tests to compare interaction scores between clusters
differentialInteraction <- function(InteractionMat,
                                    Meta,
                                    Attribute,
                                    Comparison) {
  
  # Transpose the interaction mat
  interactionDF <- as.data.frame(t(InteractionMat))
  
  # Add cell as column in meta and count mat
  interactionDF <- interactionDF %>%
    rownames_to_column(var = 'cell')
  diffIntMeta2 <- Meta %>%
    rownames_to_column(var = 'cell')
  
  # Now join meta and column
  tTestDF <- interactionDF %>%
    left_join(diffIntMeta2, by = 'cell')
  
  # Subset down to clusters we want to compare
  # Focus on 2 and 10 for now
  tTestDFFilt <- tTestDF %>%
    filter(Attribute %in% Comparison)
  
  # Make dataframe for results
  wilcoxRes <- data.frame()
  
  # For loop to perform test on all columns
  for (x in c(2:(ncol(tTestDFFilt) - ncol(Meta)))) {
    
    # As data is not normally distributed - use wilcoxon
    res <-
      wilcox.test(tTestDFFilt[, x] ~ Attribute,
                  data = tTestDFFilt,
                  paired = F)
    
    # Adjust the p-value
    FDR <-
      p.adjust(res$p.value, method = 'fdr', n = length(c(2:(
        ncol(tTestDFFilt) - 2
      ))))
    
    # Put result into dataframe
    wilcoxRes[colnames(tTestDFFilt)[x], 'Receptor-Ligand'] <-
      colnames(tTestDFFilt)[x]
    wilcoxRes[colnames(tTestDFFilt)[x], 'P-Value'] <- res$p.value
    wilcoxRes[colnames(tTestDFFilt)[x], 'FDR'] <- FDR
    
  }
}

write.csv(wilcoxRes, file = paste0(comp, '_wilcox_res.csv'))

