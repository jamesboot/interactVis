# Script for implementing, testing and developing interactVis functions

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(multcomp)

# Change working directory
setwd('/Users/jamesboot/Documents/GitHub/')

# Load functions
source('interactVis/findNeighbours.R')
source('interactVis/loadDB.R')
source('interactVis/findInteractions.R')

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
# Create new df for this
diffIntDF <- interactions$Interactions
diffIntDF$Spot1_Cluster <- NA
diffIntDF$Spot2_Cluster <- NA
for (x in 1:nrow(diffIntDF)) {
  diffIntDF$Spot1_Cluster[x] <-
    partekMetaFilt$Clusters..5comp.0.75res[partekMetaFilt$Cell.name == diffIntDF$spot1[x]]
  diffIntDF$Spot2_Cluster[x] <-
    partekMetaFilt$Clusters..5comp.0.75res[partekMetaFilt$Cell.name == diffIntDF$spot2[x]]
}
View(diffIntDF)

# Ensure Spot1_Cluster and Spot2 Cluster are factors
diffIntDF$Spot1_Cluster <- as.factor(diffIntDF$Spot1_Cluster)
diffIntDF$Spot2_Cluster <- as.factor(diffIntDF$Spot2_Cluster)

# Interactions will be annotated with spot 1 cluster
# Because spot1 is where the receptor expression is taken from (i.e. the actuator of signalling) 

# How many of each interactions per cluster are there?
interactionSum <- diffIntDF %>%
  group_by(Spot1_Cluster) %>%
  summarise(n = n())
View(interactionSum)

# Plot - How many total interactions does each cluster have, regardless of receptor-ligand?
p5a <-
  ggplot(interactionSum,
         aes(x = Spot1_Cluster, y = n)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))
p5a

# Plot - What is the average interaction score per cluster, regardless of receptor-ligand?
p5b <-
  ggplot(diffIntDF,
         aes(x = Spot1_Cluster, y = interaction_score)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))
p5b

# Which receptor ligands have the most interactions across the tissue?
recepLigSum <- diffIntDF  %>%
  group_by(spot1_complex) %>%
  summarise(n = n()) %>%
  arrange(n)

# Select the receptor ligand pair with most interactions across the tissue
coi <- tail(recepLigSum$spot1_complex, n = 1)

# Subset down to complex of interest and
# Remove rows which have 10 or fewer occurrences for each cluster
diffIntDF_filt <- diffIntDF %>%
  filter(spot1_complex == coi) %>%
  group_by(Spot1_Cluster) %>%
  filter(n() >= 10)

# Plot box plot of interaction scores for complex of interest across spot 1 clusters
p5c <-
  ggplot(diffIntDF_filt,
         aes(x = Spot1_Cluster, y = interaction_score)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ))
p5c
ggsave(
  'score_boxplot.tiff',
  plot = p5,
  width = 20,
  height = 10,
  units = 'in',
  dpi = 300
)

# Is data normally distributed?
qqnorm(diffIntDF_filt$interaction_score)
qqline(diffIntDF_filt$interaction_score)

# Descriptive stats
stats <- aggregate(interaction_score ~ Spot1_Cluster,
                   data = diffIntDF_filt,
                   function(x)
                     round(c(mean = mean(x), sd = sd(x)), 2))

# Lets compare the interaction scores in two different clusters 
# Use Welchs ANOVA as not normally distributed 
res_aov1 <- oneway.test(interaction_score ~ Spot1_Cluster,
                        data = diffIntDF_filt,
                        var.equal = FALSE)
res_aov1
res_aov2 <- aov(interaction_score ~ Spot1_Cluster,
                data = diffIntDF_filt)
summary(res_aov2)

# If there is a difference perform Tukeys post hoc test
# Tukey HSD test:
post_test <- TukeyHSD(res_aov2, conf.level = 0.99)
plot(TukeyHSD(res_aov2, conf.level=.99), las = 2)

# Extract significant results 
post_test_df <- as.data.frame(post_test$Spot1_Cluster)
post_test_df_filt <- post_test_df[post_test_df$`p adj` < 0.01, ]


