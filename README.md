# interactVis

README last updated: 06/05/2024

## :dart: Aim
- Perform receptor-ligand interaction analysis of 10X Visium spatial transcriptomics data - hence - interactVis
- Find adjacent tissue spots where there is co-expression of receptors and ligands
- Compare interaction scores of different tissue spot populations (e.g. clusters)

## :nut_and_bolt: Rationale
The processing steps are performed in a similar manner to that described by Kaufmann et al. (1).
1. Load cellphoneDB (2) database.
2. Determine the expression of each receptor and ligand in each tissue spot.
3. In the case of receptor complexes the expression value of the least expressed subunit is taken as the expression value. If any subunit is not expressed the expression value of the complex is set as 0.
4. Next interaction scores are determined for each tissue spot, and its immediate surrounding neighbours. Interaction scores are also calculated within a tissue spot. The interaction score is defined as the mean expression of receptor A and ligand B in tissue spot 1 and tissue spot 2 respectively. Receptor expression is taken from spot 1 – therefore spot 1 is the *receiver* of signalling, whilst ligand expression is always taken from spot 2; therefore spot 2 is the *sender* of signalling. This calculation processes is iterated across all tissue spots in combination with all neighbours of each tissue spot.
5. Once interaction scores are calculated each tissue spot is annotated with its identity (e.g. cluster). 
6. Interaction scores for a given cluster can then be extracted. Extracting interaction scores based on the spot 1 cluster annotation identifies the interaction scores that spots for a given cluster were *receiving*, whilst using the spot 2 cluster annotation identifies the interaction scores that spots for a given cluster were *sending*.
7. Finally, interaction scores for each receptor-ligand pair detected are compared between two identities (e.g. cluster), in terms of both what each identity is *sending* and *receiving*. Interaction scores between clusters are compared using a Wilcoxon Test OR negative bionmial GLM with log-likelihood ratio test, and p-values adjusted for multiple hypothesis testing using the Benjamini-Hochberg procedure.

## :pencil: Getting started
Implementation of the functions and analysis can be found in `Analysis_MultiSample.R`, this also contains code for visualisations of results.
### Dependencies
The following packages are needed to run scripts in this repository:
```
library(Seurat)
library(ggplot2)
library(dplyr)
library(multcomp)
library(purrr)
library(tibble)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyr)
library(lmtest)
```
### Load functions and database
```
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
```
### Load and pre-process sample data
- Load data as a Seurat object using `Load10X_Spatial()`
```
# Load spatial data
dat <- Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = 'Spatial', 
  filter.matrix = TRUE,
  to.upper = FALSE
)
```
- Pre-process sample data as per Seurat guidelines, but most importantly: `NormalizeData()`
- :exclamation: Note that the seurat object should contain the metadata / identities you wish to compare later.
### Find neighbours
This funciton uses the `tissue_positions.csv` in the sample spatial directory to make a list of all tissue spot barcodes and their neighbours.
```
# Find neighbors
neighbours <- findNeighbours('spatial/tissue_positions.csv')
```
### Calculate interactions
- This function calculates the interaction score between all neighbouring tissue spots for every receptor-ligand pair in the database.
- If you want, here you can filter the interaction scores for those 'X' standard deviations greater than the mean of all interaction scores.
- :exclamation: This step takes time to compute.
```
# Find interactions
interactions <- calculateInteractions(
  neighboursList = neighbours,
  dat = dat,
  database = database,
  filter = FALSE
)
```
### Annotate interactions
- Annotate spot 1 and spot 2 of each interaction with the metadata / identity you wish to compare groups of.
- Attribute should be the name of the metadata / identity in the Seurat object.
```
# Annotate the Interactions results with the cluster for spot 1 and spot 2
# Spot 1 is the receiver (where the receptor expression is taken from)
# Spot 2 is the sender (where the ligand expression is taken from)
diffIntDF <- annotateInteractions(
  SeuratObj = dat,
  Interactions = interactions,
  Attribute = 'Cluster'
)
```
### Summarising interactions at an identity level
- The number of interactions per cluster and between cluster can be summarised as such:
```
# How many of each interactions per cluster are there?
# And between clusters
# Find number of interactions between groups
interactionSum <- diffIntDF %>%
  group_by(Spot1_Anno, Spot2_Anno) %>%
  summarise(n = n())
write.csv(interactionSum, file = 'sample_nInt_cluster.csv')
```
- This can be visually summarised using a Chord plot
- Colours of each group can be specified using `grid.col`.
```
circos.clear()
circos.par(gap.after = 5)
chordDiagram(interactionSum, transparency = 0.5, grid.col = grid.col)
```
### Differential interaction analysis
- After the above steps we can now begin to compare the interaction scores for all receptor-ligand pairs across two groups.
- First, create a numeric matrix of interaction scores from each spot for each receptor-ligand pair. `interactionMatrix` creates two matrices: one for *receivers* (i.e. spot 1) and one for *senders* (i.e. spot 2).
- Second, create a metadata dataframe
  - Specify the Seurat object, the interaction matrix, and the attributes from the Seurat object you wish to add to the metadata dataframe.
- Finally, perform differential interaction analysis
  - Specify the:
    - Interaction matrix (InteractionMatList)
    - Metadata dataframe (MetaList)
    - Attribute name of the identities you want to compare (Attribute)
    - Two groups you want to compare (Comparison)
    - Test to be performed (Wilcoxon OR negative bionmial GLM with log-likelihood ratio test)  (Test)
```
# Create Interaction Matrix
# Function will create for both SENDERS and RECEIVERS
IntMat <- interactionMatrix(AnnoInt = diffIntDF)

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
  Comparison = c(2, 10),
  Test = 'glm'
)
```
## :page_with_curl: References

1) Kaufmann, M., Schaupp, AL., Sun, R. et al. Identification of early neurodegenerative pathways in progressive multiple sclerosis. Nat Neurosci 25, 944–955 (2022). https://doi.org/10.1038/s41593-022-01097-3
2) Efremova, M., Vento-Tormo, M., Teichmann, S.A. et al. CellPhoneDB: inferring cell–cell communication from combined expression of multi-subunit ligand–receptor complexes. Nat Protoc 15, 1484–1506 (2020). https://doi.org/10.1038/s41596-020-0292-x
