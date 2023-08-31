# Script to try calculating interaction scores as part of the interactVis package

# Set working directory
setwd('/Users/jamesboot/Documents/9.Genome Centre Files/interactVis/')

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Load spatial data
data.dir <- 'testdata/GC-TM-10271-N12-297/'
dat <- Load10X_Spatial(data.dir, 
                             filename = 'filtered_feature_bc_matrix.h5', 
                             assay = 'Spatial', 
                             filter.matrix = TRUE, 
                             to.upper = FALSE)

# Normalise data
dat <- NormalizeData(dat, 
                     assay = 'Spatial',
                     normalization.method = 'LogNormalize',
                     scale.factor = 10000)

# Load custom functions for InteractVis
source('findNeighbours.R')
source('loadDB.R')

# Find neighbors 
neighbours <- findNeighbours('testdata/GC-TM-10271-N12-297/spatial/tissue_positions.csv')

# Load database
database <- loadDB(gene.input = 'cellphonedb-data-4.0.0/data/gene_input_all.csv',
                   complex.input = 'cellphonedb-data-4.0.0/data/complex_input.csv',
                   interaction.input = 'cellphonedb-data-4.0.0/data/interaction_input.csv')

# Pre-processing before attempting to calculate interaction values ----

# Create vector containing the barcodes that are in the neighbors list and in the data (i.e. actual tissue spots)
select <- names(neighbours)[names(neighbours) %in% dat@assays$Spatial@counts@Dimnames[[2]]]

# Get the scaled data in a dataframe
scaledat <- as.data.frame(GetAssayData(object = dat, slot = "data"))

# Filter the scaled data down to all genes in the database
# Create a vector of all the genes in the database 
db_elements <- unique(as.vector(unlist(database)))

# Go through each gene and select it from scaledat, if not present insert 0s
scaledat_filt <- scaledat[0,]

for (x in 1:length(db_elements)) {
  
  if (db_elements[x] %in% row.names(scaledat)) {
    
    scaledat_filt[db_elements[x], ] <- scaledat[db_elements[x], ]
    
  } else {
    
    scaledat_filt[db_elements[x], ] <- 0
    
  }
}

# Now calculate the interaction scores ----

# Create dummy dataframe
test <- data.frame()

# Iterate through each of the barcodes in the select vector created earlier
# Store the primary barcode and the neighbour barcodes in objects
# Ensure to add barcode1 to neighbors so it checks interaction within itself
# Ensure that the neighbor barcodes are in the select vector!
for (y in 1:length(select)) {
  
  bcode1 <- select[y]
  nbours <-
    c(neighbours[[bcode1]]$neighbours[neighbours[[bcode1]]$neighbours %in% select],
      bcode1)
  print(paste0(
    'Starting barcode ',
    y,
    ' of ',
    length(select),
    '. Progress: ',
    round(((y / length(select)) * 100), digits = 1),
    '%'
  ))
  
  # Iterate through the neighbors - store barcode
  for (x in 1:length(nbours)) {
    
    bcode2 <- nbours[x]
    
    # Now for the pair of barcodes selected iterate through the database
    # Find the minimum expression value of all the genes in the complex in barcode1
    # Find the expression of all the ligands in barcode2
    for (v in 1:length(database)) {
      
      spot1_complex <- scaledat_filt[database[[v]]$Complex_Genes, bcode1]
      complex_expr <- min(spot1_complex)
      
      if (complex_expr != 0) {
        
        spot2_ligands <- data.frame(gene = database[[v]]$Partner_Genes,
                                    expr = scaledat_filt[database[[v]]$Partner_Genes, bcode2])
        
        # There may be more than 1 ligand for the complex - iterate through ligands
        # Store the barcodes, complex name & expression, ligand name & expression in dataframe
        # Append the dataframe with the new information, return to start of loop
        for (ligands in 1:nrow(spot2_ligands)) {
          
          if (spot2_ligands$expr[ligands] != 0) {
            
            new_entry <- data.frame(
              spot1 = bcode1,
              spot2 = bcode2,
              spot1_complex = names(database)[v],
              spot1_complex_expr = complex_expr,
              spot2_ligand = spot2_ligands$gene[ligands],
              spot2_ligand_expr = spot2_ligands$expr[ligands]
            )
            
            test <- rbind(test, new_entry)
            
          } else {
            next
          }
        }
      } else {
        next
      }
    }
  }
}

# Calculate interaction score
test <-
  test %>%
  rowwise() %>%
  mutate(interaction_score = mean(c_across(
    c(spot1_complex_expr, spot2_ligand_expr)
  ))) %>%
  filter(!is.infinite(interaction_score))

# Find mean and standard deviation of interaction scores
mean_intscore <- mean(test$interaction_score)
sd_intscore <- sd(test$interaction_score)

# Identify rows with interaction scores 2 SDs greater than mean
results <- test[test$interaction_score > (mean_intscore+(2*sd_intscore)), ]

# Sort by interaction score
results <- results %>%
  arrange(desc(interaction_score))

# Count occurrences of different complexes in final results
complex_counts <- results %>%
  count(spot1_complex) %>%
  arrange(desc(n))

# Extract all spots enriched for top occurring complex interaction
coi <- complex_counts$spot1_complex[12]
spots <- unique(c(results$spot1[results$spot1_complex == coi],
                  results$spot2[results$spot1_complex == coi]))

# Create dataframe to add this to meta data
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
p2a <- SpatialFeaturePlot(dat, features = database[[coi]][["Complex_Genes"]][1])
p2b <- SpatialFeaturePlot(dat, features = database[[coi]][["Complex_Genes"]][2])
p2a+p2b

p3a <- SpatialFeaturePlot(dat, features = database[[coi]][["Partner_Genes"]][1])
p3b <- SpatialFeaturePlot(dat, features = database[[coi]][["Partner_Genes"]][2])
p3c <- SpatialFeaturePlot(dat, features = database[[coi]][["Partner_Genes"]][3])
p3d <- SpatialFeaturePlot(dat, features = database[[coi]][["Partner_Genes"]][4])
p3e <- SpatialFeaturePlot(dat, features = database[[coi]][["Partner_Genes"]][5])
p3f <- SpatialFeaturePlot(dat, features = database[[coi]][["Partner_Genes"]][6])
p3g <- SpatialFeaturePlot(dat, features = database[[coi]][["Partner_Genes"]][7])
p3h <- SpatialFeaturePlot(dat, features = database[[coi]][["Partner_Genes"]][8])

p3b
 





