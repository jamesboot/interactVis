# Function to make a matrix with spot1 as columns and receptor-ligand as rows, interaction scores as values
interactionMatrix <- function(AnnoInt) {
  
  # Make a dummy list
  dumList <- list()
  
  # Create new receptor ligand column in AnnoInt
  AnnoInt <- AnnoInt %>%
    mutate(Receptor_Ligand = paste0(spot1_complex, '-', spot2_ligand))
  
  # For loop to go through each cell
  # Use dplyr to extract and summarise interaction scores for all receptor-ligands
  # Put the new dataframe for the spot into dummy list
  c <- 1
  for (bc in unique(AnnoInt$spot1)) {
    message(paste0(
      'Starting spot ',
      c,
      ' of ',
      length(unique(AnnoInt$spot1)),
      '. Progress: ',
      round((c / length(
        unique(AnnoInt$spot1)
      ) * 100), digits = 1),
      '%'
    ))
    int <- AnnoInt %>%
      filter(spot1 == bc) %>%
      group_by(Receptor_Ligand) %>%
      summarise(Mean_Int = mean(interaction_score))
    colnames(int)[2] <- bc
    dumList[[bc]] <- int
    c <- c + 1
  }
  
  # Merge list elements
  interactionMat <-
    dumList %>% reduce(full_join, by = "Receptor_Ligand") %>%
    column_to_rownames('Receptor_Ligand')
  
  # Get rid of NAs
  interactionMat[is.na(interactionMat)] <- 0
  
  return(interactionMat)
  
}