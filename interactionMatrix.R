# Function to make matrices:
# 1. RECEIVER MATRIX: Matrix with spot1 as columns, receptor-ligand as rows, interaction scores as values
# Spot 1 is the receiver (where the receptor expression is taken from)
# 2. SENDER MATRIX: Matrix with spot2 as columns, receptor-ligand as rows, interaction scores as values
# Spot 2 is the sender (where the ligand expression is taken from)

# Create function
interactionMatrix <- function(AnnoInt) {
  
  # Store the 2 matrices in a list
  IntMatrices <- list()
  
  # 1. CREATE RECEIVER MATRIX: ----
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
  
  # 2. ADD RECEIVER MATRIX TO LIST: ----
  IntMatrices$ReceiverMat <- interactionMat
  
  # 3. CREATE SENDER MATRIX: ----
  
  # Make a dummy list
  dumList <- list()
  
  # Create new receptor ligand column in AnnoInt
  AnnoInt <- AnnoInt %>%
    mutate(Receptor_Ligand = paste0(spot1_complex, '-', spot2_ligand))
  
  # For loop to go through each cell
  # Use dplyr to extract and summarise interaction scores for all receptor-ligands
  # Put the new dataframe for the spot into dummy list
  c <- 1
  for (bc in unique(AnnoInt$spot2)) {
    message(paste0(
      'Starting spot ',
      c,
      ' of ',
      length(unique(AnnoInt$spot2)),
      '. Progress: ',
      round((c / length(
        unique(AnnoInt$spot2)
      ) * 100), digits = 1),
      '%'
    ))
    int <- AnnoInt %>%
      filter(spot2 == bc) %>%
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
  
  # 4. ADD RECEIVER MATRIX TO LIST:
  IntMatrices$SenderMat <- interactionMat
  
  return(IntMatrices)
  
}