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
  message('Creating receiver matrix...')
  
  # Use dplyr to pivot longer
  tmp <- AnnoInt_Comb %>%
    dplyr::select(Receiver_bcode, Receiver_expression, R_L_Name) %>%
    pivot_wider(
      names_from = Receiver_bcode,
      values_from = Receiver_expression,
      values_fn = function(x)
        mean(x)
    ) %>%
    column_to_rownames('R_L_Name')
  
  # Convert to matrix
  interactionMat <- data.matrix(tmp)
  
  # Get rid of NAs
  interactionMat[is.na(interactionMat)] <- 0
  
  # 2. ADD RECEIVER MATRIX TO LIST: ----
  IntMatrices$ReceiverMat <- interactionMat
  
  # 3. CREATE SENDER MATRIX: ----
  message('Creating sender matrix...')
  
  # Use dplyr to pivot longer
  tmp <- AnnoInt_Comb %>%
    dplyr::select(Sender_bcode, Sender_expression, R_L_Name) %>%
    pivot_wider(
      names_from = Sender_bcode,
      values_from = Sender_expression,
      values_fn = function(x)
        mean(x)
    ) %>%
    column_to_rownames('R_L_Name')
  
  # Convert to matrix
  interactionMat <- data.matrix(tmp)
  
  # Get rid of NAs
  interactionMat[is.na(interactionMat)] <- 0
  
  # 4. ADD RECEIVER MATRIX TO LIST:
  IntMatrices$SenderMat <- interactionMat
  
  return(IntMatrices)
  
}