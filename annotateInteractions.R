# Annotate interactions with an attribute from a Seurat object

# Create function
annotateInteractions <- function(SeuratObj,
                                 Interactions,
                                 Attribute) {
  
  # Retrieve results from Interactions 
  AnnoInt <- Interactions$Interactions
  
  # Get metadata
  meta <- data.frame(barcode = rownames(SeuratObj@meta.data),
                     Attribute = SeuratObj@meta.data[Attribute])
  
  # Annotate Senders
  AnnoInt <- AnnoInt %>%
    left_join(meta, by = c('Sender_bcode' = 'barcode'))
  colnames(AnnoInt)[colnames(AnnoInt) == Attribute] <- 'Sender_Anno'
  
  # Annotate Receivers
  AnnoInt <- AnnoInt %>%
    left_join(meta, by = c('Receiver_bcode' = 'barcode'))
  colnames(AnnoInt)[colnames(AnnoInt) == Attribute] <- 'Receiver_Anno'
  
  # Ensure annotations are factors
  AnnoInt$Sender_Anno <- as.factor(AnnoInt$Sender_Anno)
  AnnoInt$Receiver_Anno <- as.factor(AnnoInt$Receiver_Anno)
  
  # Return the dataframe 
  return(AnnoInt)
}