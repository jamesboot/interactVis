# Annotate interactions with an attribute from a Seurat object

# Create function
annotateInteractions <- function(SeuratObj,
                                 Interactions,
                                 Attribute) {
  
  # Retrieve results from Interactions 
  AnnoInt <- Interactions$Interactions
  
  # Create empty columns for attribute to go in 
  AnnoInt$Spot1_Anno <- NA
  AnnoInt$Spot2_Anno <- NA
  
  # Go through all the spots and annotate
  # Report progress
  for (x in 1:nrow(AnnoInt)) {
    message(paste0(
      'Annotating interaction ',
      x,
      ' of ',
      nrow(AnnoInt),
      '. Progress: ',
      round(((
        x / nrow(AnnoInt)
      ) * 100), digits = 1),
      '%'
    ))
    AnnoInt$Spot1_Anno[x] <-
      as.numeric(as.vector(SeuratObj@meta.data[, Attribute][colnames(SeuratObj) == AnnoInt$spot1[x]]))
    AnnoInt$Spot2_Anno[x] <-
      as.numeric(as.vector(SeuratObj@meta.data[, Attribute][colnames(SeuratObj) == AnnoInt$spot2[x]]))
  }
  
  # Ensure annotations are factors
  AnnoInt$Spot1_Anno <- as.factor(AnnoInt$Spot1_Anno)
  AnnoInt$Spot2_Anno <- as.factor(AnnoInt$Spot2_Anno)
  
  # Return the dataframe 
  return(AnnoInt)
}