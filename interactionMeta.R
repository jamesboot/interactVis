# Create differential interaction meta data dataframe from Seurat object meta data
# For use in statistical testing

# Create function
createMetaData <- function(SeuratObj,
                           InteractionMatList,
                           Attributes) {
  
  # Make a list for two sets of meta data to go in
  metaList <- list()
  
  # Create RECEIVER meta data 
  metaList$ReceiverMeta <- SeuratObj@meta.data[colnames(InteractionMatList$ReceiverMat), c(Attributes)]
  colnames(metaList$ReceiverMeta) <- c(Attributes)
  
  # Create SENDER meta data
  metaList$SenderMeta <- SeuratObj@meta.data[colnames(InteractionMatList$SenderMat), c(Attributes)]
  colnames(metaList$SenderMeta) <- c(Attributes)
  
  return(metaList)
}