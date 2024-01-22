# Create differential interaction meta data dataframe from Seurat object meta data
# For use in statistical testing

# Create function
createMetaData <- function(SeuratObj,
                           InteractionMat,
                           Attributes) {
  Meta <- SeuratObj@meta.data[colnames(InteractionMat), c(Attributes)]
  colnames(Meta) <- c(Attributes)
  return(Meta)
}