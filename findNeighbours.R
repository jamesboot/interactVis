# Identify barcodes of spots adjacent to each spot from tissue_positions file output from spaceranger

# Create function
findNeighbours <- function(file) {
  
  message('Finding neighbours...')
  
  # Read in CSV
  tissue_pos <- read.csv(file)
  
  # Make function to find neighbours apply
  f <- function(x, list) {
    bcode_lkup <- x[1]
    row <- as.numeric(x[3])
    col <- as.numeric(x[4])
    
    adjacent_bcodes <-
      tissue_pos[abs(row - tissue_pos$array_row) + abs(col - tissue_pos$array_col) <= 2 &
                   abs(row - tissue_pos$array_row) + abs(col - tissue_pos$array_col) > 0,]
    
    list[[bcode_lkup]] <- list(primaryBcode = bcode_lkup,
                               neighbours = adjacent_bcodes$barcode)
    
  }
  
  # Apply function over rows (1) using apply
  output_list <- list()
  neighbours <- apply(tissue_pos, 1, FUN = f, list = output_list)
  
  # Name each element in the list by the barcode
  names(neighbours) <- tissue_pos$barcode
  
  # Return the neighbours list
  message('Done!')
  return(neighbours)
  
}
