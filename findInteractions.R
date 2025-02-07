# Calculate interaction scores of receptors-ligand pairs in database betweeen adjacent spots

# Create function
calculateInteractions <- function(neighbours,
                                  dat,
                                  database,
                                  filter,
                                  threshold) {
  
  # Generate all possible gene combinations from database
  separator <- "-"
  database2 <- lapply(database, function(x){
    vector1 <- x$Complex_Genes
    vector2 <- x$Partner_Genes
    result <- as.vector(outer(vector1, vector2, FUN = function(x, y) paste(x, y, sep = separator)))
    return(result)
  })
  
  # Create vector of valid interactions
  validInteractions <- setNames(unlist(database2, use.names=F),rep(names(database2), lengths(database2)))
  
  # Create vector containing the barcodes that are in the neighbors list and in the data (i.e. actual tissue spots)
  message('Creating vector of barcodes in neighbours list and in data object.')
  select <- names(neighbours)[names(neighbours) %in% colnames(dat)]

  # Get the scaled data in a dataframe
  message('Extracting scaled data.')
  scaledat <- as.data.frame(GetAssayData(object = dat, layer = "data"))
  
  # Filter the scaled data down to all genes in the database
  message('Filtering scaled data down to all genes in the database.')
  
  # Create a vector of all the genes in the database 
  db_elements <- unique(as.vector(unlist(database)))
  
  # Go through each gene and select it from scaledat, if not present insert 0s
  scaledat_filt <- scaledat[0,]
  message('Subsetting gene expression matrix:')
  pb <- txtProgressBar(min = 1, max = length(db_elements), style = 3)
  for (x in 1:length(db_elements)) {
    if (db_elements[x] %in% row.names(scaledat)) {
      scaledat_filt[db_elements[x], ] <- scaledat[db_elements[x], ]
    } else {
      scaledat_filt[db_elements[x], ] <- 0
    }
    setTxtProgressBar(pb, x)
  }
  close(pb)
  
  scaledat_filt$Gene <- rownames(scaledat_filt)
  
  # Create long dataframe for later
  scaledat_filt_long <- pivot_longer(scaledat_filt,
                                     cols = 1:(ncol(scaledat_filt)-1),
                                     names_to = 'barcode',
                                     values_to = 'expression')

  # Now calculate the interaction scores ----
  message('Processing barcodes...')
  
  # Define the number of CPUs to use
  numCores <- detectCores() - 1  # Use all but one core
  
  # Set up a parallel cluster
  cl <- makeCluster(numCores)
  clusterEvalQ(cl, {
    library(tidyverse)
  })
  clusterExport(cl, varlist = c('neighbours',
                                'dat',
                                'database'))
  
  # Create list to iterate over
  bcIter <- 1:length(select)

  # Iterate through each of the barcodes in the select vector created earlier
  # Store the primary barcode and the neighbour barcodes in objects
  # Ensure to add barcode1 to neighbors so it checks interaction within itself
  # Ensure that the neighbor barcodes are in the select vector!
  results <- pblapply(cl = cl, bcIter, function(y) {
    bcode1 <- select[y]
    neighbours <-
      c(neighbours[[bcode1]]$neighbours[neighbours[[bcode1]]$neighbours %in% select], bcode1)
  
    # Create a dataframe of all sender barcodes and ligand genes
    sender <- scaledat_filt_long
    colnames(sender) <- c('Sender_gene', 'Sender_bcode', 'Sender_expression')
    ligandGenes <- unlist(lapply(database, function(x)
      x[[2]]))
    sender <- sender %>%
      filter(Sender_gene %in% ligandGenes) %>%
      filter(Sender_bcode == bcode1) #%>%
    #filter(Sender_expression != 0)
    
    # Create a dataframe of all neighbour barcodes and complex genes
    receiver <- scaledat_filt_long
    colnames(receiver) <- c('Receiver_gene', 'Receiver_bcode', 'Receiver_expression')
    complexGenes <- unlist(lapply(database, function(x)
      x[[1]]))
    receiver <- receiver %>%
      filter(Receiver_gene %in% complexGenes) %>%
      filter(Receiver_bcode %in% neighbours) #%>%
    #filter(Receiver_expression != 0)
    
    # Cross the dataframes
    merged_data <- crossing(sender, receiver)
    
    # Add the interaction name to merged data
    merged_data <- merged_data %>%
      mutate(interaction_name = paste(Receiver_gene, Sender_gene, sep = '-'))

    # Filter merged data to valid interactions
    merged_data <- merged_data %>%
      filter(interaction_name %in% validInteractions)
    
    # Now for the barcodes selected iterate through the database
    
    # LAPPLY METHOD
    test <- lapply(seq_along(database), function(x) {
      complexName <- names(database)[[x]]
      complex <- database[[x]]$Complex_Genes
      ligand <- database[[x]]$Partner_Genes
      
      # Filter dataframe to these genes
      # Add R_L_Name col
      # Make a unique column which is:
      # Complex + Receiver bc + ligand + sender bc
      # This is to select the minimum expression of a complex in the receiver by grouping on unique name
      merged_data %>%
        filter(Sender_gene %in% ligand &
                 Receiver_gene %in% complex) %>%
        mutate(
          R_L_Name = complexName,
          Unique_Name = paste(R_L_Name, Receiver_bcode, #Sender_gene,
                              Sender_bcode, sep = '_')
        ) %>%
        group_by(Unique_Name) %>%
        filter(Receiver_expression == min(Receiver_expression)) %>%
        filter(Sender_expression == min(Sender_expression)) %>%
        arrange(Receiver_expression, Sender_expression) %>%
        filter(row_number() == 1) %>%
        filter(Sender_expression != 0 & Receiver_expression != 0) %>%
        ungroup() %>%
        dplyr::select(!Unique_Name)
    })
    
    # Update results
    return(test)
  })
  
  # Stop the cluster
  stopCluster(cl)

  # Combine the list into a single data.frame
  firstBind <- lapply(results, function(x){
    do.call(rbind, x)
  })
  secondBind <- do.call(rbind, firstBind)

  # Calculate interaction score
  message('Calculating interaction scores.')
  results <-
    secondBind %>%
    rowwise() %>%
    mutate(interaction_score = mean(c_across(
      c(Sender_expression, Receiver_expression)
    ))) %>%
    filter(!is.infinite(interaction_score))
  
  # Find mean and standard deviation of interaction scores
  message('Calculating mean and STDEV of interaction scores.')
  mean_intscore <- mean(results$interaction_score)
  sd_intscore <- sd(results$interaction_score)
  
  # If filter == TRUE - filter to rows with interactions 'x' STDEVs greater than mean
  if (filter == TRUE) {
    
    message('Filtering to interactions with scores 2 STDEVs greater than mean.')
    results <- results[results$interaction_score > (mean_intscore+(threshold*sd_intscore)), ]
    
    # Sort by interaction score
    results <- results %>%
      arrange(desc(interaction_score))
    
  } else if (filter == FALSE) {
    
    message('Not filtering interaction scores.')
    
    # Sort by interaction score
    results <- results %>%
      arrange(desc(interaction_score))
    
  }

  # Gather and return results
  message('Returning results.')
  
  # Count occurrences of different complexes in final results
  complex_counts <- results %>%
    count(R_L_Name) %>%
    arrange(desc(n))
  
  # Create list of results and complex counts
  interactionResults <- list(Interactions = results,
                             Complex_Summary = complex_counts)
  
  # Return interactionResults
  return(interactionResults)
  
}
