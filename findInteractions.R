# Calculate interaction scores of receptors-ligand pairs in database betweeen adjacent spots

# Create function
calculateInteractions <- function(neighboursList,
                                  dat,
                                  database,
                                  filter,
                                  threshold) {
  
  # Create vector containing the barcodes that are in the neighbors list and in the data (i.e. actual tissue spots)
  message('Creating vector of barcodes in neighbours list and in data object.')
  select <- names(neighboursList)[names(neighboursList) %in% dat@assays$Spatial@counts@Dimnames[[2]]]

  # Get the scaled data in a dataframe
  message('Extracting scaled data.')
  scaledat <- as.data.frame(GetAssayData(object = dat, slot = "data"))
  
  # Filter the scaled data down to all genes in the database
  message('Filtering scaled data down to all genes in the database.')
  
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
  message('Processing barcodes...')
  
  # Create dummy dataframe
  dummy <- data.frame()
  
  # Iterate through each of the barcodes in the select vector created earlier
  # Store the primary barcode and the neighbour barcodes in objects
  # Ensure to add barcode1 to neighbors so it checks interaction within itself
  # Ensure that the neighbor barcodes are in the select vector!
  for (y in 1:length(select)) {
  
    bcode1 <- select[y]
    nbours <-
      c(neighboursList[[bcode1]]$neighbours[neighboursList[[bcode1]]$neighbours %in% select],
        bcode1)
  
    message(paste0(
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
              
              dummy <- rbind(dummy, new_entry)
              
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
  message('Calculating interaction scores.')
  dummy <-
    dummy %>%
    rowwise() %>%
    mutate(interaction_score = mean(c_across(
      c(spot1_complex_expr, spot2_ligand_expr)
    ))) %>%
    filter(!is.infinite(interaction_score))
  
  # Find mean and standard deviation of interaction scores
  message('Calculating mean and STDEV of interaction scores.')
  mean_intscore <- mean(dummy$interaction_score)
  sd_intscore <- sd(dummy$interaction_score)
  
  # If filter == TRUE - filter to rows with interactions 'x' STDEVs greater than mean
  if (filter == TRUE) {
    
    message('Filtering to interactions with scores 2 STDEVs greater than mean.')
    dummy <- dummy[dummy$interaction_score > (mean_intscore+(threshold*sd_intscore)), ]
    
    # Sort by interaction score
    results <- dummy %>%
      arrange(desc(interaction_score))
    
  } else if (filter == FALSE) {
    
    message('Not filtering interaction scores.')
    
    # Sort by interaction score
    results <- dummy %>%
      arrange(desc(interaction_score))
    
  }

  # Gather and return results
  message('Returning results.')
  
  # Count occurrences of different complexes in final results
  complex_counts <- results %>%
    count(spot1_complex) %>%
    arrange(desc(n))
  
  # Create list of results and complex counts
  interactionResults <- list(Interactions = results,
                             Complex_Summary = complex_counts)
  
  # Return interactionResults
  return(interactionResults)
  
}
