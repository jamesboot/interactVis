# Function to perform multiple t-tests to compare interaction scores between clusters

# Create function
differentialInteraction <- function(InteractionMat,
                                    Meta,
                                    Attribute,
                                    Comparison) {
  
  # Transpose the interaction mat
  interactionDF <- as.data.frame(t(InteractionMat))
  
  # Add cell as column in meta and count mat
  interactionDF <- interactionDF %>%
    rownames_to_column(var = 'cell')
  diffIntMeta2 <- Meta %>%
    rownames_to_column(var = 'cell')
  
  # Now join meta and column
  tTestDF <- interactionDF %>%
    left_join(diffIntMeta2, by = 'cell')
  
  # Subset down to clusters we want to compare
  tTestDFFilt <- tTestDF %>%
    filter(!!as.symbol(Attribute) %in% Comparison)
  
  # Make dataframe for results
  wilcoxRes <- data.frame()
  
  # For loop to perform test on all columns
  for (x in c(2:(ncol(tTestDFFilt) - ncol(Meta)))) {
    
    # As data is not normally distributed - use wilcoxon
    res <-
      wilcox.test(tTestDFFilt[, x] ~ tTestDFFilt[, Attribute],
                  data = tTestDFFilt,
                  paired = F)
    
    # Adjust the p-value
    FDR <-
      p.adjust(res$p.value, method = 'fdr', n = length(c(2:(
        ncol(tTestDFFilt) - ncol(Meta)
      ))))
    
    # Put result into dataframe
    wilcoxRes[colnames(tTestDFFilt)[x], 'Receptor-Ligand'] <-
      colnames(tTestDFFilt)[x]
    wilcoxRes[colnames(tTestDFFilt)[x], 'P-Value'] <- res$p.value
    wilcoxRes[colnames(tTestDFFilt)[x], 'FDR'] <- FDR
    
  }
  return(wilcoxRes)
}
