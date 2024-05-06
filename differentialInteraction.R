# Function to perform multiple t-tests to compare interaction scores between clusters

# Create function
differentialInteraction <- function(InteractionMatList,
                                    MetaList,
                                    Attribute,
                                    Comparison,
                                    Test) {
  
  # Make a list for results to go in:
  diffIntRes <- list()
  
  # Perform for RECEIVERS first ----
  message('RECEIVERS first...')
  
  # Transpose the interaction mat
  interactionDF <- as.data.frame(t(InteractionMatList$ReceiverMat))
  
  # Add cell as column in meta and count mat
  interactionDF <- interactionDF %>%
    rownames_to_column(var = 'cell')
  diffIntMeta2 <- MetaList$ReceiverMeta %>%
    rownames_to_column(var = 'cell')
  
  # Now join meta and column
  tTestDF <- interactionDF %>%
    left_join(diffIntMeta2, by = 'cell')
  
  # Subset down to clusters we want to compare
  tTestDFFilt <- tTestDF %>%
    filter(!!as.symbol(Attribute) %in% Comparison)
  
  # Check if there are enough replicates of groups
  if (sum(tTestDFFilt$Cluster == Comparison[1]) > 3 &
      sum(tTestDFFilt$Cluster == Comparison[2]) > 3) {
    message('Starting statistical test for RECEIVERS.')
    
    if (Test == 'Wilcox') {
      
      # Make dataframe for results
      wilcoxRes <- data.frame()
      
      # For loop to perform test on all columns
      for (x in c(2:(ncol(tTestDFFilt) - ncol(MetaList$ReceiverMeta)))) {
        
        # As data is not normally distributed - use wilcoxon
        res <-
          wilcox.test(tTestDFFilt[, x] ~ tTestDFFilt[, Attribute],
                      data = tTestDFFilt,
                      paired = F)
        
        # Adjust the p-value
        FDR <-
          p.adjust(res$p.value, method = 'fdr', n = length(c(2:(
            ncol(tTestDFFilt) - ncol(MetaList$ReceiverMeta)
          ))))
        
        # Put result into dataframe
        wilcoxRes[colnames(tTestDFFilt)[x], 'Receptor-Ligand'] <-
          colnames(tTestDFFilt)[x]
        wilcoxRes[colnames(tTestDFFilt)[x], 'P-Value'] <- res$p.value
        wilcoxRes[colnames(tTestDFFilt)[x], 'FDR'] <- FDR
        
      }
      
      diffIntRes$ReceiverResults <- wilcoxRes
      
    } else if (Test == 'glm') {
      
      # Make dataframe for results
      glmRes <- data.frame()
      
      # For loop to perform test on all columns
      for (x in c(2:(ncol(tTestDFFilt) - ncol(MetaList$ReceiverMeta)))) {
        
        if (sum(tTestDFFilt[, x]) == 0) {
          
          next
          
        } else {
          
          # Fit null model
          fit0 <- glm.nb(tTestDFFilt[, x] ~ 1,
                         data = tTestDFFilt)
          
          # Fit proposed model
          fit1 <- glm.nb(tTestDFFilt[, x] ~ tTestDFFilt[, Attribute],
                         data = tTestDFFilt)
          
          # Compare fit of models
          lrtResult <- lrtest(fit1, fit0)
          
          # Return p-value
          p <- lrtResult$`Pr(>Chisq)`[2]
          
          # Adjust the p-value
          FDR <-
            p.adjust(p, method = 'fdr', n = length(c(2:(
              ncol(tTestDFFilt) - ncol(MetaList$ReceiverMeta)
            ))))
          
          # Put result into dataframe
          glmRes[colnames(tTestDFFilt)[x], 'Receptor-Ligand'] <-
            colnames(tTestDFFilt)[x]
          glmRes[colnames(tTestDFFilt)[x], 'P-Value'] <- p
          glmRes[colnames(tTestDFFilt)[x], 'FDR'] <- FDR
          
        }
        
      }
      
      diffIntRes$ReceiverResults <- glmRes
      
    }
    
  } else {
    message('Not enough replicates in groups for statistical test. Skipping RECEIVERS.')
  }
  
  # Now perform for SENDERS ----
  message('SENDERS second...')
  
  # Transpose the interaction mat
  interactionDF <- as.data.frame(t(InteractionMatList$SenderMat))
  
  # Add cell as column in meta and count mat
  interactionDF <- interactionDF %>%
    rownames_to_column(var = 'cell')
  diffIntMeta2 <- MetaList$SenderMeta %>%
    rownames_to_column(var = 'cell')
  
  # Now join meta and column
  tTestDF <- interactionDF %>%
    left_join(diffIntMeta2, by = 'cell')
  
  # Subset down to clusters we want to compare
  tTestDFFilt <- tTestDF %>%
    filter(!!as.symbol(Attribute) %in% Comparison)
  
  # Check if there are enough replicates of groups
  if (sum(tTestDFFilt$Cluster == Comparison[1]) > 3 &
      sum(tTestDFFilt$Cluster == Comparison[2]) > 3) {
    
    message('Starting statistical test for SENDERS.')
    
    if (Test == 'Wilcox') {
      
      # Make dataframe for results
      wilcoxRes <- data.frame()
      
      # For loop to perform test on all columns
      for (x in c(2:(ncol(tTestDFFilt) - ncol(MetaList$SenderMeta)))) {
        
        # As data is not normally distributed - use wilcoxon
        res <-
          wilcox.test(tTestDFFilt[, x] ~ tTestDFFilt[, Attribute],
                      data = tTestDFFilt,
                      paired = F)
        
        # Adjust the p-value
        FDR <-
          p.adjust(res$p.value, method = 'fdr', n = length(c(2:(
            ncol(tTestDFFilt) - ncol(MetaList$SenderMeta)
          ))))
        
        # Put result into dataframe
        wilcoxRes[colnames(tTestDFFilt)[x], 'Receptor-Ligand'] <-
          colnames(tTestDFFilt)[x]
        wilcoxRes[colnames(tTestDFFilt)[x], 'P-Value'] <- res$p.value
        wilcoxRes[colnames(tTestDFFilt)[x], 'FDR'] <- FDR
        
      }
      
      diffIntRes$SenderResults <- wilcoxRes
      
    } else if (Test == 'glm') {
      
      # Make dataframe for results
      glmRes <- data.frame()
      
      # For loop to perform test on all columns
      for (x in c(2:(ncol(tTestDFFilt) - ncol(MetaList$SenderMeta)))) {
        
        if (sum(tTestDFFilt[, x]) == 0) {
          
          next
          
        } else {
          
          # Fit null model
          fit0 <- glm.nb(tTestDFFilt[, x] ~ 1,
                         data = tTestDFFilt)
          
          # Fit proposed model
          fit1 <- glm.nb(tTestDFFilt[, x] ~ tTestDFFilt[, Attribute],
                         data = tTestDFFilt)
          
          # Compare fit of models
          lrtResult <- lrtest(fit1, fit0)
          
          # Return p-value
          p <- lrtResult$`Pr(>Chisq)`[2]
          
          # Adjust the p-value
          FDR <-
            p.adjust(p, method = 'fdr', n = length(c(2:(
              ncol(tTestDFFilt) - ncol(MetaList$SenderMeta)
            ))))
          
          # Put result into dataframe
          glmRes[colnames(tTestDFFilt)[x], 'Receptor-Ligand'] <-
            colnames(tTestDFFilt)[x]
          glmRes[colnames(tTestDFFilt)[x], 'P-Value'] <- p
          glmRes[colnames(tTestDFFilt)[x], 'FDR'] <- FDR
          
        }
        
      }
      
      diffIntRes$SenderResults <- glmRes
      
    }
    
  } else {
    message('Not enough replicates in groups for statistical test. Skipping SENDERS.')
  }
  
  # Now return list
  return(diffIntRes)
}
