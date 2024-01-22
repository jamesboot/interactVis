# Try creating DGEList from edgeR
y <- DGEList(counts = interactionMat,
             group = diffIntMeta$Cluster)

# Estimate dispersion
y <- estimateDisp(y)

# Design matrix
design <- model.matrix(~0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)

# Fit model
fit <- glmQLFit(y, design)

# Make comparisons 
myContrast <- makeContrasts(Two-Ten,
                            Two-Eight,
                            Eight-Ten,
                            Four-Twelve,
                            Seven-Twelve,
                            Four-Seven,
                            levels = design)

# For loop to go through all comparisons and perform tests
for (x in colnames(myContrast)) {
  
  message(paste0('Starting comparison: ', x))
  results <- glmQLFTest(fit, contrast = myContrast[, x])
  filtRes <- topTags(results, n = 1000, sort.by = 'PValue', p.value = 1)
  write.csv(filtRes, file = paste0(x, '_R-L_Interactions.csv'))
  
} 
