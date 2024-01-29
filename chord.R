# Script to trial making chord plots for Tom

# Packages
library(circlize)

# Set working directory
setwd('/Users/jamesboot/Documents/9.Genome Centre Files/GC-TM-10271/interactVis_trial/AllSamplesTrial/')

# Load sample results
res <- read.csv('NH12-297/NH12-297_nInt_cluster.csv')[, 2:4]

# Try plot
tiff(filename = 'Chord.tiff',
     height = 4,
     width = 4,
     units = 'in',
     res = 300)
circos.clear()
circos.par(gap.after = c(rep(5, length(unique(res[[1]])))))
chordDiagram(res, transparency = 0)
dev.off()