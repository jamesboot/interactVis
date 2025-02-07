# Script for plotting top interaction results

# Packages
library(ggplot2)
library(dplyr)

# Working dir
setwd('/Users/bootj/Documents/Milner_et_al')

# Locate files
files <- list.files(path = 'interactVis_new_outs_v2',
                    pattern = '*.csv',
                    recursive = F,
                    full.names = T)

# Create a name for each file
fileNames <- gsub('_wilcox.csv', '', basename(files))

# Create list for processing
inputList <- mapply(function(x, y) list(x, y), fileNames, files, SIMPLIFY = F)

# Loop through files and plot
test <- lapply(inputList, function(x){
  res <- read.csv(x[[2]])
  res <- res %>%
    arrange(FDR) %>%
    mutate('-log(FDR)' = -log(FDR)) %>%
    slice_head(n = 25) %>%
    mutate(`-log(FDR)` = replace(`-log(FDR)`, `-log(FDR)` == Inf, -log(1e-300)))
  
  res$Receptor.Ligand <- sub('^-', '', res$Receptor.Ligand)
  res$Receptor.Ligand <- factor(res$Receptor.Ligand,
                                levels = rev(res$Receptor.Ligand))
  
  plt <- ggplot(res, aes(x = Receptor.Ligand, y = `-log(FDR)`)) +
    geom_bar(stat = "identity", colour = 'black', fill = 'darkorange') +
    xlab('Interaction Name') +
    theme_bw() +
    coord_flip()
  
  ggsave(plot = plt,
         filename = paste0('interactVis_new_outs_v2/', x[[1]], '.pdf'),
         height = 7,
         width = 10)
    
  return(res)
})
 
# Tom wants a manual plot:
# Could you add 'Signaling_by_Opioid-OPRD1-PENK' to the AllSamples_Neuronal_SENDER bar graph?
x <- inputList[['AllSamples_Neuronal_SENDER']]
res1 <- read.csv(x[[2]])
res2 <- res1 %>%
  arrange(FDR) %>%
  mutate('-log(FDR)' = -log(FDR)) %>%
  slice_head(n = 25) %>%
  mutate(`-log(FDR)` = replace(`-log(FDR)`, `-log(FDR)` == Inf, -log(1e-300)))

res3 <- res1 %>%
  mutate('-log(FDR)' = -log(FDR)) %>%
  filter(Receptor.Ligand == 'Signaling_by_Opioid-OPRD1-PENK') %>%
  mutate(`-log(FDR)` = replace(`-log(FDR)`, `-log(FDR)` == Inf, -log(1e-300)))

res4 <- rbind(res2, res3)

res4$Receptor.Ligand <- sub('^-', '', res4$Receptor.Ligand)
res4$Receptor.Ligand <- factor(res4$Receptor.Ligand,
                              levels = rev(res4$Receptor.Ligand))

plt <- ggplot(res4, aes(x = Receptor.Ligand, y = `-log(FDR)`)) +
  geom_bar(stat = "identity", colour = 'black', fill = 'darkorange') +
  xlab('Interaction Name') +
  theme_bw() +
  coord_flip()

ggsave(plot = plt,
       filename = paste0('interactVis_new_outs_v2/', x[[1]], '_v2.pdf'),
       height = 7,
       width = 10)












