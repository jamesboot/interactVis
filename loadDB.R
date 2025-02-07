# Function to read in the database from cellphonedb ready for analysis of adjacent spots
loadDB <- function(interaction.input) {
  
  message('Loading database...')
  
  # Create blank list
  cellphonedb <- list()
  
  # We need to load the interaction input file
  # This contains the complex names and the proteins that interact with it
  interact_lookup <-
    read.csv(
      interaction.input
    )
  
  # Remove interactions with multiple '-'
  matches <- grepl("^([^\\-]*\\-){2}[^\\-]*$", interact_lookup$interactors)
  interact_lookup <- interact_lookup[!matches, ]
  
  # Filter the interact lookup to only ligand-receptors
  # Split the interactors column based on '-' character: before = ligand, after = receptor
  interact_lookup <- interact_lookup %>%
    filter(directionality == 'Ligand-Receptor') %>%
    separate(interactors, c('Ligand', 'Receptor'), sep = '-') #%>%
    #filter(is_ppi == 'True')
  
  # Iterate through rows of interact table to update database list
  for (x in 1:nrow(interact_lookup)) {
    complexGenes <- unlist(strsplit(interact_lookup$Receptor[x], split = '\\+'))
    ligandGenes <- unlist(strsplit(interact_lookup$Ligand[x], split = '\\+'))
    complexName <- paste(
      gsub(' ', '_', interact_lookup$classification[x]),
      paste(complexGenes, collapse = '_'),
      paste(ligandGenes, collapse = '_'),
      sep = '-'
    )
    cellphonedb[[complexName]] <- list(Complex_Genes = complexGenes, Partner_Genes = ligandGenes)
  }
  
  message('Done!')
  
  return(cellphonedb)
  
}
