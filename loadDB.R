# Function to read in the database from cellphonedb ready for analysis of adjacent spots
loadDB <- function(gene.input,
                   complex.input,
                   interaction.input) {
  
  message('Loading database...')
  
  # Create blank list
  cellphonedb <- list()
  
  # We need to load the gene input all file
  # This allows us to look up gene names and find their uniprot ID
  gene_lookup <-
    read.csv(
      '/Users/jamesboot/Documents/9.Genome Centre Files/interactVis/cellphonedb-data-4.0.0/data/gene_input_all.csv'
    )
  
  # We need to load the complex input file
  # This contains all info of the complexes - what proteins are involved
  complex_lookup <-
    read.csv(
      '/Users/jamesboot/Documents/9.Genome Centre Files/interactVis/cellphonedb-data-4.0.0/data/complex_input.csv'
    )
  # Make the complex name the row names
  row.names(complex_lookup) <- complex_lookup$complex_name
  
  # We need to load the interaction input file
  # This contains the complex names and the proteins that interact with it
  interact_lookup <-
    read.csv(
      '/Users/jamesboot/Documents/9.Genome Centre Files/interactVis/cellphonedb-data-4.0.0/data/interaction_input.csv'
    )
  
  # Iterate through all complexes and update database list along to way
  for (x in 1:nrow(complex_lookup)) {
    # First: given a complex - return all the gene names of the proteins in the complex
    # Get all the protein ids of proteins in the complex
    complex_name <- complex_lookup$complex_name[x]
    complex_prot <- c(
      complex_lookup[complex_name,]$uniprot_1,
      complex_lookup[complex_name,]$uniprot_2,
      complex_lookup[complex_name,]$uniprot_3,
      complex_lookup[complex_name,]$uniprot_4
    )
    # Remove blanks
    complex_prot <- complex_prot[complex_prot != '']
    # Now find the gene name
    complex_genes <-
      unique(gene_lookup$gene_name[gene_lookup$uniprot %in% complex_prot])
    
    # Second: given a complex - return the gene names of the proteins it interacts with (partner)
    # Get the protein id of the partner
    partner_prot <-
      interact_lookup$partner_b[interact_lookup$partner_a == complex_name]
    # Now find the gene names for partners
    partner_genes <-
      unique(gene_lookup$gene_name[gene_lookup$uniprot %in% partner_prot])
    
    # Add to the list (list is the database) - if length is greater than 0
    if (length(partner_genes) > 0) {
      
      cellphonedb[[complex_name]] <- list(Complex_Genes = complex_genes,
                                          Partner_Genes = partner_genes)
      
    } else {
      next
    }
  }
  
  message('Done!')
  
  return(cellphonedb)
  
}
