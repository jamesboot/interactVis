# interactVis

Scripts used to perform my own neighbourhood analysis of 10X Visium spatial transcriptomics data.

The aim of the neighbourhood analysis is to find adjacent tissue spots where there is co-expression of receptors and ligands.

loadDB.R contains a function to load the database from the cellphonedb package. The database contains complex and ligand pairs.
findNeighbours.R contains a function to identify barcodes of spots adjacent to each spot from tissue_positions file output from spaceranger.
findInteractions.R contains a function that uses the loaded database and the list of tissue spot barcodes and the associated neighbouring barcodes to calculate interaction scores for every receptor ligand pair in each tissue spot and each of its neighbouring tissue spots.
