### In house function to calculate scores of Ligand-Receptor interaction.
### Input is an AverageExpression data frame calculate with the Seurat function (avg_exp).
### Specify cell type for Ligand expression, and cell type for Receptor expression (Ligand_source & Receptor_source).
### Specify the LT of interest by providing a list of L-R vectors (list_LR_pairs).

### Output is a data frame of 2 columns (pairs & score)

LR_table <- function(avg_exp, 
                     Ligand_source,
                     Receptor_source,
                     list_LR_pairs
) {
  #create empty data frame
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  
  #iterate through all the pairs, pick Ligand of pair (pair[1]) and multiply to Receptor of pair (pair[2])
  for (pair in list_LR_pairs){
    source <- paste0("Ligand(", Ligand_source, ") -> Receptor(", Receptor_source, ")")
    name <- paste0(pair[1], "-", pair[2])
    score <- avg_exp[pair[1], Ligand_source] * avg_exp[pair[2], Receptor_source]
    
    df <- rbind(df, c(source, name, score))
  }
  
  #rename column names
  colnames(df) <- c("Source", "Ligand -> Receptor", "Score")
  
  #return results
  return(df)
}