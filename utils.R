# Shared code and library imports for all scripts
# Aleix Lafita - June 2020

library(dplyr)
library(bio3d)
library(edmcr)

######################## Parameters #########################

options(stringsAsFactors = F)

backbone = c("N", "CA", "CB", "C", "O")
atom.options = c("calpha", "backbone")

########################### Bounds ############################

# Parse the protein distance bounds constraints
protein.bounds = read.csv(
  "protein_bounds.tsv",
  sep = "\t"
) %>% mutate(
  # Distance bounds were limited to 20A
  d.u = ifelse(d.u < 20, d.u, Inf)
)


########################### Function ############################

# Apply protein bounds to unknown entries in the matrix
apply_bounds = function(distmat.df) {
  
  # Calculate the residue number differences
  distmat.df = distmat.df %>% mutate(
    resno.diff = ifelse(
      chain.x == chain.y, 
      resno.y - resno.x, 
      # Residues from two different chains, consider them as 1000 residues apart
      sign(eleno.y - eleno.x)*1000),
    # Take care of the last residue bin being representative of any difference higher
    resno.diff.bkbn = ifelse(
      abs(resno.diff) > max(protein.bounds$resno.diff), 
      sign(resno.diff)*max(protein.bounds$resno.diff), 
      resno.diff)
  )
  
  # Apply the upper and lower bounds
  distmat.df.con = merge(
    distmat.df, 
    protein.bounds, 
    by.x = c("elety.x", "elety.y", "resno.diff.bkbn"), 
    by.y = c("elety.x", "elety.y", "resno.diff"),
    all.x = T
  ) %>% mutate(
    d.lower = ifelse(is.na(d.na), d.l, d.na),
    d.upper = ifelse(is.na(d.na), d.u, d.na),
    # Set chemical bonds between atoms to fixed
    d.na = ifelse(is.na(d.na) & d.avg < 2, d.avg, d.na),
    # Set distances with very little variation to fixed too - like peptide bond plane
    d.na = ifelse(is.na(d.na) & (d.u - d.l) < 0.2 & d.avg < 5, d.avg, d.na)
  )
  
  # Make sure it is properly sorted after the merge
  distmat.df.con = distmat.df.con[order(distmat.df.con$eleno.x, distmat.df.con$eleno.y),]
  
  return(distmat.df.con)
  
}

