# Circular permutation
# Aleix Lafita - June 2020

source("utils.R")

######################## Parameters #############################

# File names
input = "e1shgA1.pdb"
output = "e1shgA1_cp33.pdb"

# Parameters for CP
index = 33 # residue index to do the circular permutant
off = 1 # length of terminal loop (offset at each side of the old termini)

# Granularity
atoms = "calpha"

# EDMC algorithm parameters
set.seed(0)
tol = 0.1

######################## Parsing #############################
message("## Parsing PDB structure...")

# Use bio3d to parse the PDB structure
data.pdb.raw = read.pdb2(input, multi = F)

# Clean the PDB
data.pdb = clean.pdb(
  data.pdb.raw,
  consecutive = F,
  force.renumber = T,
  fix.chain = T,
  fix.aa = T,
  rm.wat = T,
  rm.lig = T,
  rm.h = T
)

# Extract atoms
pdb.atoms = data.pdb$atom %>% 
  # Take the first chain only if there are more than one
  filter(chain == "A") %>%
  # Select the granularity of atoms
  filter(
    (atoms == "calpha" & elety == "CA") |
      (atoms == "backbone" & is.element(elety, backbone))
  ) %>%
  mutate(eleno = 1:length(eleno))

# Total number of atoms and residues
atomnum = nrow(pdb.atoms)
resnum = max(pdb.atoms$resno)

# Work out the indices of residues in the CP to unfold
unfold = c(1, (resnum-index-off+2):(resnum-index+off+1), resnum)

# Create a CP of the structure by changing residue numbers
pdb.atoms.cp = pdb.atoms %>%
  mutate(
    # New residue numbering
    resno = 1 + (resno + resnum - index) %% resnum,
    # Unfold the terminal linker residues by setting their coordinates to NA
    x = ifelse(is.element(resno, unfold), NA, x),
    y = ifelse(is.element(resno, unfold), NA, y),
    z = ifelse(is.element(resno, unfold), NA, z)
  )
  
pdb.atoms.cp = pdb.atoms.cp[order(pdb.atoms.cp$resno, pdb.atoms.cp$eleno),] %>%
  mutate(eleno = 1:length(eleno))

######################## Distance matrix #############################
message("## Calculating distance matrix...")

# Compute inter-atomic distances
distmat.dist = dist(pdb.atoms.cp %>% select(x, y, z), diag = T, upper = T)

# Add eleno indices to matrix
distmat.df = data.frame(
  d = matrix(as.matrix(distmat.dist), ncol = 1),
  eleno.x = 1:atomnum
) %>% mutate(
  eleno.y = as.integer(0:(length(d)-1) / atomnum) + 1
)

# Merge with atom table to include information
distmat.df = merge(distmat.df, pdb.atoms.cp, by.x = "eleno.x", by.y = "eleno")
distmat.df = merge(distmat.df, pdb.atoms.cp, by.x = "eleno.y", by.y = "eleno")

# Sort the dataframe
distmat.df = distmat.df[order(distmat.df$eleno.x, distmat.df$eleno.y),] %>%
  mutate(d.na = d)

# Apply protein bounds using the function define in utils
distmat.df.bounded = apply_bounds(distmat.df)

# Matrices of upper and lower bounds
D = matrix(distmat.df.bounded$d.na, nrow = atomnum)
L = matrix(distmat.df.bounded$d.lower, nrow = atomnum)
U = matrix(distmat.df.bounded$d.upper, nrow = atomnum)

# Use triangle inequality to restrict Inf upper bounds
U = mstUB(U)
U = ceiling(U*1000) / 1000

######################## DM Completion #############################
message("## Completing distance matrix...")

# Complete the distance matrix using EDM completion
D.edmc = edmc(D, method = "dpf", d=3, lower=L, upper=U, toler = tol)

# Extract the completed matrix
D.c = D.edmc$D

message(sprintf("##  EDMC convergence loss: %.1f", D.edmc$optval))

######################## 3D coords #############################
message("## Multidimensional scaling to 3D coords...")

# Use multidimensional scaling
fit = getConfig(D.c, d=3)
C = fit$X

######################## Save PDB file #############################

# Modify the atom table and coordinates from original
data.pdb.mds = data.pdb
data.pdb.mds$atom = pdb.atoms.cp
data.pdb.mds$xyz = matrix(t(C), nrow = 1)

write.pdb(
  data.pdb.mds,
  output
)

message(sprintf("##  Model saved to '%s'", output))
message("## Done!")

