# Domain Atrophy
# Aleix Lafita - June 2020

source("utils.R")

######################## Parameters #############################

# File names
input = "e1shgA1.pdb"
output = "e1shgA1_da8-17.pdb"

# Indices of residues to delete
dels = 8:17
off = 1 # offset of unfolded residues at each side of the deletion

# Granularity: calpha, backbone
atoms = "calpha"

# EDMC algorithm parameters
set.seed(0)
tol = 0.1

######################## Parsing #############################
message("## Parsing PDB structure...")

# Use bio3d to parse the PDB structure
data.pdb.raw = read.pdb2(input, multi = F, rm.alt = T)

# Clean the PDB
data.pdb = clean.pdb(
  data.pdb.raw,
  consecutive = F,
  force.renumber = T,
  fix.chain = T,
  fix.aa = T,
  rm.wat = T,
  rm.lig = T
)

# Atom indices to unfold due to the atrophy event
unfold = unique(c(sort(rep(dels,off)) + 1:off, sort(rep(dels,off)) - 1:off))

# Extract atoms
pdb.atoms = data.pdb$atom %>%
  # Take the first chain only if there are more than one
  filter(chain == "A") %>%
  # Remove residues to be deleted from the table
  filter(!is.element(resno, dels)) %>%
  # Select the granularity of atoms
  filter(
    (atoms == "calpha" & elety == "CA") |
      (atoms == "backbone" & is.element(elety, backbone))
  ) %>% mutate(
    eleno = 1:length(eleno),
    unfold = is.element(resno, unfold),
    resno = 1:length(resno),
    x = ifelse(!unfold, x, NA),
    y = ifelse(!unfold, y, NA),
    z = ifelse(!unfold, z, NA)
  )

# Total number of atoms and residues
resnum = length(unique(pdb.atoms$resno))
atomnum = nrow(pdb.atoms)

######################## Distance matrix #############################
message("## Calculating distance matrix...")

# Compute interatomic distances
distmat.dist = dist(pdb.atoms %>% select(x, y, z), diag = T, upper = T)

# Add eleno indices to matrix
distmat.df = data.frame(
  d = matrix(as.matrix(distmat.dist), ncol = 1),
  eleno.x = 1:atomnum
) %>% mutate(
  eleno.y = as.integer(0:(length(d)-1) / atomnum) + 1
)

# Merge with atom table to include information
distmat.df = merge(distmat.df, pdb.atoms, by.x = "eleno.x", by.y = "eleno")
distmat.df = merge(distmat.df, pdb.atoms, by.x = "eleno.y", by.y = "eleno")

# Sort the dataframe
distmat.df = distmat.df[order(distmat.df$eleno.x, distmat.df$eleno.y),] %>%
  mutate(d.na = d)

distmat.df.da = apply_bounds(distmat.df)

# Matrices of upper and lower bounds
D = matrix(distmat.df.da$d.n, nrow = atomnum)
L = matrix(distmat.df.da$d.lower, nrow = atomnum)
U = matrix(distmat.df.da$d.upper, nrow = atomnum)

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
data.pdb.mds$atom = pdb.atoms
data.pdb.mds$xyz = matrix(t(C), nrow = 1)

# If the structure is in the wrong chirality, invert the Z direction
if (!same_chirality(data.pdb, data.pdb.mds)) {
  message("##  Inverting model chirality...")
  C[,3] = -C[,3]
  data.pdb.mds$xyz = matrix(t(C), nrow = 1)
}

write.pdb(
  data.pdb.mds,
  output
)

message(sprintf("##  Model saved to '%s'", output))
message("## Done!")


