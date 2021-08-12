# Domain Swap Dimer
# Aleix Lafita - June 2020

source("utils.R")

######################## Parameters #############################

# File names
input = "e1shgA1.pdb"
output = "e1shgA1_dsd33.pdb"

# Parameters for DSD
index = 33 # residue index to do the swapping hinge loop
hinge = 1 # residues to unfold in the hinge loop

# Granularity: calpha, backbone
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
revindex = resnum - index

# The atom details of the second domain - just change the chain and eleno
pdb.atoms.2 = pdb.atoms %>%
  mutate(
    eleno = eleno + atomnum,
    chain = "B"
  )

# Bind the two domains of the dimer - assign a unique resno
pdb.atoms.dimer = rbind(pdb.atoms, pdb.atoms.2) %>%
  mutate(g.resno = ifelse(chain == "B", resnum + resno, resno))

######################## Distance matrix #############################
message("## Calculating distance matrix...")

# Compute inter-atomic distances
distmat.dist = dist(pdb.atoms.dimer %>% select(x, y, z), diag = T, upper = T)

# Add eleno indices to matrix
distmat.df = data.frame(
  d = matrix(as.matrix(distmat.dist), ncol = 1),
  eleno.x = 1:(atomnum*2)
) %>% mutate(
  eleno.y = as.integer(0:(length(d)-1) / (atomnum*2)) + 1
)

# Merge with atom table to include information
distmat.df = merge(distmat.df, pdb.atoms.dimer, by.x = "eleno.x", by.y = "eleno")
distmat.df = merge(distmat.df, pdb.atoms.dimer, by.x = "eleno.y", by.y = "eleno")

# Sort the dataframe
distmat.df = distmat.df[order(distmat.df$eleno.x, distmat.df$eleno.y),]

# DM transformation for a domain swap dimer
distmat.df.swap = distmat.df %>%
  mutate(
    d.na = d,
    # Delete inter-domain distances from non-swap - top, bottom, right and left
    d.na = ifelse(is.element(g.resno.x, (index+1):(resnum+index-1)) & is.element(g.resno.y, 1:(index-1)), NA, d.na),
    d.na = ifelse(is.element(g.resno.y, (index+1):(resnum+index-1)) & is.element(g.resno.x, 1:(index-1)), NA, d.na),
    d.na = ifelse(is.element(g.resno.x, (index+1):(resnum+index-1)) & is.element(g.resno.y, (resnum+index):(2*resnum)), NA, d.na),
    d.na = ifelse(is.element(g.resno.y, (index+1):(resnum+index-1)) & is.element(g.resno.x, (resnum+index):(2*resnum)), NA, d.na),
    # Unfold the hinge loop - residue at index and surrounding
    d.na = ifelse(is.element(g.resno.x, c((index-hinge):(index+hinge), (resnum+index-hinge):(resnum+index+hinge))) & g.resno.x != g.resno.y, NA, d.na),
    d.na = ifelse(is.element(g.resno.y, c((index-hinge):(index+hinge), (resnum+index-hinge):(resnum+index+hinge))) & g.resno.x != g.resno.y, NA, d.na)
  )

distmat.df.bounded = apply_bounds(distmat.df.swap)


# Matrices of upper and lower bounds
D = matrix(distmat.df.bounded$d.na, nrow = atomnum*2)
L = matrix(distmat.df.bounded$d.lower, nrow = atomnum*2)
U = matrix(distmat.df.bounded$d.upper, nrow = atomnum*2)

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
data.pdb.mds$atom = pdb.atoms.dimer
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


