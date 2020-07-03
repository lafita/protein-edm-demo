# Protein EDM Modelling

### Description
This repository contains code written in `R` to demonstrate a protein modelling approach based on **Euclidean Distance Matrices (EDMs)**.

It contains three `R` scripts: [cp.R](cp.R), [dsd.R](dsd.R) and [da.R](da.R), that model different structural rearrangements (circular permutation, domain swap dimerisation and domain atrophy, respectively) on the structure of an SH3 domain [e1shgA1.pdb](e1shgA1.pdb), downloaded from the [ECOD](http://prodata.swmed.edu/ecod) database.
The [utils.R](utils.R) file contains code shared by all the scripts and imports necessary libraries.

The files [ecod_domains.id](ecod_domains.id) and [protein_bounds.tsv](protein_bounds.tsv) contain the list of 20 **ECOD** domain representatives and the lookup table of distance bounds extracted from them, respectively.


### Installation

The only requirement is `R` version 3.6 with the following packages that can be installed from **CRAN**: `dplyr`, `bio3d`, and `edmcr`.
If you have problems installing the `edmcr` package, try installing `devtools` and installing it from **GitHub** using the follwing command:
```
install_github("lafita/edmcr")
```

The scripts can be run from `Rstudio` or from the command line using `Rscript` within the same directory:
```
Rscript cp.R
```

### Article

Article currently in preparation.

