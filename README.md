# Protein modelling with EDMs

[![DOI](https://zenodo.org/badge/276220857.svg)](https://zenodo.org/badge/latestdoi/276220857)

### Description
This repository contains code written in `R` to demonstrate a protein modelling approach based on **Euclidean Distance Matrices (EDMs)**.

It contains three `R` scripts: [cp.R](cp.R), [dsd.R](dsd.R) and [da.R](da.R), that model different structural rearrangements (circular permutation, domain swap dimerisation and domain atrophy, respectively) on the structure of an SH3 domain [e1shgA1.pdb](e1shgA1.pdb), downloaded from the [ECOD](http://prodata.swmed.edu/ecod) database.
The [utils.R](utils.R) file contains code shared by all the scripts and imports required packages.

The files [ecod_domains.id](ecod_domains.id) and [protein_bounds.tsv](protein_bounds.tsv) contain the list of 20 **ECOD** domain representatives and the lookup table of distance bounds extracted from them, respectively.
Running times to generate C-alpha and Backbone models of circular permutations (CP) and domain swap dimers (DSD) for each of the 20 **ECOD** domains can be found in the file [runtimes_ecod.tsv](runtimes_ecod.tsv).


### Installation

The only requirement is `R` version **3.6 with** the following packages that can be installed from **CRAN**: `dplyr`, `bio3d`, and `edmcr`.
If you have problems installing the `edmcr` package, try installing it with `devtools` from **GitHub** using:
```
install_github("lafita/edmcr")
```

Scripts can be run in `Rstudio` or from the command line using `Rscript` (within this directory):
```
Rscript cp.R
```

### Article

Article currently in preparation.

