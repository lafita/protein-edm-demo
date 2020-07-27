# Protein modelling with EDMs

### Description
This repository contains code written in `R` to demonstrate a protein modelling approach based on **Euclidean Distance Matrices (EDMs)**.

It contains three `R` scripts: [cp.R](cp.R), [dsd.R](dsd.R) and [da.R](da.R), that model different structural rearrangements (circular permutation, domain swap dimerisation and domain atrophy, respectively) on the structure of an SH3 domain [e1shgA1.pdb](e1shgA1.pdb), downloaded from the [ECOD](http://prodata.swmed.edu/ecod) database.
The [utils.R](utils.R) file contains code shared by all the scripts and imports required packages.

The files [ecod_domains.id](ecod_domains.id) and [protein_bounds.tsv](protein_bounds.tsv) contain the list of 20 **ECOD** domain representatives and the lookup table of distance bounds extracted from them, respectively.
Running times to generate C-alpha and backbone models of circular permutations (CP) and domain swap dimers (DSD) for each of the 20 **ECOD** domains can be found in the file [runtimes_ecod.tsv](runtimes_ecod.tsv).

### Installation

The only requirement is `R` version **3.6** with the following packages that can be installed from **CRAN**: `dplyr`, `bio3d`, and `edmcr`.

If you have problems installing the `edmcr` package, try installing it from **GitHub** with `devtools` using:
```
library(devtools)
install_github("lafita/edmcr")
```

### Usage

Scripts can be run in `Rstudio` or from the command line using `Rscript` (within this directory):
```
Rscript cp.R
```
This is a proof of concept only, so input parameters are hard-coded in the scripts for simplicity.
To use custom protein structures manually edit the `input` and `output` variables using the path to the new PDB file, and adjust the other options.

If you need help modelling large numbers of structures, please get in touch.

### License

The `R` source code is released under an [MIT license](LICENSE).
The data files are made available under a [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0).

### Article

>Lafita A and Bateman A. Modelling structural rearrangements in proteins using Euclidean distance matrices. F1000Research 2020, 9:728 https://doi.org/10.12688/f1000research.25235.1
