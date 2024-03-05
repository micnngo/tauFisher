# tauFisher

Circadian time prediction from genomic data. 

## System requirements

### Hardware requirements

tauFisher requires only a standard computer with enough RAM to support in-memory operations. 

### OS requirements

tauFisher has been tested on macOS and Linux:

* macOS: Monterey 12.6
* Linux: CentOS 7

### R package dependencies 

tauFisher requires R >=3.5.0 and has been tested on R 4.0.5. It also requires the following R packages: 

```
MetaCycle
utils
scales
dplyr
stringr
stats
nnet
magrittr
fda
```

## Installation 

Typical installation should take ~5 minutes. 

```
if (!requireNamespace('devtools', quietly = TRUE))
  install.packages('devtools')

devtools::install_github("micnngo/tauFisher", build_vignettes = TRUE) 
```

## Demo

Once `tauFisher` is installed, you can follow any of the three vignettes  for a step-by-step tutorial. 
We suggest the `SingleDataset` or `MultipleDatasets` vignette if you would like to simply run tauFisher with minimal coding. 

Expected run-time on a typical laptop is <10 minutes. 

There are three tutorials which you can access via: 

```
vignette("SingleDataset", package="tauFisher")
vignette("MultipleDatasets", package="tauFisher")
vignette("MultipleDatasets_ByFunctions", package="tauFisher")
```

If you only have one data set to analyze, we suggest following the `SingleDataset` vignette. 
This is best for training on part of the data and predicting the circadian time of the rest of the sample, e.g., training on two replicates and predicting on the third replicate.


If you have two data sets where you know the circadian time for one and would like to predict the circadian time of the other, we suggest following the `MultipleDatasets` vignette. 
And if you prefer executing each step of tauFisher instead of executing two functions (`train_tauFisher` and `test_tauFisher`), please follow the `MultipleDatasets_ByFunctions` vignette.

## Citation 

If you use tauFisher, please cite as:  

Paper:


Package:
[![DOI](https://zenodo.org/badge/423253625.svg)](https://zenodo.org/doi/10.5281/zenodo.10780811)

## License

This package is covered under the GNU General Public License v3.0. 
