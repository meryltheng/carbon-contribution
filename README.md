# Estimating the impact of invasive pests and diseases on ecosystem services: modelling carbon sequestration loss due to myrtle rust (*Austropuccinia psidii* exotic strains) in Australia

This repository contains the input data, code, and outputs to estimate the impact of myrtle rust on carbon sequestration in Australia. See the preprint here:


**Estimating the impact of invasive pests and diseases on ecosystem services: modelling carbon sequestration loss due to myrtle rust (*Austropuccinia psidii* exotic strains) in Australia**.
Thao P. Le, Meryl Theng, Chris M. Baker, Isobel R. Abell, Tom Kompas, Emma J. Hudgins
bioRxiv 2025.05.30.657121; doi: [10.1101/2025.05.30.657121](https://doi.org/10.1101/2025.05.30.657121)


### Folders


`Data`: raw data files

`Input`: processed data (i.e., model inputs)

`Output`: output damage estimates

`Plots`: output plots (maps)

`Python`: Python scripts to run myrtle rust impact at the species/genus level

`R`: R scripts to run spatial data preparation, spatial biomass calculations, overall impact of myrtle rust to carbon sequestration

### Code authorship

The majority of code is written by Meryl Theng.

The code in the `Python` folder is written by Thao P. Le.

### Instructions

Steps to estimate myrtle rust impact at the species/genus level are found in the `Python` folder.

To run the overall impact estimate (at the national level), run the R files in the `R` folder sequentially (from 0 to 7).

