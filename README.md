# Estimating the impact of invasive pests and diseases on ecosystem services: modelling carbon sequestration loss due to myrtle rust (*Austropuccinia psidii* exotic strains) in Australia

This repository contains the input data, code, and outputs to estimate the impact of myrtle rust on carbon sequestration in Australia. See the preprint here:


**Estimating the impact of invasive pests and diseases on ecosystem services: modelling carbon sequestration loss due to myrtle rust (*Austropuccinia psidii* exotic strains) in Australia**.
Thao P. Le, Meryl Theng, Chris M. Baker, Isobel R. Abell, Tom Kompas, Emma J. Hudgins
bioRxiv 2025.05.30.657121; doi: [10.1101/2025.05.30.657121](https://doi.org/10.1101/2025.05.30.657121)


### Repository contents


`Data/`: Raw data files, including: vegetation plot surveys, myrtle rust impacts, [spatial land use data (NVIS)](https://www.dcceew.gov.au/environment/environment-information-australia/national-vegetation-information-system/data-products#key51))

`Input/`: Processed data used as model inputs, including: biomass estimates, carbon sequestration values, downloaded [climate projections](https://www.longpaddock.qld.gov.au/qld-future-climate/data-info/tern-cmip6/)

`Output/`: Model outputs, including: estimated damage values, climate suitability (infection risk) for myrtle rust

`Plots/`: Visualisations generated from the analysis

`Python/`: Python scripts to run myrtle rust impact at the species/genus level

`R/`:R scripts for data preparation and national-level impact analysis. Tasks include: climate data processing, spatial biomass calculations, carbon sequestration impact estimation. The scripts are sequentially numbered from `0` to `7`, with each script building on the previous one. They are designed to be run in order. Running the full pipeline will generate processed inputs in the `Input/` folder and final outputs in the `Output/` folder.
**Note:** Omit the `4_Damages.R` script. This script was used in the initial version of the analysis (for the [report](https://cebra.unimelb.edu.au/__data/assets/pdf_file/0006/5054550/23C_FINAL_Report.pdf)), but has been replaced by the `6_Damages plus.R` script in the current version.


### Code authorship

The majority of code is written by Meryl Theng.

The code in the `Python` folder is written by Thao P. Le.

### Instructions

Steps to estimate myrtle rust impact at the species/genus level are found in the `Python` folder.

To run the overall impact estimate (at the national level), run the R files in the `R` folder sequentially (from 0 to 7).

