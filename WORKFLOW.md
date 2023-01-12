# WORKFLOW

This file is devoted to listing, in chronological order, the processing steps that are necessary to make the whole model to work in a proper way. In particular, the main steps are:

- download input data
- process and clean data in order to have inputs that can be used directly in the codes and are consistent between each other with respect to spatial and temporal resolution, metadata, data type/formatting
- run the model
- optimize the optimizer
- make comparisons with other products

# Input data

## Field data
From IRRmodel\TEST_SITE\TEST_SITE_\*.nc get daily
- irrigation
- rainfall
- ET
- PET
- other satellite data (THEIA, CCI active/passive, SMOS, SMAP)

From Data\Platinum_Budrio.xlsx get hourly:
- SWC
- irrigation
- rainfall
- temperature

## Radar data
Main codes:
in dir GEE:
- GEE_to_sigma0: [requires ee module and dependencies of environment google] performs download of sigma0 values from GEE, average over area of interest (in linear units, then passed to dB), normalization of orbits wrt reference orbit with incidence angle that is the closest to 40° 
- Sigma_norm (???)
- VV_VH_extraction (???)


# Model

# Optimizer

# Other products