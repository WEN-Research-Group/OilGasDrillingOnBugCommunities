# OilGasDrillingOnBugCommunities

This repository containis materials affiliated with the peer-review article **"The Legacy of Conventional Oil and Gas Development Outweighs Shale Gas Impacts on Stream Biodiversity"**, submitted to the *ACS ES&T Water*.

## Project Summary
Unconventional oil and gas development (UOGD) has expanded rapidly across the Appalachian Basin, raising concerns about its ecological effects. Although UOGD has been linked to water quality changes that may affect stream biota, a generalized relationship with stream biota remains unresolved. We evaluated how UOGD and conventional oil and gas development (COGD) influence benthic macroinvertebrate (BMI) communities across the Marcellus Shale region of Pennsylvania using one of the most comprehensive statewide BMI datasets, integrating delineated catchments, oil-gas records, and more than 6,800 BMI samples. Linear mixed-effect models and co-occurrence network analyses were used to quantify effects on BMI taxonomy, functionality, and network structure while controlling for confounding factors. We found no significant association between UOGD intensity and BMI diversity metrics, whereas COGD intensity was correlated with reduced richness, diversity, and biotic integrity—comparable in magnitude to the impact of developed land cover. Network analyses indicated altered community structures near both development types: COGD was linked to larger, more connected networks dominated by pollution-tolerant taxa, while UOGD was associated with larger but more fragmented networks. Both forms of development were tied to increases in generalist taxa and declines in specialists. Overall, UOGD exerted limited but detectable ecological effects, whereas COGD imposed broader stress on stream communities. 

## Data
BMI data, publicly available through the Pennsylvania Department of Environmental Protection (PADEP), was downloaded July 1,2024. a processed dataset containing the information necessary for this study (location, well type, and API number) was aggregated into a single dataset. Both datasets are available in the `data/` folder.

BMI samples were processed through the USGS StreamStats API to delineate catchment areas and calculate basin characteristics. These data were then geospatially joined with OGD locations to determine well presence/density within each sample catchment. Data preprocessing and quality control procedures are detailed in the Materials and Methods section of the corresponding paper, as well as Sections S1.3 and S1.5 of the supplementary information.

The final datasets used in the analysis are available in the `download_and_preprocess/processed_datasets/` folder. 

## Analysis
### 1. Linear Mixed-Effect Modeling
```bash
   mixed effect modeling.R
```
Linear mixed-effect models from the [lme4](https://www.jstatsoft.org/article/view/v067i01) package were used to estimate the effects of unconventional and conventional OGD density on taxonomic and functional metrics. Unconventional well density, conventional well density, and the % of the catchment covered by developed land were used as fixed effects. Ecoregion and acid mine drainage presence/absence were included in these models as random effects to account for the significant effect of these variables. The selection of variables considered as random effects is explored in the supplementary information. 

### 2. Co-occurrence Network Generation
```bash
   NetworkCreator_3.R
```
Co-occurrence networks were generated from a subset of BMI samples in R using the [netassoc](https://CRAN.R-project.org/package=netassoc) package, using the framework established by ([Simons 2019](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5751)). Details and justification for sample subsetting are outlined in the manuscript Materials and Methods section. For each network generated, several metrics were calculated with the [igraph](https://zenodo.org/records/18060795) package to describe network structure: network size, connectance, modularity, mean co-occurance strength, and mean pollution tolerance. 

### 3. Co-occurrence Network Analysis
```bash
   network_models.R
```
Linear mixed-effect models were used to estimate the effect of OGD on co-occurence network structure. In these models, unconventional and conventional well presence/absence were included as fixed effects, while watershed (HUC8) was considered as a random effect to account for natural differences in network structure between watersheds. 

### 4. Mapping
```bash
   map plotting.py
```
This script plots a map of Pennsylvania with an outline of the Appalachian Plateau overlain with several layers relevant to this study. These include: oil and gas well locations, BMI sampling locations and catchment areas, locations of acid-mine drainage, and level 3 ecoregion.

### 5. Supplementary Information
```bash
   SI_statistics.py
   SI_modeling_code.R
```
Statistical tests justifying the inclusion/exclusion of variables in modeling. This includes histograms for variable distributions, Kruskal and Wilcoxon tests to test significant differences between groups, regression analysis using spearmans rho for preliminary investigation, and mixed-effect modeling using only samples from wadeable streams to investigate the effect of catchment size on model results. All analysis is detailed in the Supplementary Information

### 6. Figures
Figures from the main text and SI showing the results of our analyses are stored in the `figures/` folder.

## Requirements
- R version 4.4.1 or higher
- R libraries needed: *netassoc, readxl, tidyverse, lubridate, sp, FNN, sf, ggplot2, lme4, lmerTest, broom, performance, car, reformulas*

- Python version 3.9.13 or higher
- Python libraries needed: *matplotlib, rasterio, pandas, geopandas, numpy, pyproj, shapely, requests, time, concurrent, openpyxl, seaborn, scipy, scikit_posthocs*

## Citation
If you use this code or data, please cite both [Paper DOI](https://doi.org/10.5281/zenodo.18491666) and [![Data DOI](https://zenodo.org/badge/1108703630.svg)](https://doi.org/10.5281/zenodo.18491666)

## Contact

For questions, please reach out to Ryan Olivier-Meehan at rpolivie@syr.edu or Dr. Tao Wen at twen08@syr.edu.

## Acknowledgements

This material is based upon work supported by the National Science Foundation under Grant No. OAC-2209864 to Tao Wen. The authors would also like to acknowledge the help and support from Matthew Shank (PADEP).
