# OGD-on-BMI-communities
This reposity is affiliated with the manuscript **"The Legacy of Conventional Oil and Gas Development Outweighs Shale Gas Impacts on Stream Biodiversity"** submitted to ES&T Water in Jan 2026.

## Project Summary
This study analyses the ecological impact of oil and gas development (OGD) within the Appalachian Plateau region of Pennsylvania. Using statewide benthic macroinvertebrate (BMI) samples, we quantify the relationships between OGD presence/density and metrics of community composition. By identifying wells as either conventional or unconventional (employing horizontal drilling and/or hydraulic fracturing) we investigate whether BMI communities respond differently to different extraction methods. In addition to establish metrics of community composition, we employ novel co-occurence network analysis to investigate whether alterations in community stucture are associated with OGD. Our analyses reveals that unconventional OGD exerted limited but detectable ecological effects, whereas conventional OGD imposed broader stress on stream communities. 

## Data
BMI data, publicly available through the Pennsylvania Department of Environmental Protection (PADEP), was downloaded on 7/1/2024. Authors don't have permission to share raw OGD data, information necessary for this study (location, well type, and API number) was aggregated into a single dataset. Both of these datasets are available in the data/ folder. BMI samples were processed through the USGS StreamStats API to delineate catchment areas and calculate basin characteristic and geospatially joined with OGD locations to determine well presence/density within each sample catchment. Data preproccessing and quality control steps are detailed in the manuscript Materials and Methods section as well as S1.3 and S1.5 of the supplementary information. 

The final datasets used in analysis area available in the download_and_preprocess/processed_datasets/ folder. 

## Analysis
### 1. Co-occurence Network Analysis
```bash
   NetworkCreator_3.R
```
Co-occurrence networks were generated in R using the netassoc package using methodology established by [Simons 2019](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5751)
Figures from the main text and SI showing the results of our analysis are stored in the figures/ folder. 



