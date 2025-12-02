# -*- coding: utf-8 -*-
"""
Created on Fri May 30 18:04:33 2025

@author: ryanp

Big merge of spud data from pa,ny,md,wv and oh with bmi station data that had been returned from the streamstats
api call. There are two products of this script. 
1) a dataset with each bmi station and their taxonomic metrics, ffg proportions, sampling season, stream size,
huc8, and within catchment UOGD density, COGD density, % DLC and AMD presence. This dataset is used for lmm
modeling for ffg proportions and taxonomic metrics. 

2) a long format df of taxa recorded at each sample. This contains the station id, taxa id, count (of taxa), 
,pollution tolerance of taxa, taxa ffg, amd_prescence, huc8, sampling season, stream size, ibi score, %DLC, 
U_bool, C_bool and nstations (within that huc8). this dataset is fed into the network creator
"""
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import numpy as np
from shapely import wkt
import openpyxl
#pip install openpyxl --upgrade --pre

bmi_data= pd.read_csv('C:/Users/ryanp/Downloads/RESEARCH/Data/final_stations_data.csv')
taxa = pd.read_excel("C:/Users/ryanp/Downloads/RESEARCH/Data/Macroinvertebrate_Data.xlsx",sheet_name='Taxa')
taxa_metadata = pd.read_csv('C:/Users/ryanp/Downloads/RESEARCH/Data/PADEP Aquatic Macroinvertebrate Taxa.csv')

#first applying either 2001 or 2011 NLCD data depending on data of sampling
bmi_data['Date'] = pd.to_datetime(
    bmi_data[['YEAR', 'MONTH', 'DAY']].rename(columns={
        'YEAR': 'year',
        'MONTH': 'month',
        'DAY': 'day'
    }))

one = pd.to_datetime('2001-01-01')
eleven = pd.to_datetime('2011-01-01')

bmi_data['closer_to_2011'] = (abs(bmi_data['Date'] - eleven) < abs(bmi_data['Date'] - one))

bmi_data['DEV'] = np.where(bmi_data['closer_to_2011'], bmi_data['LC11DEV'], bmi_data['LC01DEV'])

bmi_data = bmi_data[[
        'STATIONID','Date','geometry', 'HUC8', 'HUC10', 'season', 'DRNAREA',
        'US_L3NAME','SmallFrees', 'LargeFrees', 'Richness', 'richEPT', 'SHANdivers','DEV',
        'METHOD_DES','SED_DEP','data_source', 'has_aml', 'area_diff', 'point','FOREST'
       ]]

bmi_data['geometry'] = bmi_data['geometry'].apply(wkt.loads)

bmi_gdf = gpd.GeoDataFrame(
    bmi_data,
    geometry='geometry',
    crs='EPSG:4326'
)


##############combine all well data from state and enverus dataset

#ohio data available via state but merged with enverus data to determine which wells were active
OH_wells = gpd.read_file("C:/Users/ryanp/Downloads/ogwells_statewide (2).zip",layer='OGWells_statewide')
OH_activity = pd.read_csv('C:/Users/ryanp/Downloads/RESEARCH/Data/NETL_data/OH_Washington_Monroe_CountiesWellHeaders.CSV')
OH_activity2 = pd.read_csv('C:/Users/ryanp/Downloads/RESEARCH/Data/NETL_data/OH_All_Other_WellsWellHeaders.CSV')

OH_activity = OH_activity[['API14','Well Status']]
OH_activity2 = OH_activity2[['API14','Well Status']]

#merge OH well data based on API number
OH_active_wells = pd.concat([OH_activity,OH_activity2], ignore_index=True).rename(columns={'API14':'API_NO'})
OH_active_wells['API_NO'] = OH_active_wells['API_NO'].astype(str).str.strip()
OH_wells['API_NO'] = OH_wells['API_NO'].astype(str).str.strip()
OH_wells_full = pd.merge(OH_active_wells, OH_wells, on=['API_NO'], how='right')

#only keep active conventional and unconventional wells
OH_wells_full = OH_wells_full[(OH_wells_full['Well Status'] == 'ACTIVE') & (OH_wells_full['WELL_TYP'].isin(['OG_R', 'UN_R']))][['API_NO','COMP_DATE','WELL_TYP','LAT83','LONG83']]

OH_wells_full['geometry'] = OH_wells_full.apply(lambda row: Point(row['LONG83'], row['LAT83']), axis=1)
OH_spud_gdf = gpd.GeoDataFrame(OH_wells_full, geometry='geometry', crs="EPSG:4326")

#rename column names to be consistent with others(using completion date rather then spud as thats all thats given)
OH_spud_gdf['API'] = OH_spud_gdf['API_NO']
OH_spud_gdf['SPUD_DATE'] = pd.to_datetime(OH_spud_gdf['COMP_DATE'],errors='coerce').dt.strftime('%d-%m-%Y')
OH_spud_gdf['UNCONVENTIONAL'] = OH_spud_gdf['WELL_TYP'].apply(
    lambda x: 'Unconventional' if x == 'UN_R' else ('Conventional' if x == 'OG_R' else None)
)

#FINAL OH DATA
OH_spud = OH_spud_gdf[['API','geometry', 'UNCONVENTIONAL','SPUD_DATE']]



#west virginia well data entirely from enverus dataset
WV_wells = pd.read_csv('C:/Users/ryanp/Downloads/RESEARCH/Data/NETL_data/WV_All_WellsWellHeaders.CSV')

#only keeping active conventional/unconventional wells and associated data
WV_wells = WV_wells[
    (WV_wells['Well Status'] == 'ACTIVE') &
    (WV_wells['Production Type'].isin(['GAS', 'CBM', 'OIL', 'OIL & Gas']))
][[
    'API10', 'Spud Date', 'Drill Type',
    'Surface Hole Latitude (WGS84)', 'Surface Hole Longitude (WGS84)'
]]
#make point data from lat/long, standarize formatting of spud date and standardize column names
WV_wells['geometry'] = WV_wells.apply(lambda row: Point(row['Surface Hole Longitude (WGS84)'], row['Surface Hole Latitude (WGS84)']), axis=1)
WV_spud_gdf = gpd.GeoDataFrame(WV_wells, geometry='geometry', crs="EPSG:4326")
WV_spud_gdf['API'] = WV_spud_gdf['API10']
WV_spud_gdf['SPUD_DATE'] = pd.to_datetime(WV_spud_gdf['Spud Date'], errors='coerce').dt.strftime('%d-%m-%Y')
WV_spud_gdf = WV_spud_gdf.assign(UNCONVENTIONAL=WV_spud_gdf['Drill Type'].apply(
        lambda x: 'Unconventional' if x == 'H' else ('Conventional' if x == 'V' else None)
    )).dropna(subset=['UNCONVENTIONAL'])

#FINAL WV DATA
WV_spud = WV_spud_gdf[['API','geometry', 'UNCONVENTIONAL','SPUD_DATE']]



#NY data colelcted from state database
NY_wells = pd.read_csv("C:/Users/ryanp/Downloads/RESEARCH/Data/Oil__Gas____Other_Regulated_Wells__Beginning_1860_20250420.csv", on_bad_lines='skip')

#keeping only active unconventional and conventional wells
NY_active_wells = NY_wells[NY_wells['Map Symbol'].isin(['GW','OW'])][['API Well Number', 'Date Spudded', 'Surface Latitude', 'Surface Longitude']]

#standardizing column names and spud date to dateime
NY_active_wells['API'] = NY_active_wells['API Well Number']
NY_active_wells['geometry'] = NY_active_wells.apply(lambda row: Point(row['Surface Longitude'], row['Surface Latitude']), axis=1)
NY_spud_gdf = gpd.GeoDataFrame(NY_active_wells, geometry='geometry', crs="EPSG:4326")
NY_spud_gdf['SPUD_DATE'] = pd.to_datetime(NY_spud_gdf['Date Spudded'], errors='coerce').dt.strftime('%d-%m-%Y')

#all wells in NY are conventional as hydraulic fracturing is banned
NY_spud_gdf['UNCONVENTIONAL'] = 'Conventional'

#FINAL NY DATA
NY_spud = NY_spud_gdf[['API','geometry', 'UNCONVENTIONAL','SPUD_DATE']]



#PA data downloaded from state database
spud = pd.read_csv('C:/Users/ryanp/Downloads/RESEARCH/Data/Spud_External_Data.csv')
spud = spud.copy()

#creating point data from lat/long coords and spud date to datetime format
spud['LATITUDE'] = pd.to_numeric(spud['LATITUDE'], errors='coerce')
spud['LONGITUDE'] = pd.to_numeric(spud['LONGITUDE'], errors='coerce')
spud.loc[:, 'geometry'] = spud.apply(lambda row: Point(row['LONGITUDE'], row['LATITUDE']), axis=1)
spud_gdf = gpd.GeoDataFrame(spud, geometry='geometry', crs="EPSG:4326")
spud_gdf['SPUD_DATE'] = pd.to_datetime(spud_gdf['SPUD_DATE'], errors='coerce').dt.strftime('%d-%m-%Y')

#FINAL PA DATA--> only keeping active unconventional and conventional data
PA_spud = spud_gdf[
    (spud_gdf['WELL_STATUS'] == 'Active') &
    (spud_gdf['WELL_CODE_DESC'].isin(['GAS', 'OIL', 'COMB. OIL&GAS',
           'MULTIPLE WELL BORE TYPE', 'COALBED METHANE']))
][['API', 'geometry', "UNCONVENTIONAL",'SPUD_DATE']]



#MD well data downloaded from fractracker database as they previously requested it from the state
MD_wells = pd.read_csv("C:/Users/ryanp/Downloads/MD_OGwells_09292011_obtained12112023.csv")

#converting to gdf point data and spud date to datetime
MD_wells.loc[:, 'geometry'] = MD_wells.apply(lambda row: Point(row['LON_83DD'], row['LAT_83DD']), axis=1)
MD_spud_gdf = gpd.GeoDataFrame(MD_wells, geometry='geometry', crs="EPSG:4326").rename(columns={'API_':'API'})
MD_spud_gdf['SPUD_DATE'] = pd.to_datetime(MD_spud_gdf['YR_DRILLED'], errors='coerce').dt.strftime('%d-%m-%Y')
#all wells are conventional as fracking is illegal in MD
MD_spud_gdf['UNCONVENTIONAL'] = 'Conventional'

#FINAL MD data --> only considering active conventional/unconventional data (all COGD is for gas)
MD_spud = MD_spud_gdf[
    (MD_spud_gdf['ACTIVE'] == 1) &
    (MD_spud_gdf['WELL_STATU'].isin(['GAS', 'Gas Production']))
][['API','geometry','UNCONVENTIONAL','SPUD_DATE']]



#COMBINE ALL WELL DATA
ALL_spud = pd.concat([PA_spud, NY_spud, WV_spud,MD_spud, OH_spud], ignore_index=True)
#ALL_spud.to_csv('marcellus_spud.csv', index=False)

#this function uses the previously aggregated OGD gdf to get a count of how many wells (UOGD and COGD)
# !spudded prior to the sampling date! are within each sample catchment 

def get_well_counts(stations, wells, date_station='Date', date_well='SPUD_DATE', well_type='Unconventional'):
    #ran this function twice (for uogd then again for COGD)
    wells_filtered = wells[wells['UNCONVENTIONAL'] == well_type].copy()

    wells_filtered = wells_filtered.to_crs(stations.crs)
    wells_filtered = wells_filtered[[date_well, 'API', 'geometry']]

    #spatial join to see how many OGD points are in BMI polygons
    wells_in_station = gpd.sjoin(wells_filtered, stations[[date_station, 'geometry']], predicate='within')
    
    #compare well dates to sample dates and only use those where well (spud) date < sample date 
    wells_in_station = wells_in_station[ pd.to_datetime(wells_in_station[date_well], format='%d-%m-%Y', errors='coerce') <= pd.to_datetime(wells_in_station[date_station], format='%d-%m-%Y', errors='coerce') ]

    #get the number of wells for each unique sample 
    well_counts = wells_in_station.groupby('index_right').size()
    count_col = f'{well_type.lower()}'
    stations[count_col] = stations.index.map(well_counts).fillna(0).astype(int)

    return stations

#run for both UOGD and COGD and then divide by the DRNAREA to get well density per catchment
stations_with_wells = get_well_counts(bmi_gdf, ALL_spud, well_type='Unconventional')
stations_with_wells = get_well_counts(bmi_gdf, ALL_spud, well_type='Conventional')
stations_with_wells['unconventional_density'] = stations_with_wells['unconventional'] / stations_with_wells['DRNAREA']
stations_with_wells['conventional_density'] = stations_with_wells['conventional'] / stations_with_wells['DRNAREA']
#stations_with_wells.conventional_density.describe()
#stations_with_wells.unconventional_density.describe() 


#most bmi samples were labeled as either small and large and appropriate metric standardizations were applied
#while small streams generally have a catchment < 52sqmi, no clear boundary was defined. In streams that 
#could not be definitely defined so the both a small and large IBI was calculated. For simple analysis, 
#the mean of these 2 values is used if both are given

def calculate_fibi(row):
    if pd.notna(row['SmallFrees']) and pd.notna(row['LargeFrees']):
        return (row['SmallFrees'] + row['LargeFrees']) / 2
    elif pd.notna(row['SmallFrees']):
        return row['SmallFrees']
    elif pd.notna(row['LargeFrees']):
        return row['LargeFrees']
    else:
        return np.nan
    
stations_with_wells['IBI'] = stations_with_wells.apply(calculate_fibi, axis=1)
#define wadeable (small) samples
stations_with_wells['Small'] = np.where(stations_with_wells.LargeFrees.isna(), 1, 0)

#stations_with_wells.to_csv('allstationdata.csv', index=False)

#swithcing gears to get proportions of FFG per sample to add to modeling data
#need to assign PADEP labels to raw taxa count data
taxa_c = taxa.copy()
ptv_clean = taxa_metadata[['ASSESSMENT_ID','DEP_ASSESS_ID_FFG']].drop_duplicates().dropna()
ptv_clean = ptv_clean.rename(columns={'ASSESSMENT_ID': 'Taxa ID'})
taxa_w_ffg = taxa_c.merge(ptv_clean, on='Taxa ID', how='left').dropna()

#group on station to get count of each ffg at each station
ffg_counts = (
    taxa_w_ffg.groupby(['STATIONID', 'DEP_ASSESS_ID_FFG'])
    .size()
    .reset_index(name='Count')
)

#get proportions by dividing by total taxa per station
ffg_counts['proportion'] = (
    ffg_counts.groupby('STATIONID')['Count']
    .transform(lambda x: x / x.sum())
)
#pivot to wide format to merge back to station data
ffg_proportions = ffg_counts.pivot(
    index='STATIONID',
    columns='DEP_ASSESS_ID_FFG',
    values='proportion'
).fillna(0).reset_index().drop(columns=(['PI','UK']),axis=1)

#merge back to sample station data
modeling_data = stations_with_wells.merge(ffg_proportions, on='STATIONID',how='left')

#drop columns not needed for statistical analysis (LMM)
modeling_data = modeling_data.drop(
    columns=(['Date','geometry','HUC10','SmallFrees','LargeFrees','METHOD_DES','SED_DEP',
    'data_source','area_diff','point','FOREST','unconventional','conventional']),axis=1)

modeling_data.isna().sum()

#WHAT SHOULD THE DATA LOOK LIKE???
modeling_data.shape[0] #6826
modeling_data[(modeling_data.season =='spring')].shape[0] #5966
modeling_data[(modeling_data.Small ==1)&(modeling_data.season=='spring')&(modeling_data.DEV<=20)].shape[0] #4463
modeling_data[(modeling_data.Small ==1)&(modeling_data.season=='spring')&(modeling_data.DEV<=20)&(modeling_data.has_aml==0)].shape[0] #3929

'''
#getting catchment area sizes

modeling_data.DRNAREA.describe()

modeling_data[
    (modeling_data.Small ==1)&
    (modeling_data.season=='spring')&
    (modeling_data.DEV<=20)&
    (modeling_data.has_aml==0)
    ].DRNAREA.plot()
'''

#data used for mixed effect modeling
modeling_data.to_csv('BugModelingData.csv', index=False)

#merging taxa counts with station attributest and converting to long format for creating networks

#simons et al split into quintiles but the distribution of UOGD/COGD doesnt work for that, so instead we do precense/absence
stations_with_wells['U_bool'] = (stations_with_wells['unconventional'] > 0).astype(int)
stations_with_wells['C_bool'] = (stations_with_wells['conventional'] > 0).astype(int)
#keep sampling threshold from simons et al
stations_with_wells['NStations'] = (stations_with_wells.groupby('HUC8')['STATIONID'].transform('nunique'))

station_metadata = stations_with_wells[['STATIONID','HUC8','season','Small','IBI','DEV','U_bool','C_bool','has_aml','NStations']]
taxa = taxa.copy()
stations = station_metadata.copy()

taxa_w_station_features = taxa.merge(stations, on='STATIONID', how='left').dropna(subset=['HUC8'])

#fixing discrepency in taxa labeling in data pulled from stations vs taxa labeling from PADEP data dictionary
taxa_w_station_features['Taxa ID'] = taxa_w_station_features['Taxa ID'].replace(
    'Stenonema(old genus)',  
    'Stenonema')
#get taxa ffg and ptv values
ptv_clean = taxa_metadata[['ASSESSMENT_ID','DEP_ASSESS_ID_FFG','DEP_ASSESS_ID_TV']].drop_duplicates().dropna()
ptv_clean = ptv_clean.rename(columns={'ASSESSMENT_ID': 'Taxa ID'})

taxa_w_features = taxa_w_station_features.copy()
taxa_w_features = taxa_w_features.merge(ptv_clean, on='Taxa ID', how='left').dropna()

#final dataset to feed into network constructor
taxa_w_features.to_csv('NetworkAnalysisData.csv', index=False)
