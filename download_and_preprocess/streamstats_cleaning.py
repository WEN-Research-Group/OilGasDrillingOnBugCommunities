# -*- coding: utf-8 -*-
"""
Created on Wed May 28 13:11:33 2025

@author: rpolivie
"""
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
from shapely import wkt


#get the raw 1m/3m results and compare to original BMI sample catchment arearesults
raw_1m = pd.read_csv('C:/Users/rpolivie/OneDrive - Syracuse University/oil gas extraction effect on benthic macros/Final Code/download_and_preprocess/processed datasets/1m_snapped_coords_w_data.csv')
raw_3m = pd.read_csv('C:/Users/rpolivie/OneDrive - Syracuse University/oil gas extraction effect on benthic macros/Final Code/download_and_preprocess/processed datasets/3m_snapped_coords_w_data.csv')
critters = pd.read_excel('C:/Users/rpolivie/OneDrive - Syracuse University/oil gas extraction effect on benthic macros/Final Code/data/Macroinvertebrate_Data.xlsx')
pa_physio = gpd.read_file('C:/Users/rpolivie/OneDrive - Syracuse University/oil gas extraction effect on benthic macros/Final Code/data/Physiographic_Provinces_of_Pennsylvania.zip')
amd = gpd.read_file("C:/Users/rpolivie/OneDrive - Syracuse University/oil gas extraction effect on benthic macros/Final Code/data/Abandoned_Mine_Land_Inventory_-_AML_Problem_Locations.zip")
pa_eco = gpd.read_file('C:/Users/rpolivie/OneDrive - Syracuse University/oil gas extraction effect on benthic macros/Final Code/data/pa_eco.zip')



stations = critters[['STATIONID','STREAM_NAME']]
raw_3m = pd.merge(stations, raw_3m, left_index=True, right_index=True)
freestone = critters[critters['METHOD']=='Freestone']
freestone['point'] = freestone.apply(lambda row: Point(row['LONG'], row['LAT']), axis=1)
free_gdf = gpd.GeoDataFrame(freestone, geometry='point').set_crs(epsg=4326, inplace=True)
pa_eco = pa_eco.dissolve(by= 'US_L3NAME').to_crs(free_gdf.crs)
eco_stations = gpd.sjoin(free_gdf, pa_eco,predicate='intersects')
ap_stations = gpd.sjoin(eco_stations,pa_physio, how='left')
ap_stations = ap_stations[['STATIONID','STREAM_NAME','point','HUC8','HUC10','HUC12','YEAR','MONTH','DAY',"US_L3NAME",
                           'METHOD_DES','SmallFrees', 'LargeFrees','Richness','richEPT','SHANdivers','SED_DEP',
                           'CHAN_FLOW', 'COND_BANKS', 'BANK_VEG', 'GRAZING', 'RIP_VEG_ZO','AREAsqmi','PROVINCE']]

#check the differences between both 1m and 3m streamstats output and true drainage area and keep whichever closer
diff_1m = np.abs((raw_1m['DRNAREA'] - critters['AREAsqmi']) / raw_1m['DRNAREA'])
diff_3m = np.abs((raw_3m['DRNAREA'] - critters['AREAsqmi']) / raw_3m['DRNAREA'])

use_1m = diff_1m < diff_3m

result = pd.DataFrame()

result = raw_1m.copy()
result = result.where(use_1m, raw_3m)
result['data_source'] = np.where(use_1m, '1m', '3m')

full_stations = pd.merge(ap_stations,result, on='STATIONID').dropna(subset='AREAsqmi')

full_stations['area_diff'] = (full_stations['DRNAREA'] - full_stations['AREAsqmi']) / full_stations['DRNAREA'] 

trimmed_stations = full_stations[full_stations['area_diff'].abs() <= 0.2] 
 
trimmed_stations = trimmed_stations.dropna(subset='geometry') 
  
trimmed_stations['geometry'] = trimmed_stations['geometry'].astype(str).apply(wkt.loads)

trimmed_stations_gdf = gpd.GeoDataFrame(trimmed_stations, geometry='geometry').set_crs(epsg=4326, inplace=True)

final_stations = trimmed_stations_gdf[trimmed_stations_gdf['PROVINCE'].isin(['Appalachian Plateaus Province'])]
#small_frees = trimmed_stations_gdf[trimmed_stations_gdf['SmallFrees'].notna()]

#label for spring or fall sampling protocols
final_stations['season'] = final_stations['MONTH'].apply(
    lambda m: 'spring' if m in [10,11,12,1,2,3,4,5] else 'fall'
)

#label for amd presence
amd = amd.to_crs(final_stations.crs)
final_stations = final_stations.copy()

final_stations['geometry'] = final_stations['geometry'].apply(lambda g: g if g.is_valid else g.buffer(0))
joined = gpd.sjoin(amd, final_stations, how='inner', predicate='within')
stations_with_aml = joined['STATIONID'].unique()
final_stations['has_aml'] = final_stations['STATIONID'].isin(stations_with_aml).astype(int)


#BMI samples now have catchments within %20 area of original samples and ws attributes (like amd and dev LC)
#total 6826 samples (of all szns and stream sizes)
final_stations.to_csv("final_stations_data.csv")


##################################################### plotting
'''
fig, ax = plt.subplots(figsize=(10, 10))
final_stations[final_stations['has_aml'] == 0].plot(
    ax=ax, color='limegreen', markersize=20, edgecolor='black', linewidth=0.2, label='No AML'
)

amd.plot(ax=ax,markersize=1)
plt.show()



pip install matplotlib-scalebar

import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.patches import FancyArrow

#marcellus = gpd.read_file('C:/Users/rpolivie/Downloads/marcellus_extent.zip')
final_stations_gdf = gpd.GeoDataFrame(final_stations, geometry='point', crs="EPSG:4326")
states = gpd.read_file('C:/Users/rpolivie/Downloads/cb_2018_us_state_500k.zip')
huc12_gdf = gpd.read_file('C:/Users/rpolivie/OneDrive - Syracuse University/Tao DATA/PA_HUC12_clip.zip')

huc12_gdf['HUC8'] = huc12_gdf['HUC_12'].astype(str).str[:8]

# Dissolve polygons by HUC8
huc8_gdf = huc12_gdf.dissolve(by='HUC8', as_index=False)

# Optional: keep only the geometry and HUC8 column
huc8_gdf = huc8_gdf[['HUC8', 'geometry']]

# --- Define Pittsburgh coordinates only ---
pittsburgh = gpd.GeoDataFrame({
    'city': ['Pittsburgh'],
    'geometry': [Point(-79.9959, 40.4406)]
}, crs="EPSG:4326")

# --- Reproject everything to UTM (meters) ---
proj_crs = "EPSG:26917"
states_proj = states.to_crs(proj_crs)
#marcellus_proj = marcellus.to_crs(proj_crs)
stations_proj = final_stations_gdf.to_crs(proj_crs)
pittsburgh_proj = pittsburgh.to_crs(proj_crs)
pa_physio_proj = pa_physio.to_crs(proj_crs)
huc8_gdf_proj = huc8_gdf.to_crs(proj_crs)
# --- Get PA bounding box in projected coordinates ---
pa_boundary = states[states['NAME'] == 'Pennsylvania'].to_crs(proj_crs)
minx, miny, maxx, maxy = pa_boundary.total_bounds

# --- Plot ---
fig, ax = plt.subplots(figsize=(10, 10))

# Plot all state boundaries

pa_physio_proj = pa_physio_proj[pa_physio_proj.PROVINCE == 'Appalachian Plateaus Province']
# Plot Appalachian Plateau boundary
#marcellus_proj.plot(ax=ax, color='lightgrey',alpha =0.7)
huc8_gdf_proj.boundary.plot(ax=ax, color = 'blue', linewidth=0.5, alpha=0.4)
states_proj.boundary.plot(ax=ax, color='black', linewidth=1.5)
pa_physio_proj.boundary.plot(ax=ax, color='black', linewidth=0.75)



# Plot stations without AML
stations_proj[stations_proj['has_aml'] == 0].plot(
    ax=ax, color='limegreen', markersize=20, edgecolor='black', linewidth=0.2, label='No AML'
)

# Plot stations with AML
stations_proj[stations_proj['has_aml'] == 1].plot(
    ax=ax, color='red', markersize=20, edgecolor='black', linewidth=0.2, label='Has AML'
)


# Plot Pittsburgh
#pittsburgh_proj.plot(ax=ax, color='blue', markersize=150, marker='*', edgecolor='black', linewidth=1, zorder=5)
#x, y = pittsburgh_proj.geometry.x.values[0], pittsburgh_proj.geometry.y.values[0]
#ax.text(x + 15000, y + 15000, 'Pittsburgh', fontsize=12, fontweight='bold', color='blue')

# Add scale bar
scalebar = ScaleBar(1, units="m", dimension="si-length", location='lower right', scale_loc="bottom")
ax.add_artist(scalebar)

# Add north arrow
arrow_x = minx + (maxx - minx)*0.9
arrow_y = miny + (maxy - miny)*0.1
ax.annotate('N', xy=(arrow_x, arrow_y + 30000), xytext=(arrow_x, arrow_y),
            arrowprops=dict(facecolor='black', width=3, headwidth=15), ha='center', fontsize=12, fontweight='bold')

# Set limits to slightly larger than PA
x_buffer = (maxx - minx) * 0.05
y_buffer = (maxy - miny) * 0.05
ax.set_xlim(minx - x_buffer, maxx + x_buffer)
ax.set_ylim(miny - y_buffer, maxy + y_buffer)

# Formatting
plt.title('Stations in Pennsylvania Colored by AML Presence')
plt.xlabel('Easting (m)')
plt.ylabel('Northing (m)')
plt.legend()
plt.tight_layout()

#plt.savefig("AMD3.svg", format="svg")
plt.show()

###########
pa_boundary = states[states['NAME'] == 'Pennsylvania'].to_crs(proj_crs)
minx, miny, maxx, maxy = pa_boundary.total_bounds

# --- Plot ---
fig, ax = plt.subplots(figsize=(10, 10))

pa_physio_proj.boundary.plot(ax=ax, color='black', linewidth=0.75)


# Overlay Marcellus extent (light grey)
marcellus_proj.plot(ax=ax, color="lightgrey", alpha=0.25)

# Plot HUC8 boundaries
huc8_gdf_proj.boundary.plot(ax=ax, color="blue", linewidth=0.3, alpha=0.5)

# Plot Pennsylvania boundary
states_proj.boundary.plot(ax=ax, color="black", linewidth=1)

# Plot stations colored by their ecoregion
stations_proj.plot(column="US_L3NAME", ax=ax, markersize=20, legend=True, edgecolor="black", linewidth=0.2)

# --- Scale bar ---
scalebar = ScaleBar(1, units="m", dimension="si-length", location="lower right", scale_loc="bottom")
ax.add_artist(scalebar)

# --- North arrow ---
arrow_x = minx + (maxx - minx) * 0.9
arrow_y = miny + (maxy - miny) * 0.1
ax.annotate("N", xy=(arrow_x, arrow_y + 30000), xytext=(arrow_x, arrow_y),
            arrowprops=dict(facecolor="black", width=3, headwidth=12),
            ha="center", fontsize=12, fontweight="bold")

# --- Map extent ---
x_buffer = (maxx - minx) * 0.05
y_buffer = (maxy - miny) * 0.05
ax.set_xlim(minx - x_buffer, maxx + x_buffer)
ax.set_ylim(miny - y_buffer, maxy + y_buffer)

# --- Formatting ---
plt.title("Macroinvertebrate Stations by Ecoregion (US Level III)", fontsize=14, weight="bold")
plt.xlabel("Easting (m)")
plt.ylabel("Northing (m)")
plt.tight_layout()

#show and save
plt.savefig("stations_by_ecoregion.svg", format="svg", dpi=300)
plt.show()
'''