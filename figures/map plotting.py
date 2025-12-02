# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 14:21:33 2025

@author: ryanp
"""

#pip install matplotlib-scalebar
#pip install --upgrade pandas seaborn
import geopandas as gpd
import pandas as pd
from shapely import wkt
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_scalebar.scalebar import ScaleBar
from shapely.geometry import box

bmi_data = pd.read_csv('C:/Users/ryanp/Downloads/RESEARCH/Scripts/allstationdata.csv')
spud_data = pd.read_csv('C:/Users/ryanp/Downloads/RESEARCH/Data/marcellus_spud.csv')
marcellus = gpd.read_file('C:/Users/ryanp/Downloads/RESEARCH/Data/MiddleDevonianMarcellusShale.gdb').to_crs("EPSG:26917")
us_gdf = gpd.read_file("C:/Users/ryanp/Downloads/s_05mr24.zip").to_crs("EPSG:26917")
physio = gpd.read_file("C:/Users/ryanp/Downloads/RESEARCH/Data/Physiographic_Sections_of_Pennsylvania.zip").to_crs("EPSG:26917")
pa_huc = gpd.read_file("C:/Users/ryanp/Downloads/Watershed_Boundary_Dataset_HUC_8s_1534649511984479449.zip").to_crs("EPSG:26917")
pa_eco = gpd.read_file("C:/Users/ryanp/Downloads/RESEARCH/Data/pa_eco.zip").to_crs("EPSG:26917")

bmi_data['geometry'] = bmi_data['geometry'].apply(wkt.loads)
bmi_data['point'] = bmi_data['point'].apply(wkt.loads)
spud_data['geometry'] = spud_data['geometry'].apply(wkt.loads)
spud_gdf = gpd.GeoDataFrame(spud_data, geometry='geometry',crs="EPSG:4326").to_crs("EPSG:26917")
bmi_gdf = gpd.GeoDataFrame(bmi_data, geometry='geometry',crs="EPSG:4326").to_crs("EPSG:26917")
bmi_point = gpd.GeoDataFrame(bmi_data, geometry='point',crs="EPSG:4326").to_crs("EPSG:26917")
bmi_point = bmi_point[bmi_point.season=='spring']
conv_wells = spud_gdf[spud_gdf['UNCONVENTIONAL'] == 'Conventional']
unconv_wells = spud_gdf[spud_gdf['UNCONVENTIONAL'] == 'Unconventional']

physio_ap = physio[physio.PROVINCE == 'Appalachian Plateaus'].dissolve()
pa_eco = pa_eco.dissolve(by='US_L3NAME').reset_index()
pa_eco_within_ap = pa_eco[pa_eco.geometry.intersects(physio_ap.geometry.iloc[0])]


bbox_gdf= gpd.GeoDataFrame(geometry=[box(-130, 20, -60, 55)], crs="EPSG:4326").to_crs("EPSG:26917")

outside = gpd.overlay(bbox_gdf, us_gdf.dissolve(), how='difference')

bmi_gdf = bmi_gdf[bmi_gdf.season=='spring']
network_subset = bmi_gdf[(bmi_gdf.Small ==1)&(bmi_gdf.season=='spring')&(bmi_gdf.DEV<=20)&(bmi_gdf.has_aml==0)] #3929

network_subset['HUC8'] = network_subset['HUC8'].astype(str).str.zfill(8)
bmi_gdf['HUC8'] = bmi_gdf['HUC8'].astype(str).str.zfill(8)
pa_huc['huc8'] = pa_huc['huc8'].astype(str).str.zfill(8)

huc8_subset = pa_huc[pa_huc['huc8'].isin(bmi_gdf['HUC8'])]

subset_ids = set(network_subset['STATIONID'])
bmi_points_subset = bmi_point[bmi_point['STATIONID'].isin(subset_ids)]
bmi_points_non_subset = bmi_point[~bmi_point['STATIONID'].isin(subset_ids)]


pa_boundary = us_gdf[us_gdf['STATE'] == 'PA']
minx, miny, maxx, maxy = pa_boundary.total_bounds
x_buffer = (maxx - minx) * 0.1
y_buffer = (maxy - miny) * 0.1 


pa_eco_clipped = gpd.clip(pa_eco_within_ap, physio_ap)

ecoregions = [
    'Central Appalachians',
    'Ridge and Valley',
    'North Central Appalachians',
    'Northern Allegheny Plateau',
    'Erie Drift Plain',
    'Western Allegheny Plateau'
]

hex_colors = [
    '#8dd3c7',  # Central Appalachians
    '#ffffb3',  # Ridge and Valley 
    '#bebada',  # North Central Appalachians
    '#fb8072',  # Northern Allegheny Plateau
    '#80b1d3',  # Erie Drift Plain
    '#fdb462'   # Western Allegheny Plateau
]
color_map = dict(zip(ecoregions, hex_colors))
pa_eco_clipped['color'] = pa_eco_clipped['US_L3NAME'].map(color_map)


####################################

fig, ax = plt.subplots(figsize=(20, 20))

marcellus.plot(ax=ax, color='lightgrey', alpha=0.8, label='Marcellus Shale')
us_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=2)
outside.plot(ax=ax, color='grey')

huc8_subset.boundary.plot(ax=ax, edgecolor='dimgrey',linewidth=1)
# Conventional and unconventional wells
conv_wells.plot(ax=ax, color='#67a9cf', markersize=8, label='Conventional')
unconv_wells.plot(ax=ax, color='#ef8a62', markersize=8, label='Unconventional')

# Appalachian Plateau boundary
physio_ap.boundary.plot(ax=ax, edgecolor='black', linewidth=1.5)

scalebar = ScaleBar(
    dx=1,                          # 1 map unit = 1 meter
    units='m',
    dimension='si-length',
    label='Distance (km)',
    length_fraction=0.25,          
    location='upper right',        
    scale_loc='bottom',            
    color='black',                 
    box_alpha=1,                 
    pad=1.5,                       
    font_properties={'size': 14, 'weight': 'bold'}  
)
ax.add_artist(scalebar)

ax.set_facecolor('whitesmoke')
ax.set_xlim(minx - x_buffer, maxx + x_buffer)
ax.set_ylim(miny - y_buffer, maxy + y_buffer)

plt.tight_layout()
#plt.savefig("wellmap.png", format="png", dpi=300)
plt.show()

##################################### #plotting

fig, ax = plt.subplots(figsize=(20, 20))

marcellus.plot(ax=ax, color='lightgrey', alpha=0.8, label='Marcellus Shale')
outside.plot(ax=ax, color='grey', label='Outside area')


bmi_gdf.dissolve().plot(ax=ax, color='#a6bddb', edgecolor='#67a9cf', alpha=0.9, zorder=2, label='All catchments')
network_subset.dissolve().plot(ax=ax, color='blue', alpha=0.6, zorder=3, label='Network subset')

huc8_subset.boundary.plot(ax=ax, edgecolor='dimgray', linewidth=1, label='HUC8 subset')
us_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=2, label='U.S. boundary')
physio_ap.boundary.plot(ax=ax, edgecolor='black', linewidth=1.5)

bmi_points_non_subset.plot(
    ax=ax,
    color='#555555',
    marker='^',
    markersize=40,
    edgecolor='black',
    linewidth=0.5,
    zorder=4,
    label='BMI stations'
)

bmi_points_subset.plot(
    ax=ax,
    color='white',
    marker='^',
    markersize=40,
    edgecolor='black',
    linewidth=1.2,
    zorder=5,
    label='Network subset stations'
)

ax.set_facecolor('whitesmoke')
ax.set_xlim(minx - x_buffer, maxx + x_buffer)
ax.set_ylim(miny - y_buffer, maxy + y_buffer)

plt.tight_layout()
#plt.savefig("catchmentmap.png", format="png", dpi=300)
plt.show()

####################################################### plotting

fig, ax = plt.subplots(figsize=(20, 20))

marcellus.plot(ax=ax, color='lightgrey', alpha=0.8, label='Marcellus Shale')

us_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=2)

outside.plot(ax=ax, color='grey')

huc8_subset.boundary.plot(ax=ax, edgecolor='dimgrey',linewidth=1)
physio_ap.boundary.plot(ax=ax, edgecolor='black', linewidth=1.5)

bmi_point[bmi_point['has_aml'] == 0].plot(
    ax=ax, color='#66c2a5', markersize=50, marker='^',edgecolor='black', linewidth=0.2, label='No AML Present'
)

# Plot stations with AML
bmi_point[bmi_point['has_aml'] == 1].plot(
    ax=ax, color='red', markersize=50, marker='^', edgecolor='black', linewidth=0.2, label='AML Present'
)


ax.set_xlim(minx - x_buffer, maxx + x_buffer)
ax.set_ylim(miny - y_buffer, maxy + y_buffer)

plt.tight_layout()
#plt.savefig("AMLmap.png", format="png", dpi=300)
plt.show()

###################################

fig, ax = plt.subplots(figsize=(20, 20))

marcellus.plot(ax=ax, color='lightgrey', alpha=0.8, label='Marcellus Shale')
outside.plot(ax=ax, color='grey', label='Outside area')

pa_eco_clipped.plot(ax=ax, color=pa_eco_clipped['color'], edgecolor='black', linewidth=0.5)

huc8_subset.boundary.plot(ax=ax, edgecolor='dimgray', linewidth=1, label='HUC8 subset')
us_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=2, label='U.S. boundary')
physio_ap.boundary.plot(ax=ax, edgecolor='black', linewidth=1.5)

bmi_points_non_subset.plot(
    ax=ax,
    color='#555555',
    marker='^',
    markersize=40,
    edgecolor='black',
    linewidth=0.5,
    zorder=4,
    label='BMI stations'
)

bmi_points_subset.plot(
    ax=ax,
    color='white',
    marker='^',
    markersize=40,
    edgecolor='black',
    linewidth=1.2,
    zorder=5,
    label='Network subset stations'
)

scalebar = ScaleBar(
    dx=1,                          # 1 map unit = 1 meter
    units='m',
    dimension='si-length',
    label='Distance (km)',
    length_fraction=0.25,          
    location='upper right',        
    scale_loc='bottom',            
    color='black',                 
    box_alpha=1,                 
    pad=1.5,                       
    font_properties={'size': 14, 'weight': 'bold'}  
)
ax.add_artist(scalebar)


for region, color in color_map.items():
    ax.plot([], [], color=color, label=region, linewidth=10)
ax.legend(title='Ecoregions', loc='lower right', fontsize='large')

ax.set_facecolor('whitesmoke')
ax.set_xlim(minx - x_buffer, maxx + x_buffer)
ax.set_ylim(miny - y_buffer, maxy + y_buffer)

plt.tight_layout()
#plt.savefig("ecomap.png", format="png", dpi=300)
plt.show()

#######################################

minx, miny, maxx, maxy = marcellus.total_bounds
x_buffer = (maxx - minx) * 0.05 
y_buffer = (maxy - miny) * 0.05

fig, ax = plt.subplots(figsize=(20, 20))

outside.plot(ax=ax, color='grey', edgecolor='none', zorder=2)

marcellus.dissolve().plot(ax=ax, color='lightgrey', edgecolor='black', linewidth=0.5, zorder=1)

us_gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=2, zorder=3)

ax.set_xlim(minx - x_buffer, maxx + x_buffer)
ax.set_ylim(miny - y_buffer, maxy + y_buffer)

ax.set_facecolor('whitesmoke')
plt.tight_layout()
#plt.savefig("marcellusmap.png", format="png", dpi=300)

plt.show()

###########################################

variables = ['DRNAREA', 'DEV', 'conventional_density', 'unconventional_density']
titles = ['Drainage Area (DRNAREA)', 'Development (DEV)', 'Conventional Density', 'Unconventional Density']

color_map = {
    'DRNAREA': 'blue',
    'DEV': '#1b9e77',
    'conventional_density': '#67a9cf',
    'unconventional_density': '#ef8a62'
}

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.subplots_adjust(hspace=0.35, wspace=0.25)

# ---- FULL DATA (top row) ----
for i, var in enumerate(variables):
    axes[0, i].hist(
        bmi_gdf[var].dropna(),
        bins=30,
        color=color_map[var],
        edgecolor='black',
        alpha=0.8
    )
    axes[0, i].set_title(titles[i], fontsize=14, weight='bold')
    axes[0, i].set_ylabel('Count', fontsize=12)
    axes[0, i].grid(alpha=0.3)

# ---- SUBSET DATA (bottom row) ----
for i, var in enumerate(variables):
    axes[1, i].hist(
        network_subset[var].dropna(),
        bins=30,
        color=color_map[var],
        edgecolor='black',
        alpha=0.8
    )
    axes[1, i].set_xlabel(var, fontsize=12)
    axes[1, i].set_ylabel('Count', fontsize=12)
    axes[1, i].grid(alpha=0.3)

    if var in ['conventional_density', 'unconventional_density']:
        zero_count = (network_subset[var] == 0).sum()
        nonzero_count = (network_subset[var] > 0).sum()
        axes[1, i].text(
            0.7, 0.9,
            f'0: {zero_count}\n>0: {nonzero_count}',
            transform=axes[1, i].transAxes,
            fontsize=12,
            bbox=dict(facecolor='white', alpha=0.6, edgecolor='black')
        )

plt.tight_layout()
#plt.savefig("sumstats.svg", format="svg")
plt.show()