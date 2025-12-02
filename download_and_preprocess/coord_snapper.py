# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 11:43:42 2025

@author: ryanp

snaps the coordinates of PADEP samples to USGS PA streamstats streamgrid 
and gives example window of streamgrid + original sampling point + snapped point
"""

import matplotlib.pyplot as plt
import rasterio as rio
import pandas as pd
import numpy as np
from pyproj import Proj, Transformer

#file paths
raster_path = 'C:/Users/rpolivie/OneDrive - Syracuse University/oil gas extraction effect on benthic macros/Final Code/data/streamgrid.tif'
data_path = 'C:/Users/rpolivie/OneDrive - Syracuse University/oil gas extraction effect on benthic macros/Final Code/data/Macroinvertebrate_Data.xlsx'

#open the raster
with rio.open(raster_path) as src:
    transform = src.transform  
    crs_raster = src.crs  
    raster_data = src.read(1)  
    nodata = src.nodata  
    height, width = raster_data.shape  

#transforming bmi data to raster crs
critters = pd.read_excel(data_path)
coords = critters[['LAT', 'LONG', 'STATIONID']].rename(columns={'LAT': 'lat', 'LONG': 'long'})
proj_wgs84 = Proj(init="epsg:4326")  
transformer = Transformer.from_proj(proj_wgs84, crs_raster)
inverse_transformer = Transformer.from_proj(crs_raster, proj_wgs84)

#function to find the highest-value pixel within the 1 or 3 pixel radius of original bmi coords and gets those new coords
def find_best_pixel(raster, x, y, max_distance=3):
    best_x, best_y = int(round(x)), int(round(y))
    best_value = raster[best_y, best_x] if 0 <= best_x < width and 0 <= best_y < height else -np.inf  

    for dx in range(-max_distance, max_distance + 1):
        for dy in range(-max_distance, max_distance + 1):
            px, py = int(round(x + dx)), int(round(y + dy))
            if 0 <= px < width and 0 <= py < height:  
                pixel_value = raster[py, px]
                if pixel_value != nodata and pixel_value > best_value:  
                    best_value = pixel_value
                    best_x, best_y = px, py

    return best_x, best_y  

#store new coordinates
snapped_coords = []
modified_points = []

#iterate through all the original samples
for index, row in coords.iterrows():
    lat, long, stationid = row['lat'], row['long'], row['STATIONID']
    
    x, y = transformer.transform(long, lat)  
    pixel_x = (x - transform.c) / transform.a
    pixel_y = (y - transform.f) / transform.e

    best_x, best_y = find_best_pixel(raster_data, pixel_x, pixel_y)

    snapped_long, snapped_lat = inverse_transformer.transform(best_x * transform.a + transform.c, best_y * transform.e + transform.f)

    if (best_x != round(pixel_x)) or (best_y != round(pixel_y)):  
        modified_points.append((lat, long, snapped_lat, snapped_long, best_x, best_y, pixel_x, pixel_y))  

    snapped_coords.append((snapped_lat, snapped_long, stationid))

snapped_df = pd.DataFrame(snapped_coords, columns=['lat', 'long', 'STATIONID'])


#------------------------------------------------------------------------------
#the code below can be used to examine the stream raster, original sample loc and snapped sample loc
if modified_points:
    orig_lat, orig_long, snapped_lat, snapped_long, best_x, best_y, orig_pixel_x, orig_pixel_y = modified_points[83] #<-- random sample
    
    #open a 100 × 100 window around the snapped point
    window_size = 100  
    col_off = max(0, best_x - window_size // 2)
    row_off = max(0, best_y - window_size // 2)

    with rio.open(raster_path) as src:
        zoomed_raster = src.read(1, window=rio.windows.Window(col_off, row_off, window_size, window_size))

    orig_pixel_x_win = orig_pixel_x - col_off
    orig_pixel_y_win = orig_pixel_y - row_off
    snapped_pixel_x_win = best_x - col_off
    snapped_pixel_y_win = best_y - row_off

    #plot zoomed raster with original and snapped points
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(zoomed_raster, cmap="terrain", origin="upper")

    ax.plot(orig_pixel_x_win, orig_pixel_y_win, 'go', markersize=8, label="Original Point")  # Green
    ax.plot(snapped_pixel_x_win, snapped_pixel_y_win, 'ro', markersize=8, label="Snapped Point")  # Red

    ax.set_title("Zoomed View: Original vs Snapped Coordinates")
    ax.legend()
    plt.show()
else:
    print("No modified points found.")
    
#------------------------------------------------------------------------------
    
#save snapped coords (redone with 3m snapping radius for later comparison)
snapped_df.to_csv('snapped_BMI_coordinates_3m.csv', index=False)

#len(modified_points) / 15107