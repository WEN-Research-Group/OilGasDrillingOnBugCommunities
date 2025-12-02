# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 18:38:51 2025

@author: rpolivie
this script uses the USGS streamstraper API to delineate the catchment and get watershed attributes for all BMI samples.
It takes several checks when calling the API to quality control and retries if unable to reach the server
"""

import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from shapely import wkt
import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed


snapped_bug_coords = "C:/Users/rpolivie/OneDrive - Syracuse University/oil gas extraction effect on benthic macros/Final Code/download_and_preprocess/snapped_BMI_coordinates_1m.csv"
streamscrape_results = "1m_snapped_coords_w_data.csv"
safety_save = 250  #save every 250 rows

#function to query StreamStats API with retry on DRNAREA is 0 or cant connect to the server. 5 retries allowed
def scraper(lat, long, retries=5, delay=5):
    url = "https://streamstats.usgs.gov/streamstatsservices/watershed.geojson"
    params = {
        "rcode": "PA",
        "xlocation": long,
        "ylocation": lat,
        "crs": 4326,
        "includeparameters": "true",
        "includeflowtypes": "false",
        "includefeatures": "true",
        "simplify": "true"
    }

    for attempt in range(retries):
        response = requests.get(url, params=params)

        if response.status_code == 200:
            try:
                data = response.json().get("featurecollection", [])
                extracted_properties, geometry = extract_globalwatershed_features(data)

                #make sure DRNAREA exists and is not zero
                if extracted_properties and extracted_properties.get("DRNAREA", 1) != 0:
                    return extracted_properties, geometry
                else:
                    print(f"DRNAREA = 0 for lat={lat}, long={long}, retrying ({attempt+1}/{retries})...")
            except requests.exceptions.JSONDecodeError:
                print(f"JSONDecodeError for lat={lat}, long={long}, retrying ({attempt+1}/{retries})...")
        else:
            print(f"Error {response.status_code} for lat={lat}, long={long}, retrying ({attempt+1}/{retries})...")

        time.sleep(delay)

    print(f"Failed to retrieve valid data for lat={lat}, long={long} after {retries} attempts.")
    return {}, None

#get all available watershed features
def extract_globalwatershed_features(watershed_data):
    for item in watershed_data:
        if item.get("name") == "globalwatershed":
            features = item["feature"]["features"]
            if features:
                properties = features[0]["properties"]
                geometry = Polygon(features[0]["geometry"]['coordinates'][0])
                return properties, geometry
    return {}, None

coords = pd.read_csv(snapped_bug_coords)


#multithreading function
def process_row(index, lat, long):
    extracted_properties, geometry = scraper(lat, long)
    return index, extracted_properties, geometry.wkt if geometry else None

#parallel processing
results = []
with ThreadPoolExecutor(max_workers=4) as executor:
    future_to_index = {executor.submit(process_row, i, row["lat"], row["long"]): i for i, row in coords.iterrows()}
    
    for i, future in enumerate(as_completed(future_to_index)):
        index, extracted_properties, geometry_wkt = future.result()
        results.append((index, extracted_properties, geometry_wkt))

        #save every 250 rows so the whole thing doesnt fail midway (again)
        if (i + 1) % safety_save == 0:
            temp_df = pd.DataFrame(results, columns=["index", "extracted_properties", "geometry"])
            temp_df.to_csv(streamscrape_results, index=False, mode='a', header=not bool(i))

#merge results with original bmi df
results_df = pd.DataFrame(results, columns=["index", "extracted_properties", "geometry"])
coords = coords.merge(results_df, left_index=True, right_on="index").drop(columns=["index"])
properties_df = coords['extracted_properties'].apply(pd.Series)
coords = pd.concat([coords.drop(columns=['extracted_properties']), properties_df], axis=1)

coords["geometry"] = coords["geometry"].apply(lambda x: wkt.loads(x) if pd.notnull(x) else None)
coords_w_data = gpd.GeoDataFrame(coords, geometry='geometry', crs="EPSG:4326")

#save final results (repeated later using 3m snapped data)
coords_w_data.to_csv(streamscrape_results, index=False)

print("DONE!! Data saved to:", streamscrape_results)
