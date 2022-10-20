import os
import time

from numba import jit
# from osgeo import gdal

import pandas as pd
import numpy as np
import geopandas as gpd
import rioxarray as rxr

from rasterio.warp import Resampling
# import richdem as rd

# from multiprocessing import Pool

from shapely.geometry import Point

import warnings
warnings.simplefilter("ignore") 

# this version date should reflect the number
# associated with the derived basin set
version = '20220923'

# Canada Atlas Lambert projection
proj_crs = 3978

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HYSETS_DIR = os.path.join(BASE_DIR, 'source_data/HYSETS_data/')

hysets_stn_df = pd.read_csv(os.path.join(HYSETS_DIR, 'HYSETS_watershed_properties.txt'), sep=';')

hysets_stn_locs = [Point(x, y) for x, y in zip(hysets_stn_df['Centroid_Lon_deg_E'].values, hysets_stn_df['Centroid_Lat_deg_N'])]

hysets_stn_df = gpd.GeoDataFrame(hysets_stn_df, geometry=hysets_stn_locs, crs='EPSG:4326')

ext_media_path = '/media/danbot/Samsung_T5/geospatial_data/'

# porosity and permeability sources
glhymps_fpath = os.path.join(ext_media_path, 'GLHYMPS/GLHYMPS.gdb')

# land use / land cover
nalcms_fpath = os.path.join(ext_media_path, 'NALCMS_NA/NA_NALCMS_2010_v2_land_cover_30m/NA_NALCMS_2010_v2_land_cover_30m.tif')

def retrieve_raster(fpath):
    rds = rxr.open_rasterio(fpath, masked=True, mask_and_scale=True)
    crs = rds.rio.crs
    affine = rds.rio.transform()
    return rds, crs, affine

# import NALCMS raster
nalcms, nalcms_crs, nalcms_affine = retrieve_raster(nalcms_fpath)

# import GLHYMPS raster
glhymps_raster_path = '/media/danbot/Samsung_T5/geospatial_data/GLHYMPS/GLHYMPS.gdb'

derived_basin_path = os.path.join(BASE_DIR, f'processed_data/processed_basin_polygons_{version}')

dem_mosaic_path = '/media/danbot/Samsung_T5/geospatial_data/DEM_data/'
dem_mosaic_file = 'EENV_DEM_mosaic.vrt'

srdem_raster, srdem_crs, srdem_affine = retrieve_raster(os.path.join(dem_mosaic_path, dem_mosaic_file))


# def open_and_clip_raster(basin_polygon, dem_path):
#     # open dem and mask using basin polygon
#     xds = rxr.open_rasterio(dem_path, masked=True)
    
#     crs = xds.rio.crs.to_epsg()
#     if not crs:
#         crs = xds.rio.crs.to_wkt()

#     # xds = xds.rio.reproject('EPSG:3005')
#     basin_polygon = basin_polygon.to_crs(crs)
#     bounds = tuple(basin_polygon.bounds.values[0])
#     subset_raster = xds.rio.clip_box(*bounds)
#     return subset_raster.rio.clip(basin_polygon.geometry, basin_polygon.crs)

def clip_raster_to_basin(basin_polygon, raster):
    crs = raster.rio.crs.to_epsg()
    if not crs:
        crs = raster.rio.crs.to_wkt()
    
    basin_polygon = basin_polygon.to_crs(crs)
    bounds = tuple(basin_polygon.bounds.values[0])

    # print(f'raster crs = {crs}, polygon crs = {basin_polygon.crs}')
    try:
        subset_raster = raster.rio.clip_box(*bounds).copy()
        clipped_raster = subset_raster.rio.clip(basin_polygon.geometry, basin_polygon.crs, all_touched=True)
        return clipped_raster, True
    except Exception as e:
        print(e)
        return None, False


# def clip_and_process_raster_stats(stn, basin_polygon, dem_path):
#     basin_data = open_and_clip_raster(basin_polygon, dem_path)
#     # evaluate masked raster data
#     vals = basin_data.data.flatten()    
#     mean_val, min_val, max_val = np.nanmean(vals), np.nanmin(vals), np.nanmax(vals)
#     return mean_val, min_val, max_val


def get_perm_and_porosity(merged):
    merged['area_frac'] = merged['Shape_Area'] / merged['Shape_Area'].sum()
    weighted_permeability = round((merged['area_frac'] * merged['Permeability_no_permafrost']).sum(), 2)
    weighted_porosity = round((merged['area_frac'] * merged['Porosity']).sum(), 2)
    return weighted_permeability, weighted_porosity
    

def process_glhymps(polygon, fpath, proj_crs):
    # convert  polygon crs to match GLHYMPS
    wkt = 'PROJCS["World_Cylindrical_Equal_Area",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["standard_parallel_1",0.0],UNIT["Meter",1.0]]'
    # convert polygon to same crs as spatial layer (wkt above)
    poly = polygon.to_crs(wkt)
    # mask the spatial layer with the polygon when reading file
    gdf = gpd.read_file(fpath, mask=poly)
    # glhymps_gdf = glhymps_gdf.to_crs('EPSG:3857')
    merged = gpd.overlay(gdf, poly, how='intersection')
    merged = merged.to_crs(proj_crs)
    return merged


def process_basin_elevation(clipped_raster):
    # evaluate masked raster data
    values = clipped_raster.data.flatten()
    mean_val = np.nanmean(values)
    median_val = np.nanmedian(values)
    min_val = np.nanmin(values)
    max_val = np.nanmax(values)
    return mean_val, median_val, min_val, max_val


def calculate_gravelius_and_perimeter(polygon):    
    perimeter = polygon.geometry.length.values[0]
    area = polygon.geometry.area.values[0] 
    if area == 0:
        return np.nan, perimeter
    else:
        perimeter_equivalent_circle = np.sqrt(4 * np.pi * area)
        gravelius = perimeter / perimeter_equivalent_circle
    
    return gravelius, perimeter / 1000


@jit(nopython=True)
def process_slope_and_aspect(E, el_px, resolution, shape):
    # resolution = E.rio.resolution()
    # shape = E.rio.shape
    # note, distances are not meaningful in EPSG 4326
    # note, we can either do a costly reprojection of the dem
    # or just use the approximate resolution of 90x90m
    # dx, dy = 90, 90# resolution
    dx, dy = resolution
    # print(resolution)
    # print(asdfd)
    # dx, dy = 90, 90
    S, A = np.empty_like(E), np.empty_like(E)
    S[:] = np.nan # track slope (in degrees)
    A[:] = np.nan # track aspect (in degrees)
    tot_p, tot_q = 0, 0
    for i, j in el_px:
        if (i == 0) | (j == 0) | (i == shape[0]) | (j == shape[1]):
            continue
            
        E_w = E[i-1:i+2, j-1:j+2]

        if E_w.shape != (3,3):
            continue

        a = E_w[0,0]
        b = E_w[1,0]
        c = E_w[2,0]
        d = E_w[0,1]
        f = E_w[2,1]
        g = E_w[0,2]
        h = E_w[1,2]
        # skip i and j because they're already used
        k = E_w[2,2]  

        all_vals = np.array([a, b, c, d, f, g, h, k])

        val_check = np.isfinite(all_vals)

        if np.all(val_check):
            p = ((c + 2*f + k) - (a + 2*d + g)) / (8 * abs(dx))
            q = ((c + 2*b + a) - (k + 2*h + g)) / (8 * abs(dy))
            cell_slope = np.sqrt(p*p + q*q)
            S[i, j] = (180 / np.pi) * np.arctan(cell_slope)
            A[i, j] = (180.0 / np.pi) * np.arctan2(q, p)

    return S, A

def calculate_circular_mean_aspect(a):
    """
    From RavenPy:
    https://github.com/CSHS-CWRA/RavenPy/blob/1b167749cdf5984545f8f79ef7d31246418a3b54/ravenpy/utilities/analysis.py#L118
    """
    angles = a[~np.isnan(a)]
    n = len(angles)
    sine_mean = np.divide(np.sum(np.sin(np.radians(angles))), n)
    cosine_mean = np.divide(np.sum(np.cos(np.radians(angles))), n)
    vector_mean = np.arctan2(sine_mean, cosine_mean)
    degrees = np.degrees(vector_mean)
    if degrees < 0:
        return degrees + 360
    else:
        return degrees


def calculate_slope_and_aspect(raster):  
    """Calculate mean basin slope and aspect 
    according to Hill (1981).

    Args:
        clipped_raster (array): dem raster

    Returns:
        slope, aspect: scalar mean values
    """
    # print(raster.data[0])
    # print(raster.rio.crs)
    # print(asfd)
    # wkt = raster.rio.crs.to_wkt()
    # affine = raster.rio.transform()

    resolution = raster.rio.resolution()
    raster_shape = raster[0].shape

#     rdem_clipped = rd.rdarray(
#         raster.data[0], 
#         no_data=raster.rio.nodata, 
#         projection=wkt, 
#     )

#     rdem_clipped.geotransform = affine.to_gdal()
#     rdem_clipped.projection = wkt

    # # ts0 = time.time()
    # use to check slope -- works out to within 1 degree...
    # slope = rd.TerrainAttribute(rdem_clipped, attrib='slope_degrees')
    # aspect_deg = rd.TerrainAttribute(rdem_clipped, attrib='aspect')
    # # ts2 = time.time()

    el_px = np.argwhere(np.isfinite(raster.data[0]))

    # print(el_px[:2])
    # print(rdem_clipped)
    # print(asdfsd)
    S, A = process_slope_and_aspect(raster.data[0], el_px, resolution, raster_shape)

    mean_slope_deg = np.nanmean(S)
    # should be within a hundredth of a degree or so.
    # print(f'my slope: {mean_slope_deg:.4f}, rdem: {np.nanmean(slope):.4f}')
    mean_aspect_deg = calculate_circular_mean_aspect(A)
    # if(mean_aspect_deg < 0):
    #     mean_aspect_deg = 90 - mean_aspect_deg
    # elif(mean_aspect_deg > 90.0):
    #     mean_aspect_deg = 360.0 - mean_aspect_deg + 90.0
    # else:
    #     mean_aspect_deg = 90.0 - mean_aspect_deg

    # print(f'my aspect: {mean_aspect_deg:.4f}, rdem: {np.nanmean(aspect_deg):.4f}')

    return mean_slope_deg, mean_aspect_deg


def check_lulc_sum(stn, data):
    checksum = sum(list(data.values())) 
    lulc_check = 1-checksum
    if abs(lulc_check) >= 0.05:
        print(f'   ...{stn} failed checksum: {checksum}')   
    return lulc_check


def recategorize_lulc(data):    
    forest = ('Land_Use_Forest_frac', [1, 2, 3, 4, 5, 6])
    shrub = ('Land_Use_Shrubs_frac', [7, 8, 11])
    grass = ('Land_Use_Grass_frac', [9, 10, 12, 13])
    wetland = ('Land_Use_Wetland_frac', [14])
    crop = ('Land_Use_Crops_frac', [15])
    urban = ('Land_Use_Urban_frac', [16, 17])
    water = ('Land_Use_Water_frac', [18])
    snow_ice = ('Land_Use_Snow_Ice_frac', [19])
    lulc_dict = {}
    for label, p in [forest, shrub, grass, wetland, crop, urban, water, snow_ice]:
        prop_vals = round(sum([data[e] if e in data.keys() else 0 for e in p]), 2)
        lulc_dict[label] = prop_vals
    return lulc_dict
    

def get_value_proportions(data):
    # vals = data.data.flatten()
    all_vals = data.data.flatten()
    vals = all_vals[~np.isnan(all_vals)]
    n_pts = len(vals)
    unique, counts = np.unique(vals, return_counts=True)

    # create a dictionary of land cover values by coverage proportion
    prop_dict = {k: v/n_pts for k, v in zip(unique, counts)}
    prop_dict = recategorize_lulc(prop_dict)
    return prop_dict    


def process_lulc(stn, basin_polygon, nalcms_crs):
    basin_polygon = basin_polygon.to_crs(nalcms_crs)
    bounds = tuple(basin_polygon.bounds.values[0])
    clipped_raster = nalcms.rio.clip_box(*bounds)
    masked_raster = clipped_raster.rio.clip(basin_polygon.geometry, basin_polygon.crs)
    # checksum verifies proportions sum to 1
    prop_dict = get_value_proportions(masked_raster)
    lulc_check = check_lulc_sum(stn, prop_dict)
    prop_dict['lulc_check'] = lulc_check
    return pd.DataFrame(prop_dict, index=[stn])


def process_basin_characteristics(stn, reproj_raster):
    print(f'    ...processing {stn} basin characteristics')
    basin_data = {}
    # print(f'   Processing {stn}')
    basin_data['Official_ID'] = str(stn)

    # hysets_data = hysets_stn_df[hysets_stn_df['Official_ID'] == stn]

    # hysets_area = hysets_data['Drainage_Area_km2'].values[0]
    # hysets_grav = hysets_data['Gravelius'].values[0]
    # hysets_el = hysets_data['Elevation_m'].values[0]
    # hysets_aspect = hysets_data['Aspect_deg'].values[0]
    # hysets_slope = hysets_data['Slope_deg'].values[0]
    # hysets_perim = hysets_data['Perimeter'].values[0]

    mean_el, median_el, min_el, _ = process_basin_elevation(reproj_raster)
    basin_data['median_el'] = median_el
    basin_data['mean_el'] = mean_el
    
    # projected_crs = projected_polygon.crs.to_wkt()
    slope, aspect = calculate_slope_and_aspect(reproj_raster)
    basin_data['Slope_deg'] = slope
    basin_data['Aspect_deg'] = aspect

    reproj_basin = basin_polygon.to_crs(proj_crs)
    gravelius, perimeter = calculate_gravelius_and_perimeter(reproj_basin)
    basin_data['Perimeter'] = perimeter
    basin_data['Gravelius'] = gravelius

    projected_polygon = basin_polygon.to_crs(proj_crs)
    area = projected_polygon.area.values[0] / 1E6
    basin_data['Drainage_Area_km2'] = area

    # hysets_area = hysets_data['Drainage_Area_km2'].values[0]
    # hysets_grav = hysets_data['Gravelius'].values[0]
    # hysets_el = hysets_data['Elevation_m'].values[0]
    # hysets_aspect = hysets_data['Aspect_deg'].values[0]
    # hysets_slope = hysets_data['Slope_deg'].values[0]
    # hysets_perim = hysets_data['Perimeter'].values[0]

    glhymps_df = process_glhymps(basin_polygon, glhymps_fpath, proj_crs)
    weighted_permeability, weighted_porosity = get_perm_and_porosity(glhymps_df)
    basin_data['Permeability_logk_m2'] = weighted_permeability
    basin_data['Porosity_frac'] = weighted_porosity

    # process lulc
    lulc_df = process_lulc(stn, basin_polygon, nalcms_crs)
    lulc_data = lulc_df.to_dict('records')[0]
    basin_data.update(lulc_data)

    # print(f'DA: {area:.1f} ({hysets_area:.1f}) el {mean_el:.0f} ({hysets_el:.0f}) slope: {slope:.1f} ({hysets_slope:.1f}) aspect: {aspect:.1f} ({hysets_aspect:.1f}) grav: {gravelius:.1f} ({hysets_grav:.1f}) perim. {perimeter:.0f} ({hysets_perim:.0f}) ')

    return basin_data


def update_results_files(all_data, update):
    basin_characteristics = pd.DataFrame(all_data)

    if update:
        updated_results = pd.concat([results_df, basin_characteristics])
        print('updated!!!')
        print(updated_results.head())
    else:
        print('new results!!')
        updated_results = basin_characteristics
        print(updated_results.head())

    updated_results.to_csv(results_path, index=False)
    return updated_results


# 'derived' (validation) or 'baseline' (HYSETS)
which_set = 'derived'

derived_basin_polygons = list(set([e.split('_')[0] for e in os.listdir(derived_basin_path) if e.endswith(f'_{which_set}.geojson')]))

stations_to_process = [e.split('_')[0] for e in derived_basin_polygons]

t0 = time.time()

results_path = os.path.join(BASE_DIR, f'results/derived_characteristics_{version}_RB.csv')

update = False
if os.path.exists(results_path):
    results_df = pd.read_csv(results_path)
    processed_stns = results_df['Official_ID'].astype(str).values
    update = True
    stations_to_process = [e for e in stations_to_process if e not in processed_stns]

# p = Pool()

all_data = []
n = 0
t0 = time.time()

print(f'There are {len(stations_to_process)} remaining basins to be processed.')
print('')

def warp_raster(stn, clipped_raster, proj_crs, temp_dem_folder):

    raster_crs = clipped_raster.rio.crs.to_epsg()

    temp_raster_path_in = os.path.join(temp_dem_folder, f'{stn}_temp_{raster_crs}.tif')

    temp_raster_path_out = os.path.join(temp_dem_folder, f'{stn}_temp_{proj_crs}.tif')

    clipped_raster.rio.to_raster(
        temp_raster_path_in,
    )

    warp_command = f'gdalwarp -s_srs EPSG:{raster_crs} -t_srs EPSG:{proj_crs} -of gtiff {temp_raster_path_in} {temp_raster_path_out} -wo NUM_THREADS=ALL_CPUS'

    try:
        os.system(warp_command)
        return True, temp_raster_path_out
    except Exception as ex:
        print('')
        print(f'Raster reprojection failed for {stn}.')
        print('')
        print(ex)
        print('')
        return False, None


for stn in stations_to_process:
    print(f'Starting basin processing on {stn}.')
    polygon_path = os.path.join(derived_basin_path, f'{stn}_{which_set}.geojson')
    basin_polygon = gpd.read_file(polygon_path)

    clipped_raster, raster_loaded = clip_raster_to_basin(basin_polygon, srdem_raster)

    print(f'    ...raster loaded and clipped: {raster_loaded}')

    temp_dem_folder = os.path.join(BASE_DIR, 'temp/')
    if not os.path.exists(temp_dem_folder):
        os.mkdir(temp_dem_folder)

    raster_clipped, temp_raster_path = warp_raster(stn, clipped_raster, proj_crs, temp_dem_folder)

    if not raster_clipped:
        continue
    else:
        reproj_raster, warped_crs, warped_affine = retrieve_raster(temp_raster_path)

    print(f'    ...raster reprojected: {raster_clipped}')

    # reproj_raster = clipped_raster.rio.reproject(proj_crs)
        
    basin_attributes = process_basin_characteristics(stn, reproj_raster)
    all_data.append(basin_attributes)
    n += 1
    # if n % 100 == 0:
    t_end = time.time()
    unit_time = (t_end-t0) / len(all_data)
    updated_df = update_results_files(all_data, update)
    print(f'Processed basin {len(all_data)}/{len(stations_to_process)} in {t_end-t0:.0f}s ({unit_time:.2f}s/basin)')

    temp_files = os.listdir(temp_dem_folder)
    for f in temp_files:
        os.remove(os.path.join(temp_dem_folder, f))
    print('    ...temporary files removed.')
    print('')


