
import os, sys
import math
import time
import shutil
from turtle import update

import pandas as pd
import numpy as np

from contextlib import contextmanager

import rioxarray as rxr

import geopandas as gpd

from multiprocessing import Pool

# import rasterio as rio
# from rasterio.mask import mask
from rasterio.warp import Resampling

from shapely.geometry import Point, box, MultiLineString

from whitebox import WhiteboxTools

wbt = WhiteboxTools()

####################
#
# Notes:
#   -snapping using wbt doesn't work well on small basins
#       even with a buffer added to the total bounds of the
#       baseline polygon and original pour point, 
#       the snapped point can end up far away and 
#       delineate a much different basin.
#       -Could try adding custom search, or finding stream name
#
###############

t_in0 = time.time()

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HYSETS_DIR = os.path.join(BASE_DIR, 'source_data/HYSETS_data/')

hysets_stn_df = pd.read_csv(os.path.join(HYSETS_DIR, 'HYSETS_watershed_properties.txt'), sep=';')
hysets_stn_locs = [Point(x, y) for x, y in zip(hysets_stn_df['Centroid_Lon_deg_E'].values, hysets_stn_df['Centroid_Lat_deg_N'])]
hysets_stn_df = gpd.GeoDataFrame(hysets_stn_df, geometry=hysets_stn_locs, crs='EPSG:4326')
(minx, miny, maxx, maxy) = hysets_stn_df.geometry.total_bounds

version = '20220923'

# get all stations with updated WSC polygons
updated_wsc_polygons = []
updated_wsc_polygon_path = '/media/danbot/Samsung_T5/geospatial_data/WSC_data/WSC_basin_polygons_20220721/'
wsc_updated_file = os.path.join(BASE_DIR, 'setup_scripts/Prerelease_included_stations.txt')
wsc_updated_stns_df = pd.read_csv(wsc_updated_file, sep='\t', encoding='latin1')
wsc_updated_stns = wsc_updated_stns_df['STATION_NUMBER'].values

# set the projected CRS to something that will be 
# reasonably accurate for all of NA.
# Otherwise, figure out how to pick the best CRS on the fly...
# HYSETS uses Pavics which uses Canada Atlas Lambert projection
# proj_crs = 3005
proj_crs = 3978

vrt_path = os.path.join(BASE_DIR, 'processed_data/processed_dem/EENV_DEM_mosaic_4326.vrt')
vrt_path = '/media/danbot/Samsung_T5/geospatial_data/DEM_data/EENV_DEM_mosaic.vrt'
# vrt_path_reproj = f'/media/danbot/Samsung_T5/geospatial_data/DEM_data/EENV_DEM_mosaic_{proj_crs}.vrt'

# cmd = f'gdalwarp -s_srs EPSG:4326 -t_srs EPSG:{proj_crs} {vrt_path} {vrt_path_reproj}'
# os.system(cmd)
# print(asdfsd)


DEM_source = 'EENV_DEM'

# external_temp_folder = '/media/danbot/Samsung_T5/hysets_polygons/temp_rasters/'
external_temp_folder = '/home/danbot/Documents/code/hysets_basin_delineation/processed_data/temp_rasters/'

hysets_basins_path = '/media/danbot/Samsung_T5/geospatial_data/HYSETS_data/HYSETS_watershed_boundaries.zip'

output_folder = os.path.join(BASE_DIR, f'processed_data/processed_basin_polygons_{version}')

external_save_path = '/media/danbot/Samsung_T5/hysets_polygons/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

basin_df = gpd.read_file(hysets_basins_path)
basin_df = basin_df.set_crs(4326)

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

def retrieve_raster(fpath):
    rds = rxr.open_rasterio(fpath, masked=True, mask_and_scale=True)
    crs = rds.rio.crs
    affine = rds.rio.transform(recalc=False)
    resolution = [abs(r) for r in rds.rio.resolution()]
    return rds, crs, affine, resolution

_, vrt_crs, _, vrt_res = retrieve_raster(vrt_path)


def clip_raster_to_basin(basin_polygon, raster):
    crs = raster.rio.crs.to_epsg()
    if not crs:
        crs = raster.rio.crs.to_wkt()
    
    basin_polygon = basin_polygon.to_crs(crs)
    bounds = tuple(basin_polygon.bounds.values[0])

    try:
        subset_raster = raster.copy().rio.clip_box(*bounds)
        clipped_raster = subset_raster.rio.clip(basin_polygon.geometry, basin_polygon.crs, all_touched=True)
        return clipped_raster, True
    except Exception as e:
        print('Exception occurred when clipping the raster.')
        print(e)
        print('')
        return None, False


def process_basin_elevation(clipped_raster):
    # evaluate masked raster data
    values = clipped_raster.data.flatten()
    mean_val = np.nanmean(values)
    median_val = np.nanmedian(values)
    min_val = np.nanmin(values)
    max_val = np.nanmax(values)
    return mean_val, median_val, min_val, max_val


def area_comparison(hysets_area, polygon_area, derived_area, print_notes=False):
    """Compare the baseline polygon area with the published values

    Args:
        area (_type_): _description_
        area_proj (_type_): _description_
        expected_da (_type_): _description_
    """
    a1 = 100 * np.abs(polygon_area - hysets_area) / hysets_area
    a2 = 100 * np.abs(derived_area - hysets_area) / hysets_area
    if print_notes:
        if (a1 < 5) & (a2 < 5):
            print(f'E(A-A*) <= 5%: {derived_area:.2f}/{hysets_area:.2f}')
        else:
            print(f'Areas not within 5% of expected value ({hysets_area:.2f}: {polygon_area:.2f}/{derived_area:.2f}):')


def set_buffer_length(expected_area, AB_flag, flag_override):

    if AB_flag & (not flag_override):
        # if the artificial boundary flag is present,
        # return a buffer length equivalent to the expected
        # area divided by 2*dx where dx is the pixel resolution
        print('    ...HYSETS artificial boundary flag.')
        print(f'    ...expected area = {expected_area} km^2.')
        # Using 95% conf. interval for highest GC class
        # GC = 2.1 ==> c = 3.0; n = 0.6 (Sassolas-Serrayet et al 2018)
        # area_m2 = expected_area * 1E6
        square_edge_len = np.sqrt(expected_area * 1E6)
        max_len = (3.0 * expected_area**0.6)
        # print(f'max stream len = {max_len:.1f} km')
        buffer_size = int(max_len * 1E3 - 0.5 * square_edge_len)
        # print(f'    ...square edge len = {square_edge_len:.2f} m')
        print(f'    ...buffer size = {buffer_size:.2f} m')
        return buffer_size

    if expected_area < 10:
        # set the buffer to ~1km for very small basins
        buffer_length_m = 250
    elif expected_area < 100:
        # set the buffer to ~8km for small basins
        buffer_length_m = 2500
    elif expected_area < 500:
        # set the buffer to ~10km for small-medium basins
        buffer_length_m = 5000
    elif expected_area < 1000:
        # set the buffer to ~20km for medium basins
        buffer_length_m = 10000
    else:
        # set the buffer to ~50km for large basins
        buffer_length_m = 20000
    return buffer_length_m


def create_buffer_bounding_box(basin_polygon, stn_loc, expected_da, buffer_length, stn_id):
    """
    Combine the first approximation pour point with the baseline polygon
    and create a buffered bounding box to ensure the pour point falls within the clipped DEM. If not, delineations might be wrong.
    """
    assert basin_polygon.crs == stn_loc.crs == proj_crs

    combined_geom = pd.concat([basin_polygon, stn_loc])
    bbox = box(*tuple(combined_geom['geometry'].total_bounds)).buffer(buffer_length)

    baseline_bbox = gpd.GeoDataFrame(geometry=[bbox], crs=basin_polygon.crs)
    # reproject to same as vrt for clipping raster
    baseline_bbox = baseline_bbox.to_crs(vrt_crs)

    stn_loc_dd = stn_loc.to_crs(vrt_crs)
    # stn_loc_dd.to_file(os.path.join(basin_temp_folder, f'{station_id}_stn_loc.geojson'))

    ppt_check = baseline_bbox.contains(stn_loc_dd).all()
    if not ppt_check:
        print('')
        print('##############################')
        print('')
        err_msg = f'    {stn_id} pour point falls outside dem bbox.'
        print(err_msg)
        print(f'    Expected DA = {expected_da} km2')
        print('')
        print('##############################') 
        return gpd.GeoDataFrame()
    return baseline_bbox


def affine_map_vec_numba(affine, x, y):
    a, b, c, d, e, f, _, _, _ = affine
    n = x.size
    new_x = np.zeros(n, dtype=np.float64)
    new_y = np.zeros(n, dtype=np.float64)
    for i in range(n):
        new_x[i] = x[i] * a + y[i] * b + c
        new_y[i] = x[i] * d + y[i] * e + f
    return new_x, new_y


def get_nearest_point(ppt_loc, acc_rds, affine, area_threshold, expected_area):    
    acc = acc_rds.data[0]
    # max_acc = np.nanmax(acc)
    # print(f'acc max: {max_acc:.0f}')

    ppt_loc = ppt_loc.geometry.values[0]

    # mask for just the cells in the expected range
    max_cells = int((1 + area_threshold) * expected_area)
    min_cells = int((1 - area_threshold) * expected_area)

    # print(f'    Search range: {min_cells} - {max_cells} ({expected_area}) --> {np.nanmin(acc)}, {np.nanmax(acc)}')
    yi, xi = np.where((acc >= min_cells) & (acc <= max_cells))
        
    affine_tup = tuple(affine)
    
    # convert indices to coordinates
    x, y = affine_map_vec_numba(affine_tup, xi, yi)

    # convert to array of tuples
    mask_coords = np.c_[x, y]
    
    # calculate distances from station to flow accumulation points
    stn_coords = (ppt_loc.x, ppt_loc.y)

    all_pts = gpd.GeoDataFrame(geometry=[Point(x, y) for (x, y) in mask_coords], crs=proj_crs)
    # print('    max_acc: ', np.nanmax(acc))
    # diffs = np.array([np.subtract(stn_coords, c) for c in mask_coords])
    
    # dists = [np.linalg.norm(d, ord=2) for d in diffs]
    dists = all_pts.distance(ppt_loc)
    # print(f'   number of values in dists = {len(dists)}')
    # print(dists)
    if len(dists) == 0:
        return False, (None, None)
    
    min_xy = mask_coords[np.argmin(dists)]
    # min_dist = np.linalg.norm(np.subtract(min_xy, stn_coords))
    min_dist = dists.min()

    print(f'    Min. distance from pour point to threshold flow acc cell = {min_dist:.0f} m')
    # print(asdf)
    return True, min_xy


def snap_pour_point(fpath, loc, expected_da, snapped_ppt_path, snapped_shp_path):
    # I think the acc write takes time and the subsequent 
    # read fails as a result?
    acc_rds, acc_crs, acc_affine, acc_resolution = retrieve_raster(fpath)
            
    assert acc_crs == loc.crs

    dx, dy = abs(acc_resolution[0]), abs(acc_resolution[1])
    
    x, y = loc.geometry.x.values[0], loc.geometry.y.values[0]
   
    # Snap pour point to threshold accumulation cell
    # ensure the snap point is at the centre of the cell!
    x_snap, y_snap = x + 0.5 * dx, y - 0.5 * dy
    shifted = False
    distance = 0
    max_shift_distance_m = 200 # meters maximum distance to shift point
    ts = time.time()
    if expected_da < 10:
        area_thresh = 0.5
    else:
        area_thresh = 0.25

    # while (area_thresh > 0.0):
        # print(f'area thresh: {area_thresh}')
        # try:
    expected_area_m2 = expected_da * 1E6
    # print(f'    expected area is {expected_area_m} m')
    expected_cells = int(expected_area_m2 / (dx * dy))
    # print(f'    threshold cells = {expected_cells:.0f}')
    # print(f'expected area [km2]: {expected_da:.2f}')
    # print(f'expected # cells: {expected_cells:.0f}')
  
    nr_found, (x_nr, y_nr) = get_nearest_point(loc, acc_rds, acc_affine, area_thresh, expected_cells)
    
    x_snap, y_snap = x_nr + 0.5 * dx, y_nr - 0.5 * dy
    shifted = True
    distance = np.sqrt((x_nr - x)**2 + (y_nr - y)**2)

    pour_point = Point(x_snap, y_snap)
    # save the shifted location
    snapped_df = gpd.GeoDataFrame(geometry=[pour_point], crs=proj_crs)
    snapped_df.to_file(snapped_ppt_path, driver='GeoJSON')
    snapped_df.to_file(snapped_shp_path)
    print(f'    ...pour point saved, shifted from stn coords: {shifted}')

    return distance


def set_contour_interval(area, basin_polygon, temp_filled_path):
    dem = rxr.open_rasterio(temp_filled_path, masked=True)
    clipped_dem = dem.rio.clip(basin_polygon.geometry, basin_polygon.crs, all_touched=True)

    min_el, max_el = clipped_dem[0].min().item(), clipped_dem[0].max().item()
    el_range = max_el - min_el

    if el_range < 0.1:
        return 0
    elif el_range < 0.2:
        i = 0.05
    elif el_range < 1:
        i = 0.1
    elif el_range < 2:
        i =  0.2
    elif el_range < 10:
        i = 1.0
    else:
        i = float(math.ceil(el_range / 10))
    print(f'    Contour interval = {i:.1f}m for {el_range:.2f}m el. range.')
    return i


def write_contour_geojson(station_id, output_path, contours_path, contour_interval, polygon):
    # open the contour file and convert from shp to geojson
    trimmed_contour_path = os.path.join(output_path, f'{station_id}_contours.geojson')
    if contour_interval > 0:
        contours = gpd.read_file(contours_path).simplify(5)
        basin = polygon.convex_hull.buffer(100).geometry[0]
        clipped_contours = contours.clip(basin)
        
        clipped_contours.to_file(trimmed_contour_path, driver='GeoJSON')
    else:
        edf = gpd.GeoDataFrame()
        edf.to_file(trimmed_contour_path, driver='GeoJSON')


def extract_streams(basin_temp_folder, temp_flow_acc_path, temp_fdir_path, output_path, dem_res, station_id, derived_polygon_gdf, min_area=0.5):
    # generate streams based on a minimum area threshold of 
    # 0.5 km^2 (converted to pixels)
    threshold = int(min_area * 1E6 / (abs(dem_res[0]) * abs(dem_res[1])))
    
    temp_streams_raster_path = os.path.join(basin_temp_folder, f'{station_id}_streams.tif')
    temp_streams_vector_path = os.path.join(basin_temp_folder, f'{station_id}_streams.shp')
    
    with suppress_stdout():  
        wbt.extract_streams(
            temp_flow_acc_path, 
            temp_streams_raster_path, 
            threshold, 
            zero_background=False, 
        )
    with suppress_stdout():
        wbt.raster_streams_to_vector(
            temp_streams_raster_path, 
            temp_fdir_path, 
            temp_streams_vector_path, 
            esri_pntr=False, 
        )
    # convert the streams vector to geojson    
    streams_gdf = gpd.read_file(temp_streams_vector_path)
    streams_gdf = streams_gdf.set_crs(proj_crs)
    assert streams_gdf.crs == derived_polygon_gdf.crs
    
    clipped_streams = gpd.clip(streams_gdf, derived_polygon_gdf)
    # clipped_streams = clipped_streams.dissolve()
    clipped_streams.to_file(output_path, driver='GeoJSON')
    print(f'    ...Stream vectors extracted for {station_id}.')


def derive_basin_polygon(row, proj_crs, basin_temp_folder, write_file=True):
    station_id = row['Official_ID']
    gsim_flag = row['Flag_GSIM_boundaries'] 
    gsim_da = row['Drainage_Area_GSIM_km2']
    hysets_area = row['Drainage_Area_km2']

    stn_loc = gpd.GeoDataFrame(geometry=[row['geometry']], crs=f'EPSG:4326')
    stn_loc = stn_loc.to_crs(proj_crs)    
    loc_path = os.path.join(output_folder, f'{station_id}_og_ppt.geojson')

    shift_dist, derived_area_km2 = None, None

    if gsim_flag:
        hysets_area = gsim_da

    hysets_data = hysets_stn_df[hysets_stn_df['Official_ID'] == station_id]
    hysets_data = hysets_data.to_crs(proj_crs)
    AB_flag = hysets_data['Flag_Artificial_Boundaries'].values[0]   
            
    hysets_polygon = basin_df[basin_df['OfficialID'] == station_id].copy()
    hysets_polygon = hysets_polygon.to_crs(proj_crs)

    assert hysets_polygon.crs == stn_loc.crs

    hysets_polygon_area_km2 = hysets_polygon['Area'].values[0]

    prior_basin = hysets_polygon
    AB_flag_override = False
    if AB_flag:
        if station_id in wsc_updated_stns:
            region_code = station_id[:2]
            fname = f'{station_id}_DrainageBasin_BassinDeDrainage.shp'
            prior_basin_path = os.path.join(updated_wsc_polygon_path, f'{region_code}/{station_id}/{fname}')
            prior_basin = gpd.read_file(prior_basin_path)
            prior_basin = prior_basin.to_crs(hysets_polygon.crs)
            AB_flag_override = True

    buffer_length = set_buffer_length(hysets_area, AB_flag, AB_flag_override)

    print(f'   buffer length = {buffer_length}')

    baseline_bbox = create_buffer_bounding_box(prior_basin, stn_loc, hysets_area, buffer_length, station_id)

    if baseline_bbox.empty:
        print('!!!!!!!  Baseline bbox is empty')
        return None, None

    contours_path = os.path.join(basin_temp_folder, f'{station_id}_contours.geojson')
    
    if write_file:
        og_polygon_path = os.path.join(output_folder, f'{station_id}_og_polygon.geojson')

        if not os.path.exists(og_polygon_path):
            hysets_polygon.to_file(og_polygon_path)
        
        baseline_bbox_path = os.path.join(basin_temp_folder, f'{station_id}_buffered_bbox_temp.geojson')

        assert baseline_bbox.crs.to_epsg() == vrt_crs

        # baseline_bbox = baseline_bbox.to_crs(vrt_crs)
        baseline_bbox.to_file(baseline_bbox_path)
        # buff_box = baseline_bbox.to_crs(proj_crs)
        # buff_box.to_file(os.path.join(output_folder, f'{station_id}_buffered_bbox.geojson'))
        if not os.path.exists(loc_path):
            stn_loc.to_file(loc_path)  

    # write the pour point in a format compatible with WBT (shp)
    # use the reprojected pour point as we will reproject the raster next
    # pour_pt_path = os.path.join(basin_temp_folder, f'{station_id}_pour_pt.shp')        
    temp_dem_path_og_crs = os.path.join(basin_temp_folder, f'{station_id}_DEM_{vrt_crs.to_epsg()}.tif')
    # temp_dem_path_og1_crs = os.path.join(basin_temp_folder, f'{station_id}_DEM_{vrt_crs.to_epsg()}_1.tif')
    temp_dem_path_proj_crs = os.path.join(basin_temp_folder, f'{station_id}_DEM_{proj_crs}.tif')

    named_layer = f'{station_id}_buffered_bbox_temp'
    
    # cut vrt dem to bbox and reproject dem to projected crs
    clip_command = f'gdalwarp -s_srs {vrt_crs} -cutline {baseline_bbox_path} -cl {named_layer} -crop_to_cutline -multi -of gtiff {vrt_path} {temp_dem_path_og_crs} -wo NUM_THREADS=ALL_CPUS'
    
    warp_command = f'gdalwarp -s_srs {vrt_crs} -t_srs EPSG:{proj_crs} -of gtiff {temp_dem_path_og_crs} {temp_dem_path_proj_crs} -wo NUM_THREADS=ALL_CPUS'

    if not os.path.exists(temp_dem_path_proj_crs):  
        tg = time.time()
        print('    gdalwarp started...')
        with suppress_stdout():  
            os.system(clip_command)
        tg1 = time.time()
        print(f'       ...Gdalwarp (clip) completed in {tg1-tg:.2f}s')
        with suppress_stdout():  
            os.system(warp_command)
        tg2 = time.time()
        print(f'       ...Gdalwarp (reproj) completed in {tg2-tg1:.2f}s')

    dem_raster = rxr.open_rasterio(temp_dem_path_proj_crs, masked=True, default='dem')
    dem_res = dem_raster.rio.resolution()
    print('')
    print(f'   Raster loaded, clipped, and reprojected.  Starting basin delineation for {station_id}')

    temp_filled_path = os.path.join(basin_temp_folder, f'{station_id}_filled.tif')

    if not os.path.exists(temp_filled_path):    
        with suppress_stdout():
            wbt.fill_depressions(
                temp_dem_path_proj_crs, 
                temp_filled_path, 
                fix_flats=True, 
                flat_increment=None, 
                max_depth=None, 
            )
        print('   ,...depressions filled')

    temp_fdir_path = os.path.join(basin_temp_folder, f'{station_id}_fdir.tif')

    if not os.path.exists(temp_fdir_path):
        with suppress_stdout():  
            wbt.d8_pointer(
                temp_filled_path, 
                temp_fdir_path, 
                esri_pntr=False, 
            )
        print('   ,...fdir completed')

    temp_flow_acc_path = os.path.join(basin_temp_folder, f'{station_id}_acc.tif')
    if not os.path.exists(temp_flow_acc_path):
        with suppress_stdout():  
            wbt.d8_flow_accumulation(
                temp_filled_path, 
                temp_flow_acc_path, 
                out_type="cells", 
                log=False, 
                clip=False, 
                pntr=False,
                esri_pntr=False, 
            )
        print('   ....acc completed')

    snapped_ppt_path = os.path.join(output_folder, f'{station_id}_ppt_adjusted.geojson')
    snapped_shp_path = os.path.join(basin_temp_folder, f'{station_id}_ppt_adjusted.shp')
    shift_dist = None
    if (~os.path.exists(snapped_ppt_path)) | (~os.path.exists(snapped_shp_path)):
        # max_dist = 0.01 # decimal degrees (0.001 ~ 111m at the equator)
        shift_dist = snap_pour_point(temp_flow_acc_path, stn_loc, hysets_area, snapped_ppt_path, snapped_shp_path)
    
    basin_derived = False
    basin_raster_path = os.path.join(basin_temp_folder, f'{station_id}_watershed.tif')

    converted_polygon_path = os.path.join(output_folder, f'{station_id}_derived.geojson')

    if not os.path.exists(converted_polygon_path):
        # try:
        with suppress_stdout():  
            wbt.watershed(
                temp_fdir_path, 
                snapped_shp_path, 
                basin_raster_path, 
                esri_pntr=False, 
                callback=print('   ,...watershed_derived')
            ) 

        derived_basin_temp_path = os.path.join(basin_temp_folder, f'{station_id}_derived.shp')
        with suppress_stdout():
            wbt.raster_to_vector_polygons(
                basin_raster_path,
                derived_basin_temp_path,
                callback=print('    ...raster basin converted to vector polygon.'),
            )

        # open the newly created polygon .shp file
        derived_polygon_gdf = gpd.read_file(derived_basin_temp_path)
        derived_polygon_gdf = derived_polygon_gdf.set_crs(proj_crs)   
        # save polygon as .geojson file
        derived_polygon_gdf.to_file(converted_polygon_path, driver='GeoJSON')
        basin_derived = True
        print('')
        print(f'    basin polygon and stn loc files written for {station_id}.')  
        derived_area_km2 = derived_polygon_gdf.geometry.area.values[0] / 1E6
        area_comparison(hysets_area, hysets_polygon_area_km2, derived_area_km2, print_notes=True)
        
    else:
        print('Basin polygon already exists.')
        basin_derived = True
        derived_polygon_gdf = gpd.read_file(converted_polygon_path)

    stream_vector_path = os.path.join(output_folder, f'{station_id}_streams.geojson')
    if not os.path.exists(stream_vector_path):
        extract_streams(basin_temp_folder, temp_flow_acc_path, temp_fdir_path, stream_vector_path, dem_res, station_id, derived_polygon_gdf)


    contour_interval = set_contour_interval(hysets_area, hysets_polygon, temp_filled_path)

    contour_output_path = os.path.join(output_folder, f'{station_id}_contours.geojson')
    if not os.path.exists(contour_output_path):
        # wbt.contours_from_raster(
        #     temp_filled_path, 
        #     contours_path,
        #     interval=contour_interval,
        #     # base=0.0, 
        #     # smooth=3, 
        #     # tolerance=5.0, 
        # )
        tc = time.time()
        
        if basin_derived:
            command = f'gdal_contour -i {contour_interval} -e elev -a height {temp_filled_path} {contours_path}'
            if contour_interval > 0:
                try:
                    os.system(command)
                    write_contour_geojson(station_id, output_folder, contours_path, contour_interval, derived_polygon_gdf)
                except Exception as ex:
                    print('')
                    print('')
                    err = f'     ! ! Contour generation failed for {station_id}.'
                    print(err)
                    print('')
                    print(ex)
                    print('')
                    print('')
                    print('####################################')
                    print('')
                    print('')
                    return None, None
        else:
            print('### Basin file not created, contours will not be generated.')

    return shift_dist, derived_area_km2



processed_basins = sorted(list(set([e.split('_')[0] for e in os.listdir(output_folder) if e.endswith('derived.geojson')])))
processed_streams = sorted(list(set([e.split('_')[0] for e in os.listdir(output_folder) if e.endswith('streams.geojson')])))
# processed_stns = []

# print(len(basins_to_process))
# print(asdfsd)
processed_stns = processed_streams

basins_to_process = hysets_stn_df[~hysets_stn_df['Official_ID'].isin(processed_stns)].copy()

print(f'{len(processed_stns)} already processed in {output_folder}, {len(basins_to_process)} remaining.')

basins_to_process.sort_values('Drainage_Area_km2', inplace=True, ascending=True)
t_start = time.time()

bad_basins = {}

results_df = basins_to_process.copy()
# results_df = results_df[results_df['Official_ID'] == '01AD002']

def process_basin(inputs, write_files=True):
    _, row = inputs
    station_id = row['Official_ID']

    basin_temp_folder = os.path.join(external_temp_folder, f'{station_id}/')
    
    if not os.path.exists(basin_temp_folder):
        os.makedirs(basin_temp_folder)   

    shift_dist, derived_area_km2 = derive_basin_polygon(row, proj_crs, basin_temp_folder, write_files)
    basin_output_path = os.path.join(output_folder, f'{station_id}_derived.geojson')
    streams_output_path = os.path.join(output_folder, f'{station_id}_streams.geojson')

    if os.path.exists(basin_output_path) & os.path.exists(streams_output_path):
        print('Removing temp files.')
        shutil.rmtree(basin_temp_folder, ignore_errors=True)
    print('')
    print('    ________________________________________')
    # print(asdfd)
    return station_id, shift_dist, derived_area_km2


# pool = Pool(4)
# results = pool.map(process_basin, basins_to_process.iterrows())
for row in basins_to_process.iterrows():
    station_id, shift_dist, derived_area_km2 = process_basin(row)


if len(bad_basins.keys()) > 0:
    print('Issues occurred with the delineation of the following basins:')
    for k, v in bad_basins.items():
        print(f'{k}: {v}')

# results_df.to_csv(os.path.join(BASE_DIR, f'processed_data/attributes_updated_{version}.csv'))

bad_df = pd.DataFrame(bad_basins, index=range(len(bad_basins.keys())))
print(bad_df)
bad_df.to_csv(os.path.join(BASE_DIR, f'processed_data/bad_basins_{version}.csv'))



    
    

    


        

        

