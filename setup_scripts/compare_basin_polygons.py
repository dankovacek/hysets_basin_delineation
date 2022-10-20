import os
import time

import pandas as pd
import numpy as np
import geopandas as gpd

# import json

from PIL import Image #ImageOps, ImageFont, ImageDraw

# import shapely
from shapely.geometry import Point

from functools import reduce

import holoviews as hv
import hvplot.pandas

hv.extension('bokeh', logo=False)

from multiprocessing import Pool

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HYSETS_DIR = os.path.join(BASE_DIR, 'source_data/HYSETS_data/')

source_DEM = 'EENV_DEM' # EENV_DEM or USGS_3DEP

# preprocess_method = 'BURNED' # filled or (stream) burned
preprocess_method = 'FILL' # fill or (stream) burn

# area threshold or Jensen snap
# snap_method = 'ASNAP'

# target_folder = f'{source_DEM}_{preprocess_method}'

revision = '20220923'

img_width = 300
img_height = 200

# disagg_output_path = os.path.join(BASE_DIR, f'processed_data/disagg_polygons_{revision}/')

processed_polygon_path = os.path.join(BASE_DIR, f'processed_data/processed_basin_polygons_{revision}/')

ext_media_folder = '/media/danbot/Samsung_T5/hysets_polygons/'
ext_media_folder = os.path.join(BASE_DIR, f'processed_data/overlay_figs/')

fig_folder = os.path.join(ext_media_folder, f'{revision}/') 

for f in [fig_folder]:
    if not os.path.exists(f):
        os.makedirs(f)

# if not os.path.exists(disagg_output_path):
#     os.mkdir(disagg_output_path)



hysets_geojson_loc = os.path.join(HYSETS_DIR, 'Hysets_stations.geojson')

if not os.path.exists(hysets_geojson_loc):
    hysets_df = pd.read_csv(os.path.join(HYSETS_DIR, 'HYSETS_watershed_properties.txt'), sep=';')
    hysets_locs = [Point(x, y) for x, y in zip(hysets_df['Centroid_Lon_deg_E'].values, hysets_df['Centroid_Lat_deg_N'])]
    hysets_df = gpd.GeoDataFrame(hysets_df, geometry=hysets_locs, crs='EPSG:4326')
    hysets_df.to_file(hysets_geojson_loc)
else:
    hysets_df = gpd.read_file(hysets_geojson_loc)



HY_basins_path = '/media/danbot/Samsung_T5/geospatial_data/HYSETS_data/HYSETS_watershed_boundaries.zip'

HY_df = gpd.read_file(HY_basins_path)
HY_df = HY_df.set_crs(4326)
HY_df = HY_df.to_crs(3978)
HY_df['poly_area'] = (HY_df.geometry.area / 1E6).round(1)
HY_df['pub_area_km2'] = HY_df['OfficialID'].apply(lambda s: hysets_df[hysets_df['Official_ID'] == s]['Drainage_Area_km2'].values[0]).round(1)


def retrieve_data(file):
    fpath = os.path.join(processed_polygon_path, file)
    if os.path.exists(fpath):
        return gpd.read_file(fpath), False
    else:
        return gpd.GeoDataFrame(), True


def get_overlay_plot(stn):

    mapping_crs = 3857
    distance_crs = 3978
    # mapping_crs = 3978
    dem_source = 'EENV'

    hysets_file = f'{stn}_og_polygon.geojson'
    derived_file = f'{stn}_derived.geojson'
    stn_loc_file = f'{stn}_og_ppt.geojson'
    snapped_loc_file = f'{stn}_ppt_adjusted.geojson'
    contours_file = f'{stn}_contours.geojson'

    # baseline_polygon = retrieve_data(hysets_file)
    baseline_polygon, missing_file = retrieve_data(hysets_file)
    # baseline2_polygon = HY_df[HY_df['OfficialID'] == stn]['geometry']
    derived_polygon, missing_file = retrieve_data(derived_file)
    stn_loc, missing_file = retrieve_data(stn_loc_file)
    snapped_loc, missing_file = retrieve_data(snapped_loc_file)
    # contours, missing_file = retrieve_data(contours_file)  

    if missing_file:
        return None

    baseline_proj = baseline_polygon.copy().to_crs(mapping_crs)
    derived_proj = derived_polygon.copy().to_crs(mapping_crs)
    # contours_proj = contours.copy().to_crs(mapping_crs)
    
    published_area = hysets_df[hysets_df['Official_ID'] == stn]['Drainage_Area_km2'].values[0]
    
    baseline_area = baseline_polygon.copy().to_crs(distance_crs).geometry.area.values[0] / 1E6
    # baseline_a2 = baseline2_polygon.copy().to_crs(distance_crs).geometry.area.values[0] / 1E6
    derived_area = derived_polygon.copy().to_crs(distance_crs).geometry.area.values[0] / 1E6
    area_coverage = 100 * (derived_area / baseline_area)

    # assert round(baseline_area, 1) == round(baseline_a2,1)

    # print(f'{baseline_area:.1f}, {baseline_a2:.1f}, {derived_area:.1f}')

    if 100 * abs(1 - baseline_area / published_area) > 5:
        print(f'      AREA WARNING: {stn} area {baseline_area:.1f} varies from published value {published_area:.1f}')
            
    hysets_data = hysets_df[hysets_df['Official_ID'] == stn]
    area_flag = 'NOFLAG'
    if hysets_data['Flag_GSIM_boundaries'].values[0]:
        area_flag = 'GSIMFLAG'
    if hysets_data['Flag_Artificial_Boundaries'].values[0]:
        area_flag = 'ABFLAG'

    # get the common area between the validation and hysets polygons
    tp_polygon = gpd.overlay(derived_proj, baseline_proj, how='intersection').dissolve(aggfunc='sum')
    # get the portion of validation basin outside the hysets basin
    fp_polygon = derived_proj.overlay(baseline_proj, how='difference').dissolve(aggfunc='sum')
    # get the portion of hysets basin not captured by the validation polygon
    fn_polygon = baseline_proj.overlay(derived_proj, how='difference').dissolve(aggfunc='sum')

    d = {
        'tp': {'color': '#94f024', 'poly': tp_polygon},
        'fp': {'color': '#ff5444', 'poly': fp_polygon},
        'fn': {'color': '#8554ff', 'poly': fn_polygon},
        }
    
    plots = []

    # if not contours_proj.empty:   
    #     ct_plot = contours_proj.hvplot(
    #         width=img_width, height=img_height, 
    #         alpha=0.7, color='black', line_width=1)
    #     plots.append(ct_plot)

    # tile_src = 'StamenTerrainRetina'
    # tile_src = 'CartoLight'
    
    tp_exists, fp_exists = False, False
    if tp_polygon.empty:
        tpa = 0
    else:
        l = 'tp'
        ply = d[l]['poly']
        clr = d[l]['color']
        tpa = 100 * ((ply.to_crs(distance_crs).geometry.area.values[0]/1E6) / baseline_area)
        ply = ply.to_crs(mapping_crs)
        tp_plot = ply.hvplot(
            # width=img_width, height=img_height, tiles=tile_src, 
            alpha=0.7, color=clr, line_width=0)
        tp_plot.border_fill_color = None  
        tp_plot.background_fill_color = None
        plots.append(tp_plot)
        tp_exists = True
        
    if fp_polygon.empty:
        fpa = 0
    else:
        l = 'fp'
        ply = d[l]['poly']
        clr = d[l]['color']
        fpa = 100 * ((ply.to_crs(distance_crs).geometry.area.values[0]/1E6) / baseline_area)
        ply = ply.to_crs(mapping_crs)
        fp_exists = True
        if tp_exists:
            fp_plot = ply.hvplot(
                alpha=0.7, color=clr, line_width=0)
        else:
            fp_plot = ply.hvplot(
                width=img_width, height=img_height, #tiles=tile_src, 
                alpha=0.7, color=clr, line_width=0)
            fp_plot.border_fill_color = None  
            fp_plot.background_fill_color = None
        plots.append(fp_plot)
        
    if fn_polygon.empty:
        fna = 100
    else:
        l = 'fn'
        ply = d[l]['poly']
        clr = d[l]['color']
        fna = 100 * ((ply.to_crs(distance_crs).geometry.area.values[0]/1E6) / baseline_area)
        fn_polygon = ply.to_crs(mapping_crs)
        if tp_exists | fp_exists:
            fn_plot = ply.hvplot(
                # tiles=tile_src, width=img_width, height=img_height,
                alpha=0.7, color=clr, line_width=0)
        else:
            fn_plot = fn_polygon.hvplot(
                width=img_width, height=img_height, #tiles=tile_src, 
                alpha=0.7, color=clr, line_width=0)
            fn_plot.border_fill_color = None  
            fn_plot.background_fill_color = None
            
        plots.append(fn_plot)
    
    print(f'            TP: {tpa:.0f}  FP: {fpa:.0f} FN: {fna:.0f}')

    
    # create plot
    stn_loc = stn_loc.to_crs(mapping_crs)
    
    plot_title = f'{stn} {baseline_area:.0f} km2 {area_flag} TP{tpa:.0f} FP{fpa:.0f} FN{fna:.0f} {dem_source}'
    
    stn_loc_plot = stn_loc.hvplot(title=plot_title, 
                          # tiles=tile_src, width=900, height=600,
                          xaxis=None, yaxis=None, tools=[],
                          marker='star', size=300, line_width=2,
                          fill_color='deepskyblue', line_color='dodgerblue')

    stn_loc_plot.background_fill_color = None
    stn_loc_plot.border_fill_color = None  

    stn_loc_plot.opts(fontsize={
        'title': 7, 
        # 'labels': 14, 
        # 'xticks': 10, 
        # 'yticks': 10,
    })

    if snapped_loc is not None:
        snapped_loc = snapped_loc.to_crs(mapping_crs)
        og_loc_plot = snapped_loc.hvplot(
            xaxis=None, yaxis=None, tools=[],
            marker='circle', size=300, line_width=2,
            fill_color='firebrick', line_color='dodgerblue'
            )
    
        layout =  reduce(lambda x, y: x*y, plots) * stn_loc_plot * og_loc_plot
    else:
        layout =  reduce(lambda x, y: x*y, plots) * stn_loc_plot 
            
    filename = f'{stn}_{baseline_area:.2f}KM2_{area_flag}_TP{tpa:.0f}_FP{fpa:.0f}_FN{fna:.0f}.png'
    save_path = os.path.join(fig_folder, filename)
    hvplot.save(layout, save_path)


existing_files = os.listdir(fig_folder)

completed_stations = [e.split('_')[0] for e in existing_files]

t_start = time.time()
num_processed = 0
i = 0

# set up inputs for multithreading
files_to_process = list(set([e.split('_')[0] for e in os.listdir(processed_polygon_path)]))

files_to_process = np.setdiff1d(files_to_process, completed_stations)
# files_to_process = []

print(f'{len(files_to_process)} files left to process.')
# print(asfsad)
for ff in files_to_process[::-1]:
    get_overlay_plot(ff)

# with Pool() as p:
#     p.map(get_overlay_plot, files_to_process)

t_end = time.time()
print(f' Processed {len(files_to_process)} basin polygon images in {t_end-t_start:.1f}s')

# Create a mosaic image of all the overlay figures.
existing_images = os.listdir(fig_folder)

d = {}
for f in existing_images:    
    if 'Legend' not in f:
        data = f.split('_')
        stn = data[0]
        area = data[1].split('K')[0]
        flag = data[2]
        tp = data[3][2:]
        fp = data[4][2:]
        fn = data[-1].split('.')[0][2:]
        
        d[stn] = {
            'area': area,
            'flag': flag,
            'tp': int(tp),
            'fp': int(fp),
            'fn': int(fn),
            'path': os.path.join(fig_folder, f)
        }


def equal_bins(x, nbin=20):
    nlen = len(x)
    return np.interp(np.linspace(0, nlen, nbin + 1),
                     np.arange(nlen),
                     np.sort(x))


df = pd.DataFrame(d).T
df.index.name = 'stn'

df['col'] = np.nan
df['row'] = np.nan

df['tp'] = df['tp'].astype(float)
df['area'] = df['area'].astype(float)

df.reset_index(inplace=True)

def get_binned_vals(param, n_bins=20):
    bins = equal_bins(df[param], n_bins)
    max_len = 0
    for i in range(1, len(bins)):
        left = bins[i-1]
        right = bins[i]
        in_bin = (df[param] >= left) & (df[param] < right)
        # assign the column number based on the published drainage area
        df.loc[in_bin, 'col'] = i
        if in_bin.sum() > max_len:
            max_count = in_bin.sum()
    for i in range(1, len(bins)):
        left = bins[i-1]
        right = bins[i]
        # assign row values based on the true positive rate, then fp, then fn
        in_bin = (df[param] >= left) & (df[param] < right)

        conf_cols = ['tp', 'fp', 'fn']

        fdf = df.loc[in_bin, conf_cols].copy()
        
        fdf = fdf.sort_values(conf_cols, ascending=[False, True, True])

        rank_indices = fdf.index.values
        ranks = list(range(0, len(fdf)))

        fdf['rank_idx'] = rank_indices

        df.loc[rank_indices, 'row'] = max_count - ranks + 1

        # ranked_ind
    return df, bins

# if len(existing_images) < 

img_width = 300
img_height = 200

n_cols = int(0.90 * np.sqrt(len(existing_images)))
# n_cols = 10

df, eq_bins = get_binned_vals('area', n_bins=n_cols)

print(f'Setting {n_cols} bins (columns) based on {len(existing_images)} images to make roughly square.')

width = df['col'].max()
height = df['row'].max()

# print(df['row'].min())
# Create a new image where the width corresponds to the number of cols
# and the height is set to make even rows (~equal number of images per col)

collage = Image.new('RGB', (int((width)*img_width), int((height)*img_height)), (255,255,255))

# print(width, height)

for y1 in range(0, int(height)+1):
    for x1 in range(0, int(width)+1):
        paste_img_path = df[(df['col'] == x1) & (df['row']== y1)]['path'].values
  
        n_imgs = len(paste_img_path)
        if n_imgs ==  0:
            pass
        elif n_imgs == 1:    
            paste_img = Image.open(paste_img_path[0])
            collage.paste(paste_img, (int((x1-1)*img_width), int((y1-1)*img_height)))
        else:
            print(' too many images returned')
            print(paste_img_path.split('/')[-1]('_')[0])

legend_path = os.path.join(BASE_DIR, 'processed_data/overlay_figs/00-Legend.png')
legend_img = Image.open(legend_path)

collage.paste(legend_img, (int(width)*img_width, 0))


mosaic_output_fname = f'{source_DEM}_collage_{revision}.png'
mosaic_path = os.path.join(BASE_DIR, f'processed_data/overlay_figs/{mosaic_output_fname}') 
collage.save(mosaic_path)

print(f'    Created image mosaic: {mosaic_output_fname}')
print('')
print('################################')
print('')
# print('   Collecting disaggregated polygons')
# out_path = os.path.join(BASE_DIR, f'procesed_data/disaggregated_polygons_summed_{revision}.geojson')

# disagg_polygon_files = os.listdir(disagg_output_path)
# paths = [os.path.join(disagg_output_path, dp) for dp in disagg_polygon_files]

# def filter_polygons(file):
#     gdf = gpd.read_file(file)
#     stn = file.split('/')[-1].split('_')[0]
#     stn_data = hysets_df[hysets_df['Official_ID'] == stn].copy()
#     expected_area = stn_data['Drainage_Area_km2'].values[0]
#     gdf['expected_area_km2'] = expected_area
#     gdf['AB_flag'] = stn_data['Flag_Artificial_Boundaries'].values[0]
#     gdf['GSIM_flag'] = stn_data['Flag_GSIM_boundaries'].values[0]
#     gdf = gdf.explode(index_parts=False)
#     gdf['area_km2'] = gdf.geometry.area / 1E6
#     gdf['pct_area'] = gdf['area_km2'] / expected_area
#     gdf = gdf[gdf['pct_area'] > 0.05]
#     gdf.reset_index(inplace=True, drop=True)
#     return gdf

# with Pool() as p:
#     dps = p.map(filter_polygons, paths)

# output_gdf = gpd.GeoDataFrame(pd.concat(dps), crs='EPSG:3005')
# output_gdf = output_gdf.dissolve(by='polygon_type', aggfunc='sum')
# output_gdf.reset_index(inplace=True)
# output_gdf['category'] = output_gdf['polygon_type']
# output_gdf.to_file(out_path, driver='GeoJSON')