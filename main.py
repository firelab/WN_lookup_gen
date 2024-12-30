import os
from WindNinja.dem_processing import mosaic_dem_tiles, process_dem, resample_tif
from WindNinja.spatial_operations import read_boundary_file, reproject_boundary
import rioxarray as rxr
from shapely.geometry import box

def run_wind_ninja_narr_cell_parallel(narr_grid_df, dem_path):
    raster_data = rxr.open_rasterio(dem_path, masked=True)
    
    dem_bounds = raster_data.rio.bounds()
    print("DEM bounds:", dem_bounds)
    dem_bbox = box(*dem_bounds)

    # Filter polygons intersecting the DEM
    narr_grid_df = narr_grid_df[narr_grid_df.intersects(dem_bbox)]
    print(f"Polygons intersecting DEM bounds: {len(narr_grid_df)}")

    # Process polygons
    for idx, row in narr_grid_df.iterrows():
        simnum = str(idx).zfill(7)
        subpoly = row.geometry
        output_dir = os.path.join("/mnt/c/Users/dgh00/OneDrive/Desktop/WN_lookup_gen/out", simnum)
        process_dem(subpoly, 16000, raster_data, output_dir)

if __name__ == "__main__":
    
    boundary_file = 'WindNinjaData/ne_10m_admin_1_states_provinces_lakes.shp'
    dem_file_path = 'WindNinjaData/LC20_Elev_220_RS_30m.tif'
    
    resampled_dem_file_path = resample_tif(dem_file_path, 120)

    bnd = read_boundary_file(boundary_file)
    dem_crs = "EPSG:5070"
    bnd_reprojected = reproject_boundary(bnd, dem_crs)
    
    run_wind_ninja_narr_cell_parallel(bnd_reprojected, resampled_dem_file_path)
    
    mosaic_dem_tiles("/mnt/c/Users/dgh00/OneDrive/Desktop/WN_lookup_gen/out", "EPSG:32612", (120,120))