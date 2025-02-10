import os
import subprocess
import multiprocessing
from osgeo import gdal, ogr, osr

gdal.UseExceptions()

washington_tiles = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 
                    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 
                    53, 54, 55, 56, 57, 58, 59, 84, 85, 86, 87, 88, 89, 90, 115, 116, 117, 118, 119, 120, 121, 122, 
                    146, 147, 148, 149, 150, 151, 152, 153, 180, 181, 182, 183, 184, 185, 186, 187, 214, 215, 216, 
                    217, 218, 219, 220, 221, 249, 250, 251, 252, 253, 254, 255, 256, 286, 287, 288, 289, 290}

def remove_temp_files(*file_patterns):
    for file_pattern in file_patterns:
        if os.path.exists(file_pattern):
            os.remove(file_pattern)
            print(f"[INFO] Removed temporary file: {file_pattern}")

def generate_aoi_mask(aoi_raster, output_shapefile):
    print(f"[INFO] Generating AOI mask from {aoi_raster}...")
    command = f"gdal_polygonize.py '{aoi_raster}' -f 'ESRI Shapefile' '{output_shapefile}'"
    subprocess.run(command, shell=True, check=True)
    print(f"[INFO] AOI mask saved as {output_shapefile}")

def clip_and_resample_and_reproject(input_raster, output_raster, aoi_mask, aoi_raster, reprojected_raster, final_resampled_raster):
    """Clip and resample input raster to match AOI resolution and bounding box, then reproject to CONUS Albers (EPSG:5070)"""
    os.makedirs(os.path.dirname(output_raster), exist_ok=True)
    os.makedirs(os.path.dirname(reprojected_raster), exist_ok=True)
    os.makedirs(os.path.dirname(final_resampled_raster), exist_ok=True)
    print(f"[INFO] Clipping and resampling {input_raster} to match AOI {aoi_raster}...")
    options = gdal.WarpOptions(
        cutlineDSName=aoi_mask,
        cropToCutline=True,
        dstNodata=-9999
    )
    gdal.Warp(output_raster, input_raster, options=options)
    print(f"[INFO] Clipped and resampled raster saved: {output_raster}")
    print(f"[INFO] Reprojecting {output_raster} to CONUS Albers (EPSG:5070)...")
    options_reproj = gdal.WarpOptions(
        dstSRS="EPSG:5070",
        resampleAlg=gdal.GRA_Bilinear,
        dstNodata=-9999
    )
    gdal.Warp(reprojected_raster, output_raster, options=options_reproj)
    print(f"[INFO] Reprojected raster saved: {reprojected_raster}")
    print(f"[INFO] Forcing resolution to 120m x 120m for {reprojected_raster}...")
    options_force_res = gdal.WarpOptions(
        xRes=120,
        yRes=120,
        resampleAlg=gdal.GRA_Bilinear,
        dstNodata=-9999
    )
    gdal.Warp(final_resampled_raster, reprojected_raster, options=options_force_res)
    print(f"[INFO] Final raster with 120m resolution saved: {final_resampled_raster}")

def process_tile(folder, base_dir, aoi_tiles_dir, output_dir):
    """Processes a single tile - Clip, resample, and reproject"""
    folder_path = os.path.join(base_dir, folder, "dems_folder", "dem0", "momentum", "0-0-deg")
    if not os.path.exists(folder_path):
        return
    for file in os.listdir(folder_path):
        if not file.endswith("_vel.asc"):
            continue
        base_name = file.replace("_vel.asc", "")
        speed_file = os.path.join(folder_path, f"{base_name}_vel.asc")
        prj_file = os.path.join(folder_path, f"{base_name}_vel.prj")
        if not (os.path.exists(speed_file) and os.path.exists(prj_file)):
            print(f"[WARNING] Missing files in {folder_path}. Skipping.")
            return
        with open(prj_file, 'r') as prj:
            projection_wkt = prj.read()
        # Find matching AOI tile
        matching_aoi_tile = next(
            (os.path.join(aoi_tiles_dir, f) for f in os.listdir(aoi_tiles_dir) if f"utm_overlap_aoi_{folder}_" in f),
            None
        )
        if not matching_aoi_tile or not os.path.exists(matching_aoi_tile):
            print(f"[WARNING] No matching AOI tile found for tile index {folder}")
            return
        final_tile = os.path.join(output_dir, folder, f"{base_name}_clipped.tif")
        final_reproj_tile = os.path.join(output_dir, folder, f"{base_name}_reproj.tif")
        final_resampled_raster = os.path.join(output_dir, folder, f"{base_name}_120m.tif")
        aoi_mask_shp = f"/tmp/aoi_mask_{folder}.shp"
        remove_temp_files(aoi_mask_shp)
        generate_aoi_mask(matching_aoi_tile, aoi_mask_shp)
        clip_and_resample_and_reproject(speed_file, final_tile, aoi_mask_shp, matching_aoi_tile, final_reproj_tile, final_resampled_raster)
        remove_temp_files(aoi_mask_shp)

def process_tiles(base_dir, aoi_tiles_dir, output_dir):
    """Processes tiles in parallel using multiprocessing"""
    os.makedirs(output_dir, exist_ok=True)
    tile_folders = [f for f in os.listdir(base_dir) if f.isdigit() and int(f) in washington_tiles]
    num_workers = min(16, os.cpu_count())
    with multiprocessing.Pool(processes=num_workers) as pool:
        pool.starmap(process_tile, [(folder, base_dir, aoi_tiles_dir, output_dir) for folder in tile_folders])
    print("[INFO] Processing complete.")

#mosaicking
def clip_tile(tile_path):
    """Clip the outer 15% of a tile to remove hard seams."""
    ds = gdal.Open(tile_path)
    if not ds:
        print(f"[ERROR] Could not open tile: {tile_path}")
        return None
    width, height = ds.RasterXSize, ds.RasterYSize
    gt = ds.GetGeoTransform()
    new_width = int(width * 0.7)  # Keeping 70% of original width
    new_height = int(height * 0.7)  # Keeping 70% of original height
    min_x = gt[0] + (width * 0.15 * gt[1])
    max_y = gt[3] + (height * 0.15 * gt[5])
    max_x = gt[0] + (width * 0.85 * gt[1])
    min_y = gt[3] + (height * 0.85 * gt[5])
    clipped_tile_path = tile_path.replace("_reproj.tif", "_clipped.tif")
    gdal.Warp(clipped_tile_path, tile_path, outputBounds=(min_x, min_y, max_x, max_y),
              width=new_width, height=new_height, dstNodata=-9999)
    print(f"[INFO] Clipped tile saved: {clipped_tile_path}")
    return clipped_tile_path

def mosaic_tiles(output_dir, final_mosaic_path):
    """Clip tiles in parallel, smooth blending of remaining 10% overlap, and mosaic seamlessly."""
    reprojected_tiles = []
    for folder in os.listdir(output_dir):
        folder_path = os.path.join(output_dir, folder)
        if os.path.isdir(folder_path):
            for file in os.listdir(folder_path):
                if file.endswith("_reproj.tif"):
                    reprojected_tiles.append(os.path.join(folder_path, file))
    if not reprojected_tiles:
        print("[ERROR] No reprojected tiles found for mosaicking.")
        return
    print(f"[INFO] Clipping outer 15% in parallel from {len(reprojected_tiles)} tiles...")
    num_workers = min(16, os.cpu_count())
    with multiprocessing.Pool(processes=num_workers) as pool:
        clipped_tiles = pool.map(clip_tile, reprojected_tiles)
    clipped_tiles = [t for t in clipped_tiles if t is not None]
    if not clipped_tiles:
        print("[ERROR] No valid clipped tiles found for mosaicking.")
        return
    print(f"[INFO] Mosaicking {len(clipped_tiles)} clipped tiles with cubic interpolation...")
    vrt_path = os.path.join(output_dir, "temp_mosaic.vrt")
    gdal.BuildVRT(vrt_path, clipped_tiles, options=gdal.BuildVRTOptions(resampleAlg="cubic", resolution="user", xRes=120, yRes=120))
    gdal.Translate(final_mosaic_path, vrt_path, format="GTiff", options=gdal.TranslateOptions(outputType=gdal.GDT_Float32, xRes=120, yRes=120))
    print(f"[INFO] Final seamless mosaic saved: {final_mosaic_path}")
    os.remove(vrt_path)
    for clipped_tile in clipped_tiles:
        os.remove(clipped_tile)

if __name__ == "__main__":
    base_dir = "/mnt/d/CONUS"
    aoi_tiles_dir = "/mnt/d/tiles/utm_aoi_tiles"
    output_dir = "/mnt/d/tiles/windNinjax_tiles_momentum"
    # process_tiles(base_dir, aoi_tiles_dir, output_dir)
    mosaic_tiles(output_dir, os.path.join(output_dir, "final_mosaic_w_smoothing_momentum.tif"))
