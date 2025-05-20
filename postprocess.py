"""
 * Project:  WN_Lookup_Gen
 * Purpose:  Python script for converting tiles(Windninja asc and prj files) from UTM to conus albers
 * Author: Gunjan Dayani <gunjan.dayani@usda.gov>

 ******************************************************************************
 *
 * THIS SOFTWARE WAS DEVELOPED AT THE ROCKY MOUNTAIN RESEARCH STATION (RMRS)
 * MISSOULA FIRE SCIENCES LABORATORY BY EMPLOYEES OF THE FEDERAL GOVERNMENT 
 * IN THE COURSE OF THEIR OFFICIAL DUTIES. PURSUANT TO TITLE 17 SECTION 105 
 * OF THE UNITED STATES CODE, THIS SOFTWARE IS NOT SUBJECT TO COPYRIGHT 
 * PROTECTION AND IS IN THE PUBLIC DOMAIN. RMRS MISSOULA FIRE SCIENCES 
 * LABORATORY ASSUMES NO RESPONSIBILITY WHATSOEVER FOR ITS USE BY OTHER 
 * PARTIES,  AND MAKES NO GUARANTEES, EXPRESSED OR IMPLIED, ABOUT ITS QUALITY, 
 * RELIABILITY, OR ANY OTHER CHARACTERISTIC.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************

"""

import os
import subprocess
import multiprocessing
from osgeo import gdal, ogr, osr
import numpy as np
import argparse

gdal.UseExceptions()

tiles = set(range(2500))

def get_raster_bounds(raster_path):
    """Get bounding box (minx, miny, maxx, maxy) from a raster."""
    ds = gdal.Open(raster_path)
    gt = ds.GetGeoTransform()
    x_min = gt[0]
    x_max = x_min + ds.RasterXSize * gt[1]
    y_max = gt[3]
    y_min = y_max + ds.RasterYSize * gt[5]
    ds = None
    return (x_min, y_min, x_max, y_max)

def rasters_overlap(raster1_path, raster2_path):
    """Check if two rasters overlap spatially."""
    ds1 = gdal.Open(raster1_path)
    ds2 = gdal.Open(raster2_path)
    if ds1 is None or ds2 is None:
        return False

    gt1 = ds1.GetGeoTransform()
    gt2 = ds2.GetGeoTransform()

    x1_min = gt1[0]
    x1_max = x1_min + ds1.RasterXSize * gt1[1]
    y1_max = gt1[3]
    y1_min = y1_max + ds1.RasterYSize * gt1[5]

    x2_min = gt2[0]
    x2_max = x2_min + ds2.RasterXSize * gt2[1]
    y2_max = gt2[3]
    y2_min = y2_max + ds2.RasterYSize * gt2[5]

    ds1 = None
    ds2 = None

    return not (x1_max < x2_min or x1_min > x2_max or y1_min > y2_max or y1_max < y2_min)

def remove_temp_files(*file_patterns):
    for file_pattern in file_patterns:
        if os.path.exists(file_pattern):
            os.remove(file_pattern)
            print(f"[INFO] Removed temporary file: {file_pattern}")

def wind_sd_to_uv(speed_array, direction_array):
    direction_array = np.where(direction_array == 360.0, 0.0, direction_array)
    dir_rad = np.radians(direction_array)
    u_array = -speed_array * np.sin(dir_rad)
    v_array = -speed_array * np.cos(dir_rad)
    u_array = np.where((direction_array == 0.0) | (direction_array == 180.0), 0.0, u_array)
    v_array = np.where((direction_array == 90.0) | (direction_array == 270.0), 0.0, v_array)
    return u_array.astype(np.float32), v_array.astype(np.float32)

def make_tif(vel_asc, ang_asc, prj_file, output_tif):
    """
    Converts velocity and direction ASC files into a two-band GeoTIFF.
    - Band 1: U wind component
    - Band 2: V wind component
    """
    if not (os.path.exists(vel_asc) and os.path.exists(ang_asc) and os.path.exists(prj_file)):
        print(f"[WARNING] Missing required files: {vel_asc}, {ang_asc}, {prj_file}. Skipping.")
        return None

    # Read ASC metadata
    def read_asc_metadata(asc_path):
        with open(asc_path, "r") as f:
            lines = f.readlines()
        metadata = {}
        for line in lines[:6]:  # First 6 lines contain metadata
            key, value = line.strip().split()
            metadata[key.lower()] = float(value) if "." in value else int(value)
        return metadata

    metadata = read_asc_metadata(vel_asc)

    with open(prj_file, 'r') as prj:
        projection_wkt = prj.read()

    # Load ASC files and convert speed/direction to U/V components
    speed_array = np.loadtxt(vel_asc, skiprows=6)
    direction_array = np.loadtxt(ang_asc, skiprows=6)
    u_array, v_array = wind_sd_to_uv(speed_array, direction_array)

    os.makedirs(os.path.dirname(output_tif), exist_ok=True)

    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(
        output_tif,
        metadata["ncols"], metadata["nrows"],
        2,  # Two bands: U and V
        gdal.GDT_Float32,
        options=["COMPRESS=LZW"]
    )

    cellsize = metadata["cellsize"]
    geotransform = (
        metadata["xllcorner"], cellsize, 0,
        metadata["yllcorner"] + metadata["nrows"] * cellsize, 0, -cellsize
    )

    out_ds.SetGeoTransform(geotransform)
    out_ds.SetProjection(projection_wkt)

    band1 = out_ds.GetRasterBand(1)
    band1.WriteArray(u_array)
    band1.SetNoDataValue(-9999)

    band2 = out_ds.GetRasterBand(2)
    band2.WriteArray(v_array)
    band2.SetNoDataValue(-9999)

    band1, band2, out_ds = None, None, None  # Close dataset

    print(f"[INFO] Created TIF: {output_tif} (U in Band 1, V in Band 2)")

    return output_tif

def crop_to_original_tile(final_raster):
    ds = gdal.Open(final_raster)
    if ds is None:
        print(f"[WARNING] Could not open raster: {final_raster}. Skipping cropping.")
        return
    gt = ds.GetGeoTransform()
    pixel_size_x, pixel_size_y = gt[1], abs(gt[5])
    width, height = ds.RasterXSize, ds.RasterYSize
    center_x = gt[0] + (width / 2) * pixel_size_x
    center_y = gt[3] - (height / 2) * pixel_size_y
    half_extent = 45000
    min_x = center_x - half_extent
    max_x = center_x + half_extent
    min_y = center_y - half_extent
    max_y = center_y + half_extent
    ds = None
    print(f"[INFO] New bounding box: ({min_x}, {min_y}, {max_x}, {max_y})")

    gdal.Warp(final_raster, final_raster, outputBounds=(min_x, min_y, max_x, max_y), dstNodata=-9999)
    print(f"[INFO] Cropped {final_raster} to 90km x 90km centered region.")

def clip_and_resample_and_reproject(input_raster, aoi_raster, clipped_raster, reprojected_raster, final_resampled_raster):
    """
    1. Clip input_raster using aoi_raster bounds
    2. Reproject clipped raster to EPSG:5070
    3. Crop to 90km × 90km centered box
    4. Resample to 120m resolution
    """
    os.makedirs(os.path.dirname(clipped_raster), exist_ok=True)
    os.makedirs(os.path.dirname(reprojected_raster), exist_ok=True)
    os.makedirs(os.path.dirname(final_resampled_raster), exist_ok=True)

    print(f"[INFO] Clipping and resampling {input_raster} to match AOI {aoi_raster}...")

    # Step 1: Check overlap
    if not rasters_overlap(input_raster, aoi_raster):
        print(f"[WARNING] Input raster {input_raster} does not overlap AOI {aoi_raster}. Skipping.")
        return

    # Step 2: Get AOI bounds
    minx, miny, maxx, maxy = get_raster_bounds(aoi_raster)

    # Step 3: Clip to AOI bounds
    options_clip = gdal.WarpOptions(
        outputBounds=(minx, miny, maxx, maxy),
        dstNodata=-9999
    )
    gdal.Warp(clipped_raster, input_raster, options=options_clip)
    print(f"[INFO] Clipped raster saved: {clipped_raster}")

    # Step 4: Reproject to EPSG:5070
    options_reproj = gdal.WarpOptions(
        dstSRS="EPSG:5070",
        resampleAlg=gdal.GRA_Bilinear,
        dstNodata=-9999
    )
    gdal.Warp(reprojected_raster, clipped_raster, options=options_reproj)
    print(f"[INFO] Reprojected raster saved: {reprojected_raster}")

    # Step 5: Crop to 90km x 90km centered box
    print(f"[INFO] Cropping {reprojected_raster} to 90km × 90km box...")
    crop_to_original_tile(reprojected_raster)

    # Step 6: Final resample to 120m
    options_final = gdal.WarpOptions(
        xRes=120,
        yRes=120,
        resampleAlg=gdal.GRA_Bilinear,
        dstNodata=-9999
    )
    gdal.Warp(final_resampled_raster, reprojected_raster, options=options_final)
    print(f"[INFO] Final raster with 120m resolution saved: {final_resampled_raster}")

def process_tile(folder, base_dir, aoi_tiles_dir, output_dir, directions):
    """Processes a single tile - Clip, resample, and reproject."""
    for wind_dir in directions:
        folder_path = os.path.join(base_dir, folder, "dems_folder", "dem0", "momentum", wind_dir)
        if not os.path.exists(folder_path):
            print(f"[WARNING] Folder {folder}: wind direction '{wind_dir}' does not exist. Skipping.")
            continue
        
        # Look for *vel.asc files in this subfolder
        for file in os.listdir(folder_path):
            if not file.endswith("_vel.asc"):
                continue
            
            base_name = file.replace("_vel.asc", "")
            speed_file = os.path.join(folder_path, f"{base_name}_vel.asc")
            direction_file = os.path.join(folder_path, f"{base_name}_ang.asc")
            prj_file = os.path.join(folder_path, f"{base_name}_vel.prj")
            tile_wind_output_dir = os.path.join(output_dir, wind_dir, folder)
            os.makedirs(tile_wind_output_dir, exist_ok=True)
            tif_output = os.path.join(tile_wind_output_dir, f"{folder}_base.tif")
            generated_tif = make_tif(speed_file, direction_file, prj_file, tif_output)
            if not generated_tif:
                print(f"[WARNING] Could not generate TIF for {folder}, wind_dir={wind_dir}.")
                continue
            
            # Attempt to find a matching AOI tile shapefile
            matching_aoi_tile = next(
                (
                    os.path.join(aoi_tiles_dir, shp)
                    for shp in os.listdir(aoi_tiles_dir)
                    if f"utm_overlap_aoi_{folder}_" in shp
                ),
                None
            )
            if not matching_aoi_tile or not os.path.exists(matching_aoi_tile):
                print(f"[WARNING] No matching AOI tile found for folder {folder}.")
                continue
            
            clipped_tif = os.path.join(tile_wind_output_dir, f"{folder}_clipped.tif")
            reprojected_tif = os.path.join(tile_wind_output_dir, f"{folder}_reproj.tif")
            final_resampled_raster = os.path.join(tile_wind_output_dir, f"{folder}.tif")

            clip_and_resample_and_reproject(
                generated_tif,
                matching_aoi_tile,
                clipped_tif,
                reprojected_tif,
                final_resampled_raster
            )

            remove_temp_files(clipped_tif)
            remove_temp_files(reprojected_tif)
            remove_temp_files(tif_output)

def process_tiles(base_dir, aoi_tiles_dir, output_dir, directions):
    """Processes tiles in parallel using multiprocessing"""
    os.makedirs(output_dir, exist_ok=True)
    tile_folders = [f for f in os.listdir(base_dir) if f.isdigit() and int(f) in tiles]
    num_workers = 8
    with multiprocessing.Pool(processes=num_workers) as pool:
        pool.starmap(process_tile, [(folder, base_dir, aoi_tiles_dir, output_dir, directions) for folder in tile_folders])
    print("[INFO] Processing complete.")

def main():
    """
    Process WindNinja tiles based on input base directory, AOI (Area of Interest) tiles directory,
    and output directory for processed results.
    
    Arguments:
    - base_dir: Directory containing WindNinja directional outputs.
    - aoi_tiles_dir: Directory containing AOI tiles to clip/match WindNinja tiles.
    - output_dir: Directory where processed tiles will be saved.
    - directions: List of wind directions to process.
    """
    parser = argparse.ArgumentParser(description="Process WindNinja tiles for multiple directions.")
    parser.add_argument("--base_dir", required=True, help="Path to base WindNinja output directory")
    parser.add_argument("--aoi_tiles_dir", required=True, help="Path to directory containing AOI tiles")
    parser.add_argument("--output_dir", required=True, help="Path to output directory for processed files")
    args = parser.parse_args()

    base_dir = args.base_dir
    aoi_tiles_dir = args.aoi_tiles_dir
    output_dir = args.output_dir

    # List of wind directions (fixed as per WindNinja output)
    directions = [
        "0-0-deg", "22-5-deg", "45-0-deg", "67-5-deg", "90-0-deg",
        "112-5-deg", "135-0-deg", "157-5-deg", "180-0-deg", "202-5-deg",
        "225-0-deg", "247-5-deg", "270-0-deg", "292-5-deg", "315-0-deg", "337-5-deg"
    ]

    # Call your processing function (assumed to be defined elsewhere)
    process_tiles(base_dir, aoi_tiles_dir, output_dir, directions)

if __name__ == "__main__":
    main()
