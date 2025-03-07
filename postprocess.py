"""
 * Project:  WN_Lookup_Gen
 * Purpose:  Python script for converting tiles(Windninja asc and prj files) from UTM to conus albers
 * Author: Gunjan Dayani <gunjandayani015@gmail.com>

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

gdal.UseExceptions()

tiles = set(range(350))

def remove_temp_files(*file_patterns):
    for file_pattern in file_patterns:
        if os.path.exists(file_pattern):
            os.remove(file_pattern)
            print(f"[INFO] Removed temporary file: {file_pattern}")

def make_tif(vel_asc, ang_asc, prj_file, output_tif):
    """
    Converts velocity and direction ASC files into a two-band GeoTIFF.
    - Band 1: Wind speed (_vel.asc)
    - Band 2: Wind direction (_ang.asc)
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
    speed_array = np.loadtxt(vel_asc, skiprows=6)
    direction_array = np.loadtxt(ang_asc, skiprows=6)
    os.makedirs(os.path.dirname(output_tif), exist_ok=True)
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(
        output_tif,
        metadata["ncols"], metadata["nrows"],
        2,
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
    band1.WriteArray(speed_array)
    band1.SetNoDataValue(-9999)
    band2 = out_ds.GetRasterBand(2)
    band2.WriteArray(direction_array)
    band2.SetNoDataValue(-9999)
    band1, band2, out_ds = None, None, None
    print(f"[INFO] Created TIF: {output_tif}")

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
    half_extent = 40000
    min_x = center_x - half_extent
    max_x = center_x + half_extent
    min_y = center_y - half_extent
    max_y = center_y + half_extent
    ds = None

    print(f"[INFO] Cropping raster to 80km x 80km centered at ({center_x}, {center_y})")
    print(f"[INFO] New bounding box: ({min_x}, {min_y}, {max_x}, {max_y})")

    gdal.Warp(final_raster, final_raster, outputBounds=(min_x, min_y, max_x, max_y), dstNodata=-9999)
    print(f"[INFO] Cropped {final_raster} to 80km x 80km centered region.")

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

    print(f"[INFO] Cropping {reprojected_raster} to 80km * 80km and get rid of NoData on each side...")
    crop_to_original_tile(reprojected_raster)
    print(f"[INFO] Forcing resolution to 120m x 120m for {reprojected_raster}...")
    options_force_res = gdal.WarpOptions(
        xRes=120,
        yRes=120,
        resampleAlg=gdal.GRA_Bilinear,
        dstNodata=-9999
    )
    gdal.Warp(final_resampled_raster, reprojected_raster, options=options_force_res)
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

            aoi_mask_shp = f"/tmp/aoi_mask_{folder}.shp"
            remove_temp_files(aoi_mask_shp)

            # Generate a mask from the AOI tile
            generate_aoi_mask(matching_aoi_tile, aoi_mask_shp)

            clip_and_resample_and_reproject(
                generated_tif,
                clipped_tif,
                aoi_mask_shp,
                matching_aoi_tile,
                reprojected_tif,
                final_resampled_raster
            )

            remove_temp_files(aoi_mask_shp)
            remove_temp_files(clipped_tif)
            remove_temp_files(reprojected_tif)
            remove_temp_files(tif_output)

def process_tiles(base_dir, aoi_tiles_dir, output_dir, directions):
    """Processes tiles in parallel using multiprocessing"""
    os.makedirs(output_dir, exist_ok=True)
    tile_folders = [f for f in os.listdir(base_dir) if f.isdigit() and int(f) in tiles]
    num_workers = min(16, os.cpu_count())
    with multiprocessing.Pool(processes=num_workers) as pool:
        pool.starmap(process_tile, [(folder, base_dir, aoi_tiles_dir, output_dir, directions) for folder in tile_folders])
    print("[INFO] Processing complete.")

if __name__ == "__main__":
    base_dir = "/mnt/d/CONUS0"
    aoi_tiles_dir = "/mnt/d/tiles/utm_aoi_tiles"
    output_dir = "/mnt/d/processed_CONUS0"
    directions = [
        "0-0-deg", "22-5-deg", "45-0-deg", "67-5-deg", "90-0-deg",
        "112-5-deg", "135-0-deg", "157-5-deg", "180-0-deg", "202-5-deg",
        "225-0-deg", "247-5-deg", "270-0-deg", "292-5-deg", "315-0-deg", "337-5-deg"
    ]
    process_tiles(base_dir, aoi_tiles_dir, output_dir, directions)
