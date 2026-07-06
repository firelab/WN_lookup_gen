"""
 * Project:  WN_Lookup_Gen
 * Purpose:  Python script for preparing overlapped DEM tiles for WindNinja
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
import shutil
import multiprocessing as mp
import numpy as np
from osgeo import gdal
import argparse

def get_nodata_value(raster_path):
    """Extract NoData value from the first band of a raster file."""
    ds = gdal.Open(raster_path)
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    ds = None
    print(f"[INFO] Detected NoData Value: {nodata}")
    return nodata

def clean_output_directory(output_dir):
    """Delete and recreate the output directory to ensure a clean start."""
    if os.path.exists(output_dir):
        print(f"[INFO] Cleaning output directory: {output_dir}")
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

def prepare_windninja_dem(input_raster, tile_bounds, output_dir, tile_index, overlap):
    """Create a WindNinja DEM tile with the requested overlap."""
    nodata_value = get_nodata_value(input_raster)

    final_tiles_dir = os.path.join(output_dir, "final_tiles")
    os.makedirs(final_tiles_dir, exist_ok=True)

    # Compute center of the tile
    min_x, min_y, max_x, max_y = tile_bounds
    center_x, center_y = (min_x + max_x) / 2, (min_y + max_y) / 2

    # Compute tile bounds with the requested overlap percentage.
    aoi_width = (max_x - min_x) * (1 + overlap / 100)
    aoi_height = (max_y - min_y) * (1 + overlap / 100)

    aoi_bounds = (
        center_x - aoi_width / 2,
        center_y - aoi_height / 2,
        center_x + aoi_width / 2,
        center_y + aoi_height / 2,
    )
    print(f"[INFO] {overlap}% overlapped tile bounds: {aoi_bounds}")

    parent_dir = os.path.join(final_tiles_dir, str(tile_index))
    base_output_dir = os.path.join(parent_dir, "dems_folder", "dem0")
    os.makedirs(base_output_dir, exist_ok=True)
    final_resampled_path = os.path.join(base_output_dir, "dem0.tif")

    gdal.Warp(
        final_resampled_path,
        input_raster,
        outputBounds=aoi_bounds,
        xRes=30,
        yRes=30,
        dstNodata=nodata_value,
        targetAlignedPixels=True,
        resampleAlg=gdal.GRA_NearestNeighbour
    )

    ds = gdal.Open(final_resampled_path, gdal.GA_Update)
    for band_index in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_index)
        arr = band.ReadAsArray()

        # Replace rogue -9999 values and enforce correct NoData value.
        arr[(arr == -9999) | (arr == band.GetNoDataValue())] = nodata_value
        band.WriteArray(arr)
        band.SetNoDataValue(nodata_value)

    ds.FlushCache()
    ds = None

    print(f"[INFO] Final overlapped DEM tile saved: {final_resampled_path}")

def create_tiles(input_raster, output_dir, tile_size_km, overlap_percentage, max_tiles=None, num_workers=8):
    """Generate raster tiles with overlap in batches of 8 tiles at a time."""

    clean_output_directory(output_dir)
    original_tiles_dir = os.path.join(output_dir, "original_tiles")
    os.makedirs(original_tiles_dir, exist_ok=True)

    dataset = gdal.Open(input_raster)
    transform = dataset.GetGeoTransform()
    pixel_size = transform[1]
    tile_size_px = int((tile_size_km * 1000) / pixel_size)

    # Compute total number of tiles in x and y directions
    x_size, y_size = dataset.RasterXSize, dataset.RasterYSize
    nodata_value = get_nodata_value(input_raster)

    print(f"[DEBUG] Raster size: {x_size} x {y_size} pixels")
    print(f"[DEBUG] Tile size: {tile_size_px} x {tile_size_px} pixels")

    total_tiles_x = (x_size + tile_size_px - 1) // tile_size_px
    total_tiles_y = (y_size + tile_size_px - 1) // tile_size_px
    estimated_total_tiles = total_tiles_x * total_tiles_y

    print(f"[DEBUG] Estimated total tiles: {estimated_total_tiles} ({total_tiles_x} x {total_tiles_y})")

    if max_tiles is None:
        max_tiles = estimated_total_tiles

    tile_count = 0
    tile_data_list = []

    for x in range(0, x_size, tile_size_px):
        for y in range(0, y_size, tile_size_px):
            if tile_count >= max_tiles:
                break

            output_tile = os.path.join(original_tiles_dir, f"tile_{tile_count}_{x}_{y}.tif")
            gdal.Translate(output_tile, dataset, srcWin=[x, y, tile_size_px, tile_size_px])

            tile_ds = gdal.Open(output_tile, gdal.GA_Update)
            if not tile_ds:
                continue

            contains_valid_data = False
            for band_index in range(1, tile_ds.RasterCount + 1):
               band = tile_ds.GetRasterBand(band_index)
               data = band.ReadAsArray()
               data[data == -9999] = nodata_value
               if np.any(data != nodata_value):
                  contains_valid_data = True

            if not contains_valid_data:
                os.remove(output_tile)
                print(f"[INFO] Skipped tile {output_tile} (contains only NoData).")
                continue

            min_x = transform[0] + x * pixel_size
            max_x = min_x + tile_size_px * pixel_size
            max_y = transform[3] - y * abs(transform[5])
            min_y = max_y - tile_size_px * abs(transform[5])
            tile_bounds = (min_x, min_y, max_x, max_y)
            tile_data_list.append((input_raster, tile_bounds, output_dir, tile_count, overlap_percentage))
            print(f"[INFO] Tile {tile_count} queued.")
            tile_count += 1
            if len(tile_data_list) == 8 or tile_count >= max_tiles:
                with mp.Pool(processes=num_workers) as pool:
                    pool.starmap(prepare_windninja_dem, tile_data_list)
                tile_data_list = []  # Clear the list for the next batch

    print(f"[INFO] Total tiles generated: {tile_count}/{estimated_total_tiles}")

def main():
    """
    Create tiles from a large raster based on specified tile size and overlap.

    Arguments:
    - input_raster_path: Path to the large input raster file (e.g., CONUS LCP).
    - output_directory: Directory where the generated tiles will be saved.
    - tile_size: Size of each tile (in number of pixels).
    - overlap: Overlap between adjacent tiles (in number of pixels).
    - num_workers: Number of parallel workers for faster processing.
    """
    parser = argparse.ArgumentParser(description="Tile a large raster into smaller chunks with overlap.")
    parser.add_argument("--input_raster_path", required=True, help="Path to the input raster file")
    parser.add_argument("--output_directory", required=True, help="Directory to save output tiles")
    parser.add_argument("--tile_size", type=int, default=64, help="Size of each tile in km (default: 64)")
    parser.add_argument("--overlap", type=int, default=45, help="Overlap between tiles in percentage (default: 45)")
    parser.add_argument("--num_workers", type=int, default=8, help="Number of parallel workers (default: 8)")
    args = parser.parse_args()

    input_raster_path = args.input_raster_path
    output_directory = args.output_directory
    tile_size = args.tile_size
    overlap = args.overlap
    num_workers = args.num_workers

    # Call your tiling function (assumed to be defined elsewhere)
    create_tiles(input_raster_path, output_directory, tile_size, overlap, num_workers=num_workers)

if __name__ == "__main__":
    main()
