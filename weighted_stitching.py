"""
 * Project:  WN_Lookup_Gen
 * Purpose:  Python script for weighted stitching of tiles
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

 * Wind Direction Mosaic Processing Script

 * This script processes 16 wind direction mosaic grids by:
    - Identifying tiles overlapping the current resampled tile.
    - Computing tile weights based on distance from the mosaic cell center.
    - Converting wind speed and direction into U and V components.
    - Applying weights to U and V components.
    - Accumulating weighted sums of U, V, and total weights.
    - Normalizing U and V using weighted averages with Band 3.
    - Converting weighted U and V back into wind speed and direction:
        - Band 1: Wind speed
        - Band 2: Wind direction
        - Band 3: Sum of weights
"""

import os
import shutil
from osgeo import gdal, osr
import numpy as np
from multiprocessing import Pool, cpu_count
import csv
import argparse

def find_overlapping_tiles(tiles_dir, mosaic_metadata):
    """Find overlapping tiles directly from .tif files using parent folder as tile_id."""
    
    # Get all .tif files in the tiles directory
    tif_files = []
    
    for root, _, files in os.walk(tiles_dir):
        for file in files:
            if file.endswith(".tif"):
                tif_files.append(os.path.join(root, file))
    
    print(f"Found {len(tif_files)} TIF files in {tiles_dir}.")
    
    overlapping_tiles = {}

    for tif_file in tif_files:
        ds = gdal.Open(tif_file)
        if ds is None:
            print(f"[WARNING] Failed to open TIF file: {tif_file}, skipping.")
            continue

        # Extract metadata
        geo_transform = ds.GetGeoTransform()
        tile_x_size, tile_y_size = ds.RasterXSize, ds.RasterYSize
        tile_cell_size = abs(geo_transform[1])
        tile_x_min, tile_y_max = geo_transform[0], geo_transform[3]
        tile_x_max = tile_x_min + tile_x_size * tile_cell_size
        tile_y_min = tile_y_max - tile_y_size * tile_cell_size

        ds = None  # Close dataset

        # Extract the parent folder name as the tile_id
        tile_id = os.path.basename(os.path.dirname(tif_file))  # 👈 Use folder name instead of filename

        # Mosaic boundaries
        mosaic_x_min, mosaic_y_min = mosaic_metadata["xllcorner"], mosaic_metadata["yllcorner"]
        mosaic_x_max = mosaic_x_min + mosaic_metadata["ncols"] * mosaic_metadata["cellsize"]
        mosaic_y_max = mosaic_y_min + mosaic_metadata["nrows"] * mosaic_metadata["cellsize"]

        # Compute affected mosaic cells
        col_start = round((tile_x_min - mosaic_x_min) / mosaic_metadata["cellsize"])
        col_end = round((tile_x_max - mosaic_x_min) / mosaic_metadata["cellsize"])
        row_start = round((mosaic_y_max - tile_y_max) / mosaic_metadata["cellsize"])
        row_end = round((mosaic_y_max - tile_y_min) / mosaic_metadata["cellsize"])

        col_start, col_end = max(0, col_start), min(mosaic_metadata["ncols"], col_end)
        row_start, row_end = max(0, row_start), min(mosaic_metadata["nrows"], row_end)

        num_cells = (row_end - row_start) * (col_end - col_start)  # Number of affected cells
        print(f"[DEBUG] Tile {tile_id} affects {num_cells} cells.")

        for row in range(row_start, row_end):
            for col in range(col_start, col_end):
                cell_id = (row, col)
                if cell_id not in overlapping_tiles:
                    overlapping_tiles[cell_id] = []
                overlapping_tiles[cell_id].append(tile_id)

    print(f"Found {len(overlapping_tiles)} overlapping cells.")

    for cell, tiles in overlapping_tiles.items():
        if len(tiles) > 4:
            print(f"[WARNING] Cell {cell} has {len(tiles)} tiles, exceeding expected max (4)")
    
    print("Finished finding tiles for each cell.")

    return overlapping_tiles

def test_find_overlapping_tiles(input_raster, tiles_dir, output_raster):
    """ Generates a raster where each cell contains the count of overlapping tiles. """
    print(f"Loading input raster: {input_raster}")
    src_ds = gdal.Open(input_raster)
    if src_ds is None:
        raise RuntimeError(f"Failed to open input raster: {input_raster}")

    geo_transform = list(src_ds.GetGeoTransform())
    projection = src_ds.GetProjection()
    
    mosaic_metadata = {
        "ncols": src_ds.RasterXSize,
        "nrows": src_ds.RasterYSize,
        "cellsize": abs(geo_transform[1]),
        "xllcorner": geo_transform[0],
        "yllcorner": geo_transform[3] - src_ds.RasterYSize * abs(geo_transform[5]),
    }

    print(f"Mosaic Grid Info: {mosaic_metadata}")

    overlapping_tiles = find_overlapping_tiles(tiles_dir, mosaic_metadata)
    
    count_array = np.zeros((mosaic_metadata["nrows"], mosaic_metadata["ncols"]), dtype=np.uint16)

    for (row, col), tiles in overlapping_tiles.items():
        count_array[row, col] = len(tiles)

    print(f"[DEBUG] Total cells affected in mosaic: {np.count_nonzero(count_array)}")

    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(
        output_raster,
        mosaic_metadata["ncols"],
        mosaic_metadata["nrows"],
        1,
        gdal.GDT_UInt16,
        options=["COMPRESS=LZW", "BIGTIFF=YES"]
    )
    out_ds.SetGeoTransform(geo_transform)
    out_ds.SetProjection(projection)

    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(count_array)
    out_band.SetNoDataValue(0)

    out_band = None
    out_ds = None
    src_ds = None

    print(f"Overlap count raster saved: {output_raster}")

def generate_weight_matrix():
    rows, cols = 1067, 1067
    weight_matrix = np.zeros((rows, cols), dtype=np.float32)
    ignore_edge_width_pixels = 30
    fully_weighted_radius = 240

    center_x, center_y = rows // 2, cols // 2

    for i in range(rows):
        for j in range(cols):
            dist_x = abs(i - center_x)
            dist_y = abs(j - center_y)
            max_dist = max(dist_x, dist_y)
            if max_dist >= (center_x - ignore_edge_width_pixels):
                weight_matrix[i, j] = 0.0
            elif max_dist <= fully_weighted_radius:
                weight_matrix[i, j] = 1.0
            else:
                weight_matrix[i, j] = (center_x - max_dist - ignore_edge_width_pixels) / (center_x - fully_weighted_radius)

    return weight_matrix

def get_mosaic_cell_center(mosaic_metadata, row, col):
    x_center = mosaic_metadata["xllcorner"] + (col + 0.5) * mosaic_metadata["cellsize"]
    y_center = mosaic_metadata["yllcorner"] + (mosaic_metadata["nrows"] - row - 0.5) * mosaic_metadata["cellsize"]
    return x_center, y_center

def convert_to_latlon(x, y, src_epsg=5070, dest_epsg=4326):
    """Convert projected coordinates to latitude/longitude."""
    src_srs = osr.SpatialReference()
    src_srs.ImportFromEPSG(src_epsg)  # Source coordinate system (e.g., UTM)
    
    dest_srs = osr.SpatialReference()
    dest_srs.ImportFromEPSG(dest_epsg)  # WGS84 (Lat/Lon)
    
    transform = osr.CoordinateTransformation(src_srs, dest_srs)
    lon, lat, _ = transform.TransformPoint(x, y)
    
    return lon, lat

def wind_uv_to_sd(u, v):
    speed = np.sqrt(u**2 + v**2)
    inter_dir = np.degrees(np.arctan2(v, u)) - 180.0
    inter_dir = np.where(inter_dir < 0, inter_dir + 360.0, inter_dir)
    direction = (450.0 - inter_dir) % 360.0  # This is WindNinja's xy_to_n()
    return speed, direction

def process_single_tile(args):
    tile_path, mosaic_metadata, weight_matrix = args

    ds_tile = gdal.Open(tile_path)
    if ds_tile is None:
        print(f"[WARNING] Failed to open tile: {tile_path}")
        return None

    tile_id = os.path.basename(os.path.dirname(tile_path))

    tile_geo_transform = ds_tile.GetGeoTransform()
    tile_x_min, tile_y_max = tile_geo_transform[0], tile_geo_transform[3]
    tile_x_size, tile_y_size = ds_tile.RasterXSize, ds_tile.RasterYSize
    tile_cell_size = abs(tile_geo_transform[1])
    tile_x_max = tile_x_min + tile_x_size * tile_cell_size
    tile_y_min = tile_y_max - tile_y_size * tile_cell_size

    col_start = max(0, round((tile_x_min - mosaic_metadata["xllcorner"]) / mosaic_metadata["cellsize"]))
    col_end = min(mosaic_metadata["ncols"], round((tile_x_max - mosaic_metadata["xllcorner"]) / mosaic_metadata["cellsize"]))
    row_start = max(0, round((mosaic_metadata["yllcorner"] + mosaic_metadata["nrows"] * mosaic_metadata["cellsize"] - tile_y_max) / mosaic_metadata["cellsize"]))
    row_end = min(mosaic_metadata["nrows"], round((mosaic_metadata["yllcorner"] + mosaic_metadata["nrows"] * mosaic_metadata["cellsize"] - tile_y_min) / mosaic_metadata["cellsize"]))

    tile_u = ds_tile.GetRasterBand(1).ReadAsArray()
    tile_v = ds_tile.GetRasterBand(2).ReadAsArray()

    tile_u = tile_u.astype(np.float32)
    tile_v = tile_v.astype(np.float32)

    tile_x_center = (tile_x_min + tile_x_max) / 2
    tile_y_center = (tile_y_min + tile_y_max) / 2

    results = []
    for row in range(row_start, row_end):
        for col in range(col_start, col_end):
            cell_x, cell_y = get_mosaic_cell_center(mosaic_metadata, row, col)

            tile_pixel_x = int((cell_x - tile_x_min) / tile_cell_size)
            tile_pixel_y = int((tile_y_max - cell_y) / tile_cell_size)

            if 0 <= tile_pixel_x < tile_x_size and 0 <= tile_pixel_y < tile_y_size:
                u = tile_u[tile_pixel_y, tile_pixel_x]
                v = tile_v[tile_pixel_y, tile_pixel_x]

                weight = weight_matrix[tile_pixel_y, tile_pixel_x]

                if weight == 0:
                    continue

                results.append((row, col, u * weight, v * weight, weight, tile_id))

    ds_tile = None
    print(f"[INFO] Finished processing tile: {tile_path}")
    return results

def process_overlapping_tiles(tiles_dir, mosaic_metadata, output_raster):
    print("[INFO] Processing tiles using multiprocessing...")

    ds = gdal.Open(output_raster, gdal.GA_Update)
    if ds is None:
        raise RuntimeError(f"Failed to open output raster: {output_raster}")

    velocity_band = np.full((mosaic_metadata["nrows"], mosaic_metadata["ncols"]), -9999, dtype=np.float32)
    direction_band = np.full((mosaic_metadata["nrows"], mosaic_metadata["ncols"]), -9999, dtype=np.float32)
    weight_band = np.full((mosaic_metadata["nrows"], mosaic_metadata["ncols"]), -9999, dtype=np.float32)

    u_sum_band = np.zeros_like(velocity_band, dtype=np.float32)
    v_sum_band = np.zeros_like(direction_band, dtype=np.float32)
    sum_weights_band = np.zeros_like(weight_band, dtype=np.float32)

    tile_paths = [
        os.path.join(tiles_dir, tile_id, f"{tile_id}.tif")
        for tile_id in os.listdir(tiles_dir)
        if os.path.exists(os.path.join(tiles_dir, tile_id, f"{tile_id}.tif"))
    ]

    weight_matrix = generate_weight_matrix()

    args = [(path, mosaic_metadata, weight_matrix) for path in tile_paths]

    tile_contributions = {}

    with Pool(processes=30) as pool:
        for result in pool.imap_unordered(process_single_tile, args):
            if result is None:
                continue

            for row, col, u_weighted, v_weighted, weight, tile_id in result:
                if sum_weights_band[row, col] == 0:
                    u_sum_band[row, col] = 0.0
                    v_sum_band[row, col] = 0.0
                    sum_weights_band[row, col] = 0.0

                if debug:
                    if (row, col) not in tile_contributions:
                        tile_contributions[(row, col)] = 0
        
                    tile_contributions[(row, col)] += 1

                u_sum_band[row, col] += u_weighted
                v_sum_band[row, col] += v_weighted
                sum_weights_band[row, col] += weight

    print(f"[DEBUG] Before Normalization -> Sum U: {u_sum_band.min()} to {u_sum_band.max()}")
    print(f"[DEBUG] Before Normalization -> Sum V: {v_sum_band.min()} to {v_sum_band.max()}")

    valid_mask = sum_weights_band > 0
    u_avg = np.zeros_like(u_sum_band)
    v_avg = np.zeros_like(v_sum_band)
    u_avg[valid_mask] = u_sum_band[valid_mask] / sum_weights_band[valid_mask]
    v_avg[valid_mask] = v_sum_band[valid_mask] / sum_weights_band[valid_mask]

    velocity_band[valid_mask], direction_band[valid_mask] = wind_uv_to_sd(u_avg[valid_mask], v_avg[valid_mask])

    # Set NoData values
    velocity_band[~valid_mask] = -9999
    direction_band[~valid_mask] = -9999
    weight_band[valid_mask] = sum_weights_band[valid_mask]
    weight_band[~valid_mask] = -9999

    # Write to raster bands
    ds.GetRasterBand(1).WriteArray(velocity_band)
    ds.GetRasterBand(1).SetNoDataValue(-9999)
    print("[INFO] Velocity band written")

    ds.GetRasterBand(2).WriteArray(direction_band)
    ds.GetRasterBand(2).SetNoDataValue(-9999)
    print("[INFO] Direction band written")

    # ds.GetRasterBand(3).WriteArray(weight_band)
    # ds.GetRasterBand(3).SetNoDataValue(-9999)
    # print("[INFO] Weight band written")

    ds = None
    print("[INFO] Mosaic processing completed successfully!")

def process_uv_components_only(tiles_dir, mosaic_metadata, output_raster):
    print("[INFO] Processing tiles to produce U/V component bands...")

    ds = gdal.Open(output_raster, gdal.GA_Update)
    if ds is None:
        raise RuntimeError(f"Failed to open output raster: {output_raster}")

    u_band = np.full((mosaic_metadata["nrows"], mosaic_metadata["ncols"]), -9999, dtype=np.float32)
    v_band = np.full((mosaic_metadata["nrows"], mosaic_metadata["ncols"]), -9999, dtype=np.float32)
    weight_band = np.full((mosaic_metadata["nrows"], mosaic_metadata["ncols"]), -9999, dtype=np.float32)

    u_sum_band = np.zeros_like(u_band, dtype=np.float32)
    v_sum_band = np.zeros_like(v_band, dtype=np.float32)
    sum_weights_band = np.zeros_like(weight_band, dtype=np.float32)

    tile_paths = [
        os.path.join(tiles_dir, tile_id, f"{tile_id}.tif")
        for tile_id in os.listdir(tiles_dir)
        if os.path.exists(os.path.join(tiles_dir, tile_id, f"{tile_id}.tif"))
    ]

    weight_matrix = generate_weight_matrix()
    args = [(path, mosaic_metadata, weight_matrix) for path in tile_paths]

    with Pool(processes=30) as pool:
        for result in pool.imap_unordered(process_single_tile, args):
            if result is None:
                continue

            for row, col, u_weighted, v_weighted, weight, tile_id in result:
                if sum_weights_band[row, col] == 0:
                    u_sum_band[row, col] = 0.0
                    v_sum_band[row, col] = 0.0
                    sum_weights_band[row, col] = 0.0

                u_sum_band[row, col] += u_weighted
                v_sum_band[row, col] += v_weighted
                sum_weights_band[row, col] += weight

    valid_mask = sum_weights_band > 0
    u_band[valid_mask] = u_sum_band[valid_mask] / sum_weights_band[valid_mask]
    v_band[valid_mask] = v_sum_band[valid_mask] / sum_weights_band[valid_mask]

    # Set NoData values
    u_band[~valid_mask] = -9999
    v_band[~valid_mask] = -9999
    weight_band[valid_mask] = sum_weights_band[valid_mask]
    weight_band[~valid_mask] = -9999

    # Write U/V to raster bands
    ds.GetRasterBand(1).WriteArray(u_band)
    ds.GetRasterBand(1).SetNoDataValue(-9999)
    print("[INFO] U component band written")

    ds.GetRasterBand(2).WriteArray(v_band)
    ds.GetRasterBand(2).SetNoDataValue(-9999)
    print("[INFO] V component band written")

    # ds.GetRasterBand(3).WriteArray(weight_band)
    # ds.GetRasterBand(3).SetNoDataValue(-9999)
    # print("[INFO] Weight band written")

    ds = None
    print("[INFO] U/V mosaic processing completed successfully!")

def main():
    """
    Create mosaics for each wind direction by overlaying processed tiles onto a reference raster.

    Arguments:
    - base_raster: Path to the reference raster (defines output grid, projection, resolution).
    - output_directory: Directory where final mosaics will be saved.
    - tiles_base_directory: Directory containing subdirectories for each wind direction's processed tiles.
    """
    parser = argparse.ArgumentParser(description="Create WindNinja CONUS mosaics from processed tiles.")
    parser.add_argument("--base_raster", required=True, help="Path to the base reference raster (e.g., ref_raster_120m.tif)")
    parser.add_argument("--output_directory", required=True, help="Directory to save output mosaics")
    parser.add_argument("--tiles_base_directory", required=True, help="Directory containing tiles organized by wind direction")
    parser.add_argument("--uv", action="store_true", help="Generate U/V component bands instead of Speed/Direction (default is Spd/Dir)")
    args = parser.parse_args()

    base_raster = args.base_raster
    output_directory = args.output_directory
    tiles_base_directory = args.tiles_base_directory
    generate_uv = args.uv

    # Hardcoded list of wind directions (always same)
    wind_directions = [
        "0-0-deg", "22-5-deg", "45-0-deg", "67-5-deg", "90-0-deg",
        "112-5-deg", "135-0-deg", "157-5-deg", "180-0-deg", "202-5-deg",
        "225-0-deg", "247-5-deg", "270-0-deg", "292-5-deg", "315-0-deg", "337-5-deg"
    ]

    # Open base raster to extract metadata
    src_ds = gdal.Open(base_raster)
    if src_ds is None:
        raise RuntimeError(f"Failed to open base raster: {base_raster}")

    geo_transform = list(src_ds.GetGeoTransform())
    projection = src_ds.GetProjection()
    
    mosaic_metadata = {
        "ncols": src_ds.RasterXSize,
        "nrows": src_ds.RasterYSize,
        "cellsize": abs(geo_transform[1]),
        "xllcorner": geo_transform[0],
        "yllcorner": geo_transform[3] - src_ds.RasterYSize * abs(geo_transform[5]),
    }

    print(f"Mosaic Grid Info: {mosaic_metadata}")

    os.makedirs(output_directory, exist_ok=True)

    for direction in wind_directions:
        output_raster = os.path.join(output_directory, f"{direction}.tif")
        tiles_directory = os.path.join(tiles_base_directory, direction)

        # Ensure the tiles directory exists
        if not os.path.exists(tiles_directory):
            print(f"[WARNING] Tiles directory not found: {tiles_directory}, skipping {direction}.")
            continue

        print(f"[INFO] Processing direction: {direction}")

        # Create BigTIFF raster from base metadata
        driver = gdal.GetDriverByName("GTiff")
        out_ds = driver.Create(
            output_raster,
            mosaic_metadata["ncols"],
            mosaic_metadata["nrows"],
            3,  # Velocity, Direction, weight bands
            gdal.GDT_Float32,
            options=["COMPRESS=LZW", "BIGTIFF=YES"]
        )
        out_ds.SetGeoTransform(geo_transform)
        out_ds.SetProjection(projection)

        # Fill with NoData initially
        for b in range(1, 3):
            band = out_ds.GetRasterBand(b)
            band.Fill(-9999)
            band.SetNoDataValue(-9999)

        out_ds = None

        if generate_uv:
            process_uv_components_only(tiles_directory, mosaic_metadata, output_raster)
        else:
            process_overlapping_tiles(tiles_directory, mosaic_metadata, output_raster)

if __name__ == "__main__":
    main()
