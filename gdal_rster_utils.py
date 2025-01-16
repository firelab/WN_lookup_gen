import os
import numpy as np
from osgeo import gdal, osr
from concurrent.futures import ThreadPoolExecutor
import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt

def plot_tile(raster_path, output_png_path, title):
    """Plot a raster and save as PNG."""
    try:
        print(f"[DEBUG] Plotting raster: {raster_path}")
        with rasterio.open(raster_path) as raster:
            fig, ax = plt.subplots(figsize=(10, 8))
            show(raster, ax=ax, title=title)
            plt.savefig(output_png_path)
            plt.close()
        print(f"[DEBUG] Saved plot: {output_png_path}")
    except Exception as e:
        print(f"[ERROR] Failed to plot raster {raster_path}: {e}")

def replace_ocean_values(array, band_index):
    replacements = [0, 0, 0, 98, 0, 0, 0, 0]
    print(f"[DEBUG] Replacing ocean values in band {band_index}")
    return np.where(array == -9999, replacements[band_index - 1], array)

def create_tile(input_raster, output_dir, reproj_dir, x_offset, y_offset, tile_size_x, tile_size_y, tile_counter, utm_zone, buffer_percent):
    print(f"[DEBUG] Creating tile {tile_counter} at offset ({x_offset}, {y_offset})")
    x_end = x_offset + tile_size_x
    y_end = y_offset + tile_size_y

    temp_tile_path = f"/tmp/temp_tile_{tile_counter:05d}.tif"
    try:
        print(f"[DEBUG] Extracting tile to: {temp_tile_path}")
        gdal.Translate(
            temp_tile_path,
            input_raster,
            srcWin=[x_offset, y_offset, x_end - x_offset, y_end - y_offset]
        )

        dataset = gdal.Open(temp_tile_path)
        if not dataset:
            print(f"[ERROR] Failed to create tile at ({x_offset}, {y_offset}). Skipping.")
            return None

        band_count = dataset.RasterCount
        band_type = dataset.GetRasterBand(1).DataType
        replaced_tile_path = f"/tmp/replaced_tile_{tile_counter:05d}.tif"
        driver = gdal.GetDriverByName("GTiff")
        replaced_dataset = driver.Create(replaced_tile_path, dataset.RasterXSize, dataset.RasterYSize, band_count, band_type)

        print(f"[DEBUG] Replacing ocean values in tile {tile_counter}")
        for i in range(1, band_count + 1):
            band = dataset.GetRasterBand(i)
            array = band.ReadAsArray()
            replaced_array = replace_ocean_values(array, i)
            replaced_band = replaced_dataset.GetRasterBand(i)
            replaced_band.WriteArray(replaced_array)

            if band_type in [gdal.GDT_Float32, gdal.GDT_Float64]:
                replaced_band.SetNoDataValue(float(-9999))
            else:
                replaced_band.SetNoDataValue(int(-9999))

        replaced_dataset.SetGeoTransform(dataset.GetGeoTransform())
        replaced_dataset.SetProjection(dataset.GetProjection())
        replaced_dataset.FlushCache()
        replaced_dataset = None
        dataset = None

        reproj_tile_path = os.path.join(reproj_dir, f"reprojected_tile_{tile_counter:05d}.tif")
        print(f"[DEBUG] Reprojecting tile {tile_counter} to: {reproj_tile_path}")
        gdal.Warp(reproj_tile_path, replaced_tile_path, dstSRS=f"EPSG:{utm_zone}")

        clip_tile_path = os.path.join(reproj_dir, f"clipped_tile_{tile_counter:05d}.tif")
        print(f"[DEBUG] Clipping tile {tile_counter} to: {clip_tile_path}")
        raster_info = gdal.Info(reproj_tile_path, format="json")
        clip_pixels_x = int(buffer_percent * tile_size_x)
        clip_pixels_y = int(buffer_percent * tile_size_y)
        gdal.Translate(
            clip_tile_path,
            reproj_tile_path,
            srcWin=[
                clip_pixels_x,
                clip_pixels_y,
                tile_size_x - 2 * clip_pixels_x,
                tile_size_y - 2 * clip_pixels_y,
            ]
        )

        # Create a folder for graphs
        graphs_dir = os.path.join(reproj_dir, "graphs")
        os.makedirs(graphs_dir, exist_ok=True)

        plot_tile(temp_tile_path, os.path.join(graphs_dir, f"original_tile_{tile_counter:05d}.png"), "Original Tile")
        plot_tile(reproj_tile_path, os.path.join(graphs_dir, f"reprojected_tile_{tile_counter:05d}.png"), "Reprojected Tile")
        plot_tile(clip_tile_path, os.path.join(graphs_dir, f"final_clipped_tile_{tile_counter:05d}.png"), "Final Clipped Tile")

    except Exception as e:
        print(f"[ERROR] Error processing tile {tile_counter}: {e}")
    finally:
        if os.path.exists(temp_tile_path):
            os.remove(temp_tile_path)
        if os.path.exists(replaced_tile_path):
            os.remove(replaced_tile_path)

def generate_tiles(input_raster, output_dir, reproj_dir, tile_size_km=32, utm_zone=None, max_tiles=1000):
    km_to_m = 1000
    buffer_percent = 0.45

    print(f"[DEBUG] Preparing directories")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(reproj_dir, exist_ok=True)

    dataset = gdal.Open(input_raster)
    geo_transform = dataset.GetGeoTransform()
    pixel_width = geo_transform[1]
    pixel_height = abs(geo_transform[5])
    raster_x_size = dataset.RasterXSize
    raster_y_size = dataset.RasterYSize

    tile_size_m = tile_size_km * km_to_m
    tile_size_x = int(tile_size_m / pixel_width)
    tile_size_y = int(tile_size_m / pixel_height)

    print(f"[DEBUG] Tile size: {tile_size_x}x{tile_size_y} pixels")
    print(f"[DEBUG] Raster size: {raster_x_size}x{raster_y_size} pixels")

    x_offset = 0
    y_offset = 0
    tile_counter = 0

    with ThreadPoolExecutor(max_workers=8) as executor:
        tasks = []
        while y_offset < raster_y_size and tile_counter < max_tiles:
            x_offset = 0
            while x_offset < raster_x_size and tile_counter < max_tiles:
                print(f"[DEBUG] Scheduling tile {tile_counter} for creation")
                task = executor.submit(
                    create_tile,
                    input_raster,
                    output_dir,
                    reproj_dir,
                    x_offset,
                    y_offset,
                    tile_size_x,
                    tile_size_y,
                    tile_counter,
                    utm_zone,
                    buffer_percent
                )
                tasks.append(task)
                tile_counter += 1
                x_offset += tile_size_x

            y_offset += tile_size_y

        for task in tasks:
            task.result()

    print(f"[DEBUG] Generated {tile_counter} tiles.")

# Example Usage
input_raster_path = "/mnt/c/Users/dgh00/OneDrive/Desktop/CONUS2022/2023_lcp.tif"
output_tiles_dir = "/mnt/d/gdal_out/original_tiles"
reprojected_tiles_dir = "/mnt/d/gdal_out/reprojected_tiles"
generate_tiles(input_raster_path, output_tiles_dir, reprojected_tiles_dir, tile_size_km=350, utm_zone=32612, max_tiles=25)
