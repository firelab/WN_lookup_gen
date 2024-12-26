import os
from pathlib import Path
import numpy as np
from rasterio import open as rio_open
from rasterio.merge import merge
from rasterio.windows import Window
from scipy.ndimage import gaussian_filter

def km_to_pixels(km, resolution):
    """
    Converts kilometers to pixels based on the raster resolution.
    """
    return int((km * 1000) / resolution)

def tile_raster_with_overlap(input_raster, tile_size_km, overlap_km, output_folder):
    """
    Tiles a raster into smaller chunks with overlap, supporting NoData values.
    """
    with rio_open(input_raster) as src:
        resolution = abs(src.transform[0])  # Pixel size in meters
        tile_size = km_to_pixels(tile_size_km, resolution)
        overlap = km_to_pixels(overlap_km, resolution)
        
        os.makedirs(output_folder, exist_ok=True)
        nodata_value = src.nodata  # Retrieve NoData value from the source raster

        for i in range(0, src.width, tile_size - overlap):
            for j in range(0, src.height, tile_size - overlap):
                width = min(tile_size, src.width - i)
                height = min(tile_size, src.height - j)
                window = Window(i, j, width, height)
                tile_data = src.read(window=window)

                # Save the tile
                output_file = os.path.join(output_folder, f"tile_{i}_{j}.tif")
                profile = src.profile
                profile.update({
                    "width": width,
                    "height": height,
                    "transform": src.window_transform(window),
                    "nodata": nodata_value  # Preserve NoData value
                })
                with rio_open(output_file, "w", **profile) as dst:
                    dst.write(tile_data)
                print(f"Saved tile: {output_file}")

def process_tile_smoothing(input_tile, output_tile, sigma):
    """
    Smooths a raster tile using a Gaussian filter, handling NoData values.
    """
    with rio_open(input_tile) as src:
        profile = src.profile
        nodata_value = src.nodata
        data = src.read()

        # Mask NoData values
        mask = data == nodata_value
        data = np.where(mask, 0, data)  # Replace NoData with 0 for processing

        # Apply Gaussian filter to each band
        smoothed_data = np.stack([gaussian_filter(band, sigma=sigma) for band in data], axis=0)

        # Restore NoData values
        smoothed_data = np.where(mask, nodata_value, smoothed_data)

        # Save the smoothed tile
        profile.update(dtype=smoothed_data.dtype)
        with rio_open(output_tile, "w", **profile) as dst:
            dst.write(smoothed_data)
        print(f"Smoothed tile saved: {output_tile}")

def trim_overlap(input_tile, output_tile, overlap_km, resolution, is_edge_tile=False):
    """
    Trims the overlap from a raster tile. Does not trim if it's an edge tile.
    """
    overlap = km_to_pixels(overlap_km, resolution)

    with rio_open(input_tile) as src:
        profile = src.profile
        data = src.read()

        # If the tile is an edge tile, skip trimming
        if is_edge_tile:
            trimmed_data = data
        else:
            # Trim the overlap
            start_x = overlap // 2
            end_x = -overlap // 2 if overlap > 0 else None
            start_y = overlap // 2
            end_y = -overlap // 2 if overlap > 0 else None
            trimmed_data = data[:, start_y:end_y, start_x:end_x]

        # Skip invalid tiles
        if trimmed_data.shape[1] <= 0 or trimmed_data.shape[2] <= 0:
            print(f"Skipping invalid trimmed tile: {input_tile}")
            return

        # Update profile with new dimensions
        profile.update({
            "height": trimmed_data.shape[1],
            "width": trimmed_data.shape[2]
        })
        with rio_open(output_tile, "w", **profile) as dst:
            dst.write(trimmed_data)
        print(f"Trimmed tile saved: {output_tile}")

def merge_tiles(tiles_folder, output_raster):
    """
    Merges processed tiles into a single raster, preserving NoData values.
    """
    tiles = [os.path.join(tiles_folder, f) for f in os.listdir(tiles_folder) if f.endswith(".tif")]
    datasets = [rio_open(tile) for tile in tiles]
    mosaic, out_transform = merge(datasets, nodata=datasets[0].nodata)  # Use consistent NoData value
    profile = datasets[0].profile
    profile.update({
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_transform,
        "nodata": datasets[0].nodata  # Preserve NoData value
    })

    # Save the merged raster
    with rio_open(output_raster, "w", **profile) as dst:
        dst.write(mosaic)
    for ds in datasets:
        ds.close()
    print(f"Merged raster saved: {output_raster}")
