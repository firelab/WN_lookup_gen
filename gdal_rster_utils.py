import os
from osgeo import gdal
from concurrent.futures import ThreadPoolExecutor

def create_tile(input_raster, output_dir, reproj_dir, x_offset, y_offset, tile_size_x, tile_size_y, overlap_x, overlap_y, meaningful_tile_counter, utm_zone, nodata_value):
    x_end = x_offset + tile_size_x
    y_end = y_offset + tile_size_y

    temp_tile_path = f"/tmp/temp_tile_{meaningful_tile_counter:05d}.tif"
    gdal.Translate(
        temp_tile_path,
        input_raster,
        srcWin=[x_offset, y_offset, x_end - x_offset, y_end - y_offset]
    )

    temp_dataset = gdal.Open(temp_tile_path)
    temp_band = temp_dataset.GetRasterBand(1)
    stats = temp_band.GetStatistics(False, True)
    if stats[0] == stats[1] == nodata_value:
        print(f"Tile at ({x_offset}, {y_offset}) contains only NoData values. Skipping.")
        os.remove(temp_tile_path)
        return None

    os.remove(temp_tile_path)

    # Create directory structure for original tile
    parent_dir = os.path.join(output_dir, str(meaningful_tile_counter), "dems_folder", "dem0")
    os.makedirs(parent_dir, exist_ok=True)
    tile_output_path = os.path.join(parent_dir, "dem0.tif")
    gdal.Translate(
        tile_output_path,
        input_raster,
        srcWin=[x_offset, y_offset, x_end - x_offset, y_end - y_offset]
    )
    print(f"Created original tile: {tile_output_path}")

    # Create directory structure for reprojected tile
    if utm_zone:
        reproj_parent_dir = os.path.join(reproj_dir, str(meaningful_tile_counter), "dems_folder", "dem0")
        os.makedirs(reproj_parent_dir, exist_ok=True)
        reproj_output_path = os.path.join(reproj_parent_dir, "dem0.tif")
        gdal.Warp(
            reproj_output_path,
            tile_output_path,
            dstSRS=f"EPSG:{utm_zone}"
        )
        print(f"Created reprojected tile: {reproj_output_path}")

    return meaningful_tile_counter


def generate_tiles(input_raster, output_dir, reproj_dir, tile_size_km=32, overlap_km=16, utm_zone=None, max_tiles=10, nodata_value=-9999):
    km_to_m = 1000  # Convert kilometers to meters

    # Ensure output directories exist
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(reproj_dir, exist_ok=True)

    # Open the raster and get metadata
    dataset = gdal.Open(input_raster)
    geo_transform = dataset.GetGeoTransform()
    pixel_width = geo_transform[1]
    pixel_height = geo_transform[5]
    raster_x_size = dataset.RasterXSize
    raster_y_size = dataset.RasterYSize

    # Calculate tile size in pixels
    tile_size_m = tile_size_km * km_to_m
    overlap_m = overlap_km * km_to_m
    tile_size_x = int(tile_size_m / abs(pixel_width))
    tile_size_y = int(tile_size_m / abs(pixel_height))
    overlap_x = int(overlap_m / abs(pixel_width))
    overlap_y = int(overlap_m / abs(pixel_height))

    print(f"Tile size in pixels: {tile_size_x} x {tile_size_y}")
    print(f"Overlap in pixels: {overlap_x} x {overlap_y}")

    tile_counter = 0
    meaningful_tile_counter = 0
    y_offset = 0
    tasks = []

    # Parallel processing with ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=8) as executor:
        while y_offset < raster_y_size and meaningful_tile_counter < max_tiles:
            x_offset = 0
            while x_offset < raster_x_size and meaningful_tile_counter < max_tiles:
                task = executor.submit(
                    create_tile,
                    input_raster,
                    output_dir,
                    reproj_dir,
                    x_offset,
                    y_offset,
                    tile_size_x,
                    tile_size_y,
                    overlap_x,
                    overlap_y,
                    meaningful_tile_counter,
                    utm_zone,
                    nodata_value
                )
                tasks.append(task)
                meaningful_tile_counter += 1
                x_offset += tile_size_x - overlap_x
            y_offset += tile_size_y - overlap_y

        # Wait for all tasks to complete
        for task in tasks:
            result = task.result()
            if result is not None:
                tile_counter += 1

    print(f"Generated {tile_counter} meaningful tiles in total.")

# Example Usage
input_raster_path = "/mnt/c/Users/dgh00/OneDrive/Desktop/CONUS2022/2023_lcp.tif"
generate_tiles(input_raster_path,"/mnt/d/gdal_out/original_tiles_500km_250km", "/mnt/d/gdal_out/reprojected_tiles_500km_250km", tile_size_km=500, overlap_km=250, utm_zone=32612, max_tiles=50)