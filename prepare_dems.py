import os
import shutil
import numpy as np
import multiprocessing as mp
from osgeo import gdal, osr, ogr

def get_nodata_value(raster_path):
    """Extract NoData value from the first band of a raster file."""
    ds = gdal.Open(raster_path)
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    ds = None
    print(f"[INFO] Detected NoData Value: {nodata}")
    return nodata

def determine_utm_zone(lon, lat):
    """Determine the correct UTM zone from longitude (degrees) and latitude."""
    if lon < -180 or lon > 180:
        print(f"[ERROR] Invalid longitude detected: {lon}. Expected degrees, but got meters?")
        return None

    utm_zone = int((lon + 186) / 6)
    epsg_code = 32600 + utm_zone  # EPSG for UTM North
    print(f"[DEBUG] Determined UTM Zone: {epsg_code} for Lon: {lon}, Lat: {lat}")
    return epsg_code

def get_raster_bounds(raster_path):
    """Returns the bounding box of the input raster (min_x, min_y, max_x, max_y)."""
    ds = gdal.Open(raster_path)
    transform = ds.GetGeoTransform()
    width, height = ds.RasterXSize, ds.RasterYSize
    min_x = transform[0]
    max_x = min_x + width * transform[1]
    max_y = transform[3]
    min_y = max_y + height * transform[5]
    ds = None
    return min_x, min_y, max_x, max_y

def replace_ocean_values(array, band_index, nodata_value):
    """Replaces ocean values (-9999) while preserving other data."""
    replacements = [0, 0, 0, 98, 0, 0, 0, 0]
    return np.where((array == -9999) & (array != nodata_value), replacements[band_index - 1], array)

def albers_to_latlon(x, y):
    """Convert Albers Equal Area coordinates to Lat/Lon using GDAL transformation."""
    source_srs = osr.SpatialReference()
    source_srs.ImportFromEPSG(5070)  # CONUS Albers
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(4326)  # WGS 84 (Lat/Lon)
    transform = osr.CoordinateTransformation(source_srs, target_srs)
    lon, lat, _ = transform.TransformPoint(x, y)
    print(f"[DEBUG] Converted Albers to Lat/Lon: ({x}, {y}) -> ({lat}, {lon})")
    return lon, lat

def clean_output_directory(output_dir):
    """Delete and recreate the output directory to ensure a clean start."""
    if os.path.exists(output_dir):
        print(f"[INFO] Cleaning output directory: {output_dir}")
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

def find_bounding_square(utm_bounds):
    """Find the smallest square in UTM that fully contains the reprojected AOI."""
    min_x, min_y, max_x, max_y = utm_bounds
    width = max_x - min_x
    height = max_y - min_y
    square_size = max(width, height)  # Ensure AOI fits within a square

    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2

    new_min_x = center_x - square_size / 2
    new_max_x = center_x + square_size / 2
    new_min_y = center_y - square_size / 2
    new_max_y = center_y + square_size / 2

    return (new_min_x, new_min_y, new_max_x, new_max_y)

def overlap_and_reproject(input_raster, tile_bounds, output_dir, tile_index, overlap):
    """Process raster tiles with specified overlap, reproject to UTM, and ensure square bounds."""
    nodata_value = get_nodata_value(input_raster)
    raster_min_x, raster_min_y, raster_max_x, raster_max_y = get_raster_bounds(input_raster)

    temp_albers_dir = os.path.join(output_dir, "temp_albers_tiles")
    utm_overlap_dir = os.path.join(output_dir, "utm_aoi_tiles")
    final_utm_dir = os.path.join(output_dir, "final_utm_tiles")

    os.makedirs(temp_albers_dir, exist_ok=True)
    os.makedirs(utm_overlap_dir, exist_ok=True)
    os.makedirs(final_utm_dir, exist_ok=True)

    # Compute center of the tile
    min_x, min_y, max_x, max_y = tile_bounds
    center_x, center_y = (min_x + max_x) / 2, (min_y + max_y) / 2

    # Compute 350% Overlapped Bounds in Albers
    # As we are not taking whole input into memory, after deciding not tilted tile bounds(perfect square) 
    # in UTM we need to have enough data from which we can clip data. So this is kind of temporary lookup tiles
    # I have calculated this number 350 so that we have enough data to clip from.
    expanded_width = (max_x - min_x) * 4.5
    expanded_height = (max_y - min_y) * 4.5

    expanded_bounds = (
        center_x - expanded_width / 2,
        center_y - expanded_height / 2,
        center_x + expanded_width / 2,
        center_y + expanded_height / 2,
    )
    print(f"[INFO] 350% Expanded Bounds in Albers: {expanded_bounds}")

    # Compute AOI bounds based on specified overlap percentage
    aoi_width = (max_x - min_x) * (1 + overlap / 100)
    aoi_height = (max_y - min_y) * (1 + overlap / 100)

    aoi_bounds = (
        center_x - aoi_width / 2,
        center_y - aoi_height / 2,
        center_x + aoi_width / 2,
        center_y + aoi_height / 2,
    )
    print(f"[INFO] {overlap}% AOI Bounds in Albers: {aoi_bounds}")

    # Save Temporary Overlapped Albers Tiles(350% and given overlap percentage) and replace Ocean values
    # So that we don't have -9999 values (other than NoData: 32767) which is present in input landscape file. 
    temp_350_albers = os.path.join(temp_albers_dir, f"temp_350_albers_{tile_index}.tif")
    gdal.Warp(
        temp_350_albers,
        input_raster,
        outputBounds=expanded_bounds,
        dstSRS="EPSG:5070",
        dstNodata=nodata_value,
        warpOptions=['INIT_DEST=-9999'],
        resampleAlg=gdal.GRA_NearestNeighbour
    )
    ds = gdal.Open(temp_350_albers, gdal.GA_Update)
    for band_index in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_index)
        array = band.ReadAsArray()
        array = replace_ocean_values(array, band_index, nodata_value)
        band.WriteArray(array)
    ds = None

    temp_overlap_albers = os.path.join(temp_albers_dir, f"temp_overlap_albers_{tile_index}.tif")
    gdal.Warp(
        temp_overlap_albers,
        input_raster,
        outputBounds=aoi_bounds,
        dstSRS="EPSG:5070",
        dstNodata=nodata_value,
        warpOptions=['INIT_DEST=-9999'],
        resampleAlg=gdal.GRA_NearestNeighbour
    )
    ds = gdal.Open(temp_overlap_albers, gdal.GA_Update)
    for band_index in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_index)
        array = band.ReadAsArray()
        array = replace_ocean_values(array, band_index, nodata_value)
        band.WriteArray(array)
    ds = None

    # Determine UTM Zone After 350% Expansion
    min_lon, min_lat = albers_to_latlon(expanded_bounds[0], expanded_bounds[1])
    max_lon, max_lat = albers_to_latlon(expanded_bounds[2], expanded_bounds[3])
    center_lon, center_lat = (min_lon + max_lon) / 2, (min_lat + max_lat) / 2
    epsg_utm = determine_utm_zone(center_lon, center_lat)

    print("Determined UTM Zone: " + str(epsg_utm))
    
    # Reproject 350% Overlapped Tile to UTM
    temp_utm_350 = os.path.join(temp_albers_dir, f"utm_350_tile_{tile_index}_epsg{epsg_utm}.tif")
    gdal.Warp(
        temp_utm_350,
        temp_350_albers,
        dstSRS=f"EPSG:{epsg_utm}",
        dstNodata=nodata_value,
        warpOptions=['INIT_DEST=NO_DATA'],
        resampleAlg=gdal.GRA_NearestNeighbour
    )

    # Reproject given percentage overlapped AOI Tile to UTM and force 30m resolution
    temp_utm_overlap_path = os.path.join(utm_overlap_dir, f"temp_utm_overlap_aoi_{tile_index}_epsg{epsg_utm}.tif")
    gdal.Warp(
        temp_utm_overlap_path,
        temp_overlap_albers,
        dstSRS=f"EPSG:{epsg_utm}",
        dstNodata=nodata_value,
        warpOptions=['INIT_DEST=NO_DATA'],
        resampleAlg=gdal.GRA_NearestNeighbour
    )
    utm_overlap_path = os.path.join(utm_overlap_dir, f"utm_overlap_aoi_{tile_index}_epsg{epsg_utm}.tif")
    gdal.Warp(
        utm_overlap_path,
        temp_utm_overlap_path,
        xRes=30,
        yRes=30,
        dstSRS=f"EPSG:{epsg_utm}",
        dstNodata=nodata_value,
        resampleAlg=gdal.GRA_NearestNeighbour
    )
    os.remove(temp_utm_overlap_path)
    
    # Compute Final Square and not tilted dem in UTM
    utm_overlap_ds = gdal.Open(utm_overlap_path)
    transform_overlap = utm_overlap_ds.GetGeoTransform()
    width_overlap, height_overlap = utm_overlap_ds.RasterXSize, utm_overlap_ds.RasterYSize
    pixel_size_x, pixel_size_y = abs(transform_overlap[1]), abs(transform_overlap[5])

    utm_overlap_min_x = transform_overlap[0]
    utm_overlap_max_x = utm_overlap_min_x + width_overlap * pixel_size_x
    utm_overlap_max_y = transform_overlap[3]
    utm_overlap_min_y = utm_overlap_max_y - height_overlap * pixel_size_y

    final_bounds = find_bounding_square((utm_overlap_min_x, utm_overlap_min_y, utm_overlap_max_x, utm_overlap_max_y))

    # Crop to Final Square UTM Tile
    parent_dir = os.path.join(final_utm_dir, str(tile_index))
    base_output_dir = os.path.join(parent_dir, "dems_folder", "dem0")
    os.makedirs(base_output_dir, exist_ok=True)
    final_utm_path = os.path.join(base_output_dir, "dem0_temp.tif")

    gdal.Warp(
        final_utm_path,
        temp_utm_350,
        outputBounds=final_bounds,
        dstSRS=f"EPSG:{epsg_utm}",
        dstNodata=nodata_value,
        warpOptions=['INIT_DEST=NO_DATA'],
        resampleAlg=gdal.GRA_NearestNeighbour
    )
    final_resampled_utm_path = os.path.join(base_output_dir, "dem0.tif")
    gdal.Warp(
        final_resampled_utm_path,
        final_utm_path,
        xRes=30,
        yRes=30,
        dstSRS=f"EPSG:{epsg_utm}",
        dstNodata=nodata_value,
        resampleAlg=gdal.GRA_NearestNeighbour
    )

    os.remove(final_utm_path)
    os.remove(temp_350_albers)
    os.remove(temp_overlap_albers)
    os.remove(temp_utm_350)

    print(f"[INFO] Final Square UTM Tile Saved: {final_resampled_utm_path}")

def create_tiles(input_raster, output_dir, tile_size_km, overlap_percentage, max_tiles=None, num_workers=8):
    """Generate raster tiles with overlap and reproject to UTM in batches of 8 tiles at a time."""

    clean_output_directory(output_dir)
    original_tiles_dir = os.path.join(output_dir, "original_tiles")
    os.makedirs(original_tiles_dir, exist_ok=True)

    dataset = gdal.Open(input_raster)
    transform = dataset.GetGeoTransform()
    pixel_size = transform[1]
    tile_size_px = int((tile_size_km * 1000) / pixel_size)

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

                if data is not None and np.any(data != -9999):
                    contains_valid_data = True
                    data = replace_ocean_values(data, band_index, nodata_value)
                    band.WriteArray(data)

            tile_ds = None

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
            print(f"[INFO] Tile {tile_count} completed.")
            tile_count += 1

            # Process tiles every 8 iterations
            if len(tile_data_list) == 8 or tile_count >= max_tiles:
                with mp.Pool(processes=num_workers) as pool:
                    pool.starmap(overlap_and_reproject, tile_data_list)
                tile_data_list = []  # Clear the list for the next batch

    print(f"[INFO] Total tiles generated: {tile_count}/{estimated_total_tiles}")

if __name__ == "__main__":
    input_raster_path =  "/mnt/c/Users/dgh00/OneDrive/Desktop/CONUS2022/2023_lcp.tif"
    output_directory = "/mnt/d/tiles_washington"
    create_tiles(input_raster_path, output_directory, 64, 25, 350, num_workers=8)
