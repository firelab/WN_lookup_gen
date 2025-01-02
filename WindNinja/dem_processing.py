import rasterio
from shapely.geometry import box, mapping
from rasterio.merge import merge
import rioxarray as rxr
import os
import numpy as np
from scipy.ndimage import generic_filter
import shutil

def resample_tif(input_path, target_resolution):
    dem = rxr.open_rasterio(input_path, masked=True)
    x_res, y_res = dem.rio.resolution()
    scale_x = target_resolution / abs(x_res)
    scale_y = target_resolution / abs(y_res)
    resampled_dem = dem.rio.reproject(
        dem.rio.crs,
        shape=(int(dem.shape[1] / scale_y), int(dem.shape[2] / scale_x)),
    )
    output_path = os.path.join(
        os.path.dirname(input_path),
        f"{os.path.splitext(os.path.basename(input_path))[0]}_{int(target_resolution)}m.tif"
    )
    resampled_dem.rio.to_raster(output_path)
    print(f"Resampled DEM saved to: {output_path}")
    return output_path

def replace_water_and_fill_nodata(data, water_value=-9999, nodata_value=None, replacement_value=0, fill=True):
    """
    Replace water values (-9999) with a fixed value (e.g., 0) and optionally fill NoData values.
    Also, print the count of NoData values before and after processing.
    """
    # Count NoData values before processing
    nodata_count_before = np.sum(data == nodata_value) if nodata_value is not None else 0
    water_count_before = np.sum(data == water_value)
    print(f"Count of NoData values before processing: {nodata_count_before}")
    print(f"Count of water (-9999) values before processing: {water_count_before}")

    # Replace water values (-9999) with the replacement value (e.g., 0)
    water_mask = data == water_value
    data[water_mask] = replacement_value

    # Replace NoData values (if provided) with NaN for processing
    if nodata_value is not None:
        nodata_mask = data == nodata_value
        data[nodata_mask] = np.nan  # Mark NoData values as NaN

    # Fill NoData values using a moving average filter if requested
    if fill:
        def fill_function(window):
            valid_pixels = window[~np.isnan(window)]  # Exclude NaN
            return np.mean(valid_pixels) if len(valid_pixels) > 0 else np.nan

        data = generic_filter(
            data, fill_function, size=3, mode='constant', cval=np.nan
        )

    # Replace any remaining NaN values with the replacement value
    data = np.nan_to_num(data, nan=replacement_value)

    # Count NoData values after processing
    nodata_count_after = np.sum(data == nodata_value) if nodata_value is not None else 0
    print(f"Count of NoData values after processing: {nodata_count_after}")

    return data

def process_dem(subpoly, buffer_size, dem_file, output_dir, resolution):
    buffered_subpoly = subpoly.buffer(buffer_size)

    dem_bounds = dem_file.rio.bounds()
    dem_bbox = box(*dem_bounds)

    if not buffered_subpoly.intersects(dem_bbox):
        print("Polygon does NOT intersect DEM bounds. Skipping.")
        return None

    buffered_subpoly = buffered_subpoly.intersection(dem_bbox)

    buffered_geojson = mapping(buffered_subpoly)
    print("Clipping with polygon bounds:", buffered_subpoly.bounds)

    clipped_dem = dem_file.rio.clip([buffered_geojson], crs=dem_file.rio.crs, drop=True)

    reprojected_dem = clipped_dem.rio.reproject(
        "EPSG:32612",
        resolution=(resolution, resolution),
        align=True
    )

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, "processed_dem.tif")
    reprojected_dem.rio.to_raster(output_path, resolution=(resolution, resolution))

    return output_path

def mosaic_dem_tiles(input_dir, target_crs, resolution):
    dem_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.tif'):
                dem_files.append(os.path.join(root, file))

    if not dem_files:
        raise FileNotFoundError(f"No .tif files found in directory: {input_dir}")

    print(f"Found {len(dem_files)} DEM files in the directory structure.")

    src_files_to_mosaic = []
    for dem_file in dem_files:
        try:
            src = rasterio.open(dem_file)
            src_files_to_mosaic.append(src)
        except rasterio.errors.RasterioIOError:
            print(f"Error opening file: {dem_file}. Skipping...")

    if not src_files_to_mosaic:
        raise ValueError(f"No valid DEM files found in directory: {input_dir}")

    print("Merging DEM tiles...")
    mosaic, out_trans = merge(src_files_to_mosaic, res=resolution)

    out_meta = src_files_to_mosaic[0].meta.copy()

    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": target_crs
    })

    output_path = os.path.join(input_dir, "mosaicked_dem.tif")
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(mosaic)

    print(f"Saved mosaicked DEM to: {output_path}")

    return output_path
