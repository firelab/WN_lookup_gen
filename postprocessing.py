import os
import numpy as np
import subprocess
from osgeo import gdal, ogr, osr

def remove_temp_files(*file_patterns):
    """Removes temporary files matching the given patterns."""
    for file_pattern in file_patterns:
        if os.path.exists(file_pattern):
            os.remove(file_pattern)
            print(f"[INFO] Removed temporary file: {file_pattern}")

def read_asc_file(asc_path):
    """Reads an ASCII raster file and returns a NumPy array, geotransform, and projection."""
    ds = gdal.Open(asc_path)
    if not ds:
        raise ValueError(f"Could not open {asc_path}")
    
    array = ds.GetRasterBand(1).ReadAsArray()
    geotransform = ds.GetGeoTransform()
    projection = ds.GetProjection()
    ds = None
    return array, geotransform, projection

def read_prj_file(prj_path):
    """Reads a .prj file and returns the spatial reference as WKT."""
    with open(prj_path, 'r') as prj_file:
        wkt = prj_file.read()
    return wkt

def get_raster_info(raster_path):
    """Extracts bounds, resolution, and EPSG code from a raster file."""
    ds = gdal.Open(raster_path)
    if not ds:
        raise ValueError(f"Could not open {raster_path}")
    
    geotransform = ds.GetGeoTransform()
    projection = ds.GetProjection()
    width, height = ds.RasterXSize, ds.RasterYSize

    # Calculate bounding box
    min_x = geotransform[0]
    max_y = geotransform[3]
    max_x = min_x + (width * geotransform[1])
    min_y = max_y + (height * geotransform[5])  # geotransform[5] is negative

    # Extract EPSG code
    srs = osr.SpatialReference(wkt=projection)
    epsg = srs.GetAttrValue("AUTHORITY", 1) if srs.IsProjected() else None

    ds = None
    return (min_x, min_y, max_x, max_y), width, height, epsg

def write_multiband_raster(output_path, arrays, geotransform, projection):
    """Creates a multi-band raster from input arrays."""
    rows, cols = arrays[0].shape
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output_path, cols, rows, len(arrays), gdal.GDT_Float32)

    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(projection)

    for i, array in enumerate(arrays):
        band = dst_ds.GetRasterBand(i + 1)
        band.WriteArray(array)
        band.SetNoDataValue(-9999)

    dst_ds = None
    print(f"[INFO] Multi-band raster saved: {output_path}")

def generate_aoi_mask(aoi_raster, output_shapefile):
    """Generate a mask shapefile from the AOI raster using the CLI command (more reliable)."""
    print(f"[INFO] Generating AOI mask from {aoi_raster}...")
    
    command = f"gdal_polygonize.py '{aoi_raster}' -f 'ESRI Shapefile' '{output_shapefile}'"
    subprocess.run(command, shell=True, check=True)

    print(f"[INFO] AOI mask saved as {output_shapefile}")

def clip_and_resample_and_reproject(input_raster, output_raster, aoi_mask, aoi_raster, reprojected_raster):
    """Clip and resample input raster to match AOI resolution and bounding box, then reproject to CONUS Albers (EPSG:5070)"""
    os.makedirs(os.path.dirname(output_raster), exist_ok=True)
    os.makedirs(os.path.dirname(reprojected_raster), exist_ok=True)

    # Extract AOI information dynamically
    bounds, width, height, epsg = get_raster_info(aoi_raster)

    if epsg is None:
        raise ValueError(f"Could not determine EPSG code for {aoi_raster}")

    print(f"[INFO] Clipping and resampling {input_raster} to match AOI {aoi_raster}...")

    # Open input raster
    ds = gdal.Open(input_raster)
    if not ds:
        raise ValueError(f"Could not open {input_raster}")

    options = gdal.WarpOptions(
        outputBounds=bounds,
        width=width,
        height=height,
        resampleAlg=gdal.GRA_Bilinear,
        cutlineDSName=aoi_mask,
        cropToCutline=True,
        dstNodata=-9999
    )
    
    gdal.Warp(output_raster, ds, options=options)
    ds = None
    print(f"[INFO] Clipped and resampled raster saved: {output_raster}")

    print(f"[INFO] Reprojecting {output_raster} to CONUS Albers (EPSG:5070)...")

    options_reproj = gdal.WarpOptions(
        dstSRS="EPSG:5070",
        resampleAlg=gdal.GRA_Bilinear,
        dstNodata=-9999
    )
    
    gdal.Warp(reprojected_raster, output_raster, options=options_reproj)
    print(f"[INFO] Reprojected raster saved: {reprojected_raster}")

def process_tiles(base_dir, aoi_tiles_dir, output_dir):
    """Iterate over UTM tiles, generate multiband TIFFs, clip, and resample."""
    os.makedirs(output_dir, exist_ok=True)
    
    for folder in sorted(os.listdir(base_dir)):
        folder_path = os.path.join(base_dir, folder, "dems_folder", "output")
        if not os.path.exists(folder_path):
            continue

        tile_index = folder

        for file in os.listdir(folder_path):
            if file.endswith("_vel.asc"):
                base_name = file.replace("_vel.asc", "")
                speed_file = os.path.join(folder_path, f"{base_name}_vel.asc")
                direction_file = os.path.join(folder_path, f"{base_name}_ang.asc")
                cloud_file = os.path.join(folder_path, f"{base_name}_cld.asc")
                prj_file = os.path.join(folder_path, f"{base_name}_vel.prj")

                if not (os.path.exists(speed_file) and os.path.exists(direction_file) and os.path.exists(cloud_file) and os.path.exists(prj_file)):
                    print(f"[WARNING] Missing files in {folder_path}. Skipping.")
                    continue
                
                # Read ASC files
                speed_array, geotransform, _ = read_asc_file(speed_file)
                direction_array, _, _ = read_asc_file(direction_file)
                cloud_array, _, _ = read_asc_file(cloud_file)

                # Read PRJ file
                projection_wkt = read_prj_file(prj_file)

                # Output path for multiband raster
                output_raster = os.path.join(folder_path, f"{base_name}_multiband.tif")
                write_multiband_raster(output_raster, [speed_array, direction_array, cloud_array], geotransform, projection_wkt)

                # Find corresponding AOI tile by matching index
                matching_aoi_tile = next((os.path.join(aoi_tiles_dir, file) for file in os.listdir(aoi_tiles_dir) if f"utm_overlap_aoi_{tile_index}_" in file), None)

                if os.path.exists(matching_aoi_tile):
                    final_tile = os.path.join(output_dir, folder, "final_clipped_multiband.tif")
                    final_reproj_tile = os.path.join(output_dir, folder, "final_reproj_multiband.tif")
                    
                    aoi_mask_shp = f"/tmp/aoi_mask_{tile_index}.shp"
                    remove_temp_files(aoi_mask_shp)
                    generate_aoi_mask(matching_aoi_tile, aoi_mask_shp)
                    clip_and_resample_and_reproject(output_raster, final_tile, aoi_mask_shp, matching_aoi_tile, final_reproj_tile)
                    remove_temp_files(aoi_mask_shp)
                else:
                    print(f"[WARNING] No matching AOI tile found for tile index {tile_index}")
    
    print("[INFO] Postprocessing complete.")

if __name__ == "__main__":
    process_tiles("/home/gunjan/Desktop/washington", "/home/gunjan/Desktop/washington_prep_dems/utm_aoi_tiles", "/home/gunjan/Desktop/washington_prep_dems/final_clipped_tiles")
