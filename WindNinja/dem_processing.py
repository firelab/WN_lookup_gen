import rasterio
from shapely.geometry import box, mapping
from rasterio.merge import merge
import rioxarray as rxr
import os
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

def process_dem(subpoly, buffer_size, dem_file, output_dir):
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

    reprojected_dem = clipped_dem.rio.reproject("EPSG:32612")

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, "processed_dem.tif")
    reprojected_dem.rio.to_raster(output_path)

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

    for src in src_files_to_mosaic:
        src.close()

    return output_path
