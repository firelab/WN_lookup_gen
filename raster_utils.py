import os
import shutil
from shapely.geometry import box, mapping
import rioxarray as rxr
import geopandas as gpd
from rasterio.merge import merge
import rasterio
from concurrent.futures import ThreadPoolExecutor

# Raster processing utilities
def process_dem(subpoly, buffer_size, dem_file, output_dir, resolution, tile_index):
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

    # Reproject with explicit resolution and alignment
    reprojected_dem = clipped_dem.rio.reproject(
        "EPSG:32612",
        resolution=(resolution, resolution),
        align=True
    )

    # Create hierarchical folder structure
    parent_dir = os.path.join(output_dir, str(tile_index), "dems_folder", "dem0")
    os.makedirs(parent_dir, exist_ok=True)

    output_path = os.path.join(parent_dir, "dem0.tif")
    reprojected_dem.rio.to_raster(output_path, resolution=(resolution, resolution))

    print(f"Saved DEM to: {output_path}")
    return output_path

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

def mosaic_dem_tiles(input_dir, target_crs, resolution, batch_size=500):
    """
    Mosaic a large number of DEM tiles in batches to avoid memory issues.
    
    Args:
        input_dir (str): Directory containing DEM tiles.
        target_crs (str): CRS for the output mosaic.
        resolution (tuple): Resolution for the mosaic.
        batch_size (int): Number of tiles to process in each batch.
    """
    dem_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.tif'):
                dem_files.append(os.path.join(root, file))

    if not dem_files:
        raise FileNotFoundError(f"No .tif files found in directory: {input_dir}")

    print(f"Found {len(dem_files)} DEM files in the directory structure.")
    
    # Divide files into batches
    num_batches = (len(dem_files) + batch_size - 1) // batch_size
    intermediate_mosaics = []

    for batch_idx in range(num_batches):
        print(f"Processing batch {batch_idx + 1} of {num_batches}...")
        batch_files = dem_files[batch_idx * batch_size : (batch_idx + 1) * batch_size]

        # Read and merge the batch
        src_files_to_mosaic = []
        for dem_file in batch_files:
            try:
                src = rasterio.open(dem_file)
                src_files_to_mosaic.append(src)
            except rasterio.errors.RasterioIOError:
                print(f"Error opening file: {dem_file}. Skipping...")

        if not src_files_to_mosaic:
            print(f"No valid files found in batch {batch_idx + 1}. Skipping batch.")
            continue

        print("Merging batch...")
        mosaic, out_trans = merge(src_files_to_mosaic, res=resolution)

        # Update metadata for the intermediate mosaic
        out_meta = src_files_to_mosaic[0].meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": mosaic.shape[1],
            "width": mosaic.shape[2],
            "transform": out_trans,
            "crs": target_crs
        })

        # Save the intermediate mosaic
        intermediate_path = os.path.join(input_dir, f"batch_{batch_idx + 1}_mosaic.tif")
        with rasterio.open(intermediate_path, "w", **out_meta) as dest:
            dest.write(mosaic)

        print(f"Saved intermediate mosaic: {intermediate_path}")
        intermediate_mosaics.append(intermediate_path)

        # Close the source files
        for src in src_files_to_mosaic:
            src.close()

    # Merge intermediate mosaics into the final output
    print("Merging all intermediate mosaics into final output...")
    src_files_to_mosaic = [rasterio.open(fp) for fp in intermediate_mosaics]
    mosaic, out_trans = merge(src_files_to_mosaic, res=resolution)

    # Update metadata for the final mosaic
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": target_crs
    })

    final_output_path = os.path.join(input_dir, "final_mosaicked_dem.tif")
    with rasterio.open(final_output_path, "w", **out_meta) as dest:
        dest.write(mosaic)

    print(f"Saved final mosaic: {final_output_path}")

    # Close all intermediate sources
    for src in src_files_to_mosaic:
        src.close()

    # Optionally clean up intermediate mosaics
    for fp in intermediate_mosaics:
        os.remove(fp)
        print(f"Deleted intermediate mosaic: {fp}")

    return final_output_path

def load_and_prepare_boundary(boundary_file, target_crs, dem_bounds):
    """Load and reproject the boundary file, then filter by DEM bounds."""
    boundary_gdf = gpd.read_file(boundary_file, engine='pyogrio')
    print("Original boundary file CRS:", boundary_gdf.crs)

    # Filter to include only the US (excluding Alaska and Hawaii)
    boundary_gdf = boundary_gdf[(boundary_gdf['gu_a3'] == 'USA') & 
                                (boundary_gdf['name'] != "Alaska") & 
                                (boundary_gdf['name'] != "Hawaii")]
    print(f"Filtered {len(boundary_gdf)} polygons for the contiguous US.")

    # Reproject to match DEM CRS
    boundary_gdf = boundary_gdf.to_crs(target_crs)

    # Filter to include only boundaries intersecting DEM bounds
    dem_bbox = box(*dem_bounds)
    boundary_gdf = boundary_gdf[boundary_gdf.intersects(dem_bbox)]
    print(f"Filtered boundaries intersecting DEM: {len(boundary_gdf)}")

    return boundary_gdf

def read_band(band_path):
    """Read a single band raster and return its data."""
    with rasterio.open(band_path) as src:
        data = src.read(1)
        print(f"Read band from {band_path}")
        return data

def create_multiband_raster(output_path, band_paths):
    if len(band_paths) < 2:
        raise ValueError("At least two band paths are required to create a multiband raster.")
    
    # Open the first band to get metadata
    with rasterio.open(band_paths[0]) as src:
        meta = src.meta.copy()
        meta.update(count=len(band_paths))
    
    # Read all bands in parallel
    with ThreadPoolExecutor() as executor:
        bands = list(executor.map(read_band, band_paths))
    
    # Write the multiband raster
    with rasterio.open(output_path, "w", **meta) as dest:
        for idx, band_data in enumerate(bands, start=1):
            dest.write(band_data, idx)
            print(f"Added band {idx} to the multiband raster.")
    
    print(f"Multiband raster saved to: {output_path}")
