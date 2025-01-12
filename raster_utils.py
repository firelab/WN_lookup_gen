import os
import shutil
from shapely.geometry import box, mapping
import rioxarray as rxr
import geopandas as gpd
from rasterio.merge import merge
import rasterio
from concurrent.futures import ThreadPoolExecutor

def process_dem(subpoly, buffer_size, dem_data, output_dir, resolution, tile_index):
    """
    Process a multiband raster (xarray.DataArray) for a specific polygon, clipping and reprojecting all bands.

    Args:
        subpoly (shapely.geometry.Polygon): The polygon to clip.
        buffer_size (int): Buffer size in projection units.
        dem_data (xarray.DataArray): The multiband raster loaded via rioxarray.
        output_dir (str): Directory to save the output tile.
        resolution (tuple): Resolution for the output raster.
        tile_index (int): Index for naming the output tile.

    Returns:
        str: Path to the processed multiband tile.
    """
    buffered_subpoly = subpoly.buffer(buffer_size)

    # Get the bounds of the DEM data
    dem_bounds = dem_data.rio.bounds()
    dem_bbox = box(*dem_bounds)

    if not buffered_subpoly.intersects(dem_bbox):
        print("Polygon does NOT intersect DEM bounds. Skipping.")
        return None

    buffered_subpoly = buffered_subpoly.intersection(dem_bbox)
    buffered_geojson = mapping(buffered_subpoly)

    print(f"Clipping with polygon bounds: {buffered_subpoly.bounds}")

    # Clip all bands
    clipped_dem = dem_data.rio.clip([buffered_geojson], crs=dem_data.rio.crs, drop=True)

    # Reproject and save all bands
    reprojected_dem = clipped_dem.rio.reproject(
        "EPSG:32612",
        resolution=resolution,
        align=True
    )

    parent_dir = os.path.join(output_dir, str(tile_index), "dems_folder", "dem0")
    os.makedirs(parent_dir, exist_ok=True)

    output_path = os.path.join(parent_dir, "dem0.tif")
    reprojected_dem.rio.to_raster(output_path)
    print(f"Saved multiband DEM to: {output_path}")

    return output_path


def mosaic_dem_tiles(input_dir, target_crs, resolution, batch_size=500):
    """
    Mosaic multiband DEM tiles in batches.

    Args:
        input_dir (str): Directory containing DEM tiles.
        target_crs (str): CRS for the output mosaic.
        resolution (tuple): Resolution for the mosaic.
        batch_size (int): Number of tiles to process in each batch.

    Returns:
        str: Path to the final mosaicked multiband raster.
    """
    dem_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.tif'):
                dem_files.append(os.path.join(root, file))

    if not dem_files:
        raise FileNotFoundError(f"No .tif files found in directory: {input_dir}")

    print(f"Found {len(dem_files)} DEM files in the directory structure.")

    # Read metadata from the first tile
    with rasterio.open(dem_files[0]) as src:
        count = src.count

    # Merge all tiles
    print("Merging multiband DEM tiles...")
    src_files_to_mosaic = [rasterio.open(fp) for fp in dem_files]
    mosaic, out_trans = merge(src_files_to_mosaic, indexes=list(range(1, count + 1)))

    # Update metadata for the final mosaic
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": target_crs,
        "count": count
    })

    output_path = os.path.join(input_dir, "final_mosaicked_dem.tif")
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(mosaic)

    print(f"Saved final multiband mosaic: {output_path}")

    for src in src_files_to_mosaic:
        src.close()

    return output_path

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
