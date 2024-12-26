import rioxarray
import os
from shapely.geometry import mapping

def process_dem(subpoly, buffer_size, dem_file, output_dir):
    """
    Clip, buffer, and reproject DEM.
    
    Parameters:
        subpoly (shapely.geometry.Polygon): Polygon geometry representing the area of interest.
        buffer_size (float): Buffer size (in CRS units, e.g., meters or degrees).
        dem_file (xarray.Dataset): Input DEM file opened using rioxarray.
        output_dir (str): Directory where the processed DEM will be saved.
    
    Returns:
        str: Path to the processed DEM file.
    """
    # Apply buffer to the polygon
    buffered_subpoly = subpoly.buffer(buffer_size, join_style=2)  # Join_style=2 for sharp corners

    # Convert the buffered polygon to GeoJSON format for clipping
    buffered_geojson = mapping(buffered_subpoly)

    # Clip the DEM using the buffered polygon
    clipped_dem = dem_file.rio.clip([buffered_geojson], crs="EPSG:4326", drop=True)

    # Reproject the clipped DEM to UTM (EPSG:32612)
    reprojected_dem = clipped_dem.rio.reproject("EPSG:32612")

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Save the processed DEM to the output directory
    output_path = os.path.join(output_dir, "processed_dem.tif")
    reprojected_dem.rio.to_raster(output_path)

    return output_path
