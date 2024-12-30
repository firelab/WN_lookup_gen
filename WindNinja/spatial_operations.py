import geopandas as gpd
from shapely.geometry import box


def read_boundary_file(filepath):
    """Read and filter the boundary shapefile to include only US geometries."""
    bnd = gpd.read_file(filepath, engine='pyogrio')
    print("Original boundary file CRS:", bnd.crs)

    # Filter to include only the US (excluding Alaska and Hawaii)
    bnd = bnd[(bnd['gu_a3'] == 'USA') & 
              (bnd['name'] != "Alaska") & 
              (bnd['name'] != "Hawaii")]
    print(f"Filtered {len(bnd)} polygons for the contiguous US.")
    return bnd

def reproject_boundary(boundary, target_crs):
    """Reproject the boundary file."""
    return boundary.to_crs(target_crs)

def filter_boundaries(boundary, dem_bounds):
    """Filter boundaries to include only those intersecting the DEM."""
    dem_bbox = box(*dem_bounds)
    return boundary[boundary.intersects(dem_bbox)]