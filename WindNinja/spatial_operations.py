import geopandas as gpd


def read_boundary_file(filepath):
    """Read and filter boundary shapefile."""
    bnd = gpd.read_file(filepath, engine='pyogrio')
    bnd = bnd[bnd['gu_a3'] == 'USA']
    bnd = bnd[(bnd['name'] != "Alaska") & (bnd['name'] != "Hawaii")]
    return bnd


def reproject_boundary(boundary, proj4):
    """Reproject the boundary file."""
    return boundary.to_crs(proj4)