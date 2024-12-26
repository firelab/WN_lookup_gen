import matplotlib.pyplot as plt
import geopandas as gpd
import xarray as xr
from rioxarray.merge import merge_arrays


def visualize_tiles_and_buffers(narr_grid_df, dem_files):
    """Visualize the tiles and buffers on the map."""
    f, ax = plt.subplots()
    merged_raster = merge_arrays(dataarrays=dem_files, res=(240, 240), crs="EPSG:32612")
    merged_raster.rio.reproject("EPSG:5070").plot(
        levels=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 18], cmap='RdYlBu_r', ax=ax
    )
    for _, row in narr_grid_df.iterrows():
        gpd.GeoSeries(row.geometry).to_crs("EPSG:5070").plot(ax=ax, facecolor='none', edgecolor='black')
    ax.set_title('Visualization of Tiles and Buffers')
    plt.show()


def visualize_dem_and_buffers(clipped_dem, subpoly, subpolybuf):
    """Visualize the DEM in native projection with buffers."""
    f, ax = plt.subplots()
    clipped_dem.band_data[0].plot(ax=ax)
    subpoly.to_crs(clipped_dem.rio.crs).plot(ax=ax, facecolor='none', edgecolor='black')
    subpolybuf.to_crs(clipped_dem.rio.crs).plot(ax=ax, facecolor='none', edgecolor='red')
    ax.set_title("DEM in native projection with buffers")
    plt.show()


def visualize_utm_projection(clipped_dem_utm, subpoly_utm):
    """Visualize the DEM in UTM projection with buffered bounds."""
    f, ax = plt.subplots()
    clipped_dem_utm.band_data[0].plot(ax=ax)
    subpoly_utm.buffer(16000, join_style='mitre').plot(ax=ax, facecolor='none', edgecolor='black')
    ax.set_title("DEM in UTM projection with buffered bounds")
    plt.show()