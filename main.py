import os
from WindNinja.dem_processing import mosaic_dem_tiles, process_dem, resample_tif
import rioxarray as rxr
from shapely.geometry import box
import geopandas as gpd

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

class DEMTiler:
    def __init__(self, dem_path, output_dir="out", boundary_file=None):
        self.dem_path = dem_path
        self.output_dir = output_dir
        self.raster_data = rxr.open_rasterio(dem_path, masked=True)
        self.dem_bounds = self.raster_data.rio.bounds()
        self.boundary_file = boundary_file

    def process_tiles_by_fixed_size(self, tile_size, buffer_size):
        """Tile DEM based on a fixed size (e.g., 500km x 500km)."""
        print(f"Tiling DEM into fixed size tiles of {tile_size}...")
        x_min, y_min, x_max, y_max = self.dem_bounds
        x_tiles = int((x_max - x_min) // tile_size)
        y_tiles = int((y_max - y_min) // tile_size)

        for i in range(x_tiles + 1):
            for j in range(y_tiles + 1):
                tile_bbox = box(
                    x_min + i * tile_size,
                    y_min + j * tile_size,
                    min(x_min + (i + 1) * tile_size, x_max),
                    min(y_min + (j + 1) * tile_size, y_max)
                )
                tile_output_dir = os.path.join(self.output_dir, f"tile_{i}_{j}")
                process_dem(tile_bbox, buffer_size, self.raster_data, tile_output_dir)

    def process_tiles_by_boundary(self, buffer_size):
        if not self.boundary_file:
            raise ValueError("Boundary file is not provided.")
        
        boundary_gdf = gpd.read_file(self.boundary_file, engine='pyogrio')
        print("Original boundary file CRS:", boundary_gdf.crs)

        dem_crs = self.raster_data.rio.crs
        boundary_gdf = boundary_gdf.to_crs(dem_crs)

        dem_bbox = box(*self.dem_bounds)
        boundary_gdf = boundary_gdf[boundary_gdf.intersects(dem_bbox)]

        print(f"Filtered boundaries intersecting DEM: {len(boundary_gdf)}")

        for idx, row in boundary_gdf.iterrows():
            tile_id = str(idx).zfill(7)
            subpoly = row.geometry
            tile_output_dir = os.path.join(self.output_dir, tile_id)
            process_dem(subpoly, buffer_size, self.raster_data, tile_output_dir)

def main():
    dem_file_path = 'WindNinjaData/LC20_Elev_220_RS_240m.tif'
    boundary_file = 'WindNinjaData/ne_10m_admin_1_states_provinces_lakes.shp'

    # Initialize DEMTiler
    dem_tiler = DEMTiler(dem_file_path, boundary_file=boundary_file)

    # Process DEM tiles by fixed size (e.g., 500km x 500km)
    dem_tiler.process_tiles_by_fixed_size(tile_size=500000, buffer_size=100000)

    # Process DEM tiles by boundary
    #dem_tiler.process_tiles_by_boundary(buffer_size=16000)

    # Mosaic processed DEM tiles
    mosaic_dem_tiles(dem_tiler.output_dir, "EPSG:32612", (240, 240))

if __name__ == "__main__":
    main()