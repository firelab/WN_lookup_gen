import os
from raster_utils import process_dem, resample_tif, mosaic_dem_tiles, load_and_prepare_boundary
import rioxarray as rxr
from shapely.geometry import box

class DEMTiler:
    def __init__(self, dem_path, output_dir="out", boundary_file=None, resolution=240):
        self.dem_path = dem_path
        self.output_dir = output_dir
        self.raster_data = rxr.open_rasterio(dem_path, masked=True)
        self.dem_bounds = self.raster_data.rio.bounds()
        self.boundary_file = boundary_file
        self.resolution = resolution
        self.tile_counter = 0

    def process_tiles_by_hybrid(self, tile_size, buffer_size):
        """Hybrid tiling that considers both size and boundary alignment."""
        print(f"Tiling DEM into hybrid tiles of {tile_size} with boundary alignment...")
        x_min, y_min, x_max, y_max = self.dem_bounds
        x_tiles = int((x_max - x_min) // tile_size)
        y_tiles = int((y_max - y_min) // tile_size)

        for i in range(x_tiles + 1):
            for j in range(y_tiles + 1):
                print(f"Processing tile ({i}, {j})...")
                tile_bbox = box(
                    x_min + i * tile_size,
                    y_min + j * tile_size,
                    min(x_min + (i + 1) * tile_size, x_max),
                    min(y_min + (j + 1) * tile_size, y_max)
                )

                if self.boundary_file:
                    print(f"Checking boundary intersections for tile ({i}, {j})...")
                    boundary_gdf = load_and_prepare_boundary(self.boundary_file, self.raster_data.rio.crs, self.dem_bounds)
                    intersecting_boundaries = boundary_gdf[boundary_gdf.intersects(tile_bbox)]

                    if intersecting_boundaries.empty:
                        print(f"No intersecting boundaries found for tile ({i}, {j}). Skipping boundary alignment.")
                    else:
                        for idx, row in intersecting_boundaries.iterrows():
                            print(f"Processing sub-tile ({i}, {j}) for boundary {idx}...")
                            sub_tile_bbox = row.geometry.intersection(tile_bbox)
                            process_dem(sub_tile_bbox, buffer_size, self.raster_data, self.output_dir, self.resolution, self.tile_counter)
                            self.tile_counter += 1
                else:
                    print(f"No boundary file provided. Processing tile ({i}, {j}) without boundary alignment.")
                    process_dem(tile_bbox, buffer_size, self.raster_data, self.output_dir, self.resolution, self.tile_counter)
                    self.tile_counter += 1

    def process_tiles_by_boundary(self, buffer_size):
        if not self.boundary_file:
            raise ValueError("Boundary file is not provided.")
        
        print("Processing tiles by boundary...")
        boundary_gdf = load_and_prepare_boundary(self.boundary_file, self.raster_data.rio.crs, self.dem_bounds)

        for idx, row in boundary_gdf.iterrows():
            print(f"Processing boundary tile {idx}...")
            process_dem(row.geometry, buffer_size, self.raster_data, self.output_dir, self.resolution, self.tile_counter)
            self.tile_counter += 1

def main():
    dem_file_path = 'WindNinjaData/2023_lcp.tif'
    boundary_file = 'WindNinjaData/ne_10m_admin_1_states_provinces_lakes.shp'

    # Initialize DEMTiler
    dem_tiler = DEMTiler(dem_file_path, boundary_file=boundary_file, resolution=30)

    # Process DEM tiles by hybrid approach (size and boundary alignment)
    dem_tiler.process_tiles_by_hybrid(tile_size=32000, buffer_size=16000)

    # Mosaic processed DEM tiles
    # mosaic_dem_tiles(dem_tiler.output_dir, "EPSG:32612", (30, 30))

if __name__ == "__main__":
    main()