from pathlib import Path
from lab_raster_utils.raster_pipeline import (
    tile_raster_with_overlap,
    process_tile_smoothing,
    trim_overlap,
    merge_tiles,
)
import rasterio

def run_pipeline(input_raster, output_dir, tile_size_km=10, overlap_km=2, sigma=2):
    """
    Full pipeline to process a raster: tile -> smooth -> trim -> merge.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Tile raster with overlap
    tiles_dir = output_dir / "tiles"
    tiles_dir.mkdir(exist_ok=True)
    tile_raster_with_overlap(input_raster, tile_size_km, overlap_km, str(tiles_dir))

    # Step 2: Smooth each tile
    smoothed_dir = output_dir / "smoothed"
    smoothed_dir.mkdir(exist_ok=True)
    for tile in tiles_dir.glob("*.tif"):
        smoothed_tile = smoothed_dir / f"smoothed_{tile.name}"
        process_tile_smoothing(str(tile), str(smoothed_tile), sigma)

    # Step 3: Trim overlap
    trimmed_dir = output_dir / "trimmed"
    trimmed_dir.mkdir(exist_ok=True)

    with rasterio.open(input_raster) as src:
        resolution = abs(src.transform[0])  # Pixel size in meters
        width, height = src.width, src.height

    for tile in smoothed_dir.glob("*.tif"):
        # Extract coordinates from the filename
        tile_name = tile.stem  # Example: "smoothed_tile_0_0"
        try:
            _, _, x, y = tile_name.split("_")  # Updated to handle the prefix
            coords = (int(x), int(y))
        except ValueError:
            raise ValueError(f"Invalid tile filename format: {tile_name}")

        # Determine if this tile is an edge tile
        is_edge_tile = (
            coords[0] == 0 or coords[1] == 0 or  # Left or top edge
            coords[0] + tile_size_km * 1000 / resolution >= width or  # Right edge
            coords[1] + tile_size_km * 1000 / resolution >= height  # Bottom edge
        )
        trimmed_tile = trimmed_dir / f"trimmed_{tile.name}"
        trim_overlap(str(tile), str(trimmed_tile), overlap_km, resolution, is_edge_tile)

    # Step 4: Merge tiles
    merged_output = output_dir / "merged_raster.tif"
    merge_tiles(str(trimmed_dir), str(merged_output))

    print(f"Pipeline completed. Final merged raster: {merged_output}")
