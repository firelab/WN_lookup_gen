import os
import rasterio
from rasterio.windows import Window

def split_raster(input_file, output_dir, tile_size=None, tile_size_meters=None, overlap=0):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    with rasterio.open(input_file) as src:
        # Get raster dimensions and transform
        width, height = src.width, src.height
        transform = src.transform
        pixel_size_x = transform[0]  # Pixel width in x-direction (meters per pixel)
        pixel_size_y = -transform[4]  # Pixel height in y-direction (meters per pixel)

        # Calculate tile size in pixels if tile_size_meters is specified
        if tile_size_meters is not None:
            tile_size_x = int(tile_size_meters / pixel_size_x)
            tile_size_y = int(tile_size_meters / pixel_size_y)
        elif tile_size is not None:
            tile_size_x = tile_size
            tile_size_y = tile_size
        else:
            raise ValueError("Either tile_size or tile_size_meters must be specified.")

        # Loop over rows and columns to create tiles
        for i in range(0, height, tile_size_y - overlap):
            for j in range(0, width, tile_size_x - overlap):
                # Define window
                window = Window(j, i, tile_size_x, tile_size_y)

                # Adjust window if it exceeds raster dimensions
                if window.col_off + window.width > width:
                    window = window.intersect(Window(j, i, width - j, tile_size_y))
                if window.row_off + window.height > height:
                    window = window.intersect(Window(j, i, tile_size_x, height - i))

                # Read window data
                tile_data = src.read(window=window)

                # Update transform for the tile
                tile_transform = src.window_transform(window)

                # Update metadata
                tile_profile = src.profile.copy()
                tile_profile.update({
                    "height": window.height,
                    "width": window.width,
                    "transform": tile_transform
                })

                # Save the tile
                tile_filename = os.path.join(output_dir, f"tile_{i}_{j}.tif")
                with rasterio.open(tile_filename, "w", **tile_profile) as dst:
                    dst.write(tile_data)
