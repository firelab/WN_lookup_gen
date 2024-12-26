import os
from pathlib import Path
import pytest
import matplotlib.pyplot as plt
import rasterio
from lab_raster_utils.pipeline_runner import run_pipeline

@pytest.fixture
def input_raster():
    """Provide the path to the test raster in the tests folder."""
    return Path(__file__).parent / "test.tif"

# Verify NoData consistency
def test_nodata_handling(input_raster):
    with rasterio.open(input_raster) as src:
        input_nodata = src.nodata

    output_raster = "tests/output/merged_raster.tif"
    with rasterio.open(output_raster) as dst:
        output_nodata = dst.nodata

    assert input_nodata == output_nodata, "NoData value mismatch between input and output rasters"
    print(f"NoData value ({input_nodata}) preserved successfully.")

def visualize_rasters_as_grid(raster_paths, output_path, title):
    """
    Visualize multiple rasters as a grid in a single PNG.
    """
    num_rasters = len(raster_paths)
    cols = 4  # Number of columns in the grid
    rows = (num_rasters + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(15, rows * 5))
    axes = axes.flatten()  # Flatten axes for easy iteration

    for ax, raster_path in zip(axes, raster_paths):
        with rasterio.open(raster_path) as src:
            data = src.read(1)  # Read the first band
            ax.imshow(data, cmap="gray")
            ax.set_title(raster_path.stem, fontsize=10)
            ax.axis("off")

    # Turn off remaining empty subplots
    for ax in axes[len(raster_paths):]:
        ax.axis("off")

    plt.suptitle(title, fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout for title
    plt.savefig(output_path)
    plt.close()

def test_pipeline_with_visualization(input_raster):
    """
    Test the updated pipeline: tiling with overlap (in km), smoothing, trimming, and merging.
    Includes visualizations for each operation.
    """
    assert input_raster.exists(), "Input raster does not exist"

    # Output directory setup
    output_dir = Path("tests/output")
    if output_dir.exists():
        import shutil
        shutil.rmtree(output_dir)  # Clean up previous runs
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run the pipeline with tile size and overlap in kilometers
    tile_size_km = 30  # Tile size in kilometers
    overlap_km = 2     # Overlap in kilometers
    smoothing_sigma = 2  # Smoothing parameter
    run_pipeline(str(input_raster), str(output_dir), tile_size_km=tile_size_km, overlap_km=overlap_km, sigma=smoothing_sigma)

    # Visualization directory
    visualization_dir = output_dir / "visualizations"
    visualization_dir.mkdir(exist_ok=True)

    # Visualize tiles
    tiles_dir = output_dir / "tiles"
    assert tiles_dir.exists() and any(tiles_dir.glob("*.tif")), "Tiling failed"
    tile_png_path = visualization_dir / "tiles_grid.png"
    visualize_rasters_as_grid(list(tiles_dir.glob("*.tif")), tile_png_path, "Tiled Rasters with Overlap (km)")

    # Visualize smoothed tiles
    smoothed_dir = output_dir / "smoothed"
    assert smoothed_dir.exists() and any(smoothed_dir.glob("*.tif")), "Smoothing failed"
    smoothed_png_path = visualization_dir / "smoothed_grid.png"
    visualize_rasters_as_grid(list(smoothed_dir.glob("*.tif")), smoothed_png_path, "Smoothed Tiles")

    # Visualize trimmed tiles
    trimmed_dir = output_dir / "trimmed"
    assert trimmed_dir.exists() and any(trimmed_dir.glob("*.tif")), "Trimming failed"
    trimmed_png_path = visualization_dir / "trimmed_grid.png"
    visualize_rasters_as_grid(list(trimmed_dir.glob("*.tif")), trimmed_png_path, "Trimmed Tiles")

    # Visualize merged raster
    merged_output = output_dir / "merged_raster.tif"
    assert merged_output.exists(), "Merging failed"
    merged_png_path = visualization_dir / "merged.png"
    visualize_rasters_as_grid([merged_output], merged_png_path, "Merged Raster")

    # Verify input and output dimensions
    with rasterio.open(input_raster) as src:
        input_width, input_height = src.width, src.height
    with rasterio.open(merged_output) as dst:
        output_width, output_height = dst.width, dst.height

    assert (input_width, input_height) == (output_width, output_height), \
        "Merged raster dimensions do not match input raster"

    print(f"Visualizations saved in {visualization_dir}")
