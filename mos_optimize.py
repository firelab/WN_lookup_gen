"""
 * Project:  WN_Lookup_Gen
 * Purpose:  Python script for optimizing size of CONUS mosaics as a last step
 * Author: Gunjan Dayani <gunjan.dayani@usda.gov>

 ******************************************************************************
 *
 * THIS SOFTWARE WAS DEVELOPED AT THE ROCKY MOUNTAIN RESEARCH STATION (RMRS)
 * MISSOULA FIRE SCIENCES LABORATORY BY EMPLOYEES OF THE FEDERAL GOVERNMENT 
 * IN THE COURSE OF THEIR OFFICIAL DUTIES. PURSUANT TO TITLE 17 SECTION 105 
 * OF THE UNITED STATES CODE, THIS SOFTWARE IS NOT SUBJECT TO COPYRIGHT 
 * PROTECTION AND IS IN THE PUBLIC DOMAIN. RMRS MISSOULA FIRE SCIENCES 
 * LABORATORY ASSUMES NO RESPONSIBILITY WHATSOEVER FOR ITS USE BY OTHER 
 * PARTIES,  AND MAKES NO GUARANTEES, EXPRESSED OR IMPLIED, ABOUT ITS QUALITY, 
 * RELIABILITY, OR ANY OTHER CHARACTERISTIC.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************
"""

import os
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from glob import glob
import argparse
import fiona
from rasterio.features import geometry_mask

def optimize_conus_mosaics(input_dir, mask_path, output_dir):
    """
    Optimizes CONUS mosaics by masking and compressing rasters.

    Arguments:
    - input_dir: Directory containing input mosaics (.tif files).
    - mask_path: Path to the CONUS mask raster (should have 1 inside CONUS, 0 outside).
    - output_dir: Directory where optimized rasters will be saved.
    """
    os.makedirs(output_dir, exist_ok=True)
    with fiona.open(mask_path, "r") as shapefile:
        conus_shapes = [feature["geometry"] for feature in shapefile]

    # Loop through input mosaics
    for tif_path in glob(os.path.join(input_dir, "*.tif")):
        with rasterio.open(tif_path) as src:
            print(f"Processing: {os.path.basename(tif_path)}")

            band1 = src.read(1)
            band2 = src.read(2)

            # Create geometry mask for current raster
            aligned_mask = geometry_mask(
                geometries=conus_shapes,
                transform=src.transform,
                invert=True,
                out_shape=(src.height, src.width)
            ).astype(np.uint8)

            # Mask out values outside CONUS
            band1[aligned_mask != 1] = -9999
            band2[aligned_mask != 1] = -9999

            # Cast to int16
            band1 = band1.astype(np.float32)
            band2 = band2.astype(np.float32)

            # Update profile
            profile = src.profile.copy()
            profile.update({
                "count": 2,
                "dtype": "float32",
                "compress": "deflate",
                "tiled": True,
                "blockxsize": 256,
                "blockysize": 256,
                "BIGTIFF": "IF_SAFER",
                "nodata": -9999
            })

            out_path = os.path.join(output_dir, os.path.basename(tif_path))
            with rasterio.open(out_path, "w", **profile) as dst:
                dst.write(band1, 1)
                dst.write(band2, 2)

    print("All rasters masked to CONUS and cleaned (values outside set to 0).")

def main():
    parser = argparse.ArgumentParser(description="Optimize CONUS mosaics by masking and compressing rasters.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input mosaics (.tif files)")
    parser.add_argument("--mask_path", required=True, help="Path to the CONUS mask raster file")
    parser.add_argument("--output_dir", required=True, help="Directory to save processed output rasters")
    args = parser.parse_args()

    optimize_conus_mosaics(args.input_dir, args.mask_path, args.output_dir)

if __name__ == "__main__":
    main()
