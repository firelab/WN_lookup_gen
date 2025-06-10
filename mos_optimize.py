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
import xarray as xr
import argparse
import fiona
from rasterio.features import geometry_mask

def optimize_float32(input_dir, mask_path, output_dir):
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

            band1 = band1.astype(np.float32)
            band2 = band2.astype(np.float32)

            # Update profile
            profile = src.profile.copy()
            profile.update({
                "count": 2,
                "dtype": "float32",
                "compress": "deflate",
                "predictor": 2,
                "tiled": True,
                "interleave": "band",
                "blockxsize": 256,
                "blockysize": 256,
                "BIGTIFF": "IF_SAFER",
                "nodata": -9999
            })

            out_path = os.path.join(output_dir, os.path.basename(tif_path))
            with rasterio.open(out_path, "w", **profile) as dst:
                dst.write(band1, 1)
                dst.write(band2, 2)
                dst.set_band_description(1, "Wind Speed")
                dst.set_band_description(2, "Wind Direction")
                dst.update_tags(1, units="mph", long_name="Wind Speed")
                dst.update_tags(2, units="degrees", long_name="Wind Direction")

def optimize_int32(input_dir, mask_path, output_dir):
    """
    Masks rasters to CONUS extent using shapefile and converts them to int32 by scaling by 100.
    Saves the result with compression and tiling options.
    """
    os.makedirs(output_dir, exist_ok=True)

    with fiona.open(mask_path, "r") as shapefile:
        conus_shapes = [feature["geometry"] for feature in shapefile]

    for tif_path in glob(os.path.join(input_dir, "*.tif")):
        with rasterio.open(tif_path) as src:
            print(f"Processing and scaling to Int32: {os.path.basename(tif_path)}")

            band1 = src.read(1)
            band2 = src.read(2)

            aligned_mask = geometry_mask(
                geometries=conus_shapes,
                transform=src.transform,
                invert=True,
                out_shape=(src.height, src.width)
            ).astype(np.uint8)

            # Mask values outside CONUS to NaN before scaling
            band1_masked = np.where(aligned_mask == 1, band1, np.nan)
            band2_masked = np.where(aligned_mask == 1, band2, np.nan)

            # Ensure float32 before scaling
            band1_masked = band1_masked.astype(np.float32)
            band2_masked = band2_masked.astype(np.float32)

            # Scale only valid values
            band1_scaled = np.where(np.isnan(band1_masked), np.nan, band1_masked * 100)
            band2_scaled = np.where(np.isnan(band2_masked), np.nan, band2_masked * 100)

            # No need to clip direction range (int32 won't overflow)
            band1_final = np.where(np.isnan(band1_scaled), -9999, band1_scaled).astype(np.int32)
            band2_final = np.where(np.isnan(band2_scaled), -9999, band2_scaled).astype(np.int32)

            profile = src.profile.copy()
            profile.update({
                "count": 2,
                "dtype": "int32",
                "compress": "deflate",
                "predictor": 2,
                "tiled": True,
                "interleave": "band",
                "blockxsize": 256,
                "blockysize": 256,
                "BIGTIFF": "IF_SAFER",
                "nodata": -9999
            })

            out_path = os.path.join(output_dir, os.path.basename(tif_path).replace(".tif", "_int32.tif"))
            with rasterio.open(out_path, "w", **profile) as dst:
                dst.write(band1_final, 1)
                dst.write(band2_final, 2)
                dst.set_band_description(1, "Wind Speed x100")
                dst.set_band_description(2, "Wind Direction x100")
                dst.update_tags(1, units="mph x100", long_name="Wind Speed (scaled)")
                dst.update_tags(2, units="degrees x100", long_name="Wind Direction (scaled)")


def write_combined_netcdf(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    directions = [
        "0-0-deg", "22-5-deg", "45-0-deg", "67-5-deg",
        "90-0-deg", "112-5-deg", "135-0-deg", "157-5-deg",
        "180-0-deg", "202-5-deg", "225-0-deg", "247-5-deg",
        "270-0-deg", "292-5-deg", "315-0-deg", "337-5-deg"
    ]
    out_nc_path = os.path.join(output_dir, "combined_wind.nc")
    first = True

    for dir_idx, d in enumerate(directions):
        tif_path = os.path.join(input_dir, f"{d}.tif")
        with rasterio.open(tif_path) as src:
            # Read and scale
            speed = (src.read(1).astype(np.float32) * 100).astype(np.int16)
            direction = (src.read(2).astype(np.float32) * 100).astype(np.int16)

            if first:
                height, width = speed.shape
                y_coords = src.xy(0, 0)[1] - np.arange(height) * src.res[1]
                x_coords = src.xy(0, 0)[0] + np.arange(width) * src.res[0]
                x_origin = src.xy(0, 0)[0]
                y_origin = src.xy(0, 0)[1]
                x_res = src.res[0]
                y_res = src.res[1]

                # Preallocate int16 array
                ds = xr.Dataset(
                    {
                        "wind_data": (("direction", "band", "y", "x"),
                                      np.full((16, 2, height, width), -9999, dtype=np.int16))
                    },
                    coords={
                        "direction": directions,
                        "band": ["speed", "direction"],
                        "y": y_coords,
                        "x": x_coords
                    },
                    attrs={
                        "description": "Wind speed (mph) and direction (deg) mosaics (x100) from 16 directions",
                        "x_origin": x_origin,
                        "y_origin": y_origin,
                        "x_res": x_res,
                        "y_res": y_res,
                        "scale_factor": 0.01,
                        "units": "int16 (x100)"
                    }
                )
                first = False

            ds["wind_data"][dir_idx, 0, :, :] = speed
            ds["wind_data"][dir_idx, 1, :, :] = direction

    ds.to_netcdf(out_nc_path, format="NETCDF4", engine="netcdf4", encoding={
        "wind_data": {"zlib": True, "complevel": 4, "dtype": "int16"}
    })

    print(f"Saved NetCDF: {out_nc_path}")

def main():
    parser = argparse.ArgumentParser(description="Optimize CONUS mosaics by masking and compressing rasters.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input mosaics (.tif files)")
    parser.add_argument("--mask_path", required=True, help="Path to the CONUS mask raster file")
    parser.add_argument("--output_dir", required=True, help="Directory to save processed output rasters")
    args = parser.parse_args()

    optimize_int32(args.input_dir, args.mask_path, args.output_dir)
    #optimize_float32(args.input_dir, args.mask_path, args.output_dir)
    #write_combined_netcdf(args.input_dir, args.output_dir)

if __name__ == "__main__":
    main()
