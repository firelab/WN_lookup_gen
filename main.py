import os
import xarray as xr
from WindNinja.dem_processing import process_dem
from WindNinja.spatial_operations import read_boundary_file, reproject_boundary
from WindNinja.sbatch_wind_ninja_cli import create_sbatch_script, submit_sbatch_script


def run_wind_ninja_narr_cell_parallel(narr_grid_df):
    """Run WindNinja for NARR grid cells in parallel using Slurm."""
    for idx, row in narr_grid_df.iterrows():
        simnum = str(idx).zfill(7)
        subpoly = row.geometry
        output_dir = os.path.join('/mnt/fsim/WindNinjaCONUS/', simnum)
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate the DEM file for the grid cell
        dem_path = process_dem(subpoly, 16000, xr.open_dataset('WindNinjaData/LC20_Elev_220_RS_240m.tif'), output_dir)
        
        # Create sbatch script for the job
        sbatch_path = create_sbatch_script(
            output_dir=output_dir,
            simnum=simnum,
            windninja_cli_path='usr/bin/WindNinja_cli',
            elevation_file=dem_path,
            input_speed=5,
            input_direction=0
        )
        
        # Submit the sbatch script to the cluster
        submit_sbatch_script(sbatch_path)


if __name__ == "__main__":
    bndryfile = 'WindNinjaData/ne_10m_admin_1_states_provinces_lakes.shp'
    bnd = read_boundary_file(bndryfile)
    narr_proj4 = "+proj=lcc +lat_1=50 +lat_0=50 +lon_0=-107 +k_0=1 +x_0=5632642.22547 +y_0=4612545.65137 +a=6371200 +b=6371200 +units=m +no_defs"
    narr_grid_df = reproject_boundary(bnd, narr_proj4)
    run_wind_ninja_narr_cell_parallel(narr_grid_df)
