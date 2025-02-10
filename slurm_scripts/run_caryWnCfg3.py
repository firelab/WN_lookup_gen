import os
import shutil
import numpy as np
import time

pathToCurrentWindNinjaBuild = ""
dems_folder = "/output/dems_folder"
base_cli_file = "/test/base_cli.cfg"
wind_dirs = np.linspace(0.0, 337.5, 16)
nWindDirs = len(wind_dirs)

def main():
    if not os.path.exists(dems_folder):
        raise RuntimeError(f"!!! input dems_folder \"{dems_folder}\" does not exist !!!")
    if not os.path.isfile(base_cli_file):
        raise RuntimeError(f"!!! input base_cli_file \"{base_cli_file}\" does not exist !!!")
    
    dem_foldername = os.path.join(dems_folder, "dem0")
    dem_filename = f"{dem_foldername}/dem0.tif"
    
    if not os.path.exists(dem_foldername):
        raise RuntimeError(f"!!! input dem folder \"{dem_foldername}\" does not exist !!!")
    
    simTypes = ["mass", "momentum"]
    
    for simType in simTypes:
        current_simType_folder = f"{dem_foldername}/{simType}"
        if os.path.exists(current_simType_folder):
            shutil.rmtree(current_simType_folder)
        os.mkdir(current_simType_folder)
        
        for current_wind_dir in wind_dirs:
            current_wind_directory = f"{current_simType_folder}/{str(current_wind_dir).replace('.', '-')}-deg"
            os.mkdir(current_wind_directory)
            current_output_directory = current_wind_directory
            
            current_cfg_filename = f"{current_output_directory}/cli.cfg"
            shutil.copyfile(base_cli_file, current_cfg_filename)
            
            with open(current_cfg_filename, 'r') as file:
                filedata = file.read()
            
            filedata = filedata.replace("$dem_file", dem_filename)
            if simType == "momentum":
                filedata = filedata.replace("#momentum_flag              = true", "momentum_flag              = true")
            filedata = filedata.replace("$wind_dir", str(current_wind_dir))
            filedata = filedata.replace("$output_directory", current_output_directory)
            
            with open(current_cfg_filename, 'w') as file:
                file.write(filedata)
            
            wnBuildPath = pathToCurrentWindNinjaBuild
            if wnBuildPath:
                wnBuildPath += "/src/cli/"
            current_cmd = f"{wnBuildPath}WindNinja_cli {current_cfg_filename}"
            
            print(f"\nrunning \"{simType}\" simulation for dem \"dem0\" windDir \"{current_wind_dir} deg\"", flush=True)
            startTime = time.perf_counter()
            os.system(current_cmd)
            endTime = time.perf_counter()
            print(f"python os.system() timer: {endTime - startTime} seconds", flush=True)
            
            if simType == "momentum":
                foam_dirs = [d for d in os.listdir(dem_foldername) if "NINJAFOAM" in d]
                if foam_dirs:
                    shutil.rmtree(os.path.join(dem_foldername, foam_dirs[0]))
                    
if __name__ == '__main__':
    main()
