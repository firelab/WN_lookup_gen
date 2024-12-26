import os
import subprocess

def create_sbatch_script(output_dir, simnum, windninja_cli_path, elevation_file, input_speed, input_direction):
    """
    Create an sbatch script for running WindNinja.
    """
    sbatch_content = f"""#!/bin/bash
#SBATCH --job-name=windninja_{simnum}
#SBATCH --output={output_dir}/windninja_{simnum}.out
#SBATCH --error={output_dir}/windninja_{simnum}.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

{windninja_cli_path} --elevation_file {elevation_file} \\
           --initialization_method domainAverageInitialization \\
           --output_wind_height 10 \\
           --units_output_wind_height m --input_speed {input_speed} \\
           --input_speed_units mph --input_direction {input_direction} \\
           --input_wind_height 10 --units_input_wind_height m --vegetation grass \\
           --mesh_resolution 240 --units_mesh_resolution m --write_ascii_output 1 \\
           --num_threads 1
"""
    sbatch_path = os.path.join(output_dir, f"windninja_{simnum}.sbatch")
    with open(sbatch_path, "w") as f:
        f.write(sbatch_content)
    return sbatch_path


def submit_sbatch_script(sbatch_path):
    """
    Submit the sbatch script to Slurm.
    """
    try:
        subprocess.run(["sbatch", sbatch_path], check=True)
        print(f"Successfully submitted job with sbatch script: {sbatch_path}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to submit sbatch script: {sbatch_path}\nError: {e}")
