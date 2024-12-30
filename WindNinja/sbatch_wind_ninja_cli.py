import os
import subprocess

def create_sbatch_script(output_dir, simnum, windninja_sif_path, elevation_file, input_speed, input_direction):
    """
    Create an sbatch script for running WindNinja with a Singularity .sif file.
    """
    sbatch_content = f"""#!/bin/bash
#SBATCH --job-name=WindNinja_Pipeline
#SBATCH --output=WindNinja_Pipeline_%j.log
#SBATCH --error=WindNinja_Pipeline_%j.err
#SBATCH --nodes=22
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

# Variables
LOCAL_DIR="/data/windninja"
WINDNINJA_SIF="{windninja_sif_path}"
ELEVATION_FILE="{elevation_file}"
INPUT_SPEED="{input_speed}"
INPUT_DIRECTION="{input_direction}"

# Remove and recreate local directory
srun --nodes=$SLURM_NNODES --ntasks=$SLURM_NNODES --exclusive bash -c "
    rm -rf $LOCAL_DIR &&
    mkdir -p $LOCAL_DIR &&
    rsync -a $WINDNINJA_SIF $LOCAL_DIR/
"
echo ".sif file copied to all compute nodes."

# Execute WindNinja
echo "Starting WindNinja simulation..."
srun --nodes=1 --ntasks=1 --exclusive bash -c "
    singularity exec -B /data:/data $LOCAL_DIR/$(basename $WINDNINJA_SIF) windninja_cli \\
        --elevation_file $ELEVATION_FILE \\
        --initialization_method domainAverageInitialization \\
        --output_wind_height 10 \\
        --units_output_wind_height m --input_speed $INPUT_SPEED \\
        --input_speed_units mph --input_direction $INPUT_DIRECTION \\
        --input_wind_height 10 --units_input_wind_height m --vegetation grass \\
        --mesh_resolution 240 --units_mesh_resolution m --write_ascii_output 1 \\
        --num_threads 32
"
echo "WindNinja simulation completed."
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
