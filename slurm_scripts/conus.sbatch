#!/bin/bash
#SBATCH --job-name=WindNinja_Pipeline
#SBATCH --output=windninja_pipeline_%j.log
#SBATCH --error=windninja_pipeline_%j.err
#SBATCH --nodes=18
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2500

# Variables
SHARED_STORAGE="/mnt/fsim/windninja/CONUS"
WINDNINJA_SIF="/mnt/fsim/windninja/wn_latest6.sif"
export RETRY_COUNT=3
echo "Starting WindNinja Pipeline Job..."

# Remove and recreate local directory on all nodes
echo "Cleaning up old directories and copying .sif to nodes..."
srun --exclusive bash -c "rsync -av --progress ${WINDNINJA_SIF} /data/"
srun --exclusive bash -c "rsync -av --progress /home/ohpc/CONUS/run_caryWnCfg3.py /test/"
srun --exclusive bash -c "rsync -av --progress /home/ohpc/CONUS/run.sh /test/"
srun --exclusive bash -c "rsync -av --progress /home/ohpc/CONUS/base_cli.cfg /test/"
echo ".sif file copied to all compute nodes."

# Collect all folders into a task queue
task_queue=()
# for folder in "${FOLDERS[@]}"; do
# for folder in 0 118  14   151  181  19   217  24   254  287  31  38  44  53  6   88 1    119  146  152  182  2    218  249  255  288  32  39  45  54  7   89 10   12   147  153  183  20   219  25   256  289  33  4   46  55  8   9 11   120  148  16   184  21   22   250  26   29   34  40  47  56  84  90 115  121  149  17   185  214  220  251  27   290  35  41  48  57  85  116  122  15   18   186  215  221  252  28   3    36  42  49  58  86 117  13   150  180  187  216  23   253  286  30   37  43  5   59  87; do
for folder in 118 135 136 137 147 148 149 150 169 178 181 182 183 184 185 212 214 215 216 217 218 242 248 249 250 265 266 268 269 270 271 272 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 52 83 95;do
    if [ -d "/mnt/fsim/windninja/CONUS/$folder" ]; then
        task_queue+=("$folder")
    fi
done

total_jobs=${#task_queue[@]}
echo "Total folders to process: $total_jobs"

# Process base-io directories using available nodes
running_jobs=0
job_index=0

while [ $job_index -lt $total_jobs ] || [ $running_jobs -gt 0 ]; do
    if [ $running_jobs -lt $((SLURM_NNODES * 4)) ] && [ $job_index -lt $total_jobs ]; then
        echo "Launching batch of 4 jobs..."
        for i in {0..3}; do
            if [ $job_index -lt $total_jobs ]; then
                folder="${task_queue[$job_index]}"
                echo "Assigning folder $folder to task $i"
                job_index=$((job_index + 1))
                srun --exclusive -N1 -n1 --cpus-per-task=10 bash -c "/mnt/fsim/windninja/CONUS/src/slurm_run_conus.sh ${folder}"&
                ((running_jobs+=1))
            fi
        done
    else
        wait -n
        ((running_jobs-=1))
    fi
done
wait
echo "All tasks completed."
