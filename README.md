# WindNinja Preprocessing and Postprocessing Pipeline

This pipeline prepares DEMs, runs WindNinja simulations, and postprocesses the output, including smoothing and mosaicking.

---
## **Dependency requirements for python scripts**
GDAL
Multiprocessing
Subprocess
Numpy

## **Step 1: Prepare DEMs for WindNinja**
This step generates overlapped DEM tiles for WindNinja input while preserving the input raster projection.

### **Modify Paths:**
Edit `prepare_dems.py` and set:
```python
input_raster_path = <path for landscape file>
output_directory = <path for output folder>
```

### **You will get below file structure**
```
output_directory/
├── original_tiles/
│   ├── tile_0_x_y.tif
│   ├── tile_1_x_y.tif
│   └── ...
└── final_tiles/
    ├── 0/
    │   ├── dems_folder/
    │   │   ├── dem0/
    │   │   │   ├── dem0.tif
    │   │   │   └── ...
```


## **Step 2: Run WindNinja Simulations on Slurm**
The `conus.sbatch` script runs WindNinja simulations on an HPC cluster using Slurm.

### **Modify Paths**
Before running the WindNinja simulations, ensure that the paths in `conus.sbatch` are correctly set according to your environment:
- **SHARED_STORAGE**: Path where WindNinja output folder final_tiles is stored
- **WINDNINJA_SIF**: Path to the Singularity image for WindNinja.
- **Scripts and Configs**: Ensure `run_caryWnCfg3.py`, `run.sh`, and `base_cli.cfg` exist in the specified directory.

```bash
SHARED_STORAGE= <path>
WINDNINJA_SIF= <path>
```
```bash
sbatch conus.sbatch
```

### **To monitor logs:**
```bash
tail -f windninja_pipeline_<job_id>.log
```

## **Step 3: Postprocessing (Smoothing and Mosaicking)**

### **Modify Paths**
Before running this step, update the paths in `postprocess.py`:
- `base_dir`: Path to the processed WindNinja tiles.
- `output_dir`: Directory where processed outputs will be saved.

Modify these in `postprocess.py`:

```python
base_dir = <path>
output_dir = <path>
```

### **Run:**
```python
python3 postprocess.py
```

### **You will get below file structure**
```
output_directory/
├── windNinjax_tiles_momentum/
│   ├── 0/
│   │   ├── dem0_0_5_120m_reproj.tif
│   │   ├── dem0_0_5_120m_clipped.tif
│   │   ├── dem0_0_5_120m.tif
│   │   └── ...
│   ├── 1/
│   │   ├── dem0_0_5_120m_reproj.tif
│   │   ├── dem0_0_5_120m_clipped.tif
│   │   ├── dem0_0_5_120m.tif
│   │   └── ...
│   ├── final_mosaic_w_smoothing_momentum.tif
```
