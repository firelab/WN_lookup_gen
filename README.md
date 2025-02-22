# WindNinja Preprocessing and Postprocessing Pipeline

This pipeline prepares DEMs, runs WindNinja simulations, and postprocesses the output, including smoothing and mosaicking.

---
## **Dependency requirements for python scripts**
GDAL
Multiprocessing
Subprocess
Numpy

## **Step 1: Prepare DEMs for WindNinja**
This step generates UTM-projected tiles with overlap for WindNinja input.

## **Strategy to geenrate square tile in UTM projection to get rid of NoData**
<img src="https://github.com/user-attachments/assets/a61b526b-cffa-4156-a64e-70f14d6673ca" width="600">

(Note: we are cliping out above square tile from 350% overlapped tile. I am generating 350% overlapped tile as I won't have whole input landscape file of CONUS(237GB) into my memory to get data after. I am saving it as intermediatory output and deleting it after getting square tile.)

### **Modify Paths:**
Edit `prepare_dems.py` and set:
```python
input_raster_path = <path for landscape file>
output_directory = <path for output folder>
```

### **Run:**
```python
python3 prepare_dems.py
```

### **You will get below file structure**
```
output_directory/
├── original_tiles/
│   ├── tile_0_x_y.tif
│   ├── tile_1_x_y.tif
│   └── ...
├── temp_albers_tiles/
│   ├── temp_350_albers_0.tif
│   ├── temp_overlap_albers_0.tif
│   ├── temp_utm_overlap_aoi_0_epsgXXXX.tif
│   ├── temp_utm_350_tile_0_epsgXXXX.tif
│   └── ...
├── utm_aoi_tiles/
│   ├── utm_overlap_aoi_0_epsgXXXX.tif
│   ├── utm_overlap_aoi_1_epsgXXXX.tif
│   └── ...
└── final_utm_tiles/
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
- **SHARED_STORAGE**: Path where WindNinja output folder final_utm_tiles is stored
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
- `aoi_tiles_dir`: Directory containing AOI tiles. (which was generated by prepare_dems scripts(Step 1))
- `output_dir`: Directory where processed outputs will be saved.

Modify these in `postprocess.py`:

```python
base_dir = <path>
aoi_tiles_dir = <path>
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
