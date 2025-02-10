#!/bin/bash
export CPL_DEBUG=NINJAFOAM
source /opt/openfoam8/etc/bashrc
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export FOAM_USER_LIBBIN=/usr/local/lib/
OUTPUT_FOLDER="/output"
LOG_FILE="${OUTPUT_FOLDER}/simulation.log"
python3 /test/run_caryWnCfg3.py > "${LOG_FILE}" 2>&1
