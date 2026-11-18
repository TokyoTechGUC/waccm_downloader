#!/bin/bash
#$ -N waccm_download
#$ -cwd
#$ -l cpu_4=1
#$ -l h_rt=05:00:00
#$ -p -5
#$ -o logs/waccm_download_$JOB_ID.out
#$ -e logs/waccm_download_$JOB_ID.err

# ==============================================================================
# TSUBAME 4.0 SGE Job Script: WACCM Subsetting & Merging Downloader
# ==============================================================================

set -euo pipefail
mkdir -p logs

# Load system environment / Guconda Python
export PATH="/gs/bs/tga-guc-lab/dependencies/dependencies_intel/conda/envs/guconda/bin:$PATH"

# Determine script directory (handles SGE execution via SGE_O_WORKDIR)
SCRIPT_DIR="${SGE_O_WORKDIR:-/gs/bs/tga-guc-lab/users/alvin/wrf_ai/tools/waccm_downloader}"
PYTHON_SCRIPT="${SCRIPT_DIR}/download_waccm.py"

if [ ! -f "$PYTHON_SCRIPT" ]; then
    if [ -f "/gs/bs/tga-guc-lab/users/alvin/guc_mods/waccm_downloader/download_waccm.py" ]; then
        PYTHON_SCRIPT="/gs/bs/tga-guc-lab/users/alvin/guc_mods/waccm_downloader/download_waccm.py"
    else
        PYTHON_SCRIPT="/gs/bs/tga-guc-lab/users/alvin/wrf_ai/tools/waccm_downloader/download_waccm.py"
    fi
fi

# Input environment variables
NAMELIST_WPS="${NAMELIST_WPS:-./namelist.wps}"
GEO_EM="${GEO_EM:-}"
START_DATE="${START_DATE:-}"
END_DATE="${END_DATE:-}"
BUFFER="${BUFFER:-1.0}"
MOZBC_INP="${MOZBC_INP:-}"
VARS="${VARS:-all}"
KEEP_DAILY="${KEEP_DAILY:-true}"
OUTPUT_FILE="${OUTPUT_FILE:-}"
OUTPUT_DIR="${OUTPUT_DIR:-/gs/bs/tga-guc-lab/compilation/wrf-chem-boundaries/mozbc}"

echo "============================================================"
echo "TSUBAME 4.0 Job ID : ${JOB_ID:-local}"
echo "Running Node       : $(hostname)"
echo "Working Dir        : ${SGE_O_WORKDIR:-$(pwd)}"
echo "Python Location    : $(which python3)"
echo "Namelist WPS Path  : $NAMELIST_WPS"
if [ -n "$GEO_EM" ]; then echo "Geo_em NetCDF Path : $GEO_EM"; fi
if [ -n "$START_DATE" ]; then echo "Start Date        : $START_DATE"; fi
if [ -n "$END_DATE" ]; then echo "End Date          : $END_DATE"; fi
echo "Buffer             : +/- $BUFFER deg"
if [ -n "$MOZBC_INP" ]; then echo "MOZBC Input Config : $MOZBC_INP"; fi
echo "Variables (--vars) : $VARS"
echo "Output Directory   : $OUTPUT_DIR"
if [ -n "$OUTPUT_FILE" ]; then echo "Output File       : $OUTPUT_FILE"; fi
echo "============================================================"

# Build python command arguments
CMD=(python3 "$PYTHON_SCRIPT" --namelist "$NAMELIST_WPS" --buffer "$BUFFER" --output-dir "$OUTPUT_DIR" --vars "$VARS")

if [ -n "$MOZBC_INP" ]; then
    CMD+=(--mozbc-inp "$MOZBC_INP")
fi

if [ "$KEEP_DAILY" = "true" ] || [ "$KEEP_DAILY" = "1" ]; then
    CMD+=(--keep-daily)
fi

if [ -n "$GEO_EM" ]; then
    CMD+=(--geo-em "$GEO_EM")
fi

if [ -n "$START_DATE" ]; then
    CMD+=(--start-date "$START_DATE")
fi

if [ -n "$END_DATE" ]; then
    CMD+=(--end-date "$END_DATE")
fi

if [ -n "$OUTPUT_FILE" ]; then
    CMD+=(--output "$OUTPUT_FILE")
fi

# Execute downloader
"${CMD[@]}"

echo "Job finished successfully at $(date)"
