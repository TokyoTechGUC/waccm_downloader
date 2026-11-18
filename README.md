# waccm_downloader

Automated NCAR WACCM / CAM-chem THREDDS NetCDF Subset Service (NCSS) downloader and MOZBC preprocessor for WRF-Chem on TSUBAME 4.0.

## Quick Start

```bash
git clone https://github.com/TokyoTechGUC/waccm_downloader.git
cd waccm_downloader

qsub -g tga-guc-lab \
     -v GEO_EM="/path/to/geo_em.d01.nc",START_DATE="2026-07-25",END_DATE="2026-09-01",MOZBC_INP="/path/to/mozbc_d01.inp",OUTPUT_DIR="$(pwd)",KEEP_DAILY="true" \
     submit_waccm_download.tsubame.bash
```

## Features
- Exact Domain 1 bounding-box extraction from `geo_em.d01.nc` or `namelist.wps` (`±1.0°` buffer).
- Dynamic variable extraction from any `mozbc` `.inp` (`spc_map`) via `--mozbc-inp`, validated against remote OPeNDAP `.dds`.
- Multi-tier fallback: single-shot NCSS (`temporal=all`) → Chunked NCSS → per-variable OPeNDAP (`dodsC`).
- Automatic injection of MOZBC-required 1D coordinates (`hyam`, `hybm`, `P0`, `date`, `datesec`).
