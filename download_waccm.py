#!/usr/bin/env python3
"""
WACCM / CAM-chem Automated Subsetting & Merging Downloader for WRF-Chem
=======================================================================

Automated tool to download spatially subsetted WACCM global atmospheric chemistry 
datasets from UCAR GDEX / THREDDS server across a specified date range, and 
merge them into a single NetCDF file (`waccm-YYYY-MM-DD-YYYY-MM-DD.nc`) for MOZBC.

Default Output Directory:
    /gs/bs/tga-guc-lab/compilation/wrf-chem-boundaries/mozbc/

Usage:
    python3 download_waccm.py --start-date 2022-01-25 --end-date 2022-01-30 --namelist ./namelist.wps

    # Or let it read dates and domain bounds automatically from namelist.wps:
    python3 download_waccm.py --namelist ./namelist.wps

Author: GUC Lab / WRF-Chem TSUBAME Workflow
"""

import os
import sys
import math
import re
import argparse
import datetime
import subprocess
from urllib.parse import urlparse, urlencode

DEFAULT_MOZBC_DIR = "/gs/bs/tga-guc-lab/compilation/wrf-chem-boundaries/mozbc"


def parse_date_string(date_str):
    """
    Parse date string like '2022-01-25_00:00:00', '2022-01-25', '20220125' into datetime.date object.
    """
    if not date_str:
        return None
    clean = str(date_str).strip().strip("'").strip('"')
    if '_' in clean:
        clean = clean.split('_')[0]

    for fmt in ['%Y-%m-%d', '%Y/%m/%d', '%Y%m%d']:
        try:
            return datetime.datetime.strptime(clean, fmt).date()
        except ValueError:
            pass
    return None


def parse_namelist_wps(namelist_path):
    """
    Parse namelist.wps to extract domain parameters and simulation start/end dates.
    """
    if not os.path.exists(namelist_path):
        raise FileNotFoundError(f"Namelist file not found: {namelist_path}")

    params = {}
    current_section = None

    with open(namelist_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            clean_line = line.split('!')[0].split('#')[0].strip()
            if not clean_line:
                continue

            if clean_line.startswith('&'):
                current_section = clean_line[1:].strip().lower()
                continue
            elif clean_line == '/':
                current_section = None
                continue

            if current_section in ['share', 'geogrid']:
                if '=' in clean_line:
                    key, val = clean_line.split('=', 1)
                    key = key.strip().lower()
                    val = val.strip().rstrip(',')
                    params[key] = val

    def get_val(key, default=None, val_type=float):
        if key not in params:
            return default
        raw = params[key]
        tokens = [t.strip().strip("'").strip('"') for t in re.split(r'[\s,]+', raw) if t.strip()]
        if not tokens:
            return default
        try:
            return val_type(tokens[0])
        except (ValueError, TypeError):
            return default

    def get_str(key, default=""):
        if key not in params:
            return default
        raw = params[key].strip().strip("'").strip('"')
        tokens = [t.strip().strip("'").strip('"') for t in re.split(r'[\s,]+', raw) if t.strip()]
        return tokens[0] if tokens else default

    parsed = {
        'max_dom': get_val('max_dom', 1, int),
        'start_date': get_str('start_date', ''),
        'end_date': get_str('end_date', ''),
        'e_we': get_val('e_we', 100, int),
        'e_sn': get_val('e_sn', 100, int),
        'dx': get_val('dx', 9000.0, float),
        'dy': get_val('dy', 9000.0, float),
        'map_proj': get_str('map_proj', 'lambert').lower(),
        'ref_lat': get_val('ref_lat', 0.0, float),
        'ref_lon': get_val('ref_lon', 0.0, float),
        'truelat1': get_val('truelat1', 0.0, float),
        'truelat2': get_val('truelat2', 0.0, float),
        'stand_lon': get_val('stand_lon', 0.0, float)
    }

    if parsed['stand_lon'] == 0.0 and parsed['ref_lon'] != 0.0:
        parsed['stand_lon'] = parsed['ref_lon']
    if parsed['truelat2'] == 0.0:
        parsed['truelat2'] = parsed['truelat1']

    return parsed


def compute_bounds_from_namelist(wps_dict):
    """
    Calculate min/max latitude and longitude for domain 1 from namelist parameters.
    """
    map_proj = wps_dict['map_proj']
    e_we = wps_dict['e_we']
    e_sn = wps_dict['e_sn']
    dx = wps_dict['dx']
    dy = wps_dict['dy']
    ref_lat = wps_dict['ref_lat']
    ref_lon = wps_dict['ref_lon']
    stand_lon = wps_dict['stand_lon']
    truelat1 = wps_dict['truelat1']
    truelat2 = wps_dict['truelat2']

    R = 6370000.0  # WRF earth radius

    i_c = (e_we + 1.0) / 2.0
    j_c = (e_sn + 1.0) / 2.0

    corner_indices = [
        (1.0, 1.0),
        (float(e_we), 1.0),
        (1.0, float(e_sn)),
        (float(e_we), float(e_sn))
    ]

    lats, lons = [], []

    if map_proj in ['lat-lon', 'cylindrical_equidistant']:
        for i, j in corner_indices:
            lon = ref_lon + (i - i_c) * dx
            lat = ref_lat + (j - j_c) * dy
            lats.append(lat)
            lons.append(lon)

    elif map_proj == 'mercator':
        k0 = math.cos(math.radians(truelat1))
        x0 = R * k0 * math.radians(ref_lon - stand_lon)
        y0 = R * k0 * math.log(math.tan(math.pi / 4.0 + math.radians(ref_lat) / 2.0))

        for i, j in corner_indices:
            x = x0 + (i - i_c) * dx
            y = y0 + (j - j_c) * dy
            lon = stand_lon + math.degrees(x / (R * k0))
            lat = math.degrees(2.0 * math.atan(math.exp(y / (R * k0))) - math.pi / 2.0)
            lats.append(lat)
            lons.append(lon)

    elif map_proj == 'lambert':
        lat1 = math.radians(truelat1)
        lat2 = math.radians(truelat2)

        if abs(truelat1 - truelat2) < 1e-5:
            n = math.sin(lat1)
        else:
            n = math.log(math.cos(lat1) / math.cos(lat2)) / math.log(
                math.tan(math.pi / 4.0 + lat2 / 2.0) / math.tan(math.pi / 4.0 + lat1 / 2.0)
            )

        F = (math.cos(lat1) * math.pow(math.tan(math.pi / 4.0 + lat1 / 2.0), n)) / n
        r0 = (R * F) / math.pow(math.tan(math.pi / 4.0 + math.radians(ref_lat) / 2.0), n)

        theta0 = n * math.radians(ref_lon - stand_lon)
        # Position of domain center relative to pole
        x0 = r0 * math.sin(theta0)
        y0 = r0 * math.cos(theta0)

        for i, j in corner_indices:
            dx_off = (i - i_c) * dx
            dy_off = (j - j_c) * dy

            if n > 0:
                x = x0 + dx_off
                y_from_pole = y0 - dy_off
                r = math.sqrt(x * x + y_from_pole * y_from_pole)
                theta = math.atan2(x, y_from_pole)
            else:
                x = x0 + dx_off
                y_from_pole = y0 + dy_off
                r = math.sqrt(x * x + y_from_pole * y_from_pole)
                theta = math.atan2(x, -y_from_pole)

            lon = stand_lon + math.degrees(theta / n)
            lat = math.degrees(2.0 * math.atan(math.pow((R * F) / r, 1.0 / n)) - math.pi / 2.0)
            lats.append(lat)
            lons.append(lon)

    elif map_proj in ['polar', 'polar_stereographic']:
        lat1 = math.radians(truelat1)
        k0 = 1.0 + math.sin(lat1)
        r0 = R * k0 * math.tan(math.pi / 4.0 - math.radians(ref_lat) / 2.0)

        theta0 = math.radians(ref_lon - stand_lon)
        x0 = r0 * math.sin(theta0)
        y0 = -r0 * math.cos(theta0)

        for i, j in corner_indices:
            x = x0 + (i - i_c) * dx
            y = y0 + (j - j_c) * dy

            r = math.sqrt(x * x + y * y)
            theta = math.atan2(x, -y)

            lon = stand_lon + math.degrees(theta)
            lat = math.degrees(math.pi / 2.0 - 2.0 * math.atan(r / (R * k0)))
            lats.append(lat)
            lons.append(lon)
    else:
        raise ValueError(f"Unsupported map_proj: {map_proj}")

    return min(lats), max(lats), min(lons), max(lons)


def compute_bounds_from_geo_em(geo_em_path):
    """
    Directly extract min/max lat and lon from geo_em.d01.nc if available.
    """
    try:
        import xarray as xr
        ds = xr.open_dataset(geo_em_path)
        lat_min = float(ds.XLAT_M.min())
        lat_max = float(ds.XLAT_M.max())
        lon_min = float(ds.XLONG_M.min())
        lon_max = float(ds.XLONG_M.max())
        ds.close()
        return lat_min, lat_max, lon_min, lon_max
    except Exception:
        import netCDF4 as nc
        ds = nc.Dataset(geo_em_path)
        xlat = ds.variables['XLAT_M'][:]
        xlong = ds.variables['XLONG_M'][:]
        ds.close()
        return float(xlat.min()), float(xlat.max()), float(xlong.min()), float(xlong.max())


def construct_daily_url(date_obj, url_template=None):
    """
    Construct GDEX THREDDS URL for a specific date.
    """
    year = date_obj.strftime('%Y')
    date_str = date_obj.strftime('%Y-%m-%d')

    if url_template:
        url = url_template.replace('{YEAR}', year).replace('{DATE}', date_str)
        return url

    # Default UCAR GDEX WACCM 6-hourly dataset path
    default_url = (
        f"https://tds.gdex.ucar.edu/thredds/fileServer/files/g/d313006/6hourly/{year}/"
        f"f.e22.beta02.FWSD.f09_f09_mg17.cesm2_2_beta02.forecast.001.cam.h3.{date_str}-00000.nc"
    )
    return default_url


def build_ncss_url(waccm_url, north, south, west, east, var_subset="all"):
    """
    Transform GDEX THREDDS URL into NCSS query URL with bounding parameters.
    """
    parsed = urlparse(waccm_url)
    path = parsed.path

    if '/fileServer/' in path:
        ncss_path = path.replace('/fileServer/', '/ncss/grid/')
    elif '/dodsC/' in path:
        ncss_path = path.replace('/dodsC/', '/ncss/grid/')
    elif '/ncss/grid/' in path:
        ncss_path = path
    else:
        ncss_path = re.sub(r'/thredds/(fileServer|dodsC|ncss)/', '/thredds/ncss/grid/', path)
        if '/thredds/ncss/grid/' not in ncss_path:
            ncss_path = '/thredds/ncss/grid/' + path.lstrip('/')

    base_url = f"{parsed.scheme}://{parsed.netloc}{ncss_path}"

    query_params = {
        'north': f"{north:.4f}",
        'south': f"{south:.4f}",
        'west': f"{west:.4f}",
        'east': f"{east:.4f}",
        'temporal': 'all',
        'accept': 'netcdf'
    }

    if var_subset and var_subset.lower() != 'all':
        query_params['var'] = var_subset
    else:
        query_params['var'] = 'all'

    return f"{base_url}?{urlencode(query_params, safe=',')}"


def build_opendap_url(waccm_url):
    """
    Transform GDEX THREDDS fileServer/NCSS URL into an OPeNDAP (dodsC) URL.
    """
    parsed = urlparse(waccm_url)
    path = parsed.path
    if '/fileServer/' in path:
        dods_path = path.replace('/fileServer/', '/dodsC/')
    elif '/ncss/grid/' in path:
        dods_path = path.replace('/ncss/grid/', '/dodsC/')
    else:
        dods_path = re.sub(r'/thredds/(fileServer|ncss/grid|ncss)/', '/thredds/dodsC/', path)
    return f"{parsed.scheme}://{parsed.netloc}{dods_path}"


def parse_mozbc_inp(inp_paths):
    """
    Parse one or more MOZBC .inp files (e.g., mozbc_d01.inp, WACCM.inp) for any WRF-Chem
    chemical option (chem_opt) and extract all required WACCM source variables from spc_map,
    plus mandatory vertical interpolation state variables (PS, T).
    """
    if isinstance(inp_paths, str):
        inp_paths = [p.strip() for p in inp_paths.split(',') if p.strip()]

    extracted_vars = ['PS', 'T']
    seen = set(extracted_vars)

    for inp_path in inp_paths:
        if not os.path.exists(inp_path):
            raise FileNotFoundError(f"MOZBC input file not found: {inp_path}")

        with open(inp_path, 'r', encoding='utf-8', errors='ignore') as f:
            txt = f.read()

        # Remove Fortran comments (! ...) line by line
        clean_lines = [line.split('!')[0] for line in txt.splitlines()]
        clean_txt = '\n'.join(clean_lines)

        # Extract optional moz_var_suffix
        m_suf = re.search(r'moz_var_suffix\s*=\s*[\'"]([^\'"]*)[\'"]', clean_txt, re.I)
        suffix = m_suf.group(1).strip() if m_suf else ''

        # Extract spc_map block up to next namelist variable assignment or '/'
        m_spc = re.search(r'spc_map\s*=\s*(.*?)(?=\n\s*[a-zA-Z_][a-zA-Z0-9_]*\s*=|\n\s*/|$)', clean_txt, re.S | re.I)
        spc_block = m_spc.group(1) if m_spc else clean_txt

        entries = re.findall(r'[\'"]([^\'"]+)[\'"]', spc_block)
        for entry in entries:
            if '->' not in entry:
                continue
            rhs = entry.split('->', 1)[1].split(';', 1)[0]
            # Strip floating-point and scientific-notation literals (e.g., 0.1216, 1.e9, 1.0E-5)
            rhs_no_nums = re.sub(r'\b\d+(?:\.\d*)?(?:[eEdD][+-]?\d+)?\b', ' ', rhs)
            for tok in re.findall(r'[A-Za-z_][A-Za-z0-9_]*', rhs_no_nums):
                var_name = tok + suffix if (suffix and not tok.endswith(suffix)) else tok
                if var_name not in seen:
                    seen.add(var_name)
                    extracted_vars.append(var_name)

    return extracted_vars


def filter_available_waccm_vars(requested_vars, sample_daily_url):
    """
    Validate requested variables against the remote OPeNDAP DDS metadata (if reachable)
    so that any variable absent in WACCM (relying on def_missing_var=.true. in mozbc)
    does not trigger an HTTP 400 Bad Request on THREDDS NCSS.
    """
    if not requested_vars:
        return requested_vars

    dds_url = build_opendap_url(sample_daily_url) + ".dds"
    try:
        import urllib.request
        req = urllib.request.Request(dds_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            dds_text = resp.read().decode('utf-8', errors='ignore')
        avail_tokens = set(re.findall(r'\b([A-Za-z_][A-Za-z0-9_]*)\s*(?:\[|;)', dds_text))
        if len(avail_tokens) > 20:
            valid_vars = []
            missing_vars = []
            for v in requested_vars:
                if v in avail_tokens:
                    valid_vars.append(v)
                else:
                    missing_vars.append(v)
            if missing_vars:
                print(
                    f" [Notice] Omitting {len(missing_vars)} variable(s) not present in remote WACCM dataset "
                    f"(will use def_missing_var in MOZBC): {', '.join(missing_vars)}"
                )
            if valid_vars:
                return valid_vars
    except Exception as e_dds:
        print(f" [Notice] Could not query remote DDS metadata ({e_dds}); proceeding with full requested variable list.")

    return requested_vars


def download_daily_subset_chunked_ncss(daily_url, s_file, north, south, west, east, var_list, chunk_size=18):
    """
    Download a daily subset by splitting the variable list into small NCSS chunks
    (default 18 variables per request) so each request completes in ~10-15s and
    never exceeds UCAR THREDDS's 60-second reverse-proxy gateway timeout (HTTP 504),
    then merge the variable chunks into s_file.
    """
    if not var_list:
        return False

    import xarray as xr
    chunks = [var_list[i:i + chunk_size] for i in range(0, len(var_list), chunk_size)]
    print(f" [Chunked-NCSS] Splitting {len(var_list)} variables into {len(chunks)} fast NCSS chunks (<= {chunk_size} vars/chunk)...", flush=True)

    part_files = []
    try:
        for c_idx, c_vars in enumerate(chunks, 1):
            part_path = f"{s_file}.part{c_idx}.nc"
            if os.path.exists(part_path):
                try:
                    os.remove(part_path)
                except OSError:
                    pass
            c_url = build_ncss_url(daily_url, north, south, west, east, ','.join(c_vars))
            print(f"   -> Chunk {c_idx}/{len(chunks)} ({len(c_vars)} vars)...", flush=True)
            cmd = ["wget", "-q", "--tries=3", "--waitretry=3", "--timeout=55", "-O", part_path, c_url]
            subprocess.run(cmd, check=True)
            if not os.path.exists(part_path) or os.path.getsize(part_path) < 50000:
                raise RuntimeError(f"Chunk {c_idx} file too small or missing")
            part_files.append(part_path)

        ds_parts = [xr.open_dataset(pf) for pf in part_files]
        merged_day = xr.merge(ds_parts, compat='override')
        merged_day.to_netcdf(s_file)
        for d in ds_parts:
            d.close()
        return os.path.exists(s_file) and os.path.getsize(s_file) > 1000000
    finally:
        for pf in part_files:
            if os.path.exists(pf):
                try:
                    os.remove(pf)
                except OSError:
                    pass


def download_daily_subset_opendap(daily_url, s_file, north, south, west, east, var_list=None):
    """
    Fallback downloader using OPeNDAP (dodsC) + xarray when THREDDS NCSS times out (HTTP 504).
    Loads variables one-by-one (~400 KB per HTTP request) so UCAR's proxy never times out.
    """
    dods_url = build_opendap_url(daily_url)
    print(f" [Fallback] Subsetting via OPeNDAP (per-variable chunking): {dods_url}", flush=True)
    import xarray as xr
    ds = xr.open_dataset(dods_url)
    try:
        lon_min_q = west % 360.0 if float(ds.lon.max()) > 180.0 else west
        lon_max_q = east % 360.0 if float(ds.lon.max()) > 180.0 else east
        sub = ds.sel(lat=slice(south, north), lon=slice(lon_min_q, lon_max_q))
        if var_list and isinstance(var_list, list):
            keep = [v for v in var_list if v in sub.variables]
            for c in ['hyam', 'hybm', 'P0', 'date', 'datesec']:
                if c in ds.variables and c not in keep:
                    keep.append(c)
        else:
            keep = list(sub.data_vars)

        loaded_vars = {}
        for v in keep:
            if v in sub.variables:
                loaded_vars[v] = sub[v].load()
            elif v in ds.variables:
                loaded_vars[v] = ds[v].load()

        out_ds = xr.Dataset(loaded_vars, coords=sub.coords, attrs=ds.attrs)
        out_ds.to_netcdf(s_file)
    finally:
        ds.close()
    return os.path.exists(s_file) and os.path.getsize(s_file) > 500000


def merge_netcdf_files(file_list, output_path, target_vars=None):
    """
    Merge multiple NetCDF files along time dimension into a single file.
    Handles mixed variable sets gracefully (e.g., if earlier daily files were downloaded
    with var=all and later daily files were downloaded with a MOZBC variable subset).
    """
    print(f"\n[Merging] Combining {len(file_list)} daily subset files into: {output_path}")

    # Method 1: xarray
    try:
        import xarray as xr
        print(" -> Using xarray for NetCDF time concatenation...")
        ds_list = [xr.open_dataset(f) for f in file_list]

        # Determine common data variables across all daily files (and restrict to target_vars if specified)
        common_vars = [
            v for v in ds_list[0].data_vars
            if all(v in ds.data_vars for ds in ds_list)
        ]
        if target_vars:
            target_set = set(target_vars) | {'hyam', 'hybm', 'P0', 'date', 'datesec'}
            filtered_common = [v for v in common_vars if v in target_set]
            if filtered_common:
                common_vars = filtered_common

        sub_list = [ds[common_vars] for ds in ds_list]
        merged = xr.concat(sub_list, dim='time')
        if 'time' in merged.coords:
            merged = merged.sortby('time')
        merged.to_netcdf(output_path)
        for ds in ds_list:
            ds.close()
        print(f" -> Merging successful ({len(common_vars)} variables)! File size: {os.path.getsize(output_path) / (1024*1024):.2f} MB")
        return True
    except Exception as e_xr:
        print(f" [Notice] xarray merge encountered issue: {e_xr}. Trying ncrcat...")

    # Method 2: NCO ncrcat
    try:
        cmd = ['ncrcat', '-O']
        if target_vars:
            cmd += ['-v', ','.join(target_vars)]
        cmd += file_list + [output_path]
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if os.path.exists(output_path):
            print(f" -> ncrcat merging successful! File size: {os.path.getsize(output_path) / (1024*1024):.2f} MB", flush=True)
            return True
    except Exception as e_nco:
        print(f" [Notice] ncrcat not available or failed: {e_nco}. Trying CDO...", flush=True)

    # Method 3: CDO mergetime
    try:
        cmd = ['cdo', '-O', 'mergetime'] + file_list + [output_path]
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if os.path.exists(output_path):
            print(f" -> CDO merging successful! File size: {os.path.getsize(output_path) / (1024*1024):.2f} MB")
            return True
    except Exception as e_cdo:
        print(f" [Error] All NetCDF merging methods failed: {e_cdo}")

    return False


def inject_mozbc_coordinates(output_path, start_date):
    """
    Ensures merged NetCDF contains hyam, hybm, P0, date, datesec required by MOZBC.
    """
    import netCDF4 as nc
    import numpy as np

    print(f"\n[Post-Processing] Verifying MOZBC required coordinates in: {output_path}")
    if not os.path.exists(output_path):
        print(f" [Warning] Output file does not exist: {output_path}")
        return False

    ds = nc.Dataset(output_path, 'r+')

    # 1. Check and inject hyam, hybm, P0
    if 'hyam' not in ds.variables or 'hybm' not in ds.variables or 'P0' not in ds.variables:
        print(" -> Injecting 1D vertical hybrid coefficients (hyam, hybm, P0)...")
        ref_paths = [
            "/gs/bs/tga-guc-lab/users/rosey/Simulations_GMMA/202312/waccm-2023-12-08-2023-12-13.nc",
            "/gs/bs/tga-guc-lab/compilation/wrf-chem-boundaries/mozbc/waccm-2023-10-11-2023-10-13.nc"
        ]
        hyam_data, hybm_data, p0_val = None, None, 100000.0
        for rpath in ref_paths:
            if os.path.exists(rpath):
                try:
                    ref_ds = nc.Dataset(rpath, 'r')
                    if 'hyam' in ref_ds.variables and 'hybm' in ref_ds.variables:
                        hyam_data = ref_ds.variables['hyam'][:]
                        hybm_data = ref_ds.variables['hybm'][:]
                    if 'P0' in ref_ds.variables:
                        p0_val = float(ref_ds.variables['P0'][:])
                    ref_ds.close()
                    if hyam_data is not None:
                        break
                except Exception:
                    pass

        # If not found locally, fetch from remote OPeNDAP
        if hyam_data is None:
            try:
                remote_url = (
                    "https://tds.gdex.ucar.edu/thredds/dodsC/files/g/d313006/6hourly/2026/"
                    "f.e22.beta02.FWSD.f09_f09_mg17.cesm2_2_beta02.forecast.001.cam.h3.2026-07-25-00000.nc"
                )
                rem_ds = nc.Dataset(remote_url, 'r')
                hyam_data = rem_ds.variables['hyam'][:]
                hybm_data = rem_ds.variables['hybm'][:]
                if 'P0' in rem_ds.variables:
                    p0_val = float(rem_ds.variables['P0'][:])
                rem_ds.close()
            except Exception as e:
                print(f" [Warning] Could not fetch remote hyam/hybm: {e}")

        if hyam_data is not None and 'lev' in ds.dimensions:
            if 'hyam' not in ds.variables:
                v_a = ds.createVariable('hyam', 'f8', ('lev',))
                v_a.long_name = "hybrid A coefficient at layer midpoints"
                v_a[:] = hyam_data
            if 'hybm' not in ds.variables:
                v_b = ds.createVariable('hybm', 'f8', ('lev',))
                v_b.long_name = "hybrid B coefficient at layer midpoints"
                v_b[:] = hybm_data
            if 'P0' not in ds.variables:
                v_p0 = ds.createVariable('P0', 'f8')
                v_p0.long_name = "reference pressure"
                v_p0.units = "Pa"
                v_p0[:] = p0_val

    # 2. Check and inject date, datesec
    time_len = len(ds.dimensions['time']) if 'time' in ds.dimensions else 0
    if time_len > 0 and ('date' not in ds.variables or 'datesec' not in ds.variables):
        print(f" -> Injecting date and datesec for {time_len} time slices...")
        import datetime, math
        if isinstance(start_date, str):
            start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d').date()
        num_days = int(math.ceil(time_len / 4.0))
        dates = []
        datesecs = []
        for d in range(num_days):
            cur_d = start_date + datetime.timedelta(days=d)
            d_int = int(cur_d.strftime('%Y%m%d'))
            for s in [0, 21600, 43200, 64800]:
                if len(dates) < time_len:
                    dates.append(d_int)
                    datesecs.append(s)

        if 'date' not in ds.variables:
            v_date = ds.createVariable('date', 'i4', ('time',))
            v_date.long_name = "current date (YYYYMMDD)"
            v_date[:] = np.array(dates, dtype='i4')
        else:
            ds.variables['date'][:] = np.array(dates, dtype='i4')

        if 'datesec' not in ds.variables:
            v_sec = ds.createVariable('datesec', 'i4', ('time',))
            v_sec.long_name = "current seconds of current date"
            v_sec.units = "s"
            v_sec[:] = np.array(datesecs, dtype='i4')
        else:
            ds.variables['datesec'][:] = np.array(datesecs, dtype='i4')

    ds.close()
    print(" -> Coordinate verification and injection complete.")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Automated WACCM / CAM-chem Spatial Subsetting & Merging Downloader for WRF-Chem."
    )
    parser.add_argument(
        '--start-date', '-s', default=None,
        help="Simulation start date (YYYY-MM-DD). If omitted, read from namelist.wps."
    )
    parser.add_argument(
        '--end-date', '-e', default=None,
        help="Simulation end date (YYYY-MM-DD). If omitted, read from namelist.wps."
    )
    parser.add_argument(
        '--namelist', '-n', default="./namelist.wps",
        help="Path to namelist.wps file (default: ./namelist.wps)"
    )
    parser.add_argument(
        '--geo-em', '-g', default=None,
        help="Optional path to geo_em.d01.nc for exact grid bounds"
    )
    parser.add_argument(
        '--url-template', default=None,
        help="Custom URL template containing {YEAR} and {DATE}"
    )
    parser.add_argument(
        '--output-dir', '-d', default=DEFAULT_MOZBC_DIR,
        help=f"Target output directory (default: {DEFAULT_MOZBC_DIR})"
    )
    parser.add_argument(
        '--output', '-o', default=None,
        help="Output merged NetCDF filename or absolute path (default: waccm-START_DATE-END_DATE.nc in output-dir)"
    )
    parser.add_argument(
        '--buffer', '-b', type=float, default=1.0,
        help="Latitude/Longitude buffer in degrees added to domain bounds (default: 1.0)"
    )
    parser.add_argument(
        '--mozbc-inp', '-m', default=None,
        help="Optional path(s) (comma-separated) to MOZBC .inp file(s) (e.g., mozbc_d01.inp, WACCM.inp) to dynamically extract required WACCM variables for any chemical option (chem_opt)."
    )
    parser.add_argument(
        '--vars', default="all",
        help="Variable selection ('all' or comma-separated list; if combined with --mozbc-inp, extra variables in --vars are merged with MOZBC variables, default: all)"
    )
    parser.add_argument(
        '--keep-daily', action='store_true',
        help="Keep individual daily subset NetCDF files after merging"
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help="Display generated daily URLs, bounding box, and merged output path without downloading"
    )

    args = parser.parse_args()

    # Step 1: Parse namelist.wps if available
    wps_params = {}
    if os.path.exists(args.namelist):
        wps_params = parse_namelist_wps(args.namelist)

    # Step 2: Determine start and end dates
    start_date = parse_date_string(args.start_date)
    if not start_date and 'start_date' in wps_params:
        start_date = parse_date_string(wps_params['start_date'])

    end_date = parse_date_string(args.end_date)
    if not end_date and 'end_date' in wps_params:
        end_date = parse_date_string(wps_params['end_date'])

    if not start_date or not end_date:
        print("[Error] Could not determine start date and end date.")
        print("        Please provide --start-date YYYY-MM-DD --end-date YYYY-MM-DD or ensure namelist.wps contains start_date and end_date.")
        sys.exit(1)

    if start_date > end_date:
        print(f"[Error] Start date ({start_date}) cannot be after end date ({end_date}).")
        sys.exit(1)

    start_str = start_date.strftime('%Y-%m-%d')
    end_str = end_date.strftime('%Y-%m-%d')
    default_filename = f"waccm-{start_str}-{end_str}.nc"

    # Handle output destination path
    if not args.output:
        os.makedirs(args.output_dir, exist_ok=True)
        final_output_path = os.path.join(args.output_dir, default_filename)
    else:
        if os.path.isabs(args.output) or '/' in args.output:
            final_output_path = args.output
            os.makedirs(os.path.dirname(os.path.abspath(final_output_path)), exist_ok=True)
        else:
            os.makedirs(args.output_dir, exist_ok=True)
            final_output_path = os.path.join(args.output_dir, args.output)

    # Working temp directory for daily files (same folder as final output)
    work_dir = os.path.dirname(os.path.abspath(final_output_path))

    # Step 3: Calculate spatial bounds
    lat_min, lat_max, lon_min, lon_max = None, None, None, None
    source_msg = ""

    geo_candidate = args.geo_em
    if not geo_candidate and os.path.exists(args.namelist):
        dir_name = os.path.dirname(os.path.abspath(args.namelist))
        auto_geo = os.path.join(dir_name, "geo_em.d01.nc")
        if os.path.exists(auto_geo):
            geo_candidate = auto_geo

    if geo_candidate and os.path.exists(geo_candidate):
        try:
            lat_min, lat_max, lon_min, lon_max = compute_bounds_from_geo_em(geo_candidate)
            source_msg = f"Extracted directly from {geo_candidate}"
        except Exception as e:
            print(f"[Warning] Failed reading geo_em ({e}), falling back to namelist.wps parsing...")

    if lat_min is None:
        if not wps_params:
            print(f"[Error] Neither geo_em file nor namelist.wps ({args.namelist}) was found.")
            sys.exit(1)
        lat_min, lat_max, lon_min, lon_max = compute_bounds_from_namelist(wps_params)
        source_msg = f"Calculated from {args.namelist} (proj: {wps_params['map_proj']}, dx: {wps_params['dx']}m)"

    south = max(-90.0, lat_min - args.buffer)
    north = min(90.0, lat_max + args.buffer)
    west = lon_min - args.buffer
    east = lon_max + args.buffer

    # Step 3.5: Resolve flexible variable selection (--mozbc-inp and/or --vars)
    resolved_var_list = None
    var_source_desc = "all (full CAM-chem dataset)"

    if args.mozbc_inp:
        moz_vars = parse_mozbc_inp(args.mozbc_inp)
        resolved_var_list = list(moz_vars)
        var_source_desc = f"MOZBC config ({args.mozbc_inp}: {len(moz_vars)} vars)"
        if args.vars and args.vars.strip().lower() not in ('all', 'mozbc'):
            extra_vars = [v.strip() for v in args.vars.split(',') if v.strip()]
            for ev in extra_vars:
                if ev not in resolved_var_list:
                    resolved_var_list.append(ev)
            var_source_desc += f" + extra --vars ({len(resolved_var_list)} total)"
    elif args.vars and args.vars.strip().lower() != 'all':
        if os.path.exists(args.vars.strip()) and args.vars.strip().endswith('.inp'):
            resolved_var_list = parse_mozbc_inp(args.vars.strip())
            var_source_desc = f"MOZBC config ({args.vars.strip()}: {len(resolved_var_list)} vars)"
        else:
            resolved_var_list = [v.strip() for v in args.vars.split(',') if v.strip()]
            if 'PS' not in resolved_var_list:
                resolved_var_list.insert(0, 'PS')
            var_source_desc = f"Custom --vars ({len(resolved_var_list)} vars)"

    if resolved_var_list and not args.dry_run:
        sample_url = construct_daily_url(start_date, args.url_template)
        resolved_var_list = filter_available_waccm_vars(resolved_var_list, sample_url)

    effective_vars_str = ','.join(resolved_var_list) if resolved_var_list else "all"

    # Step 4: Generate daily tasks
    curr_date = start_date
    daily_tasks = []
    while curr_date <= end_date:
        d_str = curr_date.strftime('%Y-%m-%d')
        daily_url = construct_daily_url(curr_date, args.url_template)
        ncss_url = build_ncss_url(daily_url, north, south, west, east, effective_vars_str)
        sub_file = os.path.join(work_dir, f"subset_waccm_{d_str}.nc")
        daily_tasks.append((d_str, daily_url, ncss_url, sub_file))
        curr_date += datetime.timedelta(days=1)

    print("=" * 75)
    print("  WACCM / CAM-chem Subsetting & Merging Downloader (WRF-Chem)")
    print("=" * 75)
    print(f"Simulation Period     : {start_str} to {end_str} ({len(daily_tasks)} days)")
    print(f"Domain Bounds Source  : {source_msg}")
    print(f"Domain 1 Raw Bounds   : Lat [{lat_min:.4f}, {lat_max:.4f}], Lon [{lon_min:.4f}, {lon_max:.4f}]")
    print(f"Buffer Added          : +/- {args.buffer:.2f} degrees")
    print(f"Subset Bounding Box   : South: {south:.4f}, North: {north:.4f}, West: {west:.4f}, East: {east:.4f}")
    print(f"Variable Selection    : {var_source_desc}")
    print(f"Target Output Path    : {final_output_path}")
    print(f"Daily Files Count     : {len(daily_tasks)}")
    print("=" * 75)

    if args.dry_run:
        print("\n[Dry-Run] Daily Download Tasks & NCSS Query URLs:")
        for d_str, d_url, n_url, s_file in daily_tasks:
            print(f"\n  Date: {d_str} -> Target: {s_file}")
            print(f"  NCSS URL: {n_url}")
        print("\n[Dry-Run] Exiting without downloading or merging.")
        return

    # Step 5: Download each daily subset
    downloaded_files = []
    failed_dates = []

    for idx, (d_str, d_url, n_url, s_file) in enumerate(daily_tasks, 1):
        print(f"\n[{idx}/{len(daily_tasks)}] Downloading subset for {d_str}...", flush=True)

        if os.path.exists(s_file):
            if os.path.getsize(s_file) > 1000000:
                size_mb = os.path.getsize(s_file) / (1024 * 1024)
                print(f" [Skip] {s_file} already downloaded ({size_mb:.2f} MB)", flush=True)
                downloaded_files.append(s_file)
                continue
            else:
                try:
                    os.remove(s_file)
                except OSError:
                    pass

        cmd = ["wget", "--tries=1", "--timeout=45", "-O", s_file, n_url]
        try:
            res = subprocess.run(cmd, check=True)
            if res.returncode == 0 and os.path.exists(s_file) and os.path.getsize(s_file) > 1000000:
                size_mb = os.path.getsize(s_file) / (1024 * 1024)
                print(f" [Success] Saved {s_file} ({size_mb:.2f} MB)", flush=True)
                downloaded_files.append(s_file)
                continue
        except subprocess.CalledProcessError:
            pass

        if os.path.exists(s_file) and os.path.getsize(s_file) <= 1000000:
            try:
                os.remove(s_file)
            except OSError:
                pass

        if resolved_var_list and len(resolved_var_list) > 18:
            print(f" [Notice] Single-shot NCSS timed out for {d_str}; trying Chunked-NCSS (18 vars/chunk)...", flush=True)
            try:
                if download_daily_subset_chunked_ncss(d_url, s_file, north, south, west, east, resolved_var_list, chunk_size=18):
                    size_mb = os.path.getsize(s_file) / (1024 * 1024)
                    print(f" [Success] Saved via Chunked-NCSS {s_file} ({size_mb:.2f} MB)", flush=True)
                    downloaded_files.append(s_file)
                    continue
            except Exception as e_chunk:
                print(f" [Warning] Chunked-NCSS failed for {d_str} ({e_chunk}), switching to OPeNDAP...", flush=True)

        print(f" [Notice] Switching to OPeNDAP per-variable subsetting for {d_str}...", flush=True)
        try:
            if download_daily_subset_opendap(d_url, s_file, north, south, west, east, resolved_var_list):
                size_mb = os.path.getsize(s_file) / (1024 * 1024)
                print(f" [Success] Saved via OPeNDAP {s_file} ({size_mb:.2f} MB)", flush=True)
                downloaded_files.append(s_file)
                continue
        except Exception as e_dods:
            print(f" [Warning] OPeNDAP failed for {d_str} ({e_dods}), trying curl...", flush=True)

        if os.path.exists(s_file) and os.path.getsize(s_file) <= 1000000:
            try:
                os.remove(s_file)
            except OSError:
                pass

        cmd_curl = ["curl", "-sS", "-L", "--retry", "2", "--retry-delay", "3", "--max-time", "90", "-o", s_file, n_url]
        try:
            subprocess.run(cmd_curl, check=True)
            if os.path.exists(s_file) and os.path.getsize(s_file) > 1000000:
                size_mb = os.path.getsize(s_file) / (1024 * 1024)
                print(f" [Success] Saved via curl {s_file} ({size_mb:.2f} MB)", flush=True)
                downloaded_files.append(s_file)
                continue
        except Exception as e_curl:
            print(f" [Error] curl also failed for {d_str}: {e_curl}", flush=True)

        if os.path.exists(s_file) and os.path.getsize(s_file) <= 1000000:
            try:
                os.remove(s_file)
            except OSError:
                pass
        failed_dates.append(d_str)

    if failed_dates:
        print(f"\n[Warning] Failed downloads for {len(failed_dates)} date(s): {', '.join(failed_dates)}")

    # Determine contiguous prefix of downloaded files starting from start_date
    contiguous_files = []
    for d_str, _, _, s_file in daily_tasks:
        if os.path.exists(s_file) and os.path.getsize(s_file) > 1000000:
            contiguous_files.append(s_file)
        else:
            print(f"[Notice] Contiguous daily sequence stops before missing date {d_str} ({len(contiguous_files)} contiguous days available).")
            break

    # Step 6: Merge contiguous daily files into single NetCDF file (waccm-START-END.nc)
    if contiguous_files:
        success = merge_netcdf_files(contiguous_files, final_output_path, target_vars=resolved_var_list)
        if success:
            inject_mozbc_coordinates(final_output_path, start_date)
            print(f"\n============================================================")
            print(f" SUCCESS: Merged {len(contiguous_files)} contiguous WACCM days saved to:")
            print(f"          {os.path.abspath(final_output_path)}")
            print(f"============================================================")

            if not args.keep_daily and not failed_dates:
                print("Cleaning up intermediate daily subset files...")
                for f in contiguous_files:
                    try:
                        os.remove(f)
                    except OSError:
                        pass
                print("Cleanup complete.")

    if failed_dates:
        sys.exit(1)

if __name__ == "__main__":
    main()
