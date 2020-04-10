from netCDF4 import Dataset as NC
from cdo import Cdo
from datetime import datetime, timedelta
import gdal
from glob import glob
from os.path import join, split
cdo = Cdo()
reftime = "2008-1-1"

components = ["vx", "vy", "vv"]

download_dir = "data"

for glacier in ["W69.10N"]:
    f_tiffs = glob(f"{download_dir}/TSX_{glacier}_*_v02.0.tif")
    
    for f_tiff in f_tiffs:
        f_nc = f_tiff.replace(".tif", ".nc")
        print(f"Converting {f_tiff} to {f_nc}")
        # use gdal's python binging to convert GeoTiff to netCDF
        # advantage of GDAL: it gets the projection information right
        # disadvantage: the variable is named "Band1", lacks metadata
        ds = gdal.Open(f_tiff)
        ds = gdal.Translate(f_nc, ds)
        ds = None

        # This deduces the mid-point (nominal) date from the filename
        _ , _ , start_date_str, end_date_str, _, var, _ = f_nc.split("_")
        start_date = datetime.strptime(start_date_str, "%d%b%y")
        end_date = datetime.strptime(end_date_str, "%d%b%y")
        nominal_date = start_date + (end_date - start_date) / 2
        # Set the time axis
        out_file = split(f_nc)[-1]
        cdo.settaxis(nominal_date.isoformat(), input=f"""-setreftime,{reftime} -setattribute,{var}@units="m year-1" -chname,Band1,{var} {f_nc}""", output=f"{download_dir}/cdo_{out_file}", options="-f nc4 -z zip_2")

# Merge the indiviual compontents and sort by time using "mergetime"
for v in components:
    fs = glob(f"{download_dir}/cdo_TSX_{glacier}_*_{v}_v02.0.nc")
    cdo.mergetime(input=fs, output=f"TSX_{glacier}_{v}_merged.nc", options="-f nc4 -z zip_2")

# Create the final merged file
cdo.merge(input=[f"TSX_{glacier}_{v}_merged.nc" for v in components], output=f"TSX_{glacier}_2008_2020.nc", options="-f nc4 -z zip_2")
