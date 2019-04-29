#!/bin/bash

infile=$1
outfile=$2
res=$3
gdal_rasterize -l vct -a threshold -tr $res $res -a_nodata 1.0 -te -653000.0 -3384500.0 879700.0 -632600.0 -ot Float16 "$1" "$2"
ncrename -v Band1,vonmises_calving_threshold $2
ncatted -a units,vonmises_calving_threshold,o,c,"MPa" $2
