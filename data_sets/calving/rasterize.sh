#!/bin/bash

infile=$1
outfile=$2
ires=$3
ores=$4

gdal_rasterize -l vct_g${ires}m -a threshold -tr $ores $ores -a_nodata 0.2 -te -653000.0 -3384500.0 879700.0 -632600.0 -ot Float16 "$1" "$2"
ncrename -v Band1,vonmises_calving_threshold $2
ncatted -a units,vonmises_calving_threshold,o,c,"MPa" $2
