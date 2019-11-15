#!/bin/bash

for grid in 100 250 500 1000 2000 5000; do
    python create_geometry.py -g $grid pism_outletglacier_g${grid}m.nc
    ncatted  -a _FillValue,topg,d,,  -a _FillValue,usurf,d,, -a _FillValue,thk,d,,  pism_outletglacier_g${grid}m.nc
    python create_jib.py -g $grid pism_jib_g${grid}m.nc
done
