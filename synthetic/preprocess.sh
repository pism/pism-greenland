#!/bin/bash

for grid in 250 500 1000 2000 5000; do
    python create_geometry.py -g $grid pism_outletglacier_g${grid}m.nc
    ncatted -a _FillValue,thk,d,,  pism_outletglacier_g${grid}m.nc
done
