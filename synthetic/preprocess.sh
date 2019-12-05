#!/bin/bash

for grid in 100 200 250 500 1000 2000; do
    python create_jib.py -s -g $grid pism_synth_jib_g${grid}m.nc
done
