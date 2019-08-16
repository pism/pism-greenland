#!/bin/bash
# rename variable to PISM's liking

ncrename -O -v sftgif,land_ice_area_fraction_retreat $1 $1
