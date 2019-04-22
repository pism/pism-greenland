#!/bin/bash

odir=$1

mkdir -p ${odir}/io

cd ${odir}_tmp
for f in ex_*.nc; do
    echo $f
    if  [ ! -f "../${odir}/spatial/${f}" ]; then
        ncks -O -4 -L 3 ${f} ../${odir}/spatial/${f}
        extract_interface.py -t ice_ocean -o ../${odir}/io/${f} ../${odir}/spatial/${f}
    fi
done
cd ../
for f in ${odir}/state/*.nc; do
    ncks -O -4 -L 3 ${f} ${f}
done
