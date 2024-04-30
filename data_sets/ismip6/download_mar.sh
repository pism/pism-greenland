#!/bin/bash

# wget -nc --reject "index.html*" -r -np -nH --cut-dirs=1 -l40 ftp://ftp.climato.be/fettweis/MARv3.9/ISMIP6/GrIS/ERA_1958-2017 $1

for gcm in UKESM1-CM6; do
    for year in {2008..2014}; do
        wget -nc --reject "index.html*" -r -np -nH --cut-dirs=1 -l40 ftp://ftp.climato.be/fettweis/MARv3.9/ISMIP6/GrIS/${gcm}-histo_1950_2014/MARv3.9-monthly-${gcm}-histo-${year}.nc $1
    done
done

