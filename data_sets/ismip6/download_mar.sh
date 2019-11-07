#!/bin/bash

wget -nc --reject "index.html*" -r -np -nH --cut-dirs=1 -l40 ftp://ftp.climato.be/fettweis/MARv3.9/ISMIP6/GrIS/ERA_1958-2017 $1

wget -nc --reject "index.html*" -r -np -nH --cut-dirs=1 -l40 ftp://ftp.climato.be/fettweis/MARv3.9/ISMIP6/GrIS/MIROC5-rcp85_2006_2100 $1

