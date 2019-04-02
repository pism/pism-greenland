#!/bin/bash

wget -nc  --reject "index.html*" -r -np -nH --cut-dirs=1 -l40 ftp://ftp.climato.be/fettweis/MARv3.9/ISMIP6/GrIS $1

