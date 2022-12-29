#!/bin/bash

# (c) 2021-22 Andy Aschwanden

set -x -e

# Produces geothermal hot spot in variable bheatflx.

# Version 2:  instead of circular hot blob near divide, an elliptical blob
#   along the NE Greenland ice stream route; the ends of the long, narrow
#   ellipse were eye-balled to be at the original center (-32000m, -1751180m)
#   and at (103000m,-1544330) (in projection coords used by SeaRISE)
# here are relevant octave computations:
# > 0.5*(-32000+103000)
# ans =  35500  # x coord of new center
# > 0.5*(-1751180 + -1544330)
# ans = -1647755  # y coord center
# > theta = atan( (-1544330 - (-1647755)) / (103000 - 35500) )
# theta =  0.99256  # rotation angle; = 56.9 deg
# > cos(theta)
# ans =  0.54655
# > sin(theta)
# ans =  0.83743
# > a = sqrt( (103000 - 35500)^2 +  (-1544330 - (-1647755))^2 )
# a =  1.2350e+05
# > b = 50000^2 / a  #  set b so that ab=R^2 where R = 50 km is orig radius
# b =  2.0242e+04

# Version 1:  The spot is at the
# source area of the NE Greenland ice stream.  The spot has the location,
# magnitude and extent suggested by
#    M. Fahnestock, et al (2001).  High geothermal heat flow, basal melt, and 
#    the origin of rapid ice flow in central Greenland, Science vol 294, 2338--2342.
# Uses NCO (ncrename, ncap2, ncks, ncatted).
# Run preprocess.py first to generate $PISMVERSION.
# center of hot spot is  (-40 W lon, 74 deg N lat)  which is
#   (x,y) = (-32000m, -1751180m)  in projection already used in $DATANAME
# parameters: radius of spot = 50000m and heat is 970 mW m-2 from Fahnstock et al 2001

# NOTE 5/20/2014:
# Switch to EPSG:3413, use coordinates in EPSG:3413 projection
# 207000 -1630000

# # Updated 9/2020

# import numpy as np
# from pyproj import Transformer
# # Transform from Bamber to EPSG grid
# transformer = Transformer.from_crs("+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=-39 +lat_0=90 +lat_ts=71 +units=m", "EPSG:3413")
# # Transform the two points
# p1 = np.array(transformer.transform(-32000, -1751180))
# p2 = np.array(transformer.transform(103000,-1544330))
# # The center is at
# p_c = 0.5 * (p1 + p2)
# = (206921.79363337, -1630125.19206727)
# a = np.sqrt( (p2[0] - p_c[0])**2 + (p2[1] - p_c[1])**2)
# = 123133
# b = 50e3**2 / a
# = 20303



set -e  -x # exit on error

add_hotspot() {

# center:
XSPOT=206000
# negate sign because -(-1630000) does not work in bash
NEGYSPOT=1630000

# parameters for ellipse and rotation in m and heat to apply there
ASPOT=123000
BSPOT=20300

COSTHETA=0.54655
SINTHETA=0.83743

GHFSPOT=970   # from Fahnstock et al 2001; in mW m-2

ncrename -v bheatflx,bheatflxSR $OUTFILE  # keep Shapiro & Ritzwoller

# do equivalent of Matlab's:  [xx,yy] = meshgrid(x,y)
ncap2 -O -s 'zero=0.0*bheatflxSR' $OUTFILE $OUTFILE
ncap2 -O -s 'xx=zero+x' $OUTFILE $OUTFILE
ncap2 -O -s 'yy=zero+y' $OUTFILE $OUTFILE
XIROT="xi=${COSTHETA}*(xx-${XSPOT})+${SINTHETA}*(yy+${NEGYSPOT})"
ncap2 -O -s $XIROT $OUTFILE $OUTFILE
ETAROT="eta=-${SINTHETA}*(xx-${XSPOT})+${COSTHETA}*(yy+${NEGYSPOT})"
ncap2 -O -s $ETAROT $OUTFILE $OUTFILE

# filled ellipse is:   xi^2/a^2 + eta^2/b^2 < 1
ELLLEFT="eleft=(xi*xi)/(${ASPOT}*${ASPOT})"
ncap2 -O -s $ELLLEFT $OUTFILE $OUTFILE
ELLRIGHT="eright=(eta*eta)/(${BSPOT}*${BSPOT})"
ncap2 -O -s $ELLRIGHT $OUTFILE $OUTFILE
ncap2 -O -s 'hotmask=(-eleft+eright-1<0)' $OUTFILE $OUTFILE

# actually create hot spot
NEWBHEATFLX="bheatflx=hotmask*${GHFSPOT}+!hotmask*bheatflxSR"
ncap2 -O -s $NEWBHEATFLX $OUTFILE $OUTFILE

# ncap2 leaves hosed attributes; start over
ncatted -a units,bheatflx,c,c,"W m-2" $OUTFILE
ncatted -a long_name,bheatflx,c,c,"basal geothermal flux" $OUTFILE
ncatted -a propose_standard_name,bheatflx,c,c,"lithosphere_upward_heat_flux" $OUTFILE

# clear out the temporary variables and only leave additional 'bheatflxSR'
ncks -O -x -v xx,yy,xi,eta,eleft,eright,hotmask,zero,bheatflxSR $OUTFILE $OUTFILE
}


# Create a buffer that is a multiple of the grid resolution
# and works for grid resolutions up to 36km.
buffer_x=148650
buffer_y=130000
xmin=$((-638000 - $buffer_x - 468000))
ymin=$((-3349600 - $buffer_y))
xmax=$((864700 + $buffer_x))
ymax=$((-657600 + $buffer_y))


for GRID in  450; do
    OUTFILE=Geothermal_heatflux_map_v2.1_g${GRID}m.nc
    gdalwarp  -overwrite  -r average -co FORMAT=NC4 -co COMPRESS=DEFLATE -co ZLEVEL=2 -s_srs EPSG:3413 -t_srs EPSG:3413 -te $xmin $ymin $xmax $ymax -tr $GRID $GRID geothermal_heat_flow_map_10km.nc $OUTFILE
    ncrename -v Band1,bheatflx $OUTFILE 
    ncatted -a units,bheatflx,o,c,"mW m-2" -a _FillValue,bheatflx,d,, -a missing_value,bheatflx,d,, $OUTFILE
    add_hotspot
done


