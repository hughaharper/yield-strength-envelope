#!/bin/csh -f
#
#  D. Sandwell FEB 10 2010
#
unset noclobber
#
# script to call curv2rigid with a smaller grid to increase speed
#
if ( $#argv < 3) then
 echo " "
 echo "Usage: age2dsf.csh age_xy.grd N dsf_out.grd"
 echo " "
 exit 1
endif 
#
#
#  downsample the age and curv grids
#
set NDEC = `echo $2`
#echo $NDEC
set INC = `grdinfo $1 -C | cut -f8`
set NEW_INC = `gmtmath -Q $INC $NDEC MUL = `
#echo $INC $NEW_INC
#
#  downsample the age grid
#
grdsample $1 -Ga2d.grd -I$NEW_INC 
#
#  run the curv2rigid code
#
./age2dsf a2d.grd dsf.grd
#
#  upsample the output grids
#
grdsample dsf.grd -G$3 -I$INC  
#
#  clean the mess
#
rm a2d.grd dsf.grd