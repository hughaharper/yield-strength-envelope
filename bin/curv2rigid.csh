#!/bin/csh -f
#
#  D. Sandwell FEB 10 2010
#
unset noclobber
#
# script to call curv2rigid with a smaller grid to increase speed
#
if ( $#argv < 7) then
 echo " "
 echo "Usage: curv2rigid.csh age_in.grd curv_in.grd mt_in.grd N rigid_out.grd moment_out.grd yld_depth.grd"
 echo " "
 exit 1
endif 
#
#
#  downsample the age and curv grids
#
set NDEC = `echo $4`
#echo $NDEC
set INC = `grdinfo $1 -C | cut -f8`
set NEW_INC = `gmtmath -Q $INC $NDEC MUL = `
#echo $INC $NEW_INC
#
#  downsample the curv and age grids
#
grdsample $1 -Ga.grd -I$NEW_INC 
grdsample $2 -Gc.grd -I$NEW_INC 
grdsample $3 -Gt.grd -I$NEW_INC 
#
#  run the curv2rigid code
#
./curv2rigid a.grd c.grd t.grd r.grd m.grd y.grd
#
#  upsample the output grids
#
grdsample r.grd -G$5 -I$INC 
grdsample m.grd -G$6 -I$INC 
grdsample y.grd -G$7 -I$INC 
#
#  clean the mess
#
rm a.grd c.grd t.grd r.grd m.grd y.grd
