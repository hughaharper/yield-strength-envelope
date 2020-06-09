#
#  scropt to test the grd_interp interpolation code.
#
# 1) go to the curv2rigid directory and create the lookup tables and
#    link them here
ln -s ../curv2rigid/*lookup .
grdinfo *lookup
#
# 2) use curve to rigid to make much smaller grids for a reasonable area
#
# age.grd			diff.grd		rigid_curv2rigid.grd
# curv.grd		moment_curv2rigid.grd	yld_curv2rigid.grd
#
#
# 3) use grd_interp to make a new rigidity grid called  rigid_out.grd 
#
ln -s ../curv2rigid/*_sub.grd .
grd_interp rigid.grd_lookup age_sub.grd curv_sub.grd rigid_out.grd 
grdmath rigid_sub.grd rigid_out.grd SUB = diff.grd
grdinfo diff.grd
#
# 4) use grd_interp to make new yld and moment grids
#
grd_interp moment.grd_lookup age_sub.grd curv_sub.grd moment_out.grd 
grdmath moment_sub.grd moment_out.grd SUB = diff.grd
grd_interp yld.grd_lookup age_sub.grd curv_sub.grd yld_out.grd 
grdmath yld_sub.grd yld_out.grd SUB = diff.grd
