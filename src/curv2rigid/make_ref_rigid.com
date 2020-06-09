#
#  script to make a three files of reference grids for the.grd_lookup_interp program
#
# 1) run make_cr_age.m to make the starting age and curvature grids
#
# 2) run curv2rigid
#
curv2rigid age.grd_lookup curv.grd_lookup rigid.grd_lookup moment.grd_lookup yld.grd_lookup
#
# 3) run.grd_lookupedit on these three files to set the limits of the x and y axes
#
grdedit rigid.grd_lookup -R1/180/-1.e-5/1.e-5 -V
grdedit moment.grd_lookup -R1/180/-1.e-5/1.e-5 -V
grdedit yld.grd_lookup -R1/180/-1.e-5/1.e-5 -V
#
# link these files to the place where they will be used
#
