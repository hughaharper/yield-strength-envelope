#!/bin/bash
# visco case for a bigger model domain

awk '{print $1,$2}' gps.por.dat > gps.sample
#awk '{print $1,$2,$3,$4}' insar.por.dat > los.sample
awk '{print $1,$2,$3,$4}' alos_asc.por.dat > los.sample

# $1 = fault.inp

cat $1 | while read id dep obs dip rec name eqname
do 
  echo "id:" $id
  echo "locking depth" $dep
  echo "name" $name

  echo "shift the fault model x by 500"
  awk 'NR<3 {print $0}' $name > header
  awk 'NR>2 {print $1+500,$2+500,$3,$4,$5,$6,$7}' $name > shift
  cat header shift > name.tmp

# secular rate
  maxwellx2 1 -60.0 $dep $obs 9999 0 name.tmp 0 U.$id.grd V.$id.grd W.$id.grd -60.0 1.e19 0.25 1 

# viscoelastic effect from previous earthquakes
  if [ $(echo "$rec == 0" | bc) -eq 0 ]
  then
#    echo $id $dep $dip $rec name.tmp $eqname > tmp3
#    episodic 1 tmp3 0 2010 0 Ueq.grd Veq.grd Weq.grd -60 1.e19 0.25 1
    run_eq.com name.tmp $eqname $rec 1900 2010 2011 -60 1.e19 0.25 1 Ueq.$id.grd Veq.$id.grd Weq.$id.grd $dep
  else
    cp emp.grd Ueq.$id.grd
    cp emp.grd Veq.$id.grd
    cp emp.grd Weq.$id.grd
  fi

  grdmath Ueq.$id.grd U.$id.grd ADD = totU.$id.grd 
  grdmath Veq.$id.grd V.$id.grd ADD = totV.$id.grd 
  grdmath Weq.$id.grd W.$id.grd ADD = totW.$id.grd 

  echo "shift the sampling locations accordingly"
  awk '{print $1+500,$2}' gps.sample | grdtrack -GtotU.$id.grd > tmp
  grdtrack tmp -GtotV.$id.grd | awk '{print $1-500,$2,$3,$4}' > gps.$id.dat

  awk '{print $1+500,$2}' gps.sample | grdtrack -GU.$id.grd > tmp
  grdtrack tmp -GV.$id.grd | awk '{print $1-500,$2,$3,$4}' > gps.secular.$id.dat

  awk '{print $1+500,$2}' gps.sample | grdtrack -GUeq.$id.grd > tmp
  grdtrack tmp -GVeq.$id.grd | awk '{print $1-500,$2,$3,$4}' > gps.eq.$id.dat

  awk '{print $1+500,$2,$3,$4}' los.sample | grdtrack -GtotU.$id.grd > tmp
  grdtrack tmp -GtotV.$id.grd | awk '{print $1-500,$2,$3,$4,$5,$6}' > tmp2
  awk '{print $1,$2,-($3*$5+$4*$6)}' tmp2 > los.$id.dat

  awk '{print $1+500,$2,$3,$4}' los.sample | grdtrack -GU.$id.grd > tmp
  grdtrack tmp -GV.$id.grd | awk '{print $1-500,$2,$3,$4,$5,$6}' > tmp2
  awk '{print $1,$2,-($3*$5+$4*$6)}' tmp2 > los.secular.$id.dat
 
  awk '{print $1+500,$2,$3,$4}' los.sample | grdtrack -GUeq.$id.grd > tmp
  grdtrack tmp -GVeq.$id.grd | awk '{print $1-500,$2,$3,$4,$5,$6}' > tmp2
  awk '{print $1,$2,-($3*$5+$4*$6)}' tmp2 > los.eq.$id.dat
  
done

rm tmp tmp2 
# rm tmp3
#rm U.grd W.grd V.grd
#rm Ueq.grd Weq.grd Veq.grd
#rm totU.grd totV.grd totW.grd
rm gps.sample
rm los.sample
rm header shift name.tmp
