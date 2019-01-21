#!/bin/bash




# get u,v,w data from -75 to -25 latitude
#mkdir ofesdata_ACC


for i in {3205..3206}

do
  export out=`printf "%05d\n" $i`
  echo "getting files at time $i"
  ncks -v 'uvel' -d time,$i -d lat,0,500 -d lon,600,1500 http://apdrc.soest.hawaii.edu:80/dods/public_ofes/OfES/ncep_0.1_global_3day/uvel -o ofesdata_ACC_lat75-25_lon60-150/uvel${out}.nc
  #ncks -v 'vvel' -d time,$i -d lat,0,500 -d lon,600,1500 http://apdrc.soest.hawaii.edu:80/dods/public_ofes/OfES/ncep_0.1_global_3day/vvel -o ofesdata_ACC/vvel${out}.nc
  #ncks -v 'wvel' -d time,$i -d lat,0,500 -d lon,600,1500 http://apdrc.soest.hawaii.edu:80/dods/public_ofes/OfES/ncep_0.1_global_3day/wvel -o ofesdata_ACC/wvel${out}.nc
  #ncks -v 'temp' -d time,$i -d lat,0,500 -d lon,600,1500 http://apdrc.soest.hawaii.edu:80/dods/public_ofes/OfES/ncep_0.1_global_3day/temp -o ofesdata_ACC/temp${out}.nc  
done


