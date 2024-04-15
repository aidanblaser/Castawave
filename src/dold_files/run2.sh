#! /bin/bash

# compile the fortran code then execute it with the ./ option



# cd /Users/nick/Documents/Research/Oceanography/BreakingPacket/numerics/Rev$


 /usr/local/bin/gfortran -o ww_write ww_write.f
./ww_write

#source /opt/intel/bin/compilervars.sh intel64


 /usr/local/bin/gfortran -framework accelerate -O2 -std=legacy -Wno-implicit-interface -o dold dold.f
./dold


