#!/usr/bin/bash

# prepare build enviromnent
source /scratch/acarof/src/CP2K/flavoured-cptk/cp2k/tools/toolchain/install/setup

# pull last change
cd flavoured-cptk
git pull /scratch/acarof/src/CP2K/flavoured-cptk/
rsync -azvP /scratch/acarof/src/CP2K/flavoured-cptk/cp2k/arch/*  cp2k/arch/

# clean & build
cd cp2k/makefiles/
#make distclean &> make.log
make -j12 ARCH="local" VERSION="sopt" &> make.log
cd ../../..

# run cp2k
cp -r trimer-bo new
cd new 
../flavoured-cptk/cp2k/exe/local/cp2k.sopt  run_TRIMER.inp > run_TRIMER.log
rm run-nacv-1.xyz
cd ..

# check against baseline
python scripts/reg_test_TRIMER.py > new/reg_test_result.txt

# store results
var1=$(grep "git:" new/run_TRIMER.log | head -1 | awk '{ print $6 }')
mv new $(date '+%Y%m%d-%H%M-TRIMER-')$var1


