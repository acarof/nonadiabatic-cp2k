#!/usr/bin/bash

system=$1
version=$2

#MAIN_PATH=/scratch/acarof/src/CP2K/flavoured-cptk
MAIN_PATH=/scratch/sgiannini/CODE_versions/adiab_propagation/flavoured-cptk/

# prepare build enviromnent
source ${MAIN_PATH}/cp2k/tools/toolchain/install/setup

# pull last change
cd ..
if [ $version = "local" ]
then
   git pull ${MAIN_PATH}
elif [ $version = "here" ]
then
   echo "do nothing"
else
   git checkout $version
fi
#rsync -azvP ${MAIN_PATH}/cp2k/arch/*  cp2k/arch/

# clean & build
cd cp2k/makefiles/
if [ $version != "here" ]
then
   make distclean &> make.log
fi
make -j12 ARCH="local" VERSION="sopt" &> make.log
cd ../../regtest-fobsh

# run cp2k
cp -r ${system}/inputs new

cd new 
../../cp2k/exe/local/cp2k.sopt --xml
python ../scripts/cp2k_strip_unsupported_keywords.py run.inp cp2k_input.xml > run_correct.inp
../../cp2k/exe/local/cp2k.sopt  run_correct.inp > run.log
rm run-nacv-1.xyz
cd ..

# check against baseline
python scripts/reg_test_TRIMER.py $system current > new/reg_test_result.txt

# store results
var1=$(grep "git:" new/run.log | head -1 | awk '{ print $6 }')
mv new ${system}/results/$(date '+%Y%m%d-%H%M-TRIMER-')$var1


