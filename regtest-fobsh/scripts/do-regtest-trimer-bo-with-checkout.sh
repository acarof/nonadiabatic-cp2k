#!/usr/bin/bash

#you need to be inside the reg-test directory and run the script as:
# ./scripts/do-regtest-trimer-bo-with-checkout.sh trimer-bo/(system) here(version)

system=$1
version=$2

#MAIN_PATH=/scratch/sgiannini/CODE_versions/new_adiabatic_prop/flavoured-cptk
MAIN_PATH=../

# prepare build enviromnent
source ${MAIN_PATH}/cp2k/tools/toolchain/install/setup
MODULEPATH=/scratch/grudorff/modulefiles module load automake fftw libint libxc mpich scalapack 

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
#make distclean &> make.log
make -j12 ARCH="local" VERSION="sopt" &> make.log
tail make.log
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
echo "okay"

# store results
#mkdir ${system}/results to comment if it is already present 
var1=$(grep "git:" new/run.log | head -1 | awk '{ print $6 }')
mv new ${system}/results/$(date '+%Y%m%d-%H%M-TRIMER-')$var1


