#! /bin/sh

#this file must be run from the julia directory of the Lindquist-Neil-Work repository

#allow EQs and numProcs to be passed as arguments
if [ "$#" -eq 1 ]; then
  numProcs=$1
else
  numProcs=20
fi

#echo Git Branch:
#git branch | grep \*
echo Git Commit:
git rev-parse HEAD

echo

echo Evaluating $EQs equations
echo Using $numProcs processors

#(cd ~/cs_research/DistributedArrays.jl/ && git checkout timings-sparse-matrix) > /dev/null

#echo
#echo
#echo "DistributedArrays (mutable) Power Method:"
#julia --color=yes -p $numProcs -O3 power-method/DA-power-method.jl $EQs

#(cd ~/cs_research/DistributedArrays.jl/ && git checkout immutable-sparse) > /dev/null

#echo
#echo
#echo "DistributedArrays (immutable) Power Method:"
#julia --color=yes -p $numProcs -O3 power-method/DA-power-method.jl $EQs

echo
echo
echo JuliaPetra Power Method:
set -x
mpiexec -bind-to core -map-by core -n $numProcs julia-1.0.0 --color=yes --project=. -O3 power-method/power-method.jl
#mpirun -np $numProcs ~/bin/julia  --track-allocation=user --color=yes julia/power-method/power-method.jl $EQs
set +x


echo
echo
echo ePetra Power Method
set -x
mpiexec -bind-to core -map-by core -n $numProcs power-method/ePetra-PowerMethod/petra_power_method_LL
set +x

#read -p "Press enter to continue"


#flag for using tcp
