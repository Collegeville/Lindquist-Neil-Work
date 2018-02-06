#! /bin/sh

#this file must be run from the julia directory of the Lindquist-Neil-Work repository

#allow EQs and numProcs to be passed as arguments
if [ "$#" -eq 2 ]; then
  EQs=$1
  numProcs=$2
else
  EQs=16000000
  numProcs=16
fi

#echo Git Branch:
#git branch | grep \*
echo Git Commit:
git rev-parse HEAD

echo

echo Evaluating $EQs equations
echo Using $numProcs processors

echo
echo
echo JuliaPetra Power Method:
mpirun -np $numProcs ~/bin/julia --color=yes -O3 power-method/power-method.jl $EQs
#mpirun -np $numProcs ~/bin/julia  --track-allocation=user --color=yes julia/power-method/power-method.jl $EQs


echo
echo
echo ePetra Power Method
mpirun -np $numProcs power-method/ePetra-PowerMethod/petra_power_method_LL $EQs

#read -p "Press enter to continue"


#flag for using tcp
