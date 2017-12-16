#! /bin/sh

EQs=4000
numProcs=4

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
mpirun -np $numProcs ~/bin/julia --color=yes -O3 /home/nslindquist/cs_research/Lindquist-Neil-Work/julia/power-method/power-method.jl $EQs
#mpirun -np $numProcs ~/bin/julia  --track-allocation=user --color=yes /home/nslindquist/cs_research/Lindquist-Neil-Work/julia/power-method/power-method.jl $EQs


echo
echo
echo ePetra Power Method
mpirun -np $numProcs /home/nslindquist/cs_research/ePetra-PowerMethod/petra_power_method_LL $EQs

read -p "Press enter to continue"


#flag for using tcp
