#! /bin/sh

cd /home/nslindquist/cs_research/Lindquist-Neil-Work/hpcg/

arch=MPI_Intel_OMP_ZFP

if [ $# != 0 ]; then
	if [ $# != 1 ]; then
		arch=$2
	fi
	sed 3s/.\*/$1\ $1\ $1/ -i $arch/bin/hpcg.dat
fi

cd $arch/bin
xhpcg
cd ../..
