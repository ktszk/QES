#!/bin/sh
##PJM settings
#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-ss"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=96:00:00"
#GE settings
#$ -N job_name
#$ -pe smp 32
#$ -cwd
#$ -V
#$ -q all.q
#$ -S /bin/bash
#Torque settings
##PBS -j oe
##PBS -o out.o
#PBS -l nodes=1:ppn=16
sw_scf=0
sw_band=0
sw_pdos=0

ncore=16
mat=GaN

cd $PBS_O_WORKDIR
if [ $sw_scf -eq 1 ] ; then
  python QE.py -scf
  mpirun -np $ncore pw.x -nk=4 <$mat.scf>$mat.scf.out
fi
if [ $sw_band -eq 1 ] ; then
  python QE.py -dos -band
  mpirun -np $ncore pw.x -nk=4 <$mat.nscf>$mat.nscf.out
  if [ $sw_pdos -eq 1 ] ; then
    mpirun -np $ncore projwfc.x -nk=4 <$mat.prjwfc>$mat.prjwfc.out
  fi
  dos.x<$mat.dos_in>$mat.dos_in.out
  mpirun -np $ncore pw.x -nk=4 <$mat.bands>$mat.bands.out
  bands.x<$mat.bands_in>$mat.bands_in.out
  plotband.x<eband.plotband>eband.plotband.out
fi
