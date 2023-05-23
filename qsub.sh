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
sw_band=1
sw_dos=0
sw_pdos=0
sw_wan=0
sw_post_wan=0
sw_ph=0
sw_post_ph=0
sw_epw=0
sw_opt=0
sw_md=0
sw_cp=0

#nthreads*ncore=num_of_smp_core
ncore=16
nthreads=1
npool=4
nbnd=0
mat=Na

export OMP_NUM_THREADS=$nthreads
if [ $nbnd -gt 0 ] ; then
  npband=-bgrp $nbnd
fi

#cd $PBS_O_WORKDIR
if [ $sw_scf -eq 1 ] ; then
  python QE.py -scf
  mpirun -np $ncore pw.x -npool $npool $npband <$mat.scf>$mat.scf.out
fi

if [ $sw_dos -eq 1 ] ; then
  python QE.py -dos
  mpirun -np $ncore pw.x -npool $npool $npband <$mat.nscf>$mat.nscf.out
  if [ $sw_pdos -eq 1 ] ; then
    mpirun -np $ncore projwfc.x -npool $npool $npband <$mat.prjwfc>$mat.prjwfc.out
  fi
  dos.x<$mat.dos_in>$mat.dos_in.out
fi

if [ $sw_band -eq 1 ] ; then
  python QE.py -band
  mpirun -np $ncore pw.x -npool $npool $npband <$mat.bands>$mat.bands.out
  bands.x<$mat.bands_in>$mat.bands_in.out
  plotband.x<eband.plotband>eband.plotband.out
fi

if [ $sw_wan -eq 1 ] ; then
  #python QE.py -wan
  mpirun -np $ncore pw.x -npool $npool $npband <$mat.nscf>$mat.nscf.out
  wannier90.x -pp $mat
  mpirun -np $ncore pw2wannier90.x<$mat.pw2wan>$mat.pw2wan.out
  wannier90.x $mat
fi

if [ $sw_post_wan -eq 1 ] ; then
  mpirun -np $ncore postw90.x $mat
fi

if [ $sw_ph -eq 1 ] ; then
  python QE.py -ph
  mpirun -np $ncore ph.x<$mat.ph>$mat.ph.out
fi

if [ $sw_post_ph -eq 1 ] ; then
  mpirun -np $ncore q2r.x<$mat.q2r>$mat.q2r.out
  mpirun -np $ncore matdyn.x<$mat.matdyn>$mat.matdyn.out
  mpirun -np $ncore matdyn.x<$mat.freq>$mat.freq.out
fi

if [ $sw_epw -eq 1 ] ; then
  mpirun -np $ncore pw.x -npool $ncore  <$mat.nscf>$mat.nscf.out
  mpirun -np $ncore epw.x -npool $ncore  <$mat.epw>$mat.epw.out
fi

if [ $sw_opt -eq 1 ] ; then
  python QE.py -opt
  mpirun -np $ncore pw.x -npool $npool <$mat.scf>$mat.scf.out
fi

if [ $sw_md -eq 1 ] ; then
  python QE.py -md
  mpirun -np $ncore pw.x -npool $npool <$mat.md>$mat.md.out
fi

if [ $sw_cp -eq 1 ] ; then
  python QE.py -cp
  mpirun -np $ncore cp.x -npool $npool <$mat.cp>$mat.cp.out
  cppp.x <$mat.cppp>$mat.cppp.out
fi
