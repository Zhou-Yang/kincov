#!/bin/bash

#SBATCH --job-name=pka
#SBATCH --partition=gpu01
#SBATCH --gpus=1
#SBATCH --ntasks-per-node=1

module load buildenv-cuda/11.1.0_455.23.05 buildenv-gcc/8.2.0 CMake/3.19.6 MPI/3.4.2 openmpi/4.1.1 buildenv-gcc/8.2.0
source /home/zy/Amber/amber20/amber.sh

for PH in {8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14};do
    pmemd.cuda -O -i prod.implicit.mdin.pH_$PH -p protein.parm7 -c protein.equil.rst7 -x $PH/protein.md1.nc -inf $PH/protein.md1.mdinfo -o $PH/protein.md1.mdout -r $PH/protein.md1.rst7 -cpin protein.equil.cpin -cpout $PH/protein.md1.cpout -cprestrt $PH/protein.md1.cprestrt
done

