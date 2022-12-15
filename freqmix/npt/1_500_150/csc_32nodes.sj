#!/bin/sh
#SBATCH --account=HEDISINT
#SBATCH --ntasks=1080
#SBATCH --cpus-per-task=2
#SBATCH --partition=L
#SBATCH --mail-type=ALL
#SBATCH --mail-user="matsuda-nayuta775@g.ecc.u-tokyo.ac.jp"
#SBATCH --job-name=unchi

srun /project/HEDISINT/lammps/lammps/src/lmp_mpi -in in.lammps.detec.mix
