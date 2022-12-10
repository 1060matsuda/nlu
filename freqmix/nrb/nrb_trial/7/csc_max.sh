#!/bin/sh
#SBATCH --account=HEDISINT
#SBATCH --ntasks=10240
#SBATCH --cpus-per-task=2
#SBATCH --partition=L
#SBATCH --time=23:59:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user="matsuda-nayuta775@g.ecc.u-tokyo.ac.jp"
#SBATCH --job-name=nrb_trial

srun /project/HEDISINT/lammps/lammps/src/lmp_mpi -in mod_in.lammps.nrb.mix
