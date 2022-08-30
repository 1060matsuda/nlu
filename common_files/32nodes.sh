#!/bin/bash
#SBATCH --nodes=32             # Number of nodes 
#SBATCH --ntasks-per-node=4    # Number of MPI ranks per node   
#SBATCH --cpus-per-task=4    # number of HW threads per task (equal to OMP_NUM_THREADS*4)8
#SBATCH --time 23:59:59         # Walltime, format: HH:MM:SS
#SBATCH -A FUSIO_ja5HeBub
#SBATCH -p m100_fua_prod
#SBATCH --mail-type=ALL
#SBATCH --mail-user="matsuda-nayuta775@g.ecc.u-tokyo.ac.jp"
#SBATCH --job-bame="32nodes"

module load spectrum_mpi/

#module purge
#module load profile/chem-phys
#module load autoload lammps/29sep2021

mpirun -np 64 /m100_work/FUSIO_ja5HeBub/lammps/src/lmp_mpi -in in_mlmd.lammps # 32 MPI tasks, 4 GPUs per  node
