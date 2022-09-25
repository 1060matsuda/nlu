#!/bin/bash
#SBATCH --nodes=2             # Number of nodes 
#SBATCH --ntasks-per-node=32   # Number of MPI ranks per node   
#SBATCH --cpus-per-task=4    # number of HW threads per task (equal to OMP_NUM_THREADS*4)8
#SBATCH -A FUSIO_ja6HeBub
#SBATCH -p m100_fua_prod
#SBATCH --qos m100_qos_fuadbg
#SBATCH --mail-type=ALL
#SBATCH --mail-user="matsuda-nayuta775@g.ecc.u-tokyo.ac.jp"
#SBATCH --job-name=vac

module load spectrum_mpi/

#module purge
#module load profile/chem-phys
#module load autoload lammps/29sep2021

mpirun -np 1024 /m100_work/FUSIO_ja5HeBub/lammps/src/lmp_mpi -in in.lammps.detec.mix # 32 MPI tasks, 4 GPUs per  node
