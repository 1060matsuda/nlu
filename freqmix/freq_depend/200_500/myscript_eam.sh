#!/bin/bash
#SBATCH --nodes=16             # Number of nodes 
#SBATCH --ntasks-per-node=32    # Number of MPI ranks per node   
#SBATCH --cpus-per-task=4
#SBATCH --time 23:59:59         # Walltime, format: HH:MM:SS
#SBATCH -A FUSIO_ja6HeBub
#SBATCH -p m100_fua_prod
#SBATCH --mail-type=ALL
#SBATCH --mail-user="matsuda-nayuta775@g.ecc.u-tokyo.ac.jp"

module load spectrum_mpi/

#module purge
#module load profile/chem-phys
#module load autoload lammps/29sep2021

#mpirun -np 64 /m100_work/FUSIO_ja5HeBub/lammps/src/lmp_mpi -in in_default.lammps # 32 MPI tasks, 4 GPUs per  node
mpirun -np 512 /m100_work/FUSIO_ja5HeBub/lammps/src/lmp_mpi -in in.lammps.detec # 32 MPI tasks, 4 GPUs per node

#mpirun -gpu -np 128  lmp_kokkos_cuda_mpi -k on g 4  -sf kk -in in_default.lammps # 4 MPI tasks, 4 GPUs #error occurs with this line ;_;
