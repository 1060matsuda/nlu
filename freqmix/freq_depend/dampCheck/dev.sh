#!/bin/bash
#SBATCH --nodes=2             # Number of nodes 
#SBATCH --ntasks-per-node=32    # Number of MPI ranks per node   
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:2            # Number of requested gpus per node, can vary between 1 and 4
#SBATCH -A FUSIO_ja6HeBub
#SBATCH -p m100_fua_prod
#SBATCH --qos=m100_qos_fuadbg
#SBATCH --mail-type=ALL
#SBATCH --mail-user="matsuda-nayuta775@g.ecc.u-tokyo.ac.jp"
#SBATCH --time 01:59:59         # Walltime, format: HH:MM:SS
#SBATCH --job-name=multidetec

module purge
module load spectrum_mpi/
module load profile/chem-phys
module load autoload lammps/29sep2021

mpirun -np 64 /m100_work/FUSIO_ja5HeBub/lammps/src/lmp_mpi -in in.lammps.detec # 32 MPI tasks, 4 GPUs per  node
#mpirur -np 128 /m100_work/FUSIO_ja5HeBub/lammps/src/lmp_mpi -in in_default.lammps # 32 MPI tasks, 4 GPUs per node
#mpirun -np 64 lmp_gpu -sf gpu -pk gpu 2 -in strain_tutorial_dev.lammps # 32 MPI tasks, 4 GPUs per  node


#mpirun -gpu -np 64  lmp_kokkos_cuda_mpi -k on g 2  -sf kk -in strain_tutorial_dev.lammps # 4 MPI tasks, 4 GPUs #error occurs with this line ;_;

###IMPORTANT##
#lmp_kokkos and lmp_gpu outputs the cuda error. lmp_mpi is recommended.
