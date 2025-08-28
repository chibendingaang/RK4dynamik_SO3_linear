#!/bin/bash
#PBS -N alphadynamics_2D_hybrid
#PBS -l select=3:ncpus=32:mpiprocs=32
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o alphadynamics_job.log

# -- Load required modules --
# Note: A parallel-enabled HDF5 library is required.
module load mpich/4.1.2
module load hdf5-parallel

cd $PBS_O_WORKDIR

# -- Count how many ranks were allocated --
NP=$(cat $PBS_NODEFILE | wc -l)

echo "Running job on $NP MPI ranks across $PBS_NUM_NODES nodes"

# --- Launch with mpiexec ---
mpiexec -np $NP ./alphadynamics_2D_hybrid_mpi spin_trajectories_2D_parallel.h5