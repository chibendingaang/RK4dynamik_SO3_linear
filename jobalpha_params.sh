#!/bin/bash
#PBS -N alphadynamics_192_min1pt20          
#PBS -l nodes=node1:ppn=32+node2:ppn=32+node3:ppn=32         # 3 nodes, 32 cores per node = 96 ranks
#PBS -l walltime=24:00:00      
#PBS -l mem=32gb              
#PBS -j oe                     # Join stdout and stderr
#PBS -o alphadynamics_job.log      

# Load MPICH -- any other modules which may be needed--
module load mpich/4.1.2   

cd $PBS_O_WORKDIR

# -- Count how many ranks were allocated --
NP=$(cat $PBS_NODEFILE | wc -l)

echo "Running job on $NP MPI ranks across $PBS_NUM_NODES nodes"

# --- Launch with mpiexec (MPICH's launcher) ---
mpiexec -np $NP -machinefile $PBS_NODEFILE ./alphadynamics_mpi spin_trajectories_parallel.h5

