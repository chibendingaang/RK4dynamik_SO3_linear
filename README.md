# Classical Heisenberg Spin Dynamics Simulations

This project provides a C++ framework for simulating the dynamics of classical Heisenberg spin systems using a 4th-order Runge-Kutta integrator and the HDF5 data format. It begins with a serial code for a **1D spin chain** and evolves to an embarrassingly parallel MPI version for distributing multiple independent 1D simulations. The framework culminates in a high-performance simulation for a **2D spin lattice** that employs a sophisticated hybrid MPI model. This advanced version uses domain decomposition with halo exchange to parallelize a single lattice simulation across a dedicated group of processes, while simultaneously running multiple of these parallel simulations within a single job for high-throughput computing. 
Note: a previous version in Fortran is also provided as much of the simulations in initial phase was carried out in Fortran. HDF5 library is simpler to use with C++.

---

## Equations of Motion

The simulations model the time evolution of classical spin vectors, $\vec{S}$, according to the Landau-Lifshitz equation.

### 1D Chain (Generalized)
The dynamics for the 1D chain are governed by the generalized equation of motion with an asymmetric exchange parameter, $\alpha$:

$$\frac{d\vec{S}_i}{dt} = \vec{S}_i \times (\vec{S}_{i+1} + \alpha \cdot \vec{S}_{i-1})$$

### 2D Lattice (Anisotropic)
For the 2D scenario, the interaction is anisotropic. The dynamics correspond to setting $\alpha = +1$ for vertical interactions and $\alpha = -1$ for horizontal interactions. The effective field on a spin $\vec{S}_{i,j}$ is the sum of contributions from its four neighbors, resulting in the equation:

$$\frac{d\vec{S}_{i,j}}{dt} = \vec{S}_{i,j} \times \left( (\vec{S}_{i,j+1} + \vec{S}_{i,j-1}) + (\vec{S}_{i+1,j} - \vec{S}_{i-1,j}) \right)$$

---

## Code Versions

This repository contains three main simulation codes representing an evolution in parallel computing strategy:

1.  `alphadynamics_2.cpp`: A serial implementation for the 1D spin chain. Good for testing and small-scale runs.
2.  `alphadynamics2_mpi.cpp`: An embarrassingly parallel version for the 1D chain. It uses MPI to run many independent 1D simulations simultaneously.
3.  `alphadynamics_2D_hybrid_mpi.cpp`: The most advanced version for the 2D lattice. It uses a hybrid MPI model to parallelize a single simulation via domain decomposition while also running multiple of these parallel simulations concurrently. This is the primary code for large-scale production runs.

---

## Getting Started

This guide focuses on the main 2D parallel code, `alphadynamics_2D_hybrid_mpi.cpp`.

### Dependencies
* A C++17 compatible compiler
* An MPI library (e.g., MPICH, OpenMPI)
* A **parallel-enabled** HDF5 library
* The [HighFive](https://github.com/BlueBrain/HighFive) header-only library for HDF5

### Compilation
A `Makefile` is provided for convenience.

1.  **Edit the `Makefile`**: You must update the `HDF5_PATH` and `HIGHFIVE_PATH` variables to point to the correct locations on your system.
    ```makefile
    # Compiler: Use the MPI C++ wrapper
    CXX = mpic++
    CXXFLAGS = -std-c++17 -O3

    # --- Paths ---
    # IMPORTANT: You must update HDF5_PATH to point to your PARALLEL HDF5 installation.
    # The HighFive path should point to the header-only library location.
    HDF5_PATH = /path/to/parallel-hdf5
    HIGHFIVE_PATH = $(HOME)/HighFive

    # Include paths
    INCLUDES = -I$(HDF5_PATH)/include -I$(HIGHFIVE_PATH)/include

    # Libraries and search paths for PARALLEL HDF5
    LIBS = -L$(HDF5_PATH)/lib -lhdf5 -lhdf5_hl

    # --- Files ---
    # Source and target executable
    SRCS = alphadynamics_2D_hybrid_mpi.cpp
    TARGET = alphadynamics_2D_mpi

    # --- Rules ---
    all: $(TARGET)

    $(TARGET): $(SRCS)
        $(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(SRCS) $(LIBS)

    clean:
        rm -f $(TARGET)
    ```
2.  **Compile the code** by running `make` in your terminal:
    ```bash
    make
    ```

### Execution on an HPC Cluster
The code is designed to be run on an HPC cluster using a scheduler like PBS.

1.  **Edit `alphaparam.hpp`**: Configure your simulation parameters (lattice size `Lx`, `Ly`, number of runs `total_runs`, etc.).
2.  **Prepare the PBS Script**: Use the provided PBS script. It requests 3 nodes with 32 cores each (for a total of 96 processes) and launches the simulation. The C++ code will automatically divide these 96 processes into 6 groups of 16 to run 6 concurrent simulations.
    ```bash
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
    mpiexec -np $NP ./alphadynamics_2D_mpi spin_trajectories_2D_parallel.h5
    ```
3.  **Submit the Job**:
    ```bash
    qsub your_script_name.pbs
    ```
The results, including the main log file (`alphadynamics_job.log`) and the HDF5 data file (`spin_trajectories_2D_parallel.h5`), will be saved in the submission directory.