// alphadynamics_2D_hybrid_mpi.cpp
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <mpi.h>
#include <highfive/H5File.hpp>

#include "alphaparam.hpp"

// A simple logger for keeping track of what's going on.
class Logger {
public:
    static void open(const std::string& filename) {
        log_file.open(filename, std::ios_base::app);
    }

    static void log(const std::string& message) {
        if (log_file.is_open()) {
            std::time_t t = std::time(nullptr);
            log_file << "[" << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S") << "] "
                     << message << std::endl;
        }
    }
private:
    static std::ofstream log_file;
};
std::ofstream Logger::log_file;


// Our trusty 3D vector class. No changes here.
class Vector3D {
public:
    double x, y, z;
    Vector3D(double x_=0.0, double y_=0.0, double z_=0.0) : x(x_), y(y_), z(z_) {}
    Vector3D operator+(const Vector3D& rhs) const { return Vector3D(x + rhs.x, y + rhs.y, z + rhs.z); }
    Vector3D operator-(const Vector3D& rhs) const { return Vector3D(x - rhs.x, y - rhs.y, z - rhs.z); }
    Vector3D operator*(double scalar) const { return Vector3D(scalar * x, scalar * y, scalar * z); }
    Vector3D cross(const Vector3D& rhs) const { return Vector3D(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x); }
    void normalize() { double norm = std::sqrt(x*x + y*y + z*z); if (norm > 0) *this = *this * (1.0 / norm); }
};

// This class handles one entire parallel simulation.
// It knows about its own little piece of the grid and how to talk to its neighbors.
class SpinLattice {
private:
    std::vector<Vector3D> local_spins; // Our local patch of spins, plus a border for halos.
    int local_nx, local_ny;            // The size of our actual patch (e.g., L/4).
    int global_Lx, global_Ly;          // The size of the full lattice across all processes.

    // MPI info for our little group of 16 processes.
    MPI_Comm sim_comm; // The private communicator for our team.
    int rank, size;    // Our personal ID (0-15) and team size (16).
    int coords[2];     // Our (x,y) position in the 4x4 grid of processes.
    int neighbors[4];  // Who is above, below, left, and right of us.

    // A helper to find a spin in our local grid, accounting for the halo border.
    int get_local_idx(int ix, int iy) const {
        return (ix + 1) * (local_ny + 2) + (iy + 1);
    }

public:
    SpinLattice(int Lx, int Ly, MPI_Comm communicator) : global_Lx(Lx), global_Ly(Ly), sim_comm(communicator) {
        MPI_Comm_rank(sim_comm, &rank);
        MPI_Comm_size(sim_comm, &size);

        if (size != 16) {
            if (rank == 0) Logger::log("Oops! A simulation group needs exactly 16 processes to form a 4x4 grid.");
            MPI_Abort(sim_comm, 1);
        }
        
        // Let's arrange our 16 processes into a nice 4x4 grid.
        int dims[2] = {4, 4};
        int periods[2] = {1, 1}; // The grid wraps around (periodic boundaries).
        MPI_Comm cart_comm;
        MPI_Cart_create(sim_comm, 2, dims, periods, 1, &cart_comm);

        // Now, find out where we are and who our neighbors are.
        MPI_Cart_coords(cart_comm, rank, 2, coords);
        MPI_Cart_shift(cart_comm, 0, 1, &neighbors[2], &neighbors[3]); // Left, Right
        MPI_Cart_shift(cart_comm, 1, 1, &neighbors[0], &neighbors[1]); // Up, Down

        // Figure out the dimensions of our personal patch of the lattice.
        local_nx = Lx / dims[0];
        local_ny = Ly / dims[1];

        // Allocate memory for our patch plus a 1-cell border for the halos.
        local_spins.resize((local_nx + 2) * (local_ny + 2));
    }

    void randomizeUnitSpins(unsigned int seed) {
        // We add our rank to the seed to make sure every process gets a unique random sequence.
        std::mt19937 rng(seed + rank);
        std::uniform_real_distribution<double> dist_z(-1.0, 1.0);
        std::uniform_real_distribution<double> dist_phi(-M_PI, M_PI);

        for (int i = 0; i < local_nx; ++i) {
            for (int j = 0; j < local_ny; ++j) {
                double Sz = dist_z(rng);
                double phi = dist_phi(rng);
                double Sxy = std::sqrt(1.0 - Sz * Sz);
                local_spins[get_local_idx(i, j)] = Vector3D(Sxy * std::cos(phi), Sxy * std::sin(phi), Sz);
            }
        }
    }
    
    // This is where the magic happens. We swap boundary info with our neighbors.
    void halo_exchange(std::vector<Vector3D>& spins_to_exchange) {
        // To send a column of data, we need to tell MPI its shape.
        MPI_Datatype col_type;
        MPI_Type_vector(local_ny, 3, (local_ny + 2) * 3, MPI_DOUBLE, &col_type);
        MPI_Type_commit(&col_type);
        
        double* send_buf = reinterpret_cast<double*>(spins_to_exchange.data());
        
        // Swap rows with our up and down neighbors.
        MPI_Sendrecv(&send_buf[get_local_idx(0, 0)*3], local_ny * 3, MPI_DOUBLE, neighbors[0], 0, 
                     &send_buf[get_local_idx(local_nx, 0)*3], local_ny * 3, MPI_DOUBLE, neighbors[1], 0, sim_comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&send_buf[get_local_idx(local_nx-1, 0)*3], local_ny * 3, MPI_DOUBLE, neighbors[1], 1, 
                     &send_buf[get_local_idx(-1, 0)*3], local_ny * 3, MPI_DOUBLE, neighbors[0], 1, sim_comm, MPI_STATUS_IGNORE);

        // Swap columns with our left and right neighbors using our custom column type.
        MPI_Sendrecv(&send_buf[get_local_idx(0, 0)*3], 1, col_type, neighbors[2], 2, 
                     &send_buf[get_local_idx(0, local_ny)*3], 1, col_type, neighbors[3], 2, sim_comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&send_buf[get_local_idx(0, local_ny-1)*3], 1, col_type, neighbors[3], 3,
                     &send_buf[get_local_idx(0, -1)*3], 1, col_type, neighbors[2], 3, sim_comm, MPI_STATUS_IGNORE);

        MPI_Type_free(&col_type);
    }
    
    // Calculate the time derivative for a single spin using the 2D formula.
    Vector3D rk_rhs(const std::vector<Vector3D>& S, int ix, int iy) const {
        const Vector3D& S_ij = S[get_local_idx(ix, iy)];
        const Vector3D& S_up = S[get_local_idx(ix, iy + 1)];
        const Vector3D& S_down = S[get_local_idx(ix, iy - 1)];
        const Vector3D& S_right = S[get_local_idx(ix + 1, iy)];
        const Vector3D& S_left = S[get_local_idx(ix - 1, iy)];

        Vector3D H_eff = S_up + S_down + S_right - S_left;
        return S_ij.cross(H_eff);
    }

    void rk4Step(double dt) {
        std::vector<Vector3D> k1(local_spins.size()), k2(local_spins.size()), k3(local_spins.size()), k4(local_spins.size()), temp(local_spins.size());

        // For each stage of the RK4 integrator, we need fresh boundary data.
        halo_exchange(local_spins);
        for (int i = 0; i < local_nx; ++i) for (int j = 0; j < local_ny; ++j) k1[get_local_idx(i, j)] = rk_rhs(local_spins, i, j);
        for (size_t i=0; i < local_spins.size(); ++i) temp[i] = local_spins[i] + k1[i] * (0.5 * dt);

        halo_exchange(temp);
        for (int i = 0; i < local_nx; ++i) for (int j = 0; j < local_ny; ++j) k2[get_local_idx(i, j)] = rk_rhs(temp, i, j);
        for (size_t i=0; i < local_spins.size(); ++i) temp[i] = local_spins[i] + k2[i] * (0.5 * dt);
        
        halo_exchange(temp);
        for (int i = 0; i < local_nx; ++i) for (int j = 0; j < local_ny; ++j) k3[get_local_idx(i, j)] = rk_rhs(temp, i, j);
        for (size_t i=0; i < local_spins.size(); ++i) temp[i] = local_spins[i] + k3[i] * dt;
        
        halo_exchange(temp);
        for (int i = 0; i < local_nx; ++i) for (int j = 0; j < local_ny; ++j) k4[get_local_idx(i, j)] = rk_rhs(temp, i, j);
        
        // And now for the final update...
        for (int i = 0; i < local_nx; ++i) {
            for (int j = 0; j < local_ny; ++j) {
                int idx = get_local_idx(i, j);
                local_spins[idx] = local_spins[idx] + (k1[idx] + k2[idx] * 2.0 + k3[idx] * 2.0 + k4[idx]) * (dt / 6.0);
                local_spins[idx].normalize();
            }
        }
    }
    
    // Gathers up the spin data from our local patch, ready for saving.
    std::vector<std::vector<double>> get_current_state() {
        std::vector<std::vector<double>> current_state;
        current_state.reserve(local_nx * local_ny);
        for (int i = 0; i < local_nx; ++i) {
            for (int j = 0; j < local_ny; ++j) {
                const auto& spin = local_spins[get_local_idx(i, j)];
                current_state.push_back({spin.x, spin.y, spin.z});
            }
        }
        return current_state;
    }

    // This function saves our little patch of the simulation into one big HDF5 file.
    void save_trajectory(const std::vector<std::vector<std::vector<double>>>& trajectory_history,
                         unsigned int seed, int run_index, const std::string& filename) {
        using namespace HighFive;
        
        // This tells HDF5 that we're a team and will be writing to the file together.
        FileAccessProps fapl;
        fapl.add(MPIOFileAccess<MPI_COMM_WORLD, MPI_INFO_NULL>{});
        File file(filename, File::ReadWrite | File::Create | File::OpenOrCreate, fapl);

        std::string group_name = "/run_" + std::to_string(run_index);
        Group group = file.exist(group_name) ? file.getGroup(group_name) : file.createGroup(group_name);
        
        // The first process in our group (rank 0) is in charge of writing the general info.
        if (rank == 0) {
            group.createAttribute<unsigned int>("seed", DataSpace::From(seed)).write(seed);
            group.createAttribute<int>("Lx", DataSpace::From(global_Lx)).write(global_Lx);
            group.createAttribute<int>("Ly", DataSpace::From(global_Ly)).write(global_Ly);
        }

        size_t num_saves = trajectory_history.size();
        DataSpace filespace({num_saves, (size_t)global_Lx, (size_t)global_Ly, 3});
        DataSet dataset = group.createDataSet<double>("trajectory", filespace);

        // For each saved timestep, we write our local data to the correct spot in the big array.
        for (size_t t = 0; t < num_saves; ++t) {
            std::vector<double> local_data_flat;
            local_data_flat.reserve(local_nx * local_ny * 3);
            for(const auto& spin_vec : trajectory_history[t]) {
                local_data_flat.insert(local_data_flat.end(), spin_vec.begin(), spin_vec.end());
            }
            
            // We tell HDF5 exactly which slice of the big final array this process is responsible for.
            DataSpace memspace({(size_t)local_nx, (size_t)local_ny, 3});
            dataset.select({t, (size_t)coords[0] * local_nx, (size_t)coords[1] * local_ny, 0}, 
                           {(size_t)1, (size_t)local_nx, (size_t)local_ny, 3}).write_raw(local_data_flat.data(), memspace);
        }
    }
};

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    std::string out_filename = "spin_trajectories_parallel.h5";
    if (argc > 1) out_filename = argv[1];

    if (world_rank == 0) {
        Logger::open("simulation.log");
        Logger::log("Job kicked off with " + std::to_string(world_size) + " total MPI processes.");
    }
    
    // We split all our processes into teams of 16 since the nodes typically have 16/32/48/64 cores for usage
    const int processes_per_sim = 16;
    if (world_size % processes_per_sim != 0) {
        if (world_rank == 0) Logger::log("ERROR: Total number of processes should be a multiple of 16.");
        MPI_Finalize();
        return 1;
    }

    // 'color' will be the team ID. All processes with the same color are on the same team.
    int color = world_rank / processes_per_sim;
    MPI_Comm sim_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &sim_comm);

    int sim_rank, sim_size;
    MPI_Comm_rank(sim_comm, &sim_rank);
    MPI_Comm_size(sim_comm, &sim_size);
    
    // The team ID (color) also serves as the unique run number for this simulation.
    int global_run_index = color; 
    
    using namespace SimulationParams;

    if (global_run_index < total_runs) {
        // Generate a unique seed. Multiplying by world_rank makes it unique for everyone.
        unsigned int run_seed = static_cast<unsigned int>(time(NULL)) * (world_rank + 1);

        // The leader of each team (rank 0) logs the start of their run.
        if (sim_rank == 0) {
            std::stringstream ss;
            ss << "Team " << color << " (run " << global_run_index << ") is a go! We have " << sim_size 
               << " processes on our side. Seed: " << run_seed;
            Logger::log(ss.str());
        }

        SpinLattice lattice(Lx, Ly, sim_comm);
        lattice.randomizeUnitSpins(run_seed);

        std::vector<std::vector<std::vector<double>>> trajectory_history;
        trajectory_history.reserve(num_saves);
        
        for (std::size_t step = 0; step < total_steps; ++step) {
            lattice.rk4Step(dt);
            if (step % save_interval == 0) {
                trajectory_history.push_back(lattice.get_current_state());
            }
        }
        
        // We wait until every process is done computing before we start writing the file.
        MPI_Barrier(MPI_COMM_WORLD);
        if (sim_rank == 0) Logger::log("Run " + std::to_string(global_run_index) + " finished calculations. Now saving data...");
        
        lattice.save_trajectory(trajectory_history, run_seed, global_run_index, out_filename);

        if (sim_rank == 0) Logger::log("Run " + std::to_string(global_run_index) + " is all done and saved!");
    }

    MPI_Comm_free(&sim_comm);
    if (world_rank == 0) Logger::log("All teams have finished. Job complete.");
    MPI_Finalize();
    return 0;
}