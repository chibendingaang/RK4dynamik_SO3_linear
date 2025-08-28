// alphadynamics_mpi.cpp
// MPI-enabled variant of alphadynamics_2.cpp
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <ctime>

#include <highfive/H5File.hpp>

// Include the parameters
#include "alphaparam.hpp"

#ifdef USE_MPI
  #include <mpi.h>
#endif

const double PI = 3.14159265358979323846;

class Vector3D {
public:
    double x, y, z;
    Vector3D(double x_=0.0, double y_=0.0, double z_=0.0)
        : x(x_), y(y_), z(z_) {}
    Vector3D operator+(const Vector3D& rhs) const {
        return Vector3D(x + rhs.x, y + rhs.y, z + rhs.z);
    }
    Vector3D operator-(const Vector3D& rhs) const {
        return Vector3D(x - rhs.x, y - rhs.y, z - rhs.z);
    }
    Vector3D operator*(double alpha) const {
        return Vector3D(alpha * x, alpha * y, alpha * z);
    }
    Vector3D cross(const Vector3D& rhs) const {
        return Vector3D(
            y * rhs.z - z * rhs.y,
            z * rhs.x - x * rhs.z,
            x * rhs.y - y * rhs.x
        );
    }
    void normalize() {
        double norm = std::sqrt(x*x + y*y + z*z);
        if (norm > 0)
            *this = *this * (1.0 / norm);
    }
    void print() const {
        std::cout << "(" << x << ", " << y << ", " << z << ")";
    }
};

class SpinChain {
private:
    std::vector<Vector3D> spins;
    double alpha;
    std::mt19937 rng;

public:
    SpinChain(int L, double alpha_)
        : spins(L), alpha(alpha_), rng(std::random_device{}()) {}

    void randomizeUnitSpins(unsigned int seed = 0) {
        if (seed != 0) rng.seed(seed);
        std::uniform_real_distribution<double> dist_z(-1.0, 1.0);
        std::uniform_real_distribution<double> dist_phi(-PI, PI);

        for (auto& spin : spins) {
            double Sz = dist_z(rng);
            double phi = dist_phi(rng);
            double Sxy = std::sqrt(std::max(0.0, 1.0 - Sz * Sz));
            spin = Vector3D(Sxy * std::cos(phi), Sxy * std::sin(phi), Sz);
        }
    }

    Vector3D rhs(int i) const {
        int L = spins.size();
        if (i == 0)
            return spins[0].cross(spins[1]);
        if (i == L - 1)
            return spins[L - 1].cross(spins[L - 2] * alpha);
        return spins[i].cross(spins[i + 1] + spins[i - 1] * alpha);
    }

    Vector3D rk_rhs(const std::vector<Vector3D>& S, int i) const {
        int L = S.size();
        if (i == 0)
            return S[0].cross(S[1]);
        if (i == L - 1)
            return S[L - 1].cross(S[L - 2] * alpha);
        return S[i].cross(S[i + 1] + S[i - 1] * alpha);
    }

    void rk4Step(double dt) {
        int L = spins.size();
        std::vector<Vector3D> k1(L), k2(L), k3(L), k4(L), temp(L);

        for (int i = 0; i < L; ++i) k1[i] = rhs(i);
        for (int i = 0; i < L; ++i) temp[i] = spins[i] + k1[i] * (0.5 * dt);
        for (int i = 0; i < L; ++i) k2[i] = rk_rhs(temp, i);
        for (int i = 0; i < L; ++i) temp[i] = spins[i] + k2[i] * (0.5 * dt);
        for (int i = 0; i < L; ++i) k3[i] = rk_rhs(temp, i);
        for (int i = 0; i < L; ++i) temp[i] = spins[i] + k3[i] * dt;
        for (int i = 0; i < L; ++i) k4[i] = rk_rhs(temp, i);

        for (int i = 0; i < L; ++i) {
            spins[i] = spins[i] + (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i]) * (dt / 6.0);
            spins[i].normalize();
        }
    }

    void evolve(int steps, double dt) {
        for (int step = 0; step < steps; ++step) {
            rk4Step(dt);
            if (step % 5000 == 0) {
                std::cout << "Step " << step << ":\n"; printSpins(); std::cout << std::endl;
            }
        }
    }

    void printSpins() const {
        for (size_t i = 0; i < spins.size(); ++i) {
            std::cout << "S[" << i + 1 << "] = ";
            spins[i].print();
            std::cout << "\n";
        }
    }

    const std::vector<Vector3D>& getSpins() const { return spins; }
};

// Generate a seed based on global run index and rank to ensure uniqueness
unsigned int generate_unique_seed(int run_index, int rank) {
    std::time_t t = std::time(nullptr);
    // Mix time, run_index and rank with large constants (simple hash)
    unsigned int seed = static_cast<unsigned int>(t);
    seed ^= static_cast<unsigned int>(run_index * 0x9e3779b1u);
    seed ^= static_cast<unsigned int>(rank * 0x85ebca6bu);
    // Additional scramble
    seed = (seed ^ (seed << 13)) ^ (seed >> 17);
    return seed ? seed : 1u; // avoid zero seed if you prefer
}

// save_trajectory: same as before, writes group /run_<run_index> with "seed", "final_spins" and "trajectory"
void save_trajectory(const SpinChain& chain,
                     const std::vector<std::vector<std::vector<double>>>& trajectory_history,
                     int run_index,
                     unsigned int seed,
                     const std::string& filename) {
    try {
        using namespace HighFive;
        // Open or create the file. Each writer opens/closes the file when its turn comes.
        File file(filename, File::ReadWrite | File::Create | File::OpenOrCreate);

        std::string group_name = "/run_" + std::to_string(run_index);
        Group group;
        try {
            group = file.createGroup(group_name);
        } catch (const HighFive::Exception&) {
            group = file.getGroup(group_name);
        }

        // Write seed attribute (overwrite if exists)
        if (group.hasAttribute("seed")) {
            auto attr = group.getAttribute("seed");
            attr.write(seed);
        } else {
            group.createAttribute<unsigned int>("seed", DataSpace::From(seed)).write(seed);
        }

        // Save final state
        std::vector<std::vector<double>> current_state;
        for (const auto& spin : chain.getSpins()) {
            current_state.push_back({spin.x, spin.y, spin.z});
        }

        if (group.exists("final_spins")) {
            group.getDataSet("final_spins").write(current_state);
        } else {
            DataSet current_dataset = group.createDataSet<double>("final_spins", DataSpace::From(current_state));
            current_dataset.write(current_state);
        }

        // Save full trajectory if non-empty
        if (!trajectory_history.empty()) {
            if (group.exists("trajectory")) {
                group.getDataSet("trajectory").write(trajectory_history);
            } else {
                DataSet trajectory_dataset = group.createDataSet<double>("trajectory", DataSpace::From(trajectory_history));
                trajectory_dataset.write(trajectory_history);
            }
        }

        std::cout << "Rank wrote run " << run_index << " to " << filename << std::endl;
    } catch (const HighFive::Exception& e) {
        std::cerr << "Error saving to HDF5 file for run " << run_index << ": " << e.what() << std::endl;
    }
}

int main(int argc, char** argv) {
    // Optional: allow overriding filename via argv
    std::string out_filename = "spin_trajectories_parallel.h5";
    if (argc > 1) out_filename = argv[1];

    #ifdef USE_MPI
    MPI_Init(&argc, &argv);
    int rank = 0, nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    #else
    int rank = 0, nprocs = 1;
    #endif

    using namespace SimulationParams;

    if (rank == 0) {
        std::cout << "Starting MPI-enabled spin simulation\n";
        std::cout << "Requested output file: " << out_filename << "\n";
        std::cout << "Parameters: L=" << L << ", alpha=" << alpha << ", dt=" << dt
                  << ", total_steps=" << total_steps << ", save_interval=" << save_interval
                  << ", num_runs=" << num_runs << ", num_saves (reserve) = " << num_saves << std::endl;
        std::cout << "MPI ranks: " << nprocs << std::endl;
    }

    // --- Distribute runs across ranks: each global run index assigned deterministically ---
    // For each global run index in [0, num_runs), if (run % nprocs == rank) this rank executes it.
    for (int global_run = 0; global_run < num_runs; ++global_run) {
        if ((global_run % nprocs) != rank) continue; // not my run

        unsigned int run_seed = generate_unique_seed(global_run, rank);
        std::cout << "[rank " << rank << "] Starting global run " << global_run << " with seed " << run_seed << std::endl;

        SpinChain chain(L, alpha);
        chain.randomizeUnitSpins(run_seed);

        // record trajectory
        std::vector<std::vector<std::vector<double>>> trajectory_history;
        trajectory_history.reserve(num_saves);

        for (std::size_t step = 0; step < total_steps; ++step) {
            chain.rk4Step(dt);
            if (step % save_interval == 0) {
                std::vector<std::vector<double>> current_state;
                current_state.reserve((size_t)L);
                for (const auto& spin : chain.getSpins()) {
                    current_state.push_back({spin.x, spin.y, spin.z});
                }
                trajectory_history.push_back(std::move(current_state));
            }
        }

        // Save result -- serialized writes across ranks to avoid race conditions.
        #ifdef USE_MPI
        // Wait so that all ranks have finished computation before writing phase begins
        MPI_Barrier(MPI_COMM_WORLD);

        // Serialise writers by rank. Each rank writes all its assigned global_runs when its turn arrives.
        for (int writer = 0; writer < nprocs; ++writer) {
            if (writer == rank) {
                save_trajectory(chain, trajectory_history, global_run, run_seed, out_filename);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        #else
        // Non-MPI: just write
        save_trajectory(chain, trajectory_history, global_run, run_seed, out_filename);
        #endif

        std::cout << "[rank " << rank << "] Completed global run " << global_run << std::endl;
    }

    #ifdef USE_MPI
    if (rank == 0) std::cout << "All ranks completed their assigned runs.\n";
    MPI_Finalize();
    #else
    std::cout << "Single-process run complete.\n";
    #endif

    return 0;
}

