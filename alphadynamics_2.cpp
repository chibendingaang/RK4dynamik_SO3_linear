// save_alphatrajectory.cpp
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <highfive/H5File.hpp>
#include <ctime>

// Include the parameters
#include "alphaparam.hpp"

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
        // Set the seed if provided
        if (seed != 0) {
            rng.seed(seed);
        }
        
        std::uniform_real_distribution<double> dist_z(-1.0, 1.0);
        std::uniform_real_distribution<double> dist_phi(-PI, PI);

        for (auto& spin : spins) {
            double Sz = dist_z(rng);
            double phi = dist_phi(rng);
            double Sxy = std::sqrt(1.0 - Sz * Sz);
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

        // k1
        for (int i = 0; i < L; ++i)
            k1[i] = rhs(i);

        // temp = S + 0.5*dt*k1
        for (int i = 0; i < L; ++i)
            temp[i] = spins[i] + k1[i] * (0.5 * dt);

        // k2
        for (int i = 0; i < L; ++i)
            k2[i] = rk_rhs(temp, i);

        // temp = S + 0.5*dt*k2
        for (int i = 0; i < L; ++i)
            temp[i] = spins[i] + k2[i] * (0.5 * dt);

        // k3
        for (int i = 0; i < L; ++i)
            k3[i] = rk_rhs(temp, i);

        // temp = S + dt*k3
        for (int i = 0; i < L; ++i)
            temp[i] = spins[i] + k3[i] * dt;

        // k4
        for (int i = 0; i < L; ++i)
            k4[i] = rk_rhs(temp, i);

        // Final update
        for (int i = 0; i < L; ++i) {
            spins[i] = spins[i] + (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i]) * (dt / 6.0);
            spins[i].normalize();
        }
    }

    void evolve(int steps, double dt) {
        for (int step = 0; step < steps; ++step) {
            rk4Step(dt);
            
            // Print spins only every 5000 steps as requested
            if (step % 5000 == 0) {
                std::cout << "Step " << step << ":" << std::endl;
                printSpins();
                std::cout << std::endl;
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

    const std::vector<Vector3D>& getSpins() const {
        return spins;
    }
};

// Function to generate a unique seed based on run_index and current time
unsigned int generate_unique_seed(int run_index) {
    // Get current time
    std::time_t t = std::time(nullptr);
    std::tm* now = std::localtime(&t);
    
    // Create a hash combining run_index and time components
    unsigned int seed = static_cast<unsigned int>(t) + run_index; 
                        //(now->tm_year + 1900)*365 + (now->tm_mon + 1)*12 +  // year*month
                        //(now->tm_mday) *                            // date
                        //((now->tm_hour) * 3600 + (now->tm_min)*60 +            // hr*min
                        //(now->tm_sec));                              // seconds
    
    return seed;
}

// Combined function for saving trajectory data
void save_trajectory(const SpinChain& chain, 
                    const std::vector<std::vector<std::vector<double>>>& trajectory_history,
                    int run_index, unsigned int seed, const std::string& filename) {
    try {
        using namespace HighFive;

        // Open or create the file
        File file(filename, File::ReadWrite | File::Create | File::OpenOrCreate);

        // Create group name
        std::string group_name = "/run_" + std::to_string(run_index);
        
        // Create the group (will throw if it already exists)
        Group group;
        try {
            group = file.createGroup(group_name);
        } catch (const HighFive::Exception& e) {
            // If group exists, get it
            group = file.getGroup(group_name);
        }

        // Write seed as attribute
        group.createAttribute<unsigned int>("seed", DataSpace::From(seed)).write(seed);

        // Save current state (last point in trajectory)
        std::vector<std::vector<double>> current_state;
        for (const auto& spin : chain.getSpins()) {
            current_state.push_back({spin.x, spin.y, spin.z});
        }
        
        // Save final state
        DataSet current_dataset = group.createDataSet<double>("final_spins", DataSpace::From(current_state));
        current_dataset.write(current_state);
        
        // Save full trajectory if available
        if (!trajectory_history.empty()) {
            DataSet trajectory_dataset = group.createDataSet<double>("trajectory", DataSpace::From(trajectory_history));
            trajectory_dataset.write(trajectory_history);
        }
        
        std::cout << "Successfully saved run " << run_index << " to " << filename << std::endl;
    } catch (const HighFive::Exception& e) {
        std::cerr << "Error saving to HDF5 file: " << e.what() << std::endl;
    }
}

int main() {
    using namespace SimulationParams;
    
    std::cout << "Starting simulation with parameters:" << std::endl;
    std::cout << "L = " << L << ", alpha = " << alpha << ", dt = " << dt << std::endl;
    std::cout << "Total steps = " << total_steps << ", save interval = " << save_interval << std::endl;
    
    // For a single run test
    SpinChain chain(L, alpha);
    chain.randomizeUnitSpins(seed);

    std::cout << "Initial spins:" << std::endl;
    chain.printSpins();

    // For recording trajectory
    std::vector<std::vector<std::vector<double>>> trajectory_history;
    trajectory_history.reserve(num_saves);
    
    // Run simulation and save snapshots
    for (std::size_t step = 0; step < total_steps; ++step) {
        chain.rk4Step(dt);
        
        // Print spins only every 5000 steps as requested
        if (step % 5000 == 0) {
            std::cout << "Step " << step << ":" << std::endl;
            chain.printSpins();
            std::cout << std::endl;
        }
        
        // Save trajectory at intervals (still every save_interval)
        if (step % save_interval == 0) {
            std::vector<std::vector<double>> current_state;
            for (const auto& spin : chain.getSpins()) {
                current_state.push_back({spin.x, spin.y, spin.z});
            }
            trajectory_history.push_back(current_state);
        }
    }
    
    std::cout << "Final spins after evolution:" << std::endl;
    chain.printSpins();
    
    // Save the final result
    save_trajectory(chain, trajectory_history, 0, seed, "spin_trajectory_L128_alpha_min_0pt95.h5");
    
    // For multiple runs
    std::cout << "\nRunning multiple simulations..." << std::endl;
    const int num_runs = 1000; // Set this to desired number of runs
    
    for (int run = 0; run < num_runs; ++run) {
        unsigned int run_seed = generate_unique_seed(run);
        std::cout << "Run " << run << " using seed: " << run_seed << std::endl;
        
        SpinChain run_chain(L, alpha);
        run_chain.randomizeUnitSpins(run_seed);
        
        // For recording trajectory of this run
        std::vector<std::vector<std::vector<double>>> run_trajectory_history;
        run_trajectory_history.reserve(num_saves);
        
        // Run simulation and collect trajectory data
        for (std::size_t step = 0; step < total_steps; ++step) {
            run_chain.rk4Step(dt);
            
            // Save trajectory at intervals
            if (step % save_interval == 0) {
                std::vector<std::vector<double>> current_state;
                for (const auto& spin : run_chain.getSpins()) {
                    current_state.push_back({spin.x, spin.y, spin.z});
                }
                run_trajectory_history.push_back(current_state);
            }
        }
        
        // Save trajectory for this run
        save_trajectory(run_chain, run_trajectory_history, run, run_seed, "spin_trajectories_L128_alpha_min_0pt95.h5");
        
        if (run % 10 == 0 || run == num_runs - 1) {
            std::cout << "Completed " << (run + 1) << " of " << num_runs << " runs" << std::endl;
        }
    }
    
    std::cout << "Simulation complete." << std::endl;
    return 0;
}

