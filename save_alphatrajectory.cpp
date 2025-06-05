// save_trajectory.hpp
#pragma once

#include <highfive/H5File.hpp> // HDF% wrapper for C++
#include <vector>
#include <string>

using TrajectoryType = std::vector<std::vector<std::vector<double>>>; // 3D array: (num_saves, L, 3)

void saveTrajectory(const TrajectoryType& trajectory, const std::string& filename, unsigned int seed) {
    using namespace HighFive;

    // Create or open the HDF5 file
    File file(filename, File::Overwrite); // Overwrite existing file if exists

    // Create a group for this run
    Group runGroup = file.createGroup("/run_000"); 

    // Write the seed as an attribute
    runGroup.createAttribute("seed", seed);

    // Create and write the trajectory dataset
    DataSet dataset = runGroup.createDataSet<double>("spins", DataSpace::From(trajectory));
    dataset.write(trajectory);
}

void save_spinchain_to_hdf5(const SpinChain& chain, int run_index, unsigned int seed, const std::string& filename) {
    using namespace HighFive;

    File file(filename, File::ReadWrite | File::Create | File::OpenOrCreate);

    std::string group_name = "/run_" + std::to_string(run_index);
    Group group = file.createGroup(group_name);

    group.createAttribute<unsigned int>("seed", DataSpace::From(seed)).write(seed);

    std::vector<std::vector<double>> spin_data;
    for (const auto& spin : chain.getSpins()) {
        spin_data.push_back({spin.x, spin.y, spin.z});
    }

    DataSet dataset = group.createDataSet<double>("spins", DataSpace::From(spin_data));
    dataset.write(spin_data);
}

#include <iostream>
#include <vector>
#include <random>
#include <cmath>

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

    void randomizeUnitSpins() {
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

    Vector3D rk_rhs(const std::vector<Vector3D>& S, int i) const {
        int L = S.size();

        if (i == 0)
            return S[0].cross(S[1]);

        if (i == L - 1)
            return S[L - 1].cross(S[L - 2] * alpha);

        return S[i].cross(S[i + 1] + S[i - 1] * alpha);
    }

    void printSpins() const {
        for (size_t i = 0; i < spins.size(); ++i) {
            std::cout << "S[" << i + 1 << "] = ";
            spins[i].print();
            std::cout << "\n";
        }
    }
};



SpinChain chain(L, alpha);

for (int run = 0; run < 1000; ++run) {
    unsigned int seed = generate_unique_seed(run); // however you want
    chain.randomizeUnitSpins(seed);
    chain.evolve(time_steps, dt);
    
    save_spinchain_to_hdf5(chain, run, seed, "spin_chains.h5");
}


int main() {
    int L = 10;
    double alpha = 1.0;
    double dt = 0.01;

    SpinChain chain(L, alpha);
    chain.randomizeUnitSpins();

    std::cout << "Initial Spins:\n";
    chain.printSpins();

    for (int step = 0; step < 100; ++step) {
        chain.rk4Step(dt);
    }

    std::cout << "\nSpins after 100 RK4 steps:\n";
    chain.printSpins();

    return 0;
}