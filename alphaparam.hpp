// params.hpp
#pragma once

namespace SimulationParams {
    constexpr std::size_t L = 128;               // System size
    constexpr std::size_t total_steps = 100000;   // Total RK4 steps
    constexpr double dt = 0.002;                 // Timestep
    constexpr std::size_t save_interval = 100;   // Save every N steps
    constexpr double alpha = 0.95;                // Alpha parameter
    constexpr unsigned int seed = 42;            // Default random seed

    constexpr std::size_t num_saves = total_steps / save_interval; // Derived
}
