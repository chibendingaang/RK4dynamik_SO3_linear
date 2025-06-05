// process_corrxt.cpp - Improved version
#include <highfive/H5File.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <chrono>
#include "params.hpp"

using namespace HighFive;

// Type aliases to make code more readable
using SpinVector = std::vector<double>;
using SpinLattice = std::vector<SpinVector>;
using Trajectory = std::vector<SpinLattice>;
using CorrelationMatrix = std::vector<std::vector<double>>;

// Function declarations
std::string formatFilename(size_t L, const std::string& alpha_str);
std::string formatGroupname(size_t conf);
std::string formatOutputname(size_t L, size_t conf);
Trajectory applyAlphaWeighting(const Trajectory& spins, double alpha, size_t L);
CorrelationMatrix computeCorrelations(const Trajectory& mconsv, size_t L, 
                                     const std::vector<size_t>& T_array,
                                     size_t fine_res, size_t safe_dist);
void saveCorrelationsToFile(const CorrelationMatrix& Cxt, const std::string& filename,
                           size_t safe_dist);

int main(int argc, char* argv[]) {
    // Parse command line arguments similar to Python
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " begin end fine_res" << std::endl;
        return 1;
    }
    
    const size_t begin = std::stoul(argv[1]);
    const size_t end = std::stoul(argv[2]);
    const size_t fine_res = std::stoul(argv[3]);
    
    // Get parameters from params.hpp
    const size_t L = params::L;
    const std::string dtsymb = params::dtsymb;
    const double dt = params::dt;
    const double Lambda = params::Lambda;
    const double Mu = params::Mu;
    const double alpha = (Lambda - Mu) / (Lambda + Mu);
    const int epss = params::epss;
    
    // Format alpha string similar to Python code
    std::string alpha_deci;
    if (alpha < 1) {
        alpha_deci = "975";
    } else {
        int deci = static_cast<int>(100 * (alpha - static_cast<int>(alpha)));
        alpha_deci = (deci < 10) ? "0" + std::to_string(deci) : std::to_string(deci);
    }
    
    std::string alpha_str = std::to_string(static_cast<int>(alpha)) + "pt" + alpha_deci;
    
    // Format dt string
    std::string dt_str = dtsymb + "emin3";
    
    // Map epsilon values to strings
    std::map<int, std::string> epstr;
    epstr[3] = "emin3";
    epstr[4] = "emin4";
    epstr[6] = "emin6";
    epstr[8] = "emin8";
    epstr[33] = "min3";
    epstr[44] = "min4";
    
    // Create time correlation array based on L
    std::vector<size_t> T_array;
    if (L == 64 || L == 128) {
        for (size_t t = 0; t <= 800; t += 5) {
            T_array.push_back(t);
        }
    } else if (L == 192 || L == 256) {
        for (size_t t = 0; t <= 1200; t += fine_res) {
            T_array.push_back(t);
        }
    } else {
        // Default case
        for (size_t t = 0; t <= 800; t += 5) {
            T_array.push_back(t);
        }
    }
    
    const size_t safe_dist = L/2;
    
    // Define output directory structure
    std::string Cxtpath = "Cxt_series_storage/L" + std::to_string(L);
    std::string hiddensubfolder = "alpha_ne_pm1";
    std::string param = "xpa2b0"; // Equivalent to self.param in Python
    
    // Open HDF5 file
    const std::string filename = "spin_trajectories_L" + std::to_string(L) + "_alpha_" + alpha_str + ".h5";
    
    try {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        File file(filename, File::ReadOnly);
        
        // Process each configuration
        for (size_t conf = begin; conf < end; ++conf) {
            std::cout << "Processing configuration " << conf << "..." << std::endl;
            
            // Get trajectory data from HDF5
            const std::string groupname = "/run_" + std::to_string(conf) + "/trajectory";
            if (!file.exist(groupname)) {
                std::cerr << "Warning: Group " << groupname << " not found\n";
                continue;
            }
            
            // Read spin trajectory
            DataSet dataset = file.getDataSet(groupname);
            Trajectory Sp_a;
            dataset.read(Sp_a);
            std::cout << "Sp_a shape: [" << Sp_a.size() << ", " << Sp_a[0].size() << ", " << Sp_a[0][0].size() << "]" << std::endl;
            
            // Apply alpha weighting
            Trajectory mconsv = applyAlphaWeighting(Sp_a, alpha, L);
            
            // Free original data to save memory
            Trajectory().swap(Sp_a);
            
            // Compute time correlations
            CorrelationMatrix Cxt = computeCorrelations(mconsv, L, T_array, fine_res, safe_dist);
            
            // Save results in format matching Python's path structure
            std::ostringstream outname;
            outname << "./" << Cxtpath << "/" << hiddensubfolder << "/Cxt_L" << L 
                    << "_t_" << dt_str << "_jump" << fine_res << "_" 
                    << epstr[epss] << "_" << param << "_" << alpha_str << "_" 
                    << conf << "to" << (conf + 1) << "config.npy";
            
            saveCorrelationsToFile(Cxt, outname.str(), safe_dist);
            std::cout << "Finished conf " << conf << ", written to " << outname.str() << std::endl;
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        std::cout << "Processing time = " << duration.count() << " seconds" << std::endl;
        
    } catch (const Exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

// Format the HDF5 input filename
std::string formatFilename(size_t L, const std::string& alpha_str) {
    std::ostringstream os;
    os << "spin_trajectories_L" << L << "_alpha_" << alpha_str << ".h5";
    return os.str();
}

// Format the group name for a specific configuration
std::string formatGroupname(size_t conf) {
    std::ostringstream os;
    os << "/run_" << conf << "/trajectory";
    return os.str();
}

// Format the output filename
std::string formatOutputname(size_t L, size_t conf) {
    std::ostringstream os;
    os << "Cxt_L" << L << "_conf_" << conf << ".txt";
    return os.str();
}

// Apply alpha-weighting to spin vectors
Trajectory applyAlphaWeighting(const Trajectory& spins, double alpha, size_t L) {
    size_t steps = spins.size();
    Trajectory mconsv(steps, SpinLattice(L, SpinVector(3)));
    
    #pragma omp parallel for collapse(2) if(steps > 100)
    for (size_t t = 0; t < steps; ++t) {
        for (size_t i = 0; i < L; ++i) {
            double factor = std::pow(alpha, -static_cast<int>(i));
            for (size_t d = 0; d < 3; ++d) {
                mconsv[t][i][d] = spins[t][i][d] * factor;
            }
        }
    }
    
    return mconsv;
}

// Compute the time correlation matrix - matches Python implementation more closely
CorrelationMatrix computeCorrelations(const Trajectory& mconsv, size_t L, 
                                     const std::vector<size_t>& T_array,
                                     size_t fine_res, size_t safe_dist) {
    size_t steps = mconsv.size();
    size_t T = T_array.size();
    // Create matrix for full range: -safe_dist to +safe_dist
    CorrelationMatrix Cxt(T, std::vector<double>(2 * safe_dist + 1, 0.0));
    
    for (size_t ti = 0; ti < T; ++ti) {
        size_t t_shift = T_array[ti];
        
        // Skip if not enough data for this time shift
        if (t_shift >= steps) {
            std::cerr << "Warning: t_shift=" << t_shift << " exceeds available steps=" << steps << std::endl;
            continue;
        }
        
        // Calculate number of time windows (matching Python's approach)
        size_t num_windows = 0;
        if (t_shift == 0) {
            num_windows = steps / fine_res;
        } else {
            num_windows = (steps - t_shift) / fine_res;
        }
        
        if (num_windows == 0) continue;
        
        // Calculate correlation for this time shift
        #pragma omp parallel for if(num_windows > 10)
        for (size_t window = 0; window < num_windows; ++window) {
            // Local accumulator to avoid thread contention
            std::vector<double> local_corr(2 * safe_dist + 1, 0.0);
            
            // Get indices for time slices (matching Python's approach)
            size_t idx1, idx2;
            if (t_shift == 0) {
                idx1 = window * fine_res;
                idx2 = idx1;
            } else {
                idx1 = window * fine_res;
                idx2 = idx1 + t_shift;
            }
            
            // Positive spatial distances (0 to safe_dist)
            for (size_t x = 0; x <= safe_dist; ++x) {
                size_t sites = L - x;  // Number of site pairs for this distance
                for (size_t site = 0; site < sites; ++site) {
                    double dot_product = 0.0;
                    for (size_t d = 0; d < 3; ++d) {
                        dot_product += mconsv[idx1][site][d] * mconsv[idx2][site + x][d];
                    }
                    local_corr[safe_dist + x] += dot_product;
                }
            }
            
            // Negative spatial distances (-safe_dist to -1)
            for (int x_neg = 1; x_neg <= static_cast<int>(safe_dist); ++x_neg) {
                size_t x_abs = static_cast<size_t>(x_neg);
                size_t sites = L - x_abs;  // Number of site pairs for this distance
                for (size_t site = 0; site < sites; ++site) {
                    double dot_product = 0.0;
                    for (size_t d = 0; d < 3; ++d) {
                        dot_product += mconsv[idx1][site + x_abs][d] * mconsv[idx2][site][d];
                    }
                    local_corr[safe_dist - x_neg] += dot_product;
                }
            }
            
            // Combine results
            #pragma omp critical
            {
                for (size_t i = 0; i < 2 * safe_dist + 1; ++i) {
                    Cxt[ti][i] += local_corr[i];
                }
            }
        }
        
        // Normalize correlations
        for (size_t i = 0; i < 2 * safe_dist + 1; ++i) {
            int x = static_cast<int>(i) - static_cast<int>(safe_dist);
            size_t x_abs = std::abs(x);
            size_t normalization = num_windows * (L - x_abs);
            Cxt[ti][i] /= normalization;
        }
    }
    
    return Cxt;
}

// Save correlation matrix to file (similar to NumPy's save format)
void saveCorrelationsToFile(const CorrelationMatrix& Cxt, const std::string& filename,
                           size_t safe_dist) {
    // Create output directory if it doesn't exist
    size_t pos = filename.find_last_of("/\\");
    if (pos != std::string::npos) {
        std::string dir = filename.substr(0, pos);
        // Create directory if it doesn't exist (system-dependent)
        #ifdef _WIN32
        std::string cmd = "mkdir " + dir + " 2> nul";
        #else
        std::string cmd = "mkdir -p " + dir + " 2>/dev/null";
        #endif
        system(cmd.c_str());
    }
    
    std::ofstream fout(filename, std::ios::binary);
    if (!fout) {
        std::cerr << "Error: Cannot open output file " << filename << std::endl;
        return;
    }
    
    // For binary format similar to NumPy's .npy
    // This is a simplification - full NumPy compatibility would require proper headers
    size_t rows = Cxt.size();
    size_t cols = Cxt[0].size();
    
    // Write dimensions
    fout.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    fout.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
    
    // Write data
    for (size_t ti = 0; ti < rows; ++ti) {
        fout.write(reinterpret_cast<const char*>(Cxt[ti].data()), cols * sizeof(double));
    }
    
    // Alternative: save as text file with full precision
    /* 
    std::ofstream fout_txt(filename + ".txt");
    if (fout_txt) {
        fout_txt << std::setprecision(12);
        for (size_t ti = 0; ti < Cxt.size(); ++ti) {
            for (size_t x = 0; x < Cxt[ti].size(); ++x) {
                fout_txt << Cxt[ti][x] << (x == Cxt[ti].size()-1 ? "\n" : " ");
            }
        }
    }
    */
}