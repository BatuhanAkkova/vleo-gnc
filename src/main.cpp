/// @file main.cpp
/// Entry point for the GNC simulation.

#include "simulation/monte_carlo.hpp"
#include "simulation/campaign_config.hpp"
#include "utils/data_loader.hpp"
#include <filesystem>
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[]) {
    using namespace gnc;
    namespace fs = std::filesystem;

    std::cout << "=== GNC C++ Simulation ===" << std::endl;

    // ── Load environment data ────────────────────────────────────────
    std::string data_dir = "data";
    if (argc > 1) data_dir = argv[1];

    TruthEnv env;
    AeroCoeffs ac;

    // Attempt to load LUT data (will skip if files are missing)
    bool has_luts = false;
    try {
        std::string grav_path = data_dir + "/egm2008.json";
        if (fs::exists(grav_path)) {
            load_gravity_coefficients(grav_path, 8, env.C_nm, env.S_nm);
            env.init(8);
            std::cout << "[OK] Gravity coefficients loaded" << std::endl;
        } else {
            std::cout << "[WARN] Gravity file not found: " << grav_path << std::endl;
        }

        std::string dens_path = data_dir + "/density.bin";
        if (fs::exists(dens_path)) {
            env.density_lut = load_density_lut(dens_path);
            std::cout << "[OK] Density LUT loaded" << std::endl;
        } else {
            std::cout << "[WARN] Density file not found: " << dens_path << std::endl;
        }

        std::string igrf_path = data_dir + "/igrf.bin";
        if (fs::exists(igrf_path)) {
            env.igrf_lut = load_igrf_lut(igrf_path);
            std::cout << "[OK] IGRF LUT loaded" << std::endl;
        } else {
            std::cout << "[WARN] IGRF file not found: " << igrf_path << std::endl;
        }

        std::string wind_path = data_dir + "/wind.bin";
        if (fs::exists(wind_path)) {
            env.wind_lut = load_wind_lut(wind_path);
            std::cout << "[OK] Wind LUT loaded" << std::endl;
        } else {
            std::cout << "[WARN] Wind file not found: " << wind_path << std::endl;
        }

        std::string cn_path = data_dir + "/cn_coeff.json";
        std::string ct_path = data_dir + "/ct_coeff.json";
        if (fs::exists(cn_path) && fs::exists(ct_path)) {
            ac = load_aero_coeffs(cn_path, ct_path);
            std::cout << "[OK] Aero coefficients loaded" << std::endl;
        } else {
            std::cout << "[WARN] Aero files not found" << std::endl;
        }
        has_luts = true;
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Data loading failed: " << e.what() << std::endl;
        return 1;
    }

    // ── Configure simulation ─────────────────────────────────────────
    SimConfig cfg;
    cfg.dt = 0.1;          // 0.1s for stability
    cfg.t_end = 27000.0;    // 5 orbits to see everything
    
    // Initial estimation errors (large enough to see convergence)
    cfg.mekf_p0 = 0.1;     // ~18 deg
    cfg.ukf_p_pos = 1000.0; // 1km

    // Initial orbit: 300 km circular, equatorial
    double r0 = constants::R_EARTH + 300e3;
    double v0 = std::sqrt(constants::MU_EARTH / r0);
    cfg.x0.setZero();
    cfg.x0(0) = r0;                    // X position [m]
    cfg.x0(4) = v0;                    // Y velocity [m/s]
    cfg.x0(6) = 0; cfg.x0(7) = 0;     // quaternion identity
    cfg.x0(8) = 0; cfg.x0(9) = 1;
    cfg.x0(10) = 0.1;                 // initial angular velocity [rad/s]
    cfg.x0(11) = -0.1;
    cfg.x0(12) = 0.05;

    // ── Monte Carlo Scenarios Matrix ────────────────────────────────
    CampaignGridConfig grid;
    if (argc > 2) {
        std::string arg2 = argv[2];
        if (arg2.size() > 5 && arg2.substr(arg2.size() - 5) == ".json") {
            grid = load_campaign_config(arg2);
        } else {
            grid.runs_per_scenario = std::atoi(argv[2]);
        }
    }

    std::vector<MCScenario> mc_scenarios = build_scenarios_from_grid(grid, cfg);

    std::cout << "\nPreparing " << mc_scenarios.size() << " scenarios";
    if (grid.adaptive_runs) {
        std::cout << " with adaptive batches (initial=" << grid.initial_runs_per_scenario
                  << ", batch=" << grid.runs_batch_size
                  << ", max=" << grid.max_runs_per_scenario << ")";
    } else {
        std::cout << " with " << grid.runs_per_scenario << " runs each";
    }
    std::cout << "..." << std::endl;

    run_monte_carlo_campaign(mc_scenarios, env, ac, grid.base_seed);

    std::cout << "\nAll scenarios complete. Results in ./results/" << std::endl;
    return 0;
}


