#pragma once
/// @file monte_carlo.hpp
/// OpenMP-parallelized Monte Carlo runner with adaptive batching support.

#include "simulation/simulation.hpp"
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <omp.h>

namespace gnc {

struct MCConfig {
    SimConfig sim;
    int n_runs = 10;
    std::string output_dir = "mc_results";
};

struct MCScenario {
    SimConfig sim;
    double rate_deg_s = std::numeric_limits<double>::quiet_NaN();
    double rho_scale = 1.0;
    double actuator_scale = 1.0;
    double stress_torque_nm = 0.0;
    double onboard_density_error_pct = 0.0;
    int onboard_gravity_degree_truncation = 2;
    std::string output_dir;
    int n_runs = 10;
    std::string name;
    bool adaptive_runs = false;
    int initial_runs = 10;
    int batch_runs = 10;
    int max_runs = 60;
    int bootstrap_resamples = 400;
    double success_ci_half_width = 0.05;
    double recovery_ci_half_width_s = 60.0;
};

namespace detail {

struct RunArtifact {
    int run_idx = 0;
    unsigned seed = 0;
    std::string output_path;
    SimulationResult sim_result{67};
};

struct BootstrapInterval {
    double low = std::numeric_limits<double>::quiet_NaN();
    double high = std::numeric_limits<double>::quiet_NaN();

    double half_width() const {
        if (!std::isfinite(low) || !std::isfinite(high)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return 0.5 * (high - low);
    }
};

inline TruthEnv clone_env(const TruthEnv& env, double rho_scale) {
    TruthEnv local_env;
    local_env.sh_degree = env.sh_degree;
    local_env.C_nm = env.C_nm;
    local_env.S_nm = env.S_nm;
    local_env.sh_buf = SHBuf(env.sh_degree);
    local_env.density_lut = env.density_lut;
    local_env.igrf_lut = env.igrf_lut;
    local_env.wind_lut = env.wind_lut;
    local_env.rho_scale = rho_scale;
    return local_env;
}

inline void write_run_sidecar(const MCScenario& sc, const RunArtifact& artifact) {
    const auto& sim = sc.sim;
    const auto& result = artifact.sim_result;
    nlohmann::json meta = {
        {"scenario", sc.name},
        {"run_idx", artifact.run_idx},
        {"seed", artifact.seed},
        {"rate_deg_s", sc.rate_deg_s},
        {"rho_scale", sc.rho_scale},
        {"actuator_scale", sc.actuator_scale},
        {"stress_torque_nm", sc.stress_torque_nm},
        {"onboard_density_error_pct", sc.onboard_density_error_pct},
        {"onboard_gravity_degree_truncation", sc.onboard_gravity_degree_truncation},
        {"detumble_control_dt", sim.detumble_control_dt},
        {"bdot_k", sim.bdot_k},
        {"transition_control_dt", sim.transition_control_dt},
        {"nominal_control_dt", sim.nominal_control_dt},
        {"enable_early_termination", sim.enable_early_termination},
        {"acquisition_threshold_deg", sim.acquisition_threshold_deg},
        {"nominal_entry_pointing_deg", sim.nominal_entry_pointing_deg},
        {"nominal_entry_consistency_deg", sim.nominal_entry_consistency_deg},
        {"nominal_entry_rate_deg", sim.nominal_entry_rate_deg},
        {"nominal_exit_pointing_deg", sim.nominal_exit_pointing_deg},
        {"nominal_exit_consistency_deg", sim.nominal_exit_consistency_deg},
        {"nominal_exit_rate_deg", sim.nominal_exit_rate_deg},
        {"required_acquisition_time_s", sim.required_acquisition_time_s},
        {"required_slew_capture_time_s", sim.required_slew_capture_time_s},
        {"required_nominal_time_s", sim.required_nominal_time_s},
        {"required_nominal_exit_time_s", sim.required_nominal_exit_time_s},
        {"acquisition_attitude_gain", sim.acquisition_attitude_gain},
        {"acquisition_rate_gain", sim.acquisition_rate_gain},
        {"rate_capture_rate_gain", sim.rate_capture_rate_gain},
        {"slew_capture_ref_rate_deg", sim.slew_capture_ref_rate_deg},
        {"acquisition_lqr_scale", sim.acquisition_lqr_scale},
        {"nominal_attitude_gain", sim.nominal_attitude_gain},
        {"nominal_rate_gain", sim.nominal_rate_gain},
        {"nominal_integral_gain", sim.nominal_integral_gain},
        {"nominal_integral_limit_deg_s", sim.nominal_integral_limit_deg_s},
        {"acquisition_max_linear_err_deg", sim.acquisition_max_linear_err_deg},
        {"acquisition_desat_trigger_frac", sim.acquisition_desat_trigger_frac},
        {"acquisition_desat_rate_gate_deg", sim.acquisition_desat_rate_gate_deg},
        {"acquisition_exit_rate_deg", sim.acquisition_exit_rate_deg},
        {"slew_capture_entry_rate_deg", sim.slew_capture_entry_rate_deg},
        {"slew_capture_exit_rate_deg", sim.slew_capture_exit_rate_deg},
        {"precision_sensor_hold_s", sim.precision_sensor_hold_s},
        {"direct_detumble_to_fine_pointing", sim.direct_detumble_to_fine_pointing},
        {"single_acquisition_mode", sim.single_acquisition_mode},
        {"rw_max_speed", sim.rw_max_speed},
        {"max_controller_degraded_time_s", sim.max_controller_degraded_time_s},
        {"convergence_metric_version", 2},
        {"settle_pointing_threshold_deg", sim.settle_pointing_threshold_deg},
        {"settle_rate_threshold_deg_s", sim.settle_rate_threshold_deg_s},
        {"settle_pointing_rms_threshold_deg", sim.settle_pointing_rms_threshold_deg},
        {"settle_rate_rms_threshold_deg_s", sim.settle_rate_rms_threshold_deg_s},
        {"settled_success", result.settled_success},
        {"early_terminated", result.early_terminated},
        {"stop_time_s", result.stop_time_s},
        {"settle_start_s", result.settle_start_s},
        {"settle_complete_s", result.settle_complete_s},
        {"settle_pointing_rms_deg", result.settle_pointing_rms_deg},
        {"settle_pointing_max_deg", result.settle_pointing_max_deg},
        {"settle_rate_rms_deg_s", result.settle_rate_rms_deg_s},
        {"settle_rate_max_deg_s", result.settle_rate_max_deg_s},
        {"control_updates", result.control_updates},
        {"onboard_density_error_pct", sim.onboard_density_error_pct},
        {"onboard_gravity_degree_truncation", sim.onboard_gravity_degree_truncation},
        {"stress_test", sim.stress_test},
        {"stress_start_dynamic", sim.stress_start_dynamic},
        {"stress_start_s", sim.stress_start_s},
        {"stress_duration_s", sim.stress_duration_s},
        {"rows", result.log.n_rows()},
        {"columns", result.log.n_cols}
    };
    std::ofstream(artifact.output_path + ".json") << meta.dump(2);
}

inline std::vector<RunArtifact> run_scenario_batch(const MCScenario& sc,
                                                   const TruthEnv& env,
                                                   const AeroCoeffs& ac,
                                                   unsigned base_seed,
                                                   int start_run_idx,
                                                   int run_count) {
    std::vector<RunArtifact> artifacts(run_count);
    const int max_threads = omp_get_max_threads();

    #pragma omp parallel for schedule(dynamic)
    for (int batch_idx = 0; batch_idx < run_count; ++batch_idx) {
        Eigen::setNbThreads(1);
        const int run_idx = start_run_idx + batch_idx;
        const unsigned seed = base_seed + static_cast<unsigned>(run_idx);

        TruthEnv local_env = clone_env(env, sc.rho_scale);
        SimulationResult sim_result = run_simulation(sc.sim, local_env, ac, seed);

        RunArtifact artifact;
        artifact.run_idx = run_idx;
        artifact.seed = seed;
        artifact.output_path = sc.output_dir + "/run_" + std::to_string(run_idx) + ".bin";
        artifact.sim_result = std::move(sim_result);
        const std::string output_path = artifact.output_path;
        const int n_rows = artifact.sim_result.log.n_rows();
        artifact.sim_result.log.save(output_path);
        write_run_sidecar(sc, artifact);
        artifacts[batch_idx] = std::move(artifact);

        #pragma omp critical
        {
            std::cout << "[MC] [" << sc.name << "] Run " << (run_idx + 1)
                      << " -> " << output_path
                      << " (" << n_rows << " steps, "
                      << max_threads << " threads available)" << std::endl;
        }
    }

    return artifacts;
}

inline BootstrapInterval bootstrap_mean_interval(const std::vector<double>& samples,
                                                 int resamples,
                                                 unsigned seed) {
    BootstrapInterval out;
    if (samples.empty()) return out;
    if (samples.size() == 1) {
        out.low = samples.front();
        out.high = samples.front();
        return out;
    }

    std::mt19937 rng(seed);
    std::uniform_int_distribution<size_t> pick(0, samples.size() - 1);
    std::vector<double> means;
    means.reserve(resamples);
    for (int i = 0; i < resamples; ++i) {
        double accum = 0.0;
        for (size_t j = 0; j < samples.size(); ++j) {
            accum += samples[pick(rng)];
        }
        means.push_back(accum / static_cast<double>(samples.size()));
    }
    std::sort(means.begin(), means.end());
    const size_t lo_idx = static_cast<size_t>(0.025 * (means.size() - 1));
    const size_t hi_idx = static_cast<size_t>(0.975 * (means.size() - 1));
    out.low = means[lo_idx];
    out.high = means[hi_idx];
    return out;
}

inline bool scenario_has_converged(const MCScenario& sc,
                                   const std::vector<RunArtifact>& artifacts) {
    if (!sc.adaptive_runs) {
        return static_cast<int>(artifacts.size()) >= sc.n_runs;
    }
    if (static_cast<int>(artifacts.size()) < sc.initial_runs) {
        return false;
    }

    std::vector<double> success_samples;
    std::vector<double> settle_times;
    success_samples.reserve(artifacts.size());
    for (const auto& artifact : artifacts) {
        success_samples.push_back(artifact.sim_result.settled_success ? 1.0 : 0.0);
        if (artifact.sim_result.settled_success && std::isfinite(artifact.sim_result.settle_complete_s)) {
            settle_times.push_back(artifact.sim_result.settle_complete_s);
        }
    }

    const BootstrapInterval success_ci = bootstrap_mean_interval(
        success_samples, sc.bootstrap_resamples, 7001u + static_cast<unsigned>(artifacts.size()));
    const double success_hw = success_ci.half_width();
    const bool success_tight = std::isfinite(success_hw) && success_hw <= sc.success_ci_half_width;

    bool recovery_tight = true;
    if (settle_times.size() >= 2) {
        const BootstrapInterval recovery_ci = bootstrap_mean_interval(
            settle_times, sc.bootstrap_resamples, 9001u + static_cast<unsigned>(artifacts.size()));
        const double recovery_hw = recovery_ci.half_width();
        recovery_tight = std::isfinite(recovery_hw) && recovery_hw <= sc.recovery_ci_half_width_s;
    }

    return success_tight && recovery_tight;
}

inline void write_scenario_adaptive_summary(const MCScenario& sc,
                                            const std::vector<RunArtifact>& artifacts) {
    if (!sc.adaptive_runs) return;

    std::vector<double> success_samples;
    std::vector<double> settle_times;
    for (const auto& artifact : artifacts) {
        success_samples.push_back(artifact.sim_result.settled_success ? 1.0 : 0.0);
        if (artifact.sim_result.settled_success && std::isfinite(artifact.sim_result.settle_complete_s)) {
            settle_times.push_back(artifact.sim_result.settle_complete_s);
        }
    }

    const BootstrapInterval success_ci = bootstrap_mean_interval(
        success_samples, sc.bootstrap_resamples, 7001u + static_cast<unsigned>(artifacts.size()));
    const BootstrapInterval recovery_ci = bootstrap_mean_interval(
        settle_times, sc.bootstrap_resamples, 9001u + static_cast<unsigned>(artifacts.size()));

    nlohmann::json summary = {
        {"scenario", sc.name},
        {"adaptive_runs", sc.adaptive_runs},
        {"completed_runs", artifacts.size()},
        {"max_runs", sc.max_runs},
        {"success_ci_low", success_ci.low},
        {"success_ci_high", success_ci.high},
        {"success_ci_half_width", success_ci.half_width()},
        {"recovery_ci_low_s", recovery_ci.low},
        {"recovery_ci_high_s", recovery_ci.high},
        {"recovery_ci_half_width_s", recovery_ci.half_width()}
    };
    std::ofstream(sc.output_dir + "/scenario_adaptive_summary.json") << summary.dump(2);
}

} // namespace detail

/// Run Monte Carlo campaign. Each run gets seed = base_seed + i.
inline void run_monte_carlo(const MCConfig& mc, TruthEnv& env,
                            const AeroCoeffs& ac, unsigned base_seed = 42) {
    MCScenario sc;
    sc.sim = mc.sim;
    sc.rho_scale = env.rho_scale;
    sc.output_dir = mc.output_dir;
    sc.n_runs = mc.n_runs;
    sc.name = "single_scenario";
    std::filesystem::create_directories(sc.output_dir);
    (void)detail::run_scenario_batch(sc, env, ac, base_seed, 0, sc.n_runs);
}

/// Run multiple scenarios. Fixed-run campaigns use run-level parallel batches;
/// adaptive campaigns add seed batches until the bootstrap CI is tight enough.
inline void run_monte_carlo_campaign(const std::vector<MCScenario>& scenarios,
                                     const TruthEnv& env,
                                     const AeroCoeffs& ac,
                                     unsigned base_seed = 42) {
    std::cout << "[MC] Starting campaign with " << scenarios.size()
              << " scenarios using up to " << omp_get_max_threads()
              << " OpenMP threads." << std::endl;

    for (const auto& sc : scenarios) {
        std::filesystem::create_directories(sc.output_dir);
        std::vector<detail::RunArtifact> artifacts;
        int next_run_idx = 0;
        const int target_runs = sc.adaptive_runs ? sc.max_runs : sc.n_runs;

        while (next_run_idx < target_runs) {
            int batch_size = sc.adaptive_runs
                ? (next_run_idx == 0 ? sc.initial_runs : sc.batch_runs)
                : (target_runs - next_run_idx);
            batch_size = std::max(1, std::min(batch_size, target_runs - next_run_idx));

            auto batch = detail::run_scenario_batch(sc, env, ac, base_seed, next_run_idx, batch_size);
            artifacts.insert(artifacts.end(),
                             std::make_move_iterator(batch.begin()),
                             std::make_move_iterator(batch.end()));
            next_run_idx += batch_size;

            if (detail::scenario_has_converged(sc, artifacts)) {
                if (sc.adaptive_runs) {
                    std::cout << "[MC] [" << sc.name << "] adaptive stop after "
                              << artifacts.size() << " runs." << std::endl;
                } else {
                    std::cout << "[MC] [" << sc.name << "] completed "
                              << artifacts.size() << " fixed runs." << std::endl;
                }
                break;
            }
        }

        detail::write_scenario_adaptive_summary(sc, artifacts);
    }
}

} // namespace gnc









