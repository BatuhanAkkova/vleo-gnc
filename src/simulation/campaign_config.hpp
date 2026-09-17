#pragma once
/// @file campaign_config.hpp
/// JSON campaign configuration for Monte Carlo scenario grids.

#include "simulation/monte_carlo.hpp"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace gnc {

struct CampaignGridConfig {
    std::string campaign_name = "campaign";
    std::vector<double> deployment_rates_deg_s = {30.0};
    std::vector<double> density_scales = {10.0};
    std::vector<double> actuator_scales = {1.0};
    std::vector<double> stress_torque_nm_values = {0.0};
    std::vector<double> onboard_density_error_pct = {0.0};
    std::vector<int> onboard_gravity_degree_truncation = {2};
    int runs_per_scenario = 5;
    unsigned base_seed = 42;
    std::string output_root = "results";
    double dt = -1.0;
    double t_end = -1.0;
    double detumble_control_dt = -1.0;
    double transition_control_dt = -1.0;
    double nominal_control_dt = -1.0;
    bool enable_early_termination = false;
    double settle_pointing_threshold_deg = 0.5;
    double settle_rate_threshold_deg_s = 0.25;
    double settle_pointing_rms_threshold_deg = 0.25;
    double settle_rate_rms_threshold_deg_s = 0.1;
    double settle_hold_time_s = 300.0;
    double settle_confirmation_time_s = 60.0;
    bool adaptive_runs = false;
    int initial_runs_per_scenario = 10;
    int runs_batch_size = 10;
    int max_runs_per_scenario = 60;
    int bootstrap_resamples = 400;
    double success_ci_half_width = 0.05;
    double recovery_time_ci_half_width_s = 60.0;
    bool stress_test = false;
    double stress_torque_nm = 0.0;
    double stress_start_s = 0.0;
    bool stress_start_dynamic = false;
    double stress_duration_s = 0.0;
    double bdot_k = -1.0;
    double acquisition_attitude_gain = -1.0;
    double acquisition_rate_gain = -1.0;
    double rate_capture_rate_gain = -1.0;
    double slew_capture_ref_rate_deg = -1.0;
    double acquisition_lqr_scale = -1.0;
    double nominal_attitude_gain = -1.0;
    double nominal_rate_gain = -1.0;
    double nominal_integral_gain = -1.0;
    double nominal_integral_limit_deg_s = -1.0;
    double acquisition_max_linear_err_deg = -1.0;
    double acquisition_desat_trigger_frac = -1.0;
    double acquisition_desat_rate_gate_deg = -1.0;
    double acquisition_exit_rate_deg = -1.0;
    double slew_capture_entry_rate_deg = -1.0;
    double slew_capture_exit_rate_deg = -1.0;
    double required_slew_capture_time_s = -1.0;
    double nominal_entry_pointing_deg = -1.0;
    double nominal_entry_consistency_deg = -1.0;
    double nominal_entry_rate_deg = -1.0;
    double required_nominal_exit_time_s = -1.0;
    double nominal_exit_pointing_deg = -1.0;
    double nominal_exit_consistency_deg = -1.0;
    double precision_sensor_hold_s = -1.0;
    bool direct_detumble_to_fine_pointing = false;
    bool single_acquisition_mode = false;
};

inline std::vector<double> json_double_list(const nlohmann::json& j,
                                            const char* key,
                                            const std::vector<double>& fallback) {
    if (!j.contains(key)) return fallback;
    if (!j.at(key).is_array()) throw std::invalid_argument(std::string(key) + " must be an array");
    return j.at(key).get<std::vector<double>>();
}

inline std::vector<int> json_int_list(const nlohmann::json& j,
                                      const char* key,
                                      const std::vector<int>& fallback) {
    if (!j.contains(key)) return fallback;
    if (!j.at(key).is_array()) throw std::invalid_argument(std::string(key) + " must be an array");
    return j.at(key).get<std::vector<int>>();
}

inline std::vector<double> json_double_axis(const nlohmann::json& j,
                                            const char* key,
                                            double scalar_fallback) {
    if (!j.contains(key)) return {scalar_fallback};
    if (j.at(key).is_array()) return j.at(key).get<std::vector<double>>();
    if (j.at(key).is_null()) return {scalar_fallback};
    return {j.at(key).get<double>()};
}

inline CampaignGridConfig load_campaign_config(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open campaign config: " + path);
    nlohmann::json j;
    f >> j;

    CampaignGridConfig cfg;
    cfg.campaign_name = j.value("campaign_name", cfg.campaign_name);    cfg.deployment_rates_deg_s = json_double_list(j, "deployment_rates_deg_s", cfg.deployment_rates_deg_s);
    cfg.density_scales = json_double_list(j, "density_scales", cfg.density_scales);
    cfg.actuator_scales = json_double_list(j, "actuator_scales", cfg.actuator_scales);
    cfg.onboard_density_error_pct = json_double_list(j, "onboard_density_error_pct", cfg.onboard_density_error_pct);
    cfg.onboard_gravity_degree_truncation = json_int_list(j, "onboard_gravity_degree_truncation", cfg.onboard_gravity_degree_truncation);
    cfg.runs_per_scenario = j.value("runs_per_scenario", cfg.runs_per_scenario);
    cfg.base_seed = j.value("base_seed", cfg.base_seed);
    cfg.output_root = j.value("output_root", cfg.output_root);
    cfg.dt = j.value("dt", cfg.dt);
    cfg.t_end = j.value("t_end", cfg.t_end);
    cfg.detumble_control_dt = j.value("detumble_control_dt", cfg.detumble_control_dt);
    cfg.transition_control_dt = j.value("transition_control_dt", cfg.transition_control_dt);
    cfg.nominal_control_dt = j.value("nominal_control_dt", cfg.nominal_control_dt);
    cfg.enable_early_termination = j.value("enable_early_termination", cfg.enable_early_termination);
    cfg.settle_pointing_threshold_deg = j.value("settle_pointing_threshold_deg", cfg.settle_pointing_threshold_deg);
    cfg.settle_rate_threshold_deg_s = j.value("settle_rate_threshold_deg_s", cfg.settle_rate_threshold_deg_s);
    cfg.settle_pointing_rms_threshold_deg = j.value("settle_pointing_rms_threshold_deg", cfg.settle_pointing_rms_threshold_deg);
    cfg.settle_rate_rms_threshold_deg_s = j.value("settle_rate_rms_threshold_deg_s", cfg.settle_rate_rms_threshold_deg_s);
    cfg.settle_hold_time_s = j.value("settle_hold_time_s", cfg.settle_hold_time_s);
    cfg.settle_confirmation_time_s = j.value("settle_confirmation_time_s", cfg.settle_confirmation_time_s);
    cfg.adaptive_runs = j.value("adaptive_runs", cfg.adaptive_runs);
    cfg.initial_runs_per_scenario = j.value("initial_runs_per_scenario", cfg.initial_runs_per_scenario);
    cfg.runs_batch_size = j.value("runs_batch_size", cfg.runs_batch_size);
    cfg.max_runs_per_scenario = j.value("max_runs_per_scenario", cfg.max_runs_per_scenario);
    cfg.bootstrap_resamples = j.value("bootstrap_resamples", cfg.bootstrap_resamples);
    cfg.success_ci_half_width = j.value("success_ci_half_width", cfg.success_ci_half_width);
    cfg.recovery_time_ci_half_width_s = j.value("recovery_time_ci_half_width_s", cfg.recovery_time_ci_half_width_s);
    cfg.stress_test = j.value("stress_test", cfg.stress_test);
    cfg.stress_torque_nm_values = json_double_axis(j, "stress_torque_nm", cfg.stress_torque_nm);
    if (j.contains("stress_start_s") && j.at("stress_start_s").is_null()) {
        cfg.stress_start_dynamic = true;
    } else {
        cfg.stress_start_s = j.value("stress_start_s", cfg.stress_start_s);
    }
    cfg.stress_duration_s = j.value("stress_duration_s", cfg.stress_duration_s);
    cfg.bdot_k = j.value("bdot_k", cfg.bdot_k);
    cfg.acquisition_attitude_gain = j.value("acquisition_attitude_gain", cfg.acquisition_attitude_gain);
    cfg.acquisition_rate_gain = j.value("acquisition_rate_gain", cfg.acquisition_rate_gain);
    cfg.rate_capture_rate_gain = j.value("rate_capture_rate_gain", cfg.rate_capture_rate_gain);
    cfg.slew_capture_ref_rate_deg = j.value("slew_capture_ref_rate_deg", cfg.slew_capture_ref_rate_deg);
    cfg.acquisition_lqr_scale = j.value("acquisition_lqr_scale", cfg.acquisition_lqr_scale);
    cfg.nominal_attitude_gain = j.value("nominal_attitude_gain", cfg.nominal_attitude_gain);
    cfg.nominal_rate_gain = j.value("nominal_rate_gain", cfg.nominal_rate_gain);
    cfg.nominal_integral_gain = j.value("nominal_integral_gain", cfg.nominal_integral_gain);
    cfg.nominal_integral_limit_deg_s = j.value("nominal_integral_limit_deg_s", cfg.nominal_integral_limit_deg_s);
    cfg.acquisition_max_linear_err_deg = j.value("acquisition_max_linear_err_deg", cfg.acquisition_max_linear_err_deg);
    cfg.acquisition_desat_trigger_frac = j.value("acquisition_desat_trigger_frac", cfg.acquisition_desat_trigger_frac);
    cfg.acquisition_desat_rate_gate_deg = j.value("acquisition_desat_rate_gate_deg", cfg.acquisition_desat_rate_gate_deg);
    cfg.acquisition_exit_rate_deg = j.value("acquisition_exit_rate_deg", cfg.acquisition_exit_rate_deg);
    cfg.slew_capture_entry_rate_deg = j.value("slew_capture_entry_rate_deg", cfg.slew_capture_entry_rate_deg);
    cfg.slew_capture_exit_rate_deg = j.value("slew_capture_exit_rate_deg", cfg.slew_capture_exit_rate_deg);
    cfg.required_slew_capture_time_s = j.value("required_slew_capture_time_s", cfg.required_slew_capture_time_s);
    cfg.nominal_entry_pointing_deg = j.value("nominal_entry_pointing_deg", cfg.nominal_entry_pointing_deg);
    cfg.nominal_entry_consistency_deg = j.value("nominal_entry_consistency_deg", cfg.nominal_entry_consistency_deg);
    cfg.nominal_entry_rate_deg = j.value("nominal_entry_rate_deg", cfg.nominal_entry_rate_deg);
    cfg.required_nominal_exit_time_s = j.value("required_nominal_exit_time_s", cfg.required_nominal_exit_time_s);
    cfg.nominal_exit_pointing_deg = j.value("nominal_exit_pointing_deg", cfg.nominal_exit_pointing_deg);
    cfg.nominal_exit_consistency_deg = j.value("nominal_exit_consistency_deg", cfg.nominal_exit_consistency_deg);
    cfg.precision_sensor_hold_s = j.value("precision_sensor_hold_s", cfg.precision_sensor_hold_s);
    cfg.direct_detumble_to_fine_pointing = j.value("direct_detumble_to_fine_pointing", cfg.direct_detumble_to_fine_pointing);
    cfg.single_acquisition_mode = j.value("single_acquisition_mode", cfg.single_acquisition_mode);    if (cfg.deployment_rates_deg_s.empty()) throw std::invalid_argument("deployment_rates_deg_s cannot be empty");
    if (cfg.density_scales.empty()) throw std::invalid_argument("density_scales cannot be empty");
    if (cfg.actuator_scales.empty()) throw std::invalid_argument("actuator_scales cannot be empty");
    if (cfg.stress_torque_nm_values.empty()) throw std::invalid_argument("stress_torque_nm cannot be empty");
    if (cfg.onboard_density_error_pct.empty()) throw std::invalid_argument("onboard_density_error_pct cannot be empty");
    if (cfg.onboard_gravity_degree_truncation.empty()) throw std::invalid_argument("onboard_gravity_degree_truncation cannot be empty");
    if (cfg.runs_per_scenario <= 0) throw std::invalid_argument("runs_per_scenario must be positive");
    if (cfg.adaptive_runs) {
        if (cfg.initial_runs_per_scenario <= 0) throw std::invalid_argument("initial_runs_per_scenario must be positive");
        if (cfg.runs_batch_size <= 0) throw std::invalid_argument("runs_batch_size must be positive");
        if (cfg.max_runs_per_scenario < cfg.initial_runs_per_scenario) throw std::invalid_argument("max_runs_per_scenario must be >= initial_runs_per_scenario");
        if (cfg.bootstrap_resamples <= 0) throw std::invalid_argument("bootstrap_resamples must be positive");
    }
    return cfg;
}

inline std::vector<MCScenario> build_scenarios_from_grid(const CampaignGridConfig& grid,
                                                         const SimConfig& base) {
    std::vector<MCScenario> scenarios;
        for (double rate_deg_s : grid.deployment_rates_deg_s) {
            for (double density : grid.density_scales) {
                for (double actuator_scale : grid.actuator_scales) {
                    for (double stress_torque : grid.stress_torque_nm_values) {
                        for (double density_error : grid.onboard_density_error_pct) {
                            for (int gravity_degree : grid.onboard_gravity_degree_truncation) {
                                SimConfig sim = base;
                                if (grid.dt > 0.0) sim.dt = grid.dt;
                                if (grid.t_end > 0.0) sim.t_end = grid.t_end;
                                sim.detumble_control_dt = grid.detumble_control_dt;
                                sim.transition_control_dt = grid.transition_control_dt;
                                sim.nominal_control_dt = grid.nominal_control_dt;
                                sim.enable_early_termination = grid.enable_early_termination;
                                sim.settle_pointing_threshold_deg = grid.settle_pointing_threshold_deg;
                                sim.settle_rate_threshold_deg_s = grid.settle_rate_threshold_deg_s;
                                sim.settle_pointing_rms_threshold_deg = grid.settle_pointing_rms_threshold_deg;
                                sim.settle_rate_rms_threshold_deg_s = grid.settle_rate_rms_threshold_deg_s;
                                sim.settle_hold_time_s = grid.settle_hold_time_s;
                                sim.settle_confirmation_time_s = grid.settle_confirmation_time_s;
                                sim.rw_max_torque = base.rw_max_torque * actuator_scale;
                                sim.bdot_max = base.bdot_max * actuator_scale;
                                sim.desat_max = base.desat_max * actuator_scale;
                                sim.stress_test = grid.stress_test;
                                sim.stress_torque_nm = stress_torque;
                                sim.stress_start_s = grid.stress_start_s;
                                sim.stress_start_dynamic = grid.stress_start_dynamic;
                                sim.stress_duration_s = grid.stress_duration_s;
                                sim.onboard_density_error_pct = density_error;
                                sim.onboard_gravity_degree_truncation = gravity_degree;
                                if (grid.bdot_k > 0.0) sim.bdot_k = grid.bdot_k;
                                if (grid.acquisition_attitude_gain > 0.0) sim.acquisition_attitude_gain = grid.acquisition_attitude_gain;
                                if (grid.acquisition_rate_gain > 0.0) sim.acquisition_rate_gain = grid.acquisition_rate_gain;
                                if (grid.rate_capture_rate_gain > 0.0) sim.rate_capture_rate_gain = grid.rate_capture_rate_gain;
                                if (grid.slew_capture_ref_rate_deg > 0.0) sim.slew_capture_ref_rate_deg = grid.slew_capture_ref_rate_deg;
                                if (grid.acquisition_lqr_scale > 0.0) sim.acquisition_lqr_scale = grid.acquisition_lqr_scale;
                                if (grid.nominal_attitude_gain > 0.0) sim.nominal_attitude_gain = grid.nominal_attitude_gain;
                                if (grid.nominal_rate_gain > 0.0) sim.nominal_rate_gain = grid.nominal_rate_gain;
                                if (grid.nominal_integral_gain > 0.0) sim.nominal_integral_gain = grid.nominal_integral_gain;
                                if (grid.nominal_integral_limit_deg_s > 0.0) sim.nominal_integral_limit_deg_s = grid.nominal_integral_limit_deg_s;
                                if (grid.acquisition_max_linear_err_deg > 0.0) sim.acquisition_max_linear_err_deg = grid.acquisition_max_linear_err_deg;
                                if (grid.acquisition_desat_trigger_frac > 0.0) sim.acquisition_desat_trigger_frac = grid.acquisition_desat_trigger_frac;
                                if (grid.acquisition_desat_rate_gate_deg > 0.0) sim.acquisition_desat_rate_gate_deg = grid.acquisition_desat_rate_gate_deg;
                                if (grid.acquisition_exit_rate_deg > 0.0) sim.acquisition_exit_rate_deg = grid.acquisition_exit_rate_deg;
                                if (grid.slew_capture_entry_rate_deg > 0.0) sim.slew_capture_entry_rate_deg = grid.slew_capture_entry_rate_deg;
                                if (grid.slew_capture_exit_rate_deg > 0.0) sim.slew_capture_exit_rate_deg = grid.slew_capture_exit_rate_deg;
                                if (grid.required_slew_capture_time_s > 0.0) sim.required_slew_capture_time_s = grid.required_slew_capture_time_s;
                                if (grid.nominal_entry_pointing_deg > 0.0) sim.nominal_entry_pointing_deg = grid.nominal_entry_pointing_deg;
                                if (grid.nominal_entry_consistency_deg > 0.0) sim.nominal_entry_consistency_deg = grid.nominal_entry_consistency_deg;
                                if (grid.nominal_entry_rate_deg > 0.0) sim.nominal_entry_rate_deg = grid.nominal_entry_rate_deg;
                                if (grid.required_nominal_exit_time_s > 0.0) sim.required_nominal_exit_time_s = grid.required_nominal_exit_time_s;
                                if (grid.nominal_exit_pointing_deg > 0.0) sim.nominal_exit_pointing_deg = grid.nominal_exit_pointing_deg;
                                if (grid.nominal_exit_consistency_deg > 0.0) sim.nominal_exit_consistency_deg = grid.nominal_exit_consistency_deg;
                                if (grid.precision_sensor_hold_s > 0.0) sim.precision_sensor_hold_s = grid.precision_sensor_hold_s;
                                sim.direct_detumble_to_fine_pointing = grid.direct_detumble_to_fine_pointing;
                                sim.single_acquisition_mode = grid.single_acquisition_mode;

                                const double rate_rad_s = rate_deg_s * 3.14159265358979323846 / 180.0;
                                sim.x0.segment<3>(10) = Vec3(1.0, -1.0, 0.5).normalized() * rate_rad_s;

                                std::ostringstream name;
                                name << "rate" << static_cast<int>(rate_deg_s)
                                     << "_rho" << density
                                     << "_act" << actuator_scale;
                                if (grid.stress_test) name << "_stress" << stress_torque;
                                if (density_error != 0.0) name << "_densErr" << density_error;
                                if (gravity_degree != 2) name << "_grav" << gravity_degree;

                                MCScenario sc;
                                sc.sim = sim;
                                sc.rate_deg_s = rate_deg_s;
                                sc.rho_scale = density;
                                sc.actuator_scale = actuator_scale;
                                sc.stress_torque_nm = stress_torque;
                                sc.onboard_density_error_pct = density_error;
                                sc.onboard_gravity_degree_truncation = gravity_degree;
                                sc.n_runs = grid.runs_per_scenario;
                                sc.name = name.str();
                                sc.output_dir = grid.output_root + "/" + sc.name;
                                sc.adaptive_runs = grid.adaptive_runs;
                                sc.initial_runs = grid.initial_runs_per_scenario;
                                sc.batch_runs = grid.runs_batch_size;
                                sc.max_runs = grid.max_runs_per_scenario;
                                sc.bootstrap_resamples = grid.bootstrap_resamples;
                                sc.success_ci_half_width = grid.success_ci_half_width;
                                sc.recovery_ci_half_width_s = grid.recovery_time_ci_half_width_s;
                                scenarios.push_back(sc);
                            }
                        }
                    }
                }
            }
        }
    return scenarios;
}

} // namespace gnc







