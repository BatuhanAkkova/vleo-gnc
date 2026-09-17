#pragma once
/// @file simulation.hpp
/// Main simulation loop: truth propagation + GNC estimation + control.

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/integrator.hpp"
#include "dynamics/truth_dynamics.hpp"
#include "dynamics/gnc_dynamics.hpp"
#include "gnc/sensors.hpp"
#include "gnc/actuators.hpp"
#include "gnc/controller.hpp"
#include "gnc/mekf.hpp"
#include "gnc/mode_manager.hpp"
#include "gnc/switching_analysis.hpp"
#include "gnc/ukf.hpp"
#include "environment/gravity.hpp"
#include "environment/magnetic.hpp"
#include "utils/logger.hpp"
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <iostream>
#include <stdexcept>

namespace gnc {

struct SimConfig {
    double dt       = 0.1;
    double t_end    = 5400.0;
    double detumble_control_dt = -1.0;
    double transition_control_dt = -1.0;
    double nominal_control_dt  = -1.0;
    bool enable_early_termination = false;
    double settle_pointing_threshold_deg = 0.5;
    double settle_rate_threshold_deg_s = 0.25;
    double settle_pointing_rms_threshold_deg = 0.25;
    double settle_rate_rms_threshold_deg_s = 0.1;
    double settle_hold_time_s = 300.0;
    double settle_confirmation_time_s = 60.0;
    double mjd0     = 60388.5;
    double gmst0    = 0.0;
    double w_earth  = constants::W_EARTH;
    int    sh_degree = 8;

    Vec13 x0 = Vec13::Zero();

    double q_weight   = 6.0;
    double w_weight   = 18.0;
    double u_weight   = 1800.0;
    double bdot_k     = 8e5;
    double bdot_max   = 0.3;
    double desat_k    = 2e5;
    double desat_max  = 0.3;
    bool detumble_rw_assist_enable = true;
    double detumble_rw_assist_entry_rate_deg = 12.0;
    double detumble_rw_assist_exit_rate_deg = 4.0;
    double detumble_rw_rate_gain = 0.05;
    double detumble_rw_torque_scale = 1.0;
    double detumble_rw_desat_blend = 0.35;
    double rw_max_torque = 0.002;
    double rw_max_speed = 150.0;
    double acquisition_threshold_deg = 5.0;
    double acquisition_max_linear_err_deg = 20.0;
    double nominal_entry_pointing_deg = 15.0;
    double nominal_entry_consistency_deg = 15.0;
    double nominal_entry_rate_deg = 0.8;
    double required_nominal_time_s = 8.0;
    double required_nominal_exit_time_s = 2.0;
    double nominal_exit_pointing_deg = 20.0;
    double nominal_exit_consistency_deg = 20.0;
    double slew_capture_entry_rate_deg = 1.5;
    double slew_capture_exit_rate_deg = 3.0;
    double nominal_exit_rate_deg = 2.0;
    double acquisition_exit_rate_deg = 10.0;
    double required_acquisition_time_s = 20.0;
    double required_slew_capture_time_s = 4.0;

    double min_detumble_time = 10.0;
    bool direct_detumble_to_fine_pointing = false;
    bool single_acquisition_mode = false;
    double disturbance_torque_bound_nm = 2e-6;
    double max_controller_degraded_time_s = 20.0;
    bool stress_test = false;
    double stress_torque_nm = 0.0;
    double stress_start_s = 0.0;
    bool stress_start_dynamic = false;
    double stress_duration_s = 0.0;
    double onboard_density_error_pct = 0.0;
    int onboard_gravity_degree_truncation = 2;

    double gyro_rrw = 1e-6;
    double gyro_arw = 2.76e-4;
    double gyro_bias_limit = 0.5;
    double mag_noise = 0.2e-6;
    double mag_bias  = 50e-9;
    double gps_pos_noise = 10.0;
    double gps_vel_noise = 0.1;
    double st_noise_deg  = 0.005;
    double st_max_rate   = 2.0;
    double precision_sensor_hold_s = 4.0;
    double acquisition_desat_trigger_frac = 0.9;
    double acquisition_desat_rate_gate_deg = 1.5;
    double acquisition_attitude_gain = 0.0035;
    double rate_capture_rate_gain = 0.03;
    double slew_capture_ref_rate_deg = 0.8;
    double acquisition_rate_gain = 0.03;
    double acquisition_lqr_scale = 0.1;
    double nominal_attitude_gain = 0.003;
    double nominal_rate_gain = 0.04;
    double nominal_integral_gain = 0.0002;
    double nominal_integral_limit_deg_s = 600.0;

    double mekf_p0    = 0.01;
    double mekf_q_att = 1e-6;
    double mekf_q_bias= 1e-8;
    double mekf_r_mag = 1e-6;
    double mekf_r_st  = 1e-8;



    double ukf_p_pos  = 100.0;
    double ukf_p_vel  = 0.01;
    double ukf_p_rho  = 0.01;
    double ukf_q_diag = 1e-2;
    double ukf_r_pos  = 225.0;
    double ukf_r_vel  = 0.01;
    double ukf_alpha  = 1.0;
    double ukf_beta   = 2.0;
    double ukf_kappa  = 0.0;

    std::string output_path = "results.bin";
};

struct SimulationResult {
    Logger log;
    bool settled_success = false;
    bool early_terminated = false;
    double stop_time_s = 0.0;
    double settle_start_s = std::numeric_limits<double>::quiet_NaN();
    double settle_complete_s = std::numeric_limits<double>::quiet_NaN();
    double settle_pointing_rms_deg = std::numeric_limits<double>::quiet_NaN();
    double settle_pointing_max_deg = std::numeric_limits<double>::quiet_NaN();
    double settle_rate_rms_deg_s = std::numeric_limits<double>::quiet_NaN();
    double settle_rate_max_deg_s = std::numeric_limits<double>::quiet_NaN();
    int control_updates = 0;

    explicit SimulationResult(int n_cols) : log(n_cols) {}
};

inline double normalized_control_period(double candidate_dt, double fallback_dt) {
    if (candidate_dt <= 0.0) return fallback_dt;
    return std::max(candidate_dt, fallback_dt);
}

inline double quaternion_angle_deg(const Vec4& q_a, const Vec4& q_b) {
    const double dot = std::clamp(std::abs(q_a.dot(q_b)), 0.0, 1.0);
    return 2.0 * std::acos(dot) * 180.0 / M_PI;
}

inline Vec3 wheel_speed_safe_torque_command(const Vec3& body_torque, const Vec3& wheel_speeds,
                                            double dt, double wheel_inertia, double max_wheel_speed,
                                            bool& limited) {
    if (wheel_inertia <= 0.0 || max_wheel_speed <= 0.0 || dt <= 0.0) return body_torque;
    Vec3 safe = body_torque;
    for (int i = 0; i < 3; ++i) {
        const double next_speed = wheel_speeds(i) - safe(i) * dt / wheel_inertia;
        if (next_speed > max_wheel_speed) {
            safe(i) = (wheel_speeds(i) - max_wheel_speed) * wheel_inertia / dt;
            limited = true;
        } else if (next_speed < -max_wheel_speed) {
            safe(i) = (wheel_speeds(i) + max_wheel_speed) * wheel_inertia / dt;
            limited = true;
        }
    }
    return safe;
}

inline SimulationResult run_simulation(const SimConfig& cfg, TruthEnv& env,
                                       const AeroCoeffs& ac, unsigned seed = 42) {
#ifndef GNC_USE_OSQP
#endif
    std::mt19937 rng(seed);
    const double nan = std::numeric_limits<double>::quiet_NaN();

    Vec13 x = cfg.x0;

    Gyroscope gyro(cfg.gyro_rrw, cfg.gyro_arw, 0.0, cfg.dt);
    {
        const double bl = cfg.gyro_bias_limit * M_PI / 180.0;
        std::uniform_real_distribution<double> ud(-bl, bl);
        gyro.bias = Vec3(ud(rng), ud(rng), ud(rng));
    }
    Magnetometer mag(cfg.mag_noise, cfg.mag_bias, &rng);
    GPS gps(cfg.gps_pos_noise, cfg.gps_vel_noise);
    StarTracker st(cfg.st_noise_deg, cfg.st_max_rate);

    ReactionWheel rw[3];
    for (auto& w : rw) w = ReactionWheel(cfg.rw_max_torque, cfg.rw_max_speed, 1e-4, 0.0);

    BDotController bdot_ctrl(cfg.bdot_k, cfg.bdot_max);
    LQRController lqr_ctrl(constants::inertia(), cfg.q_weight, cfg.w_weight, cfg.u_weight);

    using M6 = Eigen::Matrix<double,6,6>;
    M6 P0 = M6::Identity() * cfg.mekf_p0;
    M6 Q_proc = M6::Zero();
    Q_proc.block<3,3>(0,0) = Mat3::Identity() * cfg.mekf_q_att;
    Q_proc.block<3,3>(3,3) = Mat3::Identity() * cfg.mekf_q_bias;
    Mat3 R_meas = Mat3::Identity() * cfg.mekf_r_mag;
    Mat3 R_st = Mat3::Identity() * cfg.mekf_r_st;

    Vec4 q_init = x.segment<4>(6);
    q_init = quat_normalize(q_init + Vec4(0.01, -0.01, 0.02, 0.0).normalized() * 0.1);
    MEKF mekf(q_init, P0, Q_proc, R_meas, R_st);

    UVec7 ukf_x0 = UVec7::Zero();
    ukf_x0.head<3>() = x.segment<3>(0) + Vec3(100, -200, 50);
    ukf_x0.segment<3>(3) = x.segment<3>(3) + Vec3(0.1, -0.1, 0.05);
    ukf_x0(6) = 0.0;

    UMat7 ukf_P0 = UMat7::Zero();
    ukf_P0.block<3,3>(0,0) = Mat3::Identity() * cfg.ukf_p_pos;
    ukf_P0.block<3,3>(3,3) = Mat3::Identity() * cfg.ukf_p_vel;
    ukf_P0(6,6) = cfg.ukf_p_rho;
    UMat7 ukf_Q = UMat7::Identity() * cfg.ukf_q_diag;
    UMat6 ukf_R = UMat6::Zero();
    ukf_R.block<3,3>(0,0) = Mat3::Identity() * cfg.ukf_r_pos;
    ukf_R.block<3,3>(3,3) = Mat3::Identity() * cfg.ukf_r_vel;

    UKF ukf(ukf_x0, ukf_P0, ukf_Q, ukf_R,
            cfg.ukf_alpha, cfg.ukf_beta, cfg.ukf_kappa);
    ukf.onboard_density_scale = 1.0 + cfg.onboard_density_error_pct / 100.0;
    ukf.onboard_gravity_degree = cfg.onboard_gravity_degree_truncation;
    ukf.onboard_C_nm = env.C_nm;
    ukf.onboard_S_nm = env.S_nm;
    if (cfg.onboard_gravity_degree_truncation > 1 && env.C_nm.rows() > 0 && env.S_nm.rows() > 0) {
        ukf.onboard_sh_buf = SHBuf(cfg.onboard_gravity_degree_truncation);
    }

    constexpr int LOG_COLS = 107;
    SimulationResult result(LOG_COLS);

    Vec13 k1, k2, k3, k4, x_tmp;

    const int n_steps = static_cast<int>(cfg.t_end / cfg.dt);
    double t = 0.0;
    Vec3 ctrl_torque = Vec3::Zero();
    Vec3 m_mtq_cmd = Vec3::Zero();
    ModeManager mode_manager({
        cfg.acquisition_threshold_deg * M_PI / 180.0,
        cfg.slew_capture_entry_rate_deg * M_PI / 180.0,
        cfg.nominal_entry_rate_deg * M_PI / 180.0,
        cfg.nominal_exit_rate_deg * M_PI / 180.0,
        cfg.required_acquisition_time_s,
        cfg.required_slew_capture_time_s,
        cfg.required_nominal_time_s,
        cfg.min_detumble_time,
        cfg.acquisition_exit_rate_deg * M_PI / 180.0,
        cfg.slew_capture_exit_rate_deg * M_PI / 180.0,
        cfg.direct_detumble_to_fine_pointing,
        cfg.single_acquisition_mode,
        cfg.required_nominal_exit_time_s
    });
    const DwellTimeInputs dwell_inputs{
        cfg.nominal_entry_rate_deg * M_PI / 180.0,
        cfg.nominal_exit_rate_deg * M_PI / 180.0,
        cfg.rw_max_torque,
        cfg.disturbance_torque_bound_nm,
        constants::inertia()
    };
    double active_stress_start_s = cfg.stress_start_s;

    const double detumble_control_period_s = normalized_control_period(cfg.detumble_control_dt, cfg.dt);
    const double transition_control_period_s = normalized_control_period(
        cfg.transition_control_dt > 0.0 ? cfg.transition_control_dt : cfg.nominal_control_dt,
        cfg.dt);
    const double nominal_control_period_s = normalized_control_period(cfg.nominal_control_dt, cfg.dt);
    double time_since_control_update_s = std::numeric_limits<double>::infinity();
    double controller_degraded_time_s = 0.0;
    double precision_sensor_hold_timer_s = 0.0;
    double attitude_ready_hold_timer_s = 0.0;
    double acquisition_star_hold_timer_s = 0.0;

    const Vec4 q_target(0, 0, 0, 1);
    const Vec3 w_target = Vec3::Zero();
    double settled_candidate_start_s = nan;
    bool settled_window_confirmed = false;
    double settle_pointing_sq_sum = 0.0;
    double settle_rate_sq_sum = 0.0;
    double settle_pointing_max = 0.0;
    double settle_rate_max = 0.0;
    int settle_sample_count = 0;

    Vec3 step_acquisition_raw_torque = Vec3::Zero();
    double step_nominal_raw_torque_utilization = nan;
    double step_nominal_torque_saturated = 0.0;
    double step_nominal_raw_torque_norm = nan;
    double acquisition_entry_rate_deg_s = nan;
    Vec3 nominal_attitude_error_int = Vec3::Zero();
    Vec3 step_nominal_integral_state = Vec3::Zero();
    Vec3 step_nominal_integral_torque = Vec3::Zero();
    Vec3 step_nominal_state_torque = Vec3::Zero();
    double step_nominal_integral_freeze = 0.0;
    double step_nominal_integral_update = 0.0;

    for (int step = 0; step < n_steps; ++step) {
        const double gmst = cfg.gmst0 + cfg.w_earth * t;
        const double mjd  = cfg.mjd0 + t / 86400.0;

        Vec3 omega_meas = gyro.update(x.segment<3>(10), rng);
        Vec3 B_eci = dipole_B(x.segment<3>(0), gmst);
        Mat3 Rq = quat_to_dcm(x.segment<4>(6));
        Vec3 B_body_true = Rq * B_eci;
        Vec3 B_meas = mag.update(B_body_true, rng);
        Vec3 r_gps, v_gps;
        gps.update(x.segment<3>(0), x.segment<3>(3), r_gps, v_gps, rng);

        Vec4 q_st;
        const bool st_valid = st.read(x.segment<4>(6), x.segment<3>(10), q_st, rng);
        precision_sensor_hold_timer_s = st_valid ? cfg.precision_sensor_hold_s : std::max(0.0, precision_sensor_hold_timer_s - cfg.dt);
            const Vec3 beta_est = mekf.beta;
        const Vec4 q_est_for_control = mekf.q;
        const Vec3 w_corr = omega_meas - beta_est;
        const double w_norm = w_corr.norm();
        const double est_pointing_error_deg = quaternion_angle_deg(q_est_for_control, q_target);
        const double precision_pointing_error_deg = st_valid ? quaternion_angle_deg(q_st, q_target) : std::numeric_limits<double>::infinity();
        const double estimate_precision_consistency_deg = st_valid ? quaternion_angle_deg(q_est_for_control, q_st) : std::numeric_limits<double>::infinity();
        acquisition_star_hold_timer_s = std::max(0.0, acquisition_star_hold_timer_s - cfg.dt);
        const bool acquisition_star_hold_active =
            mode_manager.is_acquisition() && acquisition_star_hold_timer_s > 0.0 && st_valid;
        const bool attitude_entry_ready = st_valid &&
            precision_pointing_error_deg <= cfg.nominal_entry_pointing_deg &&
            estimate_precision_consistency_deg <= cfg.nominal_entry_consistency_deg;
        const bool attitude_exit_ready = st_valid &&
            precision_pointing_error_deg <= cfg.nominal_exit_pointing_deg &&
            estimate_precision_consistency_deg <= cfg.nominal_exit_consistency_deg;
        const bool attitude_ready = mode_manager.is_nominal() ? attitude_exit_ready : attitude_entry_ready;
        attitude_ready_hold_timer_s = attitude_ready ? cfg.precision_sensor_hold_s : std::max(0.0, attitude_ready_hold_timer_s - cfg.dt);
        const bool attitude_ready_latched = attitude_ready || attitude_ready_hold_timer_s > 0.0;
        const bool precision_sensor_ready = st_valid || precision_sensor_hold_timer_s > 0.0;
        const bool estimator_healthy = true;
        const bool controller_ready = controller_degraded_time_s < cfg.max_controller_degraded_time_s;
        const bool nominal_exit_rate_violation = mode_manager.is_nominal() &&
            w_norm > cfg.nominal_exit_rate_deg * M_PI / 180.0;
        const bool nominal_exit_precision_lost = mode_manager.is_nominal() && !precision_sensor_ready;
        const bool nominal_exit_controller_unready = mode_manager.is_nominal() && !controller_ready;
        const bool nominal_exit_estimator_unhealthy = mode_manager.is_nominal() && !estimator_healthy;
        const bool nominal_exit_attitude_unready = mode_manager.is_nominal() && !attitude_ready_latched;
        const bool nominal_exit_requested =
            nominal_exit_rate_violation ||
            nominal_exit_precision_lost ||
            nominal_exit_controller_unready ||
            nominal_exit_estimator_unhealthy ||
            nominal_exit_attitude_unready;
        mode_manager.update(w_norm, precision_sensor_ready, controller_ready, estimator_healthy, attitude_ready_latched, cfg.dt, t);

        if (mode_manager.switched_to_rate_capture()) {
            acquisition_entry_rate_deg_s = w_norm * 180.0 / M_PI;
        }
        if (!mode_manager.is_nominal()) {
            nominal_attitude_error_int.setZero();
        }
        if (mode_manager.switched_to_acquisition() || mode_manager.switched_to_nominal()) {
            ukf.P = ukf_P0;
            if (cfg.stress_test && cfg.stress_start_dynamic) {
                std::uniform_real_distribution<double> stress_offset(-5.0, 5.0);
                active_stress_start_s = t + stress_offset(rng);
            }
        }

        const double active_control_period_s = 
            (mode_manager.is_detumble() ? detumble_control_period_s
               : (mode_manager.is_transition() ? transition_control_period_s : nominal_control_period_s));
        const bool control_update_due = mode_manager.switched() ||
            time_since_control_update_s >= active_control_period_s - 1e-12;

        step_acquisition_raw_torque.setZero();
        step_nominal_raw_torque_utilization = nan;
        step_nominal_torque_saturated = 0.0;
        step_nominal_raw_torque_norm = nan;
        step_nominal_integral_state.setZero();
        step_nominal_integral_torque.setZero();
        step_nominal_state_torque.setZero();
        step_nominal_integral_freeze = 0.0;
        step_nominal_integral_update = 0.0;

        if (control_update_due) {
            result.control_updates += 1;

            if (mode_manager.is_detumble()) {
                const Vec3 h_rw(rw[0].speed * rw[0].inertia,
                                rw[1].speed * rw[1].inertia,
                                rw[2].speed * rw[2].inertia);
                const Vec3 wheel_speeds(rw[0].speed, rw[1].speed, rw[2].speed);
                const Vec3 bdot_m_cmd = bdot_ctrl.compute(B_meas, w_corr);
                const bool detumble_rw_assist_active =
                    cfg.detumble_rw_assist_enable &&
                    w_corr.norm() >= cfg.detumble_rw_assist_exit_rate_deg * M_PI / 180.0 &&
                    (mode_manager.time_in_mode_s >= cfg.min_detumble_time ||
                     w_corr.norm() >= cfg.detumble_rw_assist_entry_rate_deg * M_PI / 180.0);
                if (detumble_rw_assist_active) {
                    const Vec3 rw_rate_cmd = rate_capture_control_raw(
                        w_corr, Vec3::Zero(), h_rw, constants::inertia(), cfg.detumble_rw_rate_gain);
                    bool wheel_limited = false;
                    ctrl_torque = wheel_speed_safe_torque_command(
                        saturate_torque_preserve_direction(
                            rw_rate_cmd,
                            std::max(0.0, cfg.detumble_rw_torque_scale) * cfg.rw_max_torque),
                        wheel_speeds, cfg.dt, rw[0].inertia, cfg.rw_max_speed, wheel_limited);
                    const Vec3 rw_unload_cmd = momentum_desaturation(h_rw, B_meas, cfg.desat_k);
                    m_mtq_cmd = limit_dipole(
                        bdot_m_cmd + cfg.detumble_rw_desat_blend * rw_unload_cmd,
                        std::max(cfg.bdot_max, cfg.desat_max));
                } else {
                    m_mtq_cmd = bdot_m_cmd;
                    ctrl_torque = Vec3::Zero();
                }
            } else {
                const Vec3 w_est_bias_corr = omega_meas - beta_est;
                const Vec3 h_rw(rw[0].speed * rw[0].inertia,
                                rw[1].speed * rw[1].inertia,
                                rw[2].speed * rw[2].inertia);
                const Vec3 wheel_speeds(rw[0].speed, rw[1].speed, rw[2].speed);
                const Vec3 explicit_m_cmd = momentum_desaturation(h_rw, B_meas, cfg.desat_k);
                const Vec3 explicit_m_limited = limit_dipole(explicit_m_cmd, cfg.desat_max);
                const double wheel_speed_ratio = cfg.rw_max_speed > 0.0 ? wheel_speeds.cwiseAbs().maxCoeff() / cfg.rw_max_speed : 0.0;
                const bool wheel_desat_needed = wheel_speed_ratio > cfg.acquisition_desat_trigger_frac;

                                    const bool low_rate_for_desat =
                        w_est_bias_corr.norm() <= cfg.acquisition_desat_rate_gate_deg * M_PI / 180.0;
                    const Eigen::Matrix<double,3,6> acquisition_K =
                        (cfg.acquisition_lqr_scale > 0.0)
                            ? (cfg.acquisition_lqr_scale * lqr_ctrl.K)
                            : lqr_ctrl.K;
                    const Eigen::Matrix<double,3,6> nominal_K =
                        (cfg.single_acquisition_mode && cfg.acquisition_lqr_scale > 0.0)
                            ? acquisition_K
                            : lqr_ctrl.K;
                    const double nominal_max_linear_err_rad =
                        (cfg.single_acquisition_mode && cfg.acquisition_max_linear_err_deg > 0.0)
                            ? (cfg.acquisition_max_linear_err_deg * M_PI / 180.0)
                            : (5.0 * M_PI / 180.0);
                    Vec3 commanded_torque = Vec3::Zero();
                    if (mode_manager.is_acquisition()) {
                        step_acquisition_raw_torque = lqr_control_raw(
                            q_est_for_control, w_est_bias_corr, q_target, w_target,
                            acquisition_K, lqr_ctrl.t_ff,
                            cfg.acquisition_max_linear_err_deg * M_PI / 180.0);
                        commanded_torque = saturate_torque_preserve_direction(
                            step_acquisition_raw_torque, cfg.rw_max_torque);
                    } else {
                        Vec3 nominal_integral_candidate = nominal_attitude_error_int;
                        const double nominal_integral_limit_rad_s =
                            cfg.nominal_integral_limit_deg_s > 0.0
                                ? (cfg.nominal_integral_limit_deg_s * M_PI / 180.0)
                                : 0.0;
                        const Vec3 nominal_theta = quaternion_error_rotvec(q_target, q_est_for_control);
                        if (precision_sensor_ready && estimator_healthy) {
                            nominal_integral_candidate += nominal_theta * active_control_period_s;
                            const double int_norm = nominal_integral_candidate.norm();
                            if (nominal_integral_limit_rad_s > 0.0 && int_norm > nominal_integral_limit_rad_s) {
                                nominal_integral_candidate *= (nominal_integral_limit_rad_s / int_norm);
                            }
                        }

                        auto nominal_state_torque = [&](const Vec3& theta_limited) {
                            const Vec3 rate_error = w_est_bias_corr - w_target;
                            const Vec3 coriolis_comp = w_est_bias_corr.cross(constants::inertia() * w_est_bias_corr + h_rw);
                            return -cfg.nominal_attitude_gain * theta_limited
                                 - cfg.nominal_rate_gain * rate_error
                                 + coriolis_comp
                                 - lqr_ctrl.t_ff;
                        };

                        auto limited_theta = nominal_theta;
                        const double theta_norm = limited_theta.norm();
                        if (nominal_max_linear_err_rad > 0.0 && theta_norm > nominal_max_linear_err_rad) {
                            limited_theta *= (nominal_max_linear_err_rad / theta_norm);
                        }

                        auto nominal_raw_from_integral = [&](const Vec3& integral_state) {
                            return nominal_geometric_pdi_raw(
                                q_est_for_control, w_est_bias_corr, q_target, w_target,
                                h_rw, constants::inertia(), integral_state,
                                cfg.nominal_attitude_gain, cfg.nominal_rate_gain,
                                cfg.nominal_integral_gain, nominal_max_linear_err_rad,
                                lqr_ctrl.t_ff);
                        };

                        Vec3 raw_nominal_torque = nominal_raw_from_integral(nominal_integral_candidate);
                        Vec3 saturated_nominal_torque = saturate_torque_preserve_direction(
                            raw_nominal_torque, cfg.rw_max_torque);
                        const bool nominal_saturated =
                            (raw_nominal_torque - saturated_nominal_torque).norm() > 1e-12;
                        if (nominal_saturated) {
                            step_nominal_integral_freeze = 1.0;
                            raw_nominal_torque = nominal_raw_from_integral(nominal_attitude_error_int);
                            saturated_nominal_torque = saturate_torque_preserve_direction(
                                raw_nominal_torque, cfg.rw_max_torque);
                        } else {
                            nominal_attitude_error_int = nominal_integral_candidate;
                            step_nominal_integral_update = 1.0;
                        }

                        step_nominal_integral_state = nominal_attitude_error_int;
                        step_nominal_state_torque = nominal_state_torque(limited_theta);
                        step_nominal_integral_torque = -cfg.nominal_integral_gain * nominal_attitude_error_int;
                        step_nominal_raw_torque_norm = raw_nominal_torque.norm();
                        if (cfg.rw_max_torque > 0.0) {
                            step_nominal_raw_torque_utilization =
                                raw_nominal_torque.cwiseAbs().maxCoeff() / cfg.rw_max_torque;
                        }
                        step_nominal_torque_saturated = nominal_saturated ? 1.0 : 0.0;
                        commanded_torque = saturated_nominal_torque;
                    }
                    bool wheel_limited = false;
                    ctrl_torque = wheel_speed_safe_torque_command(
                        commanded_torque, wheel_speeds, cfg.dt, rw[0].inertia,
                        cfg.rw_max_speed, wheel_limited);
                    if (wheel_desat_needed && low_rate_for_desat) {
                        const double excess = std::clamp((wheel_speed_ratio - cfg.acquisition_desat_trigger_frac) /
                            std::max(1e-6, 1.0 - cfg.acquisition_desat_trigger_frac), 0.0, 1.0);
                        m_mtq_cmd = explicit_m_limited * excess;
                    } else {
                        m_mtq_cmd = Vec3::Zero();
                    }
                    controller_degraded_time_s = std::max(0.0, controller_degraded_time_s - 2.0 * cfg.dt);
                
            }

            time_since_control_update_s = 0.0;
        } else {
            time_since_control_update_s += cfg.dt;
        }

        Vec3 t_rw_actual = Vec3::Zero();
        for (int i = 0; i < 3; ++i) {
            t_rw_actual(i) = rw[i].update(-ctrl_torque(i), cfg.dt, rng);
        }
        if (cfg.stress_test && t >= active_stress_start_s && t < active_stress_start_s + cfg.stress_duration_s) {
            t_rw_actual.x() += cfg.stress_torque_nm;
        }

        mekf.predict(omega_meas, cfg.dt);
        mekf.update(B_meas, dipole_B(ukf.x.head<3>(), gmst));
        if (mode_manager.is_nominal() && st_valid) {
            mekf.update_star_tracker(q_st);
        }

        const Vec4 q_est_for_orbit = mekf.q;
        ukf.predict(cfg.dt, gmst, q_est_for_orbit, ac);
        UVec6 z_gps;
        z_gps.head<3>() = r_gps;
        z_gps.tail<3>() = v_gps;
        ukf.update(z_gps);

        const double rho_true = nrlmsise_density(x.head<3>(), gmst, env.density_lut) * env.rho_scale;
        const double rho_est = std::exp(ukf.x(6)) * exp_density(ukf.x.head<3>());
        const Vec4 q_est_log = mekf.q;
        const Vec3 beta_log = mekf.beta;
        const Mat6 P_log = mekf.P;
        const DwellTimeResult dwell = compute_dwell_time_margin(dwell_inputs, mode_manager.time_in_mode_s);
        const double filter_health = 1.0;

        double row[LOG_COLS];
        row[0] = t;
        for (int i = 0; i < 13; ++i) row[1+i] = x(i);
        for (int i = 0; i < 4; ++i) row[14+i] = q_est_log(i);
        for (int i = 0; i < 3; ++i) row[18+i] = (omega_meas - beta_log)(i);
        for (int i = 0; i < 3; ++i) row[21+i] = ukf.x(i);
        for (int i = 0; i < 3; ++i) row[24+i] = ukf.x(3+i);
        for (int i = 0; i < 3; ++i) row[27+i] = beta_log(i);
        for (int i = 0; i < 6; ++i) row[30+i] = P_log(i,i);
        for (int i = 0; i < 7; ++i) row[36+i] = ukf.P(i,i);
        row[43] = rho_true;
        row[44] = ukf.x(6);
        row[45] = rho_est;
        for (int i = 0; i < 3; ++i) row[46+i] = rw[i].speed;
        for (int i = 0; i < 3; ++i) row[49+i] = t_rw_actual(i);
        for (int i = 0; i < 3; ++i) row[52+i] = m_mtq_cmd(i);
        row[55] = mode_manager.mode_value();
        row[56] = mode_manager.time_in_mode_s;
        row[57] = mode_manager.valid_dwell_s;
        row[58] = mode_manager.switch_event_value();
        row[59] = mode_manager.last_switch_rate_rad_s;
        row[60] = dwell.required_dwell_s;
        row[61] = dwell.margin_s;
        row[66] = 0.0;
        row[67] = nan;
        row[68] = nan;
        row[69] = filter_health;
        row[70] = 1.0;
        row[71] = (st_valid ? 1.0 : 0.0);
        row[75] = st_valid ? 1.0 : 0.0;
        row[76] = controller_degraded_time_s;
        row[77] = 0.0;
        row[78] = step_acquisition_raw_torque.x();
        row[79] = step_acquisition_raw_torque.y();
        row[80] = step_acquisition_raw_torque.z();
        row[81] = mode_manager.is_acquisition() && std::isfinite(acquisition_entry_rate_deg_s)
            ? (x.segment<3>(10).norm() * 180.0 / M_PI - acquisition_entry_rate_deg_s)
            : 0.0;
        row[82] = acquisition_star_hold_timer_s;
        row[83] = 0.0;
        row[84] = mode_manager.transition_submode_value();
        row[85] = st_valid ? precision_pointing_error_deg : nan;
        row[86] = st_valid ? estimate_precision_consistency_deg : nan;
        row[87] = nominal_exit_rate_violation ? 1.0 : 0.0;
        row[88] = nominal_exit_precision_lost ? 1.0 : 0.0;
        row[89] = nominal_exit_controller_unready ? 1.0 : 0.0;
        row[90] = nominal_exit_estimator_unhealthy ? 1.0 : 0.0;
        row[91] = nominal_exit_attitude_unready ? 1.0 : 0.0;
        row[92] = nominal_exit_requested ? 1.0 : 0.0;
        row[93] = step_nominal_raw_torque_utilization;
        row[94] = step_nominal_torque_saturated;
        row[95] = step_nominal_raw_torque_norm;
        row[96] = step_nominal_integral_state.x();
        row[97] = step_nominal_integral_state.y();
        row[98] = step_nominal_integral_state.z();
        row[99] = step_nominal_integral_torque.x();
        row[100] = step_nominal_integral_torque.y();
        row[101] = step_nominal_integral_torque.z();
        row[102] = step_nominal_state_torque.x();
        row[103] = step_nominal_state_torque.y();
        row[104] = step_nominal_state_torque.z();
        row[105] = step_nominal_integral_freeze;
        row[106] = step_nominal_integral_update;
        result.log.log(row, LOG_COLS);

        const double true_pointing_error_deg = quaternion_angle_deg(x.segment<4>(6), q_target);
        const double true_rate_deg_s = x.segment<3>(10).norm() * 180.0 / M_PI;
        const bool within_convergence_max_bounds =
            mode_manager.is_nominal() &&
            true_pointing_error_deg <= cfg.settle_pointing_threshold_deg &&
            true_rate_deg_s <= cfg.settle_rate_threshold_deg_s;

        if (within_convergence_max_bounds) {
            if (!std::isfinite(settled_candidate_start_s)) {
                settled_candidate_start_s = t;
                settle_pointing_sq_sum = 0.0;
                settle_rate_sq_sum = 0.0;
                settle_pointing_max = 0.0;
                settle_rate_max = 0.0;
                settle_sample_count = 0;
            }

            settle_pointing_sq_sum += true_pointing_error_deg * true_pointing_error_deg;
            settle_rate_sq_sum += true_rate_deg_s * true_rate_deg_s;
            settle_pointing_max = std::max(settle_pointing_max, true_pointing_error_deg);
            settle_rate_max = std::max(settle_rate_max, true_rate_deg_s);
            settle_sample_count += 1;

            const double satisfied_for_s = t - settled_candidate_start_s + cfg.dt;
            const double settle_required_s = cfg.settle_hold_time_s + cfg.settle_confirmation_time_s;
            const double pointing_rms_deg = std::sqrt(settle_pointing_sq_sum / std::max(settle_sample_count, 1));
            const double rate_rms_deg_s = std::sqrt(settle_rate_sq_sum / std::max(settle_sample_count, 1));
            const bool within_convergence_rms_bounds =
                pointing_rms_deg <= cfg.settle_pointing_rms_threshold_deg &&
                rate_rms_deg_s <= cfg.settle_rate_rms_threshold_deg_s;

            if (!settled_window_confirmed && satisfied_for_s >= settle_required_s && within_convergence_rms_bounds) {
                settled_window_confirmed = true;
                result.settled_success = true;
                result.settle_start_s = settled_candidate_start_s;
                result.settle_complete_s = t;
                result.settle_pointing_rms_deg = pointing_rms_deg;
                result.settle_pointing_max_deg = settle_pointing_max;
                result.settle_rate_rms_deg_s = rate_rms_deg_s;
                result.settle_rate_max_deg_s = settle_rate_max;
            }
            if (cfg.enable_early_termination && settled_window_confirmed) {
                result.early_terminated = true;
                result.stop_time_s = t;
                break;
            }
        } else {
            settled_candidate_start_s = nan;
            settle_pointing_sq_sum = 0.0;
            settle_rate_sq_sum = 0.0;
            settle_pointing_max = 0.0;
            settle_rate_max = 0.0;
            settle_sample_count = 0;
        }

        const Vec3 h_rw_vec(rw[0].speed * rw[0].inertia,
                            rw[1].speed * rw[1].inertia,
                            rw[2].speed * rw[2].inertia);
        auto rhs = [&](const Vec13& s, Vec13& d) {
            truth_rhs(s, gmst, mjd, t_rw_actual, m_mtq_cmd, h_rw_vec, env, rng, d);
        };

        rhs(x, k1);
        x_tmp = x + 0.5 * cfg.dt * k1;
        rhs(x_tmp, k2);
        x_tmp = x + 0.5 * cfg.dt * k2;
        rhs(x_tmp, k3);
        x_tmp = x + cfg.dt * k3;
        rhs(x_tmp, k4);

        x += (cfg.dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        x.segment<4>(6) = quat_normalize(x.segment<4>(6));
        t += cfg.dt;
        result.stop_time_s = t;
    }

    return result;
}

} // namespace gnc















