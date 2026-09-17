#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <filesystem>
#include <fstream>

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "gnc/sensors.hpp"
#include "gnc/actuators.hpp"
#include "gnc/controller.hpp"
#include "gnc/mode_manager.hpp"
#include "gnc/switching_analysis.hpp"
#include "gnc/mekf.hpp"
#include "simulation/campaign_config.hpp"

using namespace gnc;

// =====================================================================
//  Sensors
// =====================================================================

TEST(Sensors, GyroscopeBiased) {
    std::mt19937 rng(42);
    Gyroscope gyro(1e-5, 1e-4, 0.0, 0.1);
    Vec3 true_w(0.01, -0.02, 0.03);
    Vec3 meas = gyro.update(true_w, rng);
    // Should be close to true_w (noise << signal)
    EXPECT_NEAR(meas.x(), true_w.x(), 0.01);
    EXPECT_NEAR(meas.y(), true_w.y(), 0.01);
    EXPECT_NEAR(meas.z(), true_w.z(), 0.01);
}

TEST(Sensors, MagnetometerBias) {
    std::mt19937 rng(42);
    Magnetometer mag(10e-9, 50e-9, &rng);
    Vec3 true_B(2e-5, -1e-5, 3e-5);
    Vec3 meas = mag.update(true_B, rng);
    // Should be close to true_B (bias + noise << signal)
    EXPECT_NEAR(meas.x(), true_B.x(), 1e-4);
    EXPECT_NEAR(meas.y(), true_B.y(), 1e-4);
    EXPECT_NEAR(meas.z(), true_B.z(), 1e-4);
}

TEST(Sensors, SunSensorEclipse) {
    std::mt19937 rng(42);
    SunSensor ss;
    Vec3 dark = Vec3::Zero();
    Vec3 out;
    EXPECT_FALSE(ss.update(dark, out, rng));
}

TEST(Sensors, GPSNoise) {
    std::mt19937 rng(42);
    GPS gps(10.0, 0.1);
    Vec3 pos(6678137.0, 0.0, 0.0);
    Vec3 vel(0.0, 7725.0, 0.0);
    Vec3 mp, mv;
    gps.update(pos, vel, mp, mv, rng);
    EXPECT_NEAR(mp.x(), pos.x(), 100.0);
    EXPECT_NEAR(mv.y(), vel.y(), 1.0);
}

TEST(Sensors, StarTrackerValidity) {
    std::mt19937 rng(42);
    StarTracker st(0.005, 1.0);  // max rate 1 deg/s
    Vec4 q(0, 0, 0, 1);
    Vec3 slow_w(0.001, 0.001, 0.001);  // well below threshold
    Vec3 fast_w(0.1, 0.1, 0.1);        // 5.7 deg/s total

    Vec4 qm;
    EXPECT_TRUE(st.read(q, slow_w, qm, rng));
    EXPECT_NEAR(qm.squaredNorm(), 1.0, 1e-10);

    EXPECT_FALSE(st.read(q, fast_w, qm, rng));
}

// =====================================================================
//  Actuators
// =====================================================================

TEST(Actuators, RWReaction) {
    std::mt19937 rng(42);
    ReactionWheel rw(0.1, 5000, 1e-4, 0.0);  // no jitter
    double body_torque = rw.update(0.05, 0.1, rng);
    EXPECT_NEAR(body_torque, -0.05, 1e-10);
    // Wheel speed should have changed
    EXPECT_NE(rw.speed, 0.0);
}

TEST(Actuators, RWSaturation) {
    std::mt19937 rng(42);
    ReactionWheel rw(0.01, 5000, 1e-4, 0.0);
    double body_torque = rw.update(1.0, 0.1, rng);  // cmd > max
    EXPECT_NEAR(body_torque, -0.01, 1e-10);
}

TEST(Actuators, MtqClamp) {
    Magnetorquer mtq(0.1);
    double cmd[3] = {0.5, -0.05, 0.2};
    double out[3];
    mtq.get_dipole(cmd, out);
    EXPECT_NEAR(out[0], 0.1, 1e-15);
    EXPECT_NEAR(out[1], -0.05, 1e-15);
    EXPECT_NEAR(out[2], 0.1, 1e-15);
}

// =====================================================================
//  Controller
// =====================================================================

TEST(Controller, BDotFirstStep) {
    BDotController bdot(1e6, 0.1);
    Vec3 B(2e-5, -1e-5, 3e-5);
    Vec3 m = bdot.compute(B, 0.1);
    EXPECT_EQ(m.norm(), 0.0);  // first step returns zero
}

TEST(Controller, BDotNonzero) {
    BDotController bdot(1e6, 0.1);
    Vec3 B1(2e-5, -1e-5, 3e-5);
    Vec3 B2(2.1e-5, -0.9e-5, 2.8e-5);
    bdot.compute(B1, 0.1);
    Vec3 m = bdot.compute(B2, 0.1);
    EXPECT_GT(m.norm(), 0.0);
}

TEST(Controller, LQRAtTarget) {
    Eigen::Matrix<double,3,6> K;
    K << 1,0,0, 0.1,0,0,
         0,1,0, 0,0.1,0,
         0,0,1, 0,0,0.1;
    Vec4 q(0,0,0,1);  // at target
    Vec3 w = Vec3::Zero();
    Vec3 u = lqr_control(q, w, q, w, K, Vec3::Zero(), 0.0);
    EXPECT_NEAR(u.norm(), 0.0, 1e-10);
}

TEST(Controller, LQRClipping) {
    Eigen::Matrix<double,3,6> K = Eigen::Matrix<double,3,6>::Zero();
    K.block<3,3>(0,0) = Mat3::Identity(); // Torque = -1.0 * theta
    
    // Create a large error: 10 degrees along X axis
    double angle_10_deg = 10.0 * M_PI / 180.0;
    Vec4 q_target(0, 0, 0, 1);
    Vec4 q_est = quat_from_rotvec(Vec3(angle_10_deg, 0, 0));
    
    Vec3 w = Vec3::Zero();
    Vec3 u = lqr_control(q_est, w, q_target, w, K, Vec3::Zero(), 0.0);
    
    // Resulting theta should be clipped to 5 degrees
    // u = -K * [theta; 0] = -Identity * theta = -theta
    double angle_5_deg = 5.0 * M_PI / 180.0;
    EXPECT_NEAR(std::abs(u.x()), angle_5_deg, 1e-6);
    EXPECT_NEAR(u.y(), 0.0, 1e-6);
    EXPECT_NEAR(u.z(), 0.0, 1e-6);
}

TEST(Controller, LQRUsesGeodesicErrorAtLargeAngle) {
    Eigen::Matrix<double,3,6> K = Eigen::Matrix<double,3,6>::Zero();
    K.block<3,3>(0,0) = Mat3::Identity();
    const Vec4 q_target(0, 0, 0, 1);
    const double angle_rad = 170.0 * M_PI / 180.0;
    const Vec4 q_est = quat_from_rotvec(Vec3(angle_rad, 0.0, 0.0));

    const Vec3 u = lqr_control_raw(
        q_est, Vec3::Zero(), q_target, Vec3::Zero(),
        K, Vec3::Zero(), 3.0);

    EXPECT_NEAR(std::abs(u.x()), angle_rad, 1e-6);
    EXPECT_NEAR(u.y(), 0.0, 1e-9);
    EXPECT_NEAR(u.z(), 0.0, 1e-9);
}


TEST(Controller, MomentumDesat) {
    Vec3 h(0.01, 0, 0);
    Vec3 B(0, 2e-5, 0);
    Vec3 m = momentum_desaturation(h, B, 1e-3);
    EXPECT_GT(m.norm(), 0.0);
}

TEST(Controller, AcquisitionGeometricPDCommandsCorrectDirection) {
    const Vec4 q_target(0.0, 0.0, 0.0, 1.0);
    const Vec4 q_est = quat_from_rotvec(Vec3(30.0 * M_PI / 180.0, 0.0, 0.0));
    const Vec3 u = acquisition_geometric_pd_control(
        q_est, Vec3::Zero(), q_target, Vec3::Zero(),
        Vec3::Zero(), constants::inertia(), 0.002, 0.004, 0.02, 20.0 * M_PI / 180.0);
    EXPECT_LT(u.x(), 0.0);
    EXPECT_NEAR(u.y(), 0.0, 1e-9);
    EXPECT_NEAR(u.z(), 0.0, 1e-9);
    EXPECT_LE(std::abs(u.x()), 0.002 + 1e-12);
}

TEST(Controller, SlewCaptureCommandsBoundedDesiredRateAtLargeAngle) {
    const Vec4 q_target(0.0, 0.0, 0.0, 1.0);
    const Vec4 q_est = quat_from_rotvec(Vec3(170.0 * M_PI / 180.0, 0.0, 0.0));
    const Vec3 u = slew_capture_geometric_raw(
        q_est, Vec3::Zero(), q_target, Vec3::Zero(),
        Vec3(0.0, 0.0, -0.003), constants::inertia(),
        0.0035, 0.05, 0.4 * M_PI / 180.0, 20.0 * M_PI / 180.0);
    EXPECT_LT(u.x(), 0.0);
    EXPECT_NEAR(u.y(), 0.0, 1e-9);
    EXPECT_NEAR(u.z(), 0.0, 1e-9);
    EXPECT_LT(std::abs(u.x()), 0.001);
}

TEST(Controller, NominalGeometricPDICommandsCorrectDirection) {
    const Vec4 q_target(0.0, 0.0, 0.0, 1.0);
    const Vec4 q_est = quat_from_rotvec(Vec3(3.0 * M_PI / 180.0, 0.0, 0.0));
    const Vec3 u = nominal_geometric_pdi_control(
        q_est, Vec3::Zero(), q_target, Vec3::Zero(),
        Vec3::Zero(), constants::inertia(), Vec3::Zero(),
        0.002, 0.003, 0.04, 0.0002, 10.0 * M_PI / 180.0);
    EXPECT_LT(u.x(), 0.0);
    EXPECT_NEAR(u.y(), 0.0, 1e-9);
    EXPECT_NEAR(u.z(), 0.0, 1e-9);
    EXPECT_LE(std::abs(u.x()), 0.002 + 1e-12);
}

// =====================================================================
//  LQR Controller with CARE
// =====================================================================

TEST(Controller, LQRControllerCARE) {
    LQRController ctrl(constants::inertia(), 1.0, 1.0, 1.0);
    // K should be 3×6 with all positive elements on diagonal
    EXPECT_GT(ctrl.K(0,0), 0.0);
    EXPECT_GT(ctrl.K(1,1), 0.0);
    EXPECT_GT(ctrl.K(2,2), 0.0);
}

TEST(ModeManager, RequiresValidDwellBeforeNominal) {
    ModeManager manager({0.10, 0.08, 0.07, 0.15, 0.30, 0.20, 0.20, 1.0, 0.14, 0.12});

    manager.update(0.09, true, 0.1, 0.5);
    EXPECT_TRUE(manager.is_detumble());
    EXPECT_DOUBLE_EQ(manager.valid_dwell_s, 0.0);

    manager.update(0.09, true, 0.1, 1.0);
    manager.update(0.09, true, 0.1, 1.1);
    EXPECT_TRUE(manager.is_detumble());

    manager.update(0.09, true, 0.1, 1.2);
    EXPECT_TRUE(manager.is_rate_capture());
    EXPECT_TRUE(manager.switched_to_rate_capture());

    manager.update(0.065, true, true, true, false, 0.1, 1.3);
    manager.update(0.065, true, true, true, false, 0.1, 1.4);
    manager.update(0.065, true, true, true, false, 0.1, 1.5);
    EXPECT_TRUE(manager.is_transition());
    EXPECT_FALSE(manager.is_nominal());

    manager.update(0.065, true, true, true, true, 0.1, 1.6);
    manager.update(0.065, true, true, true, true, 0.1, 1.7);
    manager.update(0.065, true, true, true, true, 0.1, 1.8);
    EXPECT_TRUE(manager.is_nominal());
}
TEST(ModeManager, HysteresisPreventsChatterInsideBand) {
    ModeManager manager({0.10, 0.08, 0.07, 0.15, 0.10, 0.20, 0.20, 0.0, 0.14, 0.12});

    manager.update(0.09, true, 0.1, 0.0);
    EXPECT_TRUE(manager.is_rate_capture());

    manager.update(0.095, true, 0.1, 0.1);
    EXPECT_TRUE(manager.is_rate_capture());
    EXPECT_FALSE(manager.switched());
}

TEST(ModeManager, SingleAcquisitionTransitionsDirectlyToNominal) {
    ModeManager manager({0.10, 0.08, 0.07, 0.15, 0.10, 0.20, 0.20, 0.0, 0.14, 0.12, false, true});

    manager.update(0.09, true, 0.1, 0.0);
    ASSERT_TRUE(manager.is_rate_capture());

    manager.update(0.065, true, true, true, false, 0.1, 0.1);
    EXPECT_TRUE(manager.is_rate_capture());
    EXPECT_FALSE(manager.switched());

    manager.update(0.065, true, true, true, true, 0.1, 0.2);
    EXPECT_TRUE(manager.is_rate_capture());
    EXPECT_FALSE(manager.switched());

    manager.update(0.065, true, true, true, true, 0.1, 0.3);
    EXPECT_TRUE(manager.is_nominal());
    EXPECT_TRUE(manager.switched_to_nominal());
}


TEST(ModeManager, AcquisitionUsesDedicatedExitThreshold) {
    ModeManager manager({0.10, 0.08, 0.07, 0.15, 0.10, 0.20, 0.20, 0.0, 0.14, 0.12});

    manager.update(0.09, true, 0.1, 0.0);
    ASSERT_TRUE(manager.is_rate_capture());

    manager.update(0.12, true, true, true, false, 0.1, 0.1);
    EXPECT_TRUE(manager.is_rate_capture());
    EXPECT_FALSE(manager.switched());

    manager.update(0.15, true, true, true, false, 0.1, 0.2);
    EXPECT_TRUE(manager.is_detumble());
    EXPECT_TRUE(manager.switched_to_detumble());
}

TEST(ModeManager, ExitsNominalAboveHighThreshold) {
    ModeManager manager({0.10, 0.08, 0.07, 0.15, 0.10, 0.10, 0.10, 0.0, 0.14, 0.12, false, false, 0.2});

    manager.update(0.09, true, 0.1, 0.0);
    ASSERT_TRUE(manager.is_rate_capture());
    manager.update(0.075, true, 0.1, 0.1);
    ASSERT_TRUE(manager.is_slew_capture());
    manager.update(0.065, true, true, true, true, 0.1, 0.2);
    ASSERT_TRUE(manager.is_nominal());

    manager.update(0.16, true, true, true, true, 0.1, 0.3);
    EXPECT_TRUE(manager.is_nominal());
    EXPECT_FALSE(manager.switched());

    manager.update(0.16, true, true, true, true, 0.1, 0.4);
    EXPECT_TRUE(manager.is_detumble());
    EXPECT_TRUE(manager.switched_to_detumble());
    EXPECT_LT(manager.switch_event_value(), 0.0);
}

TEST(ModeManager, SupportsDirectDetumbleToFinePointing) {
    ModeManager manager({0.10, 0.08, 0.07, 0.15, 0.30, 0.20, 0.0, 0.0, 0.14, 0.10, true, false, 0.2});

    manager.update(0.06, true, true, true, true, 0.1, 0.0);
    EXPECT_TRUE(manager.is_nominal());
    EXPECT_TRUE(manager.switched_to_nominal());

    manager.update(0.20, true, true, true, true, 0.1, 0.1);
    EXPECT_TRUE(manager.is_nominal());
    EXPECT_FALSE(manager.switched());

    manager.update(0.20, true, true, true, true, 0.1, 0.2);
    EXPECT_TRUE(manager.is_detumble());
    EXPECT_TRUE(manager.switched_to_detumble());
}


TEST(CampaignConfig, BuildsArchitectureGrid) {
    CampaignGridConfig grid;
    grid.deployment_rates_deg_s = {10.0, 20.0};
    grid.density_scales = {1.0};
    grid.actuator_scales = {0.5};
    grid.runs_per_scenario = 2;
    grid.output_root = "C:/tmp/gnc_cpp_test_results";

    SimConfig base;
    base.x0.setZero();
    auto scenarios = build_scenarios_from_grid(grid, base);
    EXPECT_EQ(scenarios.size(), 4u);
    EXPECT_EQ(scenarios.front().n_runs, 2);
    EXPECT_NEAR(scenarios.front().sim.rw_max_torque, base.rw_max_torque * 0.5, 1e-15);
}

TEST(CampaignConfig, ExpandsStressAndMismatchAxes) {
    CampaignGridConfig grid;
    grid.deployment_rates_deg_s = {50.0};
    grid.density_scales = {1.0};
    grid.actuator_scales = {1.0};
    grid.stress_test = true;
    grid.stress_start_dynamic = true;
    grid.stress_torque_nm_values = {0.001, 0.002};
    grid.onboard_density_error_pct = {0.0, 25.0};
    grid.onboard_gravity_degree_truncation = {2, 8};

    SimConfig base;
    auto scenarios = build_scenarios_from_grid(grid, base);
    EXPECT_EQ(scenarios.size(), 8u);
    EXPECT_TRUE(scenarios.front().sim.stress_test);
    EXPECT_TRUE(scenarios.front().sim.stress_start_dynamic);
    EXPECT_NEAR(scenarios.front().sim.stress_torque_nm, 0.001, 1e-15);
}

TEST(SwitchingAnalysis, PositiveDwellRequirementAndMargin) {
    DwellTimeInputs in;
    in.omega_lo_rad_s = 0.1;
    in.omega_hi_rad_s = 0.2;
    in.max_control_torque_nm = 0.001;
    in.disturbance_torque_bound_nm = 0.0;
    in.inertia = Mat3::Identity() * 0.01;

    const auto result = compute_dwell_time_margin(in, 2.0);
    EXPECT_GT(result.required_dwell_s, 0.0);
    EXPECT_GT(result.margin_s, 0.0);
}



// =====================================================================
//  MEKF
// =====================================================================

TEST(MEKF, PredictPreservesNorm) {
    Vec4 q(0, 0, 0, 1);
    Vec3 beta = Vec3::Zero();
    Mat6 P = Mat6::Identity() * 0.01;
    Mat6 Q = Mat6::Identity() * 1e-6;

    Vec3 omega(0.01, -0.02, 0.03);
    mekf_predict(q, beta, omega, P, Q, 0.1);

    EXPECT_NEAR(q.squaredNorm(), 1.0, 1e-10);
    // P should remain symmetric
    EXPECT_NEAR((P - P.transpose()).norm(), 0.0, 1e-14);
}

TEST(MEKF, UpdateCorrects) {
    Vec4 q(0, 0, 0, 1);
    Vec3 beta = Vec3::Zero();
    Mat6 P = Mat6::Identity() * 0.01;
    Mat3 R_meas = Mat3::Identity() * 1e-6;

    // z_meas ≈ z_pred → small correction
    Vec3 z_ref(1, 0, 0);  // ECI reference
    Vec3 z_meas = z_ref;   // already aligned (identity quat)

    Vec4 q_before = q;
    mekf_update(q, beta, P, z_meas, z_ref, R_meas);

    // Should barely change since measurement matches prediction
    double dq = (q - q_before).norm();
    EXPECT_LT(dq, 0.01);
    EXPECT_NEAR(q.squaredNorm(), 1.0, 1e-10);
}

TEST(MEKF, StarTrackerUpdate) {
    Vec4 q(0, 0, 0, 1);
    Vec3 beta = Vec3::Zero();
    Mat6 P = Mat6::Identity() * 0.01;
    Mat3 R_st = Mat3::Identity() * 1e-8;

    // Exact measurement → should converge
    Vec4 q_meas(0, 0, 0, 1);
    mekf_update_star(q, beta, P, q_meas, R_st);

    EXPECT_NEAR(q.squaredNorm(), 1.0, 1e-10);
    // Covariance should shrink
    EXPECT_LT(P.trace(), 6.0 * 0.01);
}

// =====================================================================
//  UKF
// =====================================================================

#include "gnc/ukf.hpp"

TEST(UKF, WeightsSum) {
    UKFWeights w;
    w.compute(1e-3, 2.0, 0.0);
    double sum_wm = 0.0;
    for (int i = 0; i < UKF_NS; ++i) sum_wm += w.wm[i];
    EXPECT_NEAR(sum_wm, 1.0, 1e-10);
}

TEST(UKF, SigmaPointSymmetry) {
    UVec7 x = UVec7::Zero();
    x.head<3>() = Vec3(6678137.0, 0.0, 0.0);
    x.segment<3>(3) = Vec3(0.0, 7725.0, 0.0);
    x(6) = 0.0;

    UMat7 P = UMat7::Zero();
    P.block<3,3>(0,0) = Mat3::Identity() * 100.0;   // 10 m σ
    P.block<3,3>(3,3) = Mat3::Identity() * 0.01;     // 0.1 m/s σ
    P(6,6) = 0.01;

    UKFWeights w;
    w.compute(1.0, 2.0, 0.0);   // alpha=1 for orbit UKF

    SigmaMatrix sigmas;
    sigma_points(x, P, w.lambda, sigmas);

    // Weighted mean should recover x exactly (linear operation)
    UVec7 x_mean = UVec7::Zero();
    for (int i = 0; i < UKF_NS; ++i)
        x_mean += w.wm[i] * sigmas.row(i).transpose();
    EXPECT_NEAR((x_mean - x).norm(), 0.0, 1e-6);
}

TEST(UKF, PredictShiftsPosition) {
    UVec7 x = UVec7::Zero();
    x.head<3>() = Vec3(6678137.0, 0.0, 0.0);
    x.segment<3>(3) = Vec3(0.0, 7725.0, 0.0);

    UMat7 P = UMat7::Zero();
    P.block<3,3>(0,0) = Mat3::Identity() * 100.0;
    P.block<3,3>(3,3) = Mat3::Identity() * 0.01;
    P(6,6) = 0.01;
    UMat7 Q = UMat7::Identity() * 1e-6;
    UMat6 R = UMat6::Identity() * 10.0;

    AeroCoeffs ac;
    ac.cn = Eigen::MatrixXd::Zero(1, 3);
    ac.ct = Eigen::MatrixXd::Zero(1, 3);

    Vec4 q(0, 0, 0, 1);
    UKF ukf(x, P, Q, R, 1.0, 2.0, 0.0);
    UVec7 x_before = ukf.x;
    ukf.predict(0.1, q, ac);

    // After 0.1s: position shifts by ~v*dt ≈ 772.5 m
    double dr = (ukf.x.head<3>() - x_before.head<3>()).norm();
    EXPECT_GT(dr, 100.0);
    EXPECT_LT(dr, 1500.0);
}

TEST(UKF, UpdateShrinksCov) {
    UVec7 x = UVec7::Zero();
    x.head<3>() = Vec3(6678137.0, 0.0, 0.0);
    x.segment<3>(3) = Vec3(0.0, 7725.0, 0.0);

    UMat7 P = UMat7::Zero();
    P.block<3,3>(0,0) = Mat3::Identity() * 1000.0;    // 31 m σ
    P.block<3,3>(3,3) = Mat3::Identity() * 1.0;       // 1 m/s σ
    P(6,6) = 0.01;
    UMat7 Q = UMat7::Identity() * 1e-6;
    UMat6 R = UMat6::Identity() * 10.0;  // GPS ~3 m

    AeroCoeffs ac;
    ac.cn = Eigen::MatrixXd::Zero(1, 3);
    ac.ct = Eigen::MatrixXd::Zero(1, 3);

    Vec4 q(0, 0, 0, 1);
    UKF ukf(x, P, Q, R, 1.0, 2.0, 0.0);
    ukf.predict(0.1, q, ac);

    double trace_before = ukf.P.trace();

    // GPS measurement = predicted state (perfect measurement)
    UVec6 z;
    z.head<3>() = ukf.x.head<3>();
    z.tail<3>() = ukf.x.segment<3>(3);
    ukf.update(z);

    EXPECT_LT(ukf.P.trace(), trace_before);
}

TEST(UKF, OnboardGravityDegreeChangesAcceleration) {
    Eigen::MatrixXd C_nm = Eigen::MatrixXd::Zero(9, 9);
    Eigen::MatrixXd S_nm = Eigen::MatrixXd::Zero(9, 9);
    C_nm(2, 0) = -4.84165371736e-4;
    C_nm(4, 0) = 5.0e-7;
    C_nm(8, 0) = -2.0e-7;

    const Vec3 r(constants::R_EARTH + 450e3, 120e3, 310e3);
    const double gmst = 0.37;

    SHBuf sh2(2);
    SHBuf sh4(4);
    SHBuf sh8(8);
    const Vec3 a2 = onboard_gravity_accel(r, gmst, 2, C_nm, S_nm, sh2);
    const Vec3 a4 = onboard_gravity_accel(r, gmst, 4, C_nm, S_nm, sh4);
    const Vec3 a8 = onboard_gravity_accel(r, gmst, 8, C_nm, S_nm, sh8);

    EXPECT_GT((a2 - a4).norm(), 1e-10);
    EXPECT_GT((a4 - a8).norm(), 1e-10);
}

