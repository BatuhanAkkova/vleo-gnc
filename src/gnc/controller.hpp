#pragma once
/// @file controller.hpp
/// Control laws: LQR, B-dot, momentum desaturation, dipole limiting.

#include "math/types.hpp"
#include "math/quaternion.hpp"
#include <cmath>
#include <algorithm>

namespace gnc {

inline Vec3 quaternion_error_rotvec(const Vec4& q_target, const Vec4& q_est) {
    Vec4 q_inv = quat_conjugate(q_target);
    Vec4 q_err = quat_product(q_inv, q_est);
    if (q_err.w() < 0.0) q_err = -q_err;

    const Vec3 q_vec(q_err.x(), q_err.y(), q_err.z());
    const double sin_half = q_vec.norm();
    if (sin_half <= 1e-9) return Vec3::Zero();

    const double angle = 2.0 * std::atan2(sin_half, std::clamp(q_err.w(), -1.0, 1.0));
    return q_vec * (angle / sin_half);
}

inline Vec3 saturate_torque_preserve_direction(const Vec3& u, double max_torque) {
    if (max_torque <= 0.0) return u;
    const double max_abs = u.cwiseAbs().maxCoeff();
    if (max_abs <= max_torque) return u;
    return u * (max_torque / max_abs);
}

// ── LQR core ────────────────────────

inline Vec3 lqr_control_raw(const Vec4& q_est, const Vec3& w_est,
                            const Vec4& q_target, const Vec3& w_target,
                            const Eigen::Matrix<double,3,6>& K,
                            const Vec3& t_ff, double max_linear_err_rad) {
    Vec3 theta = quaternion_error_rotvec(q_target, q_est);

    const double err_magnitude = theta.norm();
    if (max_linear_err_rad > 0.0 && err_magnitude > max_linear_err_rad) {
        theta *= (max_linear_err_rad / err_magnitude);
    }

    Vec3 w_err = w_est - w_target;

    Eigen::Matrix<double,6,1> x;
    x.head<3>() = theta;
    x.tail<3>() = w_err;

    return -(K * x) - t_ff;
}

inline Vec3 lqr_control_with_cap(const Vec4& q_est, const Vec3& w_est,
                                 const Vec4& q_target, const Vec3& w_target,
                                 const Eigen::Matrix<double,3,6>& K,
                                 const Vec3& t_ff, double max_torque,
                                 double max_linear_err_rad) {
    const Vec3 u = lqr_control_raw(q_est, w_est, q_target, w_target, K, t_ff,
                                   max_linear_err_rad);
    return saturate_torque_preserve_direction(u, max_torque);
}

inline Vec3 lqr_control(const Vec4& q_est, const Vec3& w_est,
                        const Vec4& q_target, const Vec3& w_target,
                        const Eigen::Matrix<double,3,6>& K,
                        const Vec3& t_ff, double max_torque) {
    return lqr_control_with_cap(q_est, w_est, q_target, w_target, K, t_ff,
                                max_torque, 5.0 * M_PI / 180.0);
}



inline Vec3 rate_capture_control_raw(const Vec3& w_est, const Vec3& w_target,
                                     const Vec3& h_rw, const Mat3& inertia,
                                     double rate_gain) {
    const Vec3 rate_error = w_est - w_target;
    const Vec3 coriolis_comp = w_est.cross(inertia * w_est + h_rw);
    return -rate_gain * rate_error + coriolis_comp;
}

inline Vec3 slew_capture_geometric_raw(const Vec4& q_est, const Vec3& w_est,
                                       const Vec4& q_target, const Vec3& w_target,
                                       const Vec3& h_rw, const Mat3& inertia,
                                       double attitude_gain, double rate_gain,
                                       double ref_rate_limit_rad_s,
                                       double max_attitude_error_rad) {
    (void)attitude_gain;
    Vec4 q_inv = quat_conjugate(q_target);
    Vec4 q_err = quat_product(q_inv, q_est);
    if (q_err.w() < 0.0) q_err = -q_err;

    const Vec3 q_vec(q_err.x(), q_err.y(), q_err.z());
    const double sin_half = q_vec.norm();
    Vec3 rotvec = Vec3::Zero();
    if (sin_half > 1e-9) {
        const double angle = 2.0 * std::atan2(sin_half, std::clamp(q_err.w(), -1.0, 1.0));
        rotvec = q_vec * (angle / sin_half);
    }

    Vec3 limited_rotvec = rotvec;
    const double rotvec_norm = limited_rotvec.norm();
    const double effective_max_error = max_attitude_error_rad > 1e-9 ? max_attitude_error_rad : rotvec_norm;
    if (effective_max_error > 1e-9 && rotvec_norm > effective_max_error) {
        limited_rotvec *= (effective_max_error / rotvec_norm);
    }

    Vec3 w_ref = w_target;
    const double limited_norm = limited_rotvec.norm();
    if (limited_norm > 1e-9 && ref_rate_limit_rad_s > 0.0 && effective_max_error > 1e-9) {
        const double slew_fraction = std::clamp(limited_norm / effective_max_error, 0.0, 1.0);
        w_ref -= ref_rate_limit_rad_s * slew_fraction * (limited_rotvec / limited_norm);
    }

    const Vec3 rate_error = w_est - w_ref;
    // During large-angle acquisition, direct rate damping is more robust than
    // compensating wheel/body gyroscopic coupling from a still-misaligned frame.
    return -rate_gain * rate_error;
}

inline Vec3 acquisition_geometric_pd_raw(const Vec4& q_est, const Vec3& w_est,
                                         const Vec4& q_target, const Vec3& w_target,
                                         const Vec3& h_rw, const Mat3& inertia,
                                         double attitude_gain, double rate_gain,
                                         double max_linear_err_rad) {
    Vec4 q_inv = quat_conjugate(q_target);
    Vec4 q_err = quat_product(q_inv, q_est);
    if (q_err.w() < 0.0) q_err = -q_err;

    const Vec3 q_vec(q_err.x(), q_err.y(), q_err.z());
    const double sin_half = q_vec.norm();
    Vec3 rotvec = Vec3::Zero();
    double angle = 0.0;
    if (sin_half > 1e-9) {
        angle = 2.0 * std::atan2(sin_half, std::clamp(q_err.w(), -1.0, 1.0));
        rotvec = q_vec * (angle / sin_half);
    }

    double attitude_scale = 1.0;
    if (max_linear_err_rad > 0.0 && angle > max_linear_err_rad) {
        attitude_scale = max_linear_err_rad / angle;
    }

    const Vec3 rate_error = w_est - w_target;
    const Vec3 coriolis_comp = w_est.cross(inertia * w_est + h_rw);
    return -attitude_gain * attitude_scale * rotvec - rate_gain * rate_error + coriolis_comp;
}

inline Vec3 acquisition_geometric_pd_control(const Vec4& q_est, const Vec3& w_est,
                                             const Vec4& q_target, const Vec3& w_target,
                                             const Vec3& h_rw, const Mat3& inertia,
                                             double max_torque, double attitude_gain,
                                             double rate_gain, double max_linear_err_rad) {
    const Vec3 u = acquisition_geometric_pd_raw(
        q_est, w_est, q_target, w_target, h_rw, inertia,
        attitude_gain, rate_gain, max_linear_err_rad);
    return saturate_torque_preserve_direction(u, max_torque);
}

inline Vec3 nominal_geometric_pdi_raw(const Vec4& q_est, const Vec3& w_est,
                                      const Vec4& q_target, const Vec3& w_target,
                                      const Vec3& h_rw, const Mat3& inertia,
                                      const Vec3& attitude_integral,
                                      double attitude_gain, double rate_gain,
                                      double integral_gain, double max_linear_err_rad,
                                      const Vec3& t_ff = Vec3::Zero()) {
    Vec3 rotvec = quaternion_error_rotvec(q_target, q_est);
    const double rotvec_norm = rotvec.norm();
    if (max_linear_err_rad > 0.0 && rotvec_norm > max_linear_err_rad) {
        rotvec *= (max_linear_err_rad / rotvec_norm);
    }

    const Vec3 rate_error = w_est - w_target;
    const Vec3 coriolis_comp = w_est.cross(inertia * w_est + h_rw);
    return -attitude_gain * rotvec
         - rate_gain * rate_error
         - integral_gain * attitude_integral
         + coriolis_comp
         - t_ff;
}

inline Vec3 nominal_geometric_pdi_control(const Vec4& q_est, const Vec3& w_est,
                                          const Vec4& q_target, const Vec3& w_target,
                                          const Vec3& h_rw, const Mat3& inertia,
                                          const Vec3& attitude_integral,
                                          double max_torque, double attitude_gain,
                                          double rate_gain, double integral_gain,
                                          double max_linear_err_rad,
                                          const Vec3& t_ff = Vec3::Zero()) {
    const Vec3 u = nominal_geometric_pdi_raw(
        q_est, w_est, q_target, w_target, h_rw, inertia, attitude_integral,
        attitude_gain, rate_gain, integral_gain, max_linear_err_rad, t_ff);
    return saturate_torque_preserve_direction(u, max_torque);
}

inline Vec3 momentum_desaturation(const Vec3& h_sys, const Vec3& B_body,
                                  double k_desat) {
    double B_norm2 = B_body.squaredNorm();
    if (B_norm2 < 1e-12) return Vec3::Zero();
    return k_desat * h_sys.cross(B_body) / B_norm2;
}

// ── B-dot control ────────────────────────────────────────────────────

inline Vec3 bdot_control(const Vec3& B_body, const Vec3& B_prev,
                         double dt, double k_bdot) {
    if (dt < 1e-6) return Vec3::Zero();
    Vec3 B_dot = (B_body - B_prev) / dt;
    return -k_bdot * B_dot;
}

// ── Dipole limiter ───────────────────────────────────────────────────

inline Vec3 limit_dipole(const Vec3& m_cmd, double max_dipole) {
    if (max_dipole <= 0.0) return m_cmd;
    double m_norm = m_cmd.norm();
    if (m_norm > max_dipole) return m_cmd * (max_dipole / m_norm);
    return m_cmd;
}

// ── BDot Controller class ────────────────────────────────────────────

struct BDotController {
    double k_bdot;
    double max_dipole;
    Vec3 B_prev;
    bool first = true;

    BDotController(double k = 1e6, double md = 0.1)
        : k_bdot(k), max_dipole(md), B_prev(Vec3::Zero()) {}

    Vec3 compute(const Vec3& B_body, double dt) {
        if (first) { B_prev = B_body; first = false; return Vec3::Zero(); }
        Vec3 m_cmd = bdot_control(B_body, B_prev, dt, k_bdot);
        Vec3 m_mtq = limit_dipole(m_cmd, max_dipole);
        B_prev = B_body;
        return m_mtq;
    }

    Vec3 compute(const Vec3& B_body, const Vec3& w_body) {
        Vec3 B_dot = -w_body.cross(B_body);
        Vec3 m_cmd = -k_bdot * B_dot;
        return limit_dipole(m_cmd, max_dipole);
    }
};

// ── LQR Controller class (solves CARE at construction) ───────────────
// Requires Eigen unsupported CARE module — implemented in controller.cpp.

struct LQRController {
    Mat3 I_sc;
    Eigen::Matrix<double,3,6> K;
    Vec3 t_ff = Vec3::Zero();

    /// Construct and solve CARE immediately.
    LQRController(const Mat3& inertia, double q_w, double w_w, double u_w);

    /// Returns (u_rw, m_mtq).
    std::pair<Vec3, Vec3> compute(const Vec4& q_est, const Vec3& w_est,
                                  const Vec4& q_target, const Vec3& w_target,
                                  const Vec3* B_body = nullptr,
                                  const Vec3* h_sys = nullptr,
                                  const Vec3* t_ff_in = nullptr,
                                  double max_torque = 0.0,
                                  double k_desat = 1e-3,
                                  double max_dip = 0.1) const;
};

} // namespace gnc
