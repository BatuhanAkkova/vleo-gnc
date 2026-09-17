#pragma once
/// @file aerodynamics.hpp
/// CLL panel-method aerodynamic forces and torques.

#include "math/types.hpp"
#include "math/constants.hpp"
#include <cmath>
#include <random>
#include <vector>
#include <string>

namespace gnc {

// ── Polynomial coefficient storage ───────────────────────────────────

struct AeroCoeffs {
    Eigen::MatrixXd cn;  // Nx3 matrix: [val, power_s, power_a]
    Eigen::MatrixXd ct;  // Nx3 matrix: [val, power_s, power_a]
};

/// Evaluate 2D polynomial sum(c[i] * s^p_s[i] * a^p_a[i]).
inline double eval_poly(double s_norm, double angle_norm,
                        const Eigen::MatrixXd& coeffs) {
    double res = 0.0;
    for (int i = 0; i < coeffs.rows(); ++i) {
        res += coeffs(i,0) * std::pow(s_norm, coeffs(i,1))
                            * std::pow(angle_norm, coeffs(i,2));
    }
    return res;
}

// ── CLL Monte Carlo (truth dynamics) ─────────────────────────────────

inline void cll_force_torque(const Vec3& v_rel_body, double rho,
                             std::mt19937& rng,
                             Vec3& f_out, Vec3& t_out) {
    f_out.setZero(); t_out.setZero();
    double v_i = v_rel_body.norm();
    if (v_i < 1e-9) return;

    constexpr double kB = constants::BOLTZMANN;
    constexpr double Tw = constants::T_W;
    constexpr double m_O = constants::MOLAR_MASS_O / constants::AVOGADRO;
    double v_th = std::sqrt(2.0 * kB * Tw / m_O);
    double alpha_t = constants::SIGMA_T * (2.0 - constants::SIGMA_T);
    double q_dyn = 0.5 * rho * v_i * v_i;
    constexpr int M = constants::NUM_PARTICLES;

    std::uniform_real_distribution<double> udist(0.0, 1.0);
    const auto& norms = constants::normals();
    const auto& cents = constants::centers();

    for (int i = 0; i < constants::N_PLATES; ++i) {
        Vec3 n = norms.row(i);
        double cos_theta = -n.dot(v_rel_body) / v_i;
        if (cos_theta <= 0.0) continue;

        double v_in_n = v_i * cos_theta;
        double v_in_t = std::sqrt(std::max(0.0, v_i*v_i - v_in_n*v_in_n));
        double sum_pn = 0.0, sum_pt = 0.0;

        for (int j = 0; j < M; ++j) {
            double k0 = std::max(udist(rng), 1e-15);
            double k1 = udist(rng);
            double k2 = udist(rng);
            double k3 = std::max(udist(rng), 1e-15);

            double term = -constants::ALPHA_N * std::log(k0)
                + (1.0 - constants::ALPHA_N) * (v_in_n/v_th)*(v_in_n/v_th)
                + 2.0 * std::cos(2.0*M_PI*k1)
                  * std::sqrt(std::max(0.0, -constants::ALPHA_N*(1.0-constants::ALPHA_N)*std::log(k0)))
                  * (v_in_n/v_th);
            double v_rn = v_th * std::sqrt(std::max(0.0, term));
            double v_rt = (1.0-constants::SIGMA_T)*v_in_t
                + v_th * std::sqrt(std::max(0.0, alpha_t))
                  * std::cos(2.0*M_PI*k2) * std::sqrt(-std::log(k3));

            sum_pn += v_in_n + v_rn;
            sum_pt += v_in_t - v_rt;
        }

        double cn = 2.0 * (sum_pn / M) / v_i;
        double ct = 2.0 * (sum_pt / M) / v_i;

        // Normal / tangential decomposition
        double dot_in_n = -v_rel_body.dot(n);
        Vec3 fn_vec = dot_in_n * n;
        Vec3 ft_vec = -v_rel_body - fn_vec;
        double t_mag = ft_vec.norm();
        Vec3 t_hat = (t_mag > 1e-9) ? (ft_vec / t_mag).eval() : Vec3::Zero();

        Vec3 fp = constants::AREAS[i] * q_dyn * (cn * (-n) + ct * t_hat);
        f_out += fp;

        Vec3 lever = cents.row(i);  // CoM at origin
        t_out += lever.cross(fp);
    }
}

// ── CLL interpolation (GNC dynamics) ─────────────────────────────────

inline void cll_interp(const Vec3& v_rel_body, double rho,
                       const AeroCoeffs& ac,
                       Vec3& f_out, Vec3& t_out) {
    f_out.setZero(); t_out.setZero();
    double v_i = v_rel_body.norm();
    if (v_i < 1e-9) return;

    constexpr double kB = constants::BOLTZMANN;
    constexpr double Tw = constants::T_W;
    constexpr double m_O = constants::MOLAR_MASS_O / constants::AVOGADRO;
    double v_th = std::sqrt(2.0 * kB * Tw / m_O);
    double s = v_i / v_th;
    double s_norm = (s - 5.0) / (20.0 - 5.0);
    double q_dyn = 0.5 * rho * v_i * v_i;

    const auto& norms = constants::normals();
    const auto& cents = constants::centers();

    for (int i = 0; i < constants::N_PLATES; ++i) {
        Vec3 n = norms.row(i);
        double cos_theta = -n.dot(v_rel_body) / v_i;
        if (cos_theta <= 0.0) continue;

        double angle = std::acos(std::clamp(cos_theta, -1.0, 1.0)) * (180.0/M_PI);
        double angle_norm = angle / 90.0;
        double cn = eval_poly(s_norm, angle_norm, ac.cn);
        double ct = eval_poly(s_norm, angle_norm, ac.ct);
        if (angle < 1e-6) ct = 0.0;
        if (angle > 89.999) { cn = 0.0; ct = 0.0; }

        double dot_in_n = -v_rel_body.dot(n);
        Vec3 fn_vec = dot_in_n * n;
        Vec3 ft_vec = -v_rel_body - fn_vec;
        double t_mag = ft_vec.norm();
        Vec3 t_hat = (t_mag > 1e-9) ? (ft_vec / t_mag).eval() : Vec3::Zero();

        Vec3 fp = constants::AREAS[i] * q_dyn * (cn * (-n) + ct * t_hat);
        f_out += fp;

        Vec3 lever = cents.row(i);
        t_out += lever.cross(fp);
    }
}

// ── Load GSI coefficients from JSON ──────────────────────────────────
AeroCoeffs load_aero_coeffs(const std::string& cn_path, const std::string& ct_path);

} // namespace gnc
