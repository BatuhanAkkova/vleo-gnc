#pragma once
/// @file gravity.hpp
/// Gravity models: J2, spherical harmonics, gravity gradient, third body, sun/moon position.

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/frames.hpp"
#include <cmath>

namespace gnc {

// ── J2 Gravitational Acceleration (ECI) ──────────────────────────────
inline Vec3 j2_accel(const Vec3& r) {
    constexpr double MU = constants::MU_EARTH;
    constexpr double RE = constants::R_EARTH;
    constexpr double J2 = constants::J2_EARTH;

    double x = r.x(), y = r.y(), z = r.z();
    double r_mag = r.norm();
    double r2 = r_mag * r_mag;
    double r3 = r2 * r_mag;
    double r5 = r2 * r3;

    double pre = 1.5 * J2 * MU * RE * RE / r5;
    double z_sq_r2 = (z * z) / r2;

    return Vec3(
        -(MU / r3) * x + pre * x * (5.0 * z_sq_r2 - 1.0),
        -(MU / r3) * y + pre * y * (5.0 * z_sq_r2 - 1.0),
        -(MU / r3) * z + pre * z * (5.0 * z_sq_r2 - 3.0)
    );
}

// ── Spherical Harmonics ──────────────────────────────────────────────

/// Preallocated buffers for spherical harmonics.
struct SHBuf {
    int n_max;
    Eigen::MatrixXd P;          // (n_max+2) x (n_max+2)
    Eigen::VectorXd ratio_pow;  // (n_max+2)
    Eigen::VectorXd cm, sm;     // (n_max+1)

    explicit SHBuf(int nmax) : n_max(nmax),
        P(Eigen::MatrixXd::Zero(nmax + 2, nmax + 2)),
        ratio_pow(Eigen::VectorXd::Zero(nmax + 2)),
        cm(Eigen::VectorXd::Zero(nmax + 1)),
        sm(Eigen::VectorXd::Zero(nmax + 1)) {}
};

/// Zero-allocation spherical-harmonics gravitational acceleration (ECI).
inline void spherical_harmonics_accel(const Vec3& r_eci, int max_degree,
                                      const Eigen::MatrixXd& C_nm,
                                      const Eigen::MatrixXd& S_nm,
                                      double gmst, SHBuf& buf, Vec3& a_eci_out) {
    constexpr double mu = constants::MU_EARTH;
    constexpr double Re = constants::R_EARTH;

    double cos_g = std::cos(gmst), sin_g = std::sin(gmst);
    double rx_bf =  r_eci.x() * cos_g + r_eci.y() * sin_g;
    double ry_bf = -r_eci.x() * sin_g + r_eci.y() * cos_g;
    double rz_bf =  r_eci.z();

    double r_sq = rx_bf*rx_bf + ry_bf*ry_bf + rz_bf*rz_bf;
    double r = std::sqrt(r_sq);
    double rho_sq = rx_bf*rx_bf + ry_bf*ry_bf;
    double rho = std::sqrt(rho_sq);
    if (rho < 1e-9) rho = 1e-9;

    double sin_phi = rz_bf / r;
    double cos_phi = rho / r;
    if (std::abs(cos_phi) < 1e-10) cos_phi = 1e-10;

    int n_max = max_degree;
    auto& P = buf.P;

    // Legendre recursion
    P(0,0) = 1.0;
    P(1,0) = std::sqrt(3.0) * sin_phi;
    P(1,1) = std::sqrt(3.0) * cos_phi;
    for (int n = 2; n <= n_max; ++n) {
        P(n,n) = std::sqrt((2*n+1.0)/(2.0*n)) * cos_phi * P(n-1,n-1);
        for (int m = 0; m < n; ++m) {
            double anm = std::sqrt(((2.0*n-1)*(2.0*n+1)) / ((n-m)*(double)(n+m)));
            double bnm = std::sqrt(((2.0*n+1)*(n+m-1.0)*(n-m-1.0)) /
                                   ((double)(n-m)*(n+m)*(2.0*n-3)));
            P(n,m) = anm * sin_phi * P(n-1,m) - bnm * P(n-2,m);
        }
    }

    // ratio_pow
    auto& rp = buf.ratio_pow;
    rp(0) = 1.0;
    double ratio_r = Re / r;
    for (int n = 1; n <= n_max; ++n) rp(n) = rp(n-1) * ratio_r;

    // Longitude trig recursion
    double cos_lam = rx_bf / rho, sin_lam = ry_bf / rho;
    auto& cm = buf.cm;
    auto& sm = buf.sm;
    cm(0) = 1.0; sm(0) = 0.0;
    for (int m = 1; m <= n_max; ++m) {
        cm(m) = cm(m-1)*cos_lam - sm(m-1)*sin_lam;
        sm(m) = sm(m-1)*cos_lam + cm(m-1)*sin_lam;
    }

    double sum_r = 0.0, sum_lat = 0.0, sum_lon = 0.0;
    for (int n = 2; n <= n_max; ++n) {
        double rn_term = rp(n);
        for (int m = 0; m <= n; ++m) {
            double C = C_nm(n,m), S = S_nm(n,m);
            double trig = C * cm(m) + S * sm(m);
            double d_trig = m * (S * cm(m) - C * sm(m));
            double gamma = P(n,m) * trig;
            sum_r += -(n + 1) * gamma * rn_term;

            double factor, p_prev;
            if (n == m) { factor = 0.0; p_prev = 0.0; }
            else {
                factor = std::sqrt(((2.0*n+1)*(n*n - m*m)) / (2.0*n-1));
                p_prev = P(n-1,m);
            }
            double der_P = (n * sin_phi * P(n,m) - factor * p_prev) / cos_phi;
            sum_lat += der_P * trig * rn_term;
            sum_lon += P(n,m) * d_trig * rn_term;
        }
    }

    double pre = mu / r_sq;
    double fr = pre * sum_r;
    double fphi = pre * sum_lat;
    double flam = pre * (1.0 / cos_phi) * sum_lon;

    double ax_bf = (fr - mu/r_sq)*cos_phi*cos_lam - fphi*sin_phi*cos_lam - flam*sin_lam;
    double ay_bf = (fr - mu/r_sq)*cos_phi*sin_lam - fphi*sin_phi*sin_lam + flam*cos_lam;
    double az_bf = (fr - mu/r_sq)*sin_phi          + fphi*cos_phi;

    a_eci_out.x() =  ax_bf*cos_g - ay_bf*sin_g;
    a_eci_out.y() =  ax_bf*sin_g + ay_bf*cos_g;
    a_eci_out.z() =  az_bf;
}

// ── Gravity Gradient Torque ──────────────────────────────────────────

inline Vec3 gravity_gradient_torque(const Vec3& r, const Vec4& q) {
    double r_mag = r.norm();
    if (r_mag < 1e-6) return Vec3::Zero();

    Vec3 u_eci = r / r_mag;
    Mat3 R = quat_to_dcm(q);
    Vec3 u_body = R * u_eci;

    const Mat3& I = constants::inertia();
    Vec3 h = I * u_body;

    double coeff = 3.0 * constants::MU_EARTH / (r_mag * r_mag * r_mag);
    return coeff * u_body.cross(h);
}

// ── Sun Position (low-precision analytical, ECI) ─────────────────────

inline void sun_pos_inplace(double mjd, Vec3& out) {
    double T = (mjd - 51544.5) / 36525.0;
    double M = (357.5277233 + 35999.05034 * T) * (M_PI / 180.0);
    double L = 280.4606184 + 36000.77005361 * T;
    double L_rad = (L + 1.914666471*std::sin(M) + 0.019994643*std::sin(2*M)) * (M_PI/180.0);
    double eps = (23.439291 - 0.0130042 * T) * (M_PI / 180.0);
    double r_AU = 1.000140612 - 0.016708617*std::cos(M) - 0.000139589*std::cos(2*M);
    double r_m = r_AU * constants::AU;

    out.x() = r_m * std::cos(L_rad);
    out.y() = r_m * std::cos(eps) * std::sin(L_rad);
    out.z() = r_m * std::sin(eps) * std::sin(L_rad);
}

inline Vec3 sun_pos(double mjd) { Vec3 v; sun_pos_inplace(mjd, v); return v; }

// ── Moon Position (low-precision analytical, ECI) ────────────────────

inline void moon_pos_inplace(double mjd, Vec3& out) {
    double T = (mjd - 51544.5) / 36525.0;
    double L0 = (218.316 + 481267.881 * T) * (M_PI/180.0);
    double l  = (134.963 + 477198.867 * T) * (M_PI/180.0);
    double lp = (357.528 + 35999.050  * T) * (M_PI/180.0);
    double F  = ( 93.272 + 483202.017 * T) * (M_PI/180.0);
    double D  = (297.850 + 445267.111 * T) * (M_PI/180.0);

    double d_lam = 6.289*std::sin(l) - 1.274*std::sin(l-2*D) + 0.658*std::sin(2*D)
                 + 0.214*std::sin(2*l) - 0.186*std::sin(lp) - 0.114*std::sin(2*F)
                 - 0.059*std::sin(2*l-2*D);
    double lam = L0 + d_lam * (M_PI/180.0);

    double d_beta = 5.128*std::sin(F) + 0.280*std::sin(l+F)
                  + 0.278*std::sin(l-F) + 0.173*std::sin(2*D-F);
    double beta = d_beta * (M_PI/180.0);

    double dist_km = 385000.56 - 20905.355*std::cos(l) - 3699.111*std::cos(2*D-l)
                   - 2955.968*std::cos(2*D) - 569.925*std::cos(2*l);
    double r_m = dist_km * 1000.0;

    double eps = (23.439291 - 0.0130042*T) * (M_PI/180.0);

    double cb = std::cos(beta), sb = std::sin(beta);
    double cl = std::cos(lam),  sl = std::sin(lam);
    double ce = std::cos(eps),  se = std::sin(eps);

    out.x() = r_m * cb * cl;
    out.y() = r_m * (ce*cb*sl - se*sb);
    out.z() = r_m * (se*cb*sl + ce*sb);
}

inline Vec3 moon_pos(double mjd) { Vec3 v; moon_pos_inplace(mjd, v); return v; }

// ── Third Body Acceleration ──────────────────────────────────────────

inline void third_body_accel(const Vec3& r_eci, double mjd,
                             Vec3& r_moon, Vec3& a_out) {
    double r_mag = r_eci.norm();
    if (r_mag < 1.0) { a_out.setZero(); return; }

    moon_pos_inplace(mjd, r_moon);

    // Moon term
    Vec3 dm = r_moon - r_eci;
    double d_mag_m = dm.norm();
    double r_mag_m = r_moon.norm();
    double f_moon   = constants::MU_MOON / (d_mag_m * d_mag_m * d_mag_m);
    double f_moon_r = constants::MU_MOON / (r_mag_m * r_mag_m * r_mag_m);

    a_out = f_moon*dm - f_moon_r*r_moon;
}

// ── Load EGM2008 Gravity Coefficients from JSON ──────────────────────
// Implemented in gravity.cpp (requires nlohmann/json)
void load_gravity_coefficients(const std::string& egm_path, int n_max,
                               Eigen::MatrixXd& C_nm, Eigen::MatrixXd& S_nm);

} // namespace gnc
