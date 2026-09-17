#pragma once
/// @file truth_dynamics.hpp
/// High-fidelity truth dynamics RHS for the 13-state vector.

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "environment/gravity.hpp"
#include "environment/density.hpp"
#include "environment/magnetic.hpp"
#include "environment/wind.hpp"
#include "environment/aerodynamics.hpp"
#include "environment/srp.hpp"
#include <random>

namespace gnc {

/// All preloaded environment data needed by the truth RHS.
struct TruthEnv {
    int sh_degree = 0;
    Eigen::MatrixXd C_nm, S_nm;
    SHBuf sh_buf{0};
    DensityLUT density_lut;
    IgrfLUT igrf_lut;
    WindLUT wind_lut;
    double rho_scale = 1.0;

    void init(int deg) {
        sh_degree = deg;
        sh_buf = SHBuf(deg);
    }
};

/// Compute truth dynamics RHS (13-state derivative).
inline void truth_rhs(const Vec13& state, double gmst, double mjd,
                      const Vec3& ctrl_torque, const Vec3& m_mtq,
                      const Vec3& h_rw,
                      TruthEnv& env, std::mt19937& rng,
                      Vec13& dst) {
    // Decompose state
    Vec3 r = state.segment<3>(0);
    Vec3 v = state.segment<3>(3);
    Vec4 q = state.segment<4>(6);
    Vec3 w = state.segment<3>(10);

    // ── Gravity (SH + third body) ────────────────────────────────────
    Vec3 a_grav;
    spherical_harmonics_accel(r, env.sh_degree, env.C_nm, env.S_nm,
                              gmst, env.sh_buf, a_grav);
    Vec3 t_gg = gravity_gradient_torque(r, q);
    Vec3 r_moon, a_third;
    third_body_accel(r, mjd, r_moon, a_third);
    Vec3 r_sun = sun_pos(mjd);

    // ── Density + Wind ───────────────────────────────────────────────
    double rho = nrlmsise_density(r, gmst, env.density_lut) * env.rho_scale;
    Vec3 v_rel_eci = wind_hwm14(r, v, gmst, env.wind_lut);

    // ── DCM (body from ECI) ──────────────────────────────────────────
    Mat3 R = quat_to_dcm(q);
    Vec3 v_rel_body = R * v_rel_eci;

    // ── Aero (CLL Monte Carlo) ───────────────────────────────────────
    Vec3 f_aero, t_aero;
    cll_force_torque(v_rel_body, rho, rng, f_aero, t_aero);

    // ── SRP ──────────────────────────────────────────────────────────
    Vec3 dr_sun_eci = r_sun - r;
    Vec3 r_sun_body = R * dr_sun_eci;
    double sf = shadow_factor(r, r_sun);
    Vec3 f_srp, t_srp;
    srp_force_torque(r_sun_body, sf, f_srp, t_srp);

    // ── Magnetic ─────────────────────────────────────────────────────
    Vec3 B_eci = igrf_B(r, gmst, env.igrf_lut);
    Vec3 B_body = R * B_eci;
    Vec3 m_total = constants::m_body() + m_mtq;
    Vec3 t_mag = B_torque(m_total, B_body);

    // ── Translational: dv ────────────────────────────────────────────
    Vec3 f_total_body = f_aero + f_srp;
    Vec3 f_total_eci = R.transpose() * f_total_body;

    dst.segment<3>(0) = v;
    dst.segment<3>(3) = a_grav + a_third + f_total_eci / constants::MASS;

    // ── Quaternion kinematics ────────────────────────────────────────
    double qx = q.x(), qy = q.y(), qz = q.z(), qw = q.w();
    double wx = w.x(), wy = w.y(), wz = w.z();
    dst(6)  =  0.5 * ( qw*wx + qy*wz - qz*wy);
    dst(7)  =  0.5 * ( qw*wy + qz*wx - qx*wz);
    dst(8)  =  0.5 * ( qw*wz + qx*wy - qy*wx);
    dst(9)  = -0.5 * ( qx*wx + qy*wy + qz*wz);

    // ── Rotational dynamics: dw = I^1(T - wxH_sys) ─────────────────────
    const Mat3& I = constants::inertia();
    const Mat3& I_inv = constants::inv_inertia();
    Vec3 L = I * w + h_rw;
    Vec3 t_total = t_aero + t_srp + t_gg + t_mag + ctrl_torque - w.cross(L);
    dst.segment<3>(10) = I_inv * t_total;
}

} // namespace gnc
