#pragma once
/// @file gnc_dynamics.hpp
/// Reduced-fidelity GNC dynamics RHS for the 13-state vector.

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "environment/gravity.hpp"
#include "environment/density.hpp"
#include "environment/magnetic.hpp"
#include "environment/wind.hpp"
#include "environment/aerodynamics.hpp"
#include "environment/srp.hpp"

namespace gnc {

/// Compute GNC (reduced-fidelity) dynamics RHS.
inline void gnc_rhs(const Vec13& state, double gmst, double mjd,
                    const Vec3& ctrl_torque, const AeroCoeffs& ac,
                    Vec13& dst) {
    Vec3 r = state.segment<3>(0);
    Vec3 v = state.segment<3>(3);
    Vec4 q = state.segment<4>(6);
    Vec3 w = state.segment<3>(10);

    // ── J2 Gravity ───────────────────────────────────────────────────
    Vec3 a_grav = j2_accel(r);
    Vec3 t_gg = gravity_gradient_torque(r, q);

    // ── Density/Wind (simple models) ─────────────────────────────────
    double rho = exp_density(r);
    Vec3 v_rel_eci = wind_simple(r, v);

    // ── DCM ──────────────────────────────────────────────────────────
    Mat3 R = quat_to_dcm(q);
    Vec3 v_rel_body = R * v_rel_eci;

    // ── Aero (interp) ────────────────────────────────────────────────
    Vec3 f_aero, t_aero;
    cll_interp(v_rel_body, rho, ac, f_aero, t_aero);

    // ── SRP (simplified: shadow = 0.5) ───────────────────────────────
    Vec3 r_sun_eci = sun_pos(mjd);
    Vec3 r_sun_body = R * (r_sun_eci - r);
    Vec3 f_srp, t_srp;
    srp_force_torque(r_sun_body, 0.5, f_srp, t_srp);

    // ── Magnetic (dipole) ────────────────────────────────────────────
    Vec3 B_eci = dipole_B(r, gmst);
    Vec3 B_body = R * B_eci;
    Vec3 t_mag = B_torque(constants::m_body(), B_body);

    // ── Translational ────────────────────────────────────────────────
    Vec3 f_eci = R.transpose() * f_aero;
    dst.segment<3>(0) = v;
    dst.segment<3>(3) = a_grav + f_eci / constants::MASS;

    // ── Quaternion kinematics ────────────────────────────────────────
    double qx = q.x(), qy = q.y(), qz = q.z(), qw = q.w();
    double wx = w.x(), wy = w.y(), wz = w.z();
    dst(6)  =  0.5 * ( qw*wx + qy*wz - qz*wy);
    dst(7)  =  0.5 * ( qw*wy + qz*wx - qx*wz);
    dst(8)  =  0.5 * ( qw*wz + qx*wy - qy*wx);
    dst(9)  = -0.5 * ( qx*wx + qy*wy + qz*wz);

    // ── Rotational dynamics ──────────────────────────────────────────
    const Mat3& I = constants::inertia();
    const Mat3& I_inv = constants::inv_inertia();
    Vec3 L = I * w;
    Vec3 t_total = t_aero + t_gg + t_mag + ctrl_torque - w.cross(L);
    dst.segment<3>(10) = I_inv * t_total;
}

} // namespace gnc
