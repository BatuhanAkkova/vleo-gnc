#pragma once
/// @file srp.hpp
/// Solar Radiation Pressure force/torque and shadow model.

#include "math/types.hpp"
#include "math/constants.hpp"
#include <cmath>

namespace gnc {

// ── Shadow factor ────────────────────────────────────────────────────

inline double shadow_factor(const Vec3& r_body, const Vec3& r_sun) {
    double d_sun  = r_sun.norm();
    double d_body = r_body.norm();

    double cos_theta = r_body.dot(r_sun) / (d_body * d_sun);
    // (void)sin_theta; // not used in current implementation

    // Apparent radii
    double a_sun   = std::asin(constants::R_SUN  / d_sun);
    double a_earth = std::asin(constants::R_EARTH / d_body);

    // Separation angle
    Vec3 u_sun   = (r_sun - r_body).normalized();
    Vec3 u_earth = (-r_body) / d_body;
    double sep_angle = std::acos(std::clamp(u_sun.dot(u_earth), -1.0, 1.0));

    // Full sun
    if (sep_angle >= a_earth + a_sun) return 1.0;
    // Full umbra
    if (sep_angle <= a_earth - a_sun) return 0.0;
    // Penumbra (linear)
    double x = (sep_angle - (a_earth - a_sun)) / (2.0 * a_sun);
    return std::clamp(x, 0.0, 1.0);
}

// ── SRP force / torque ───────────────────────────────────────────────

inline void srp_force_torque(const Vec3& r_sun_body, double sf,
                             Vec3& force_out, Vec3& torque_out) {
    force_out.setZero();
    torque_out.setZero();
    if (sf < 1e-9) return;

    double s_mag = r_sun_body.norm();
    if (s_mag < 1e-9) return;

    Vec3 sun_hat = r_sun_body / s_mag;
    double P = constants::P_SUN * sf;

    const auto& norms  = constants::normals();
    const auto& levs   = constants::levers();

    for (int i = 0; i < constants::N_PLATES; ++i) {
        Vec3 n = norms.row(i);
        double cos_theta = n.dot(sun_hat);
        if (cos_theta <= 0.0) continue;

        double rs = constants::RS[i];
        double rd = constants::RD[i];
        double a  = constants::AREAS[i];

        double f_coeff = -P * a * cos_theta;
        double c_s = 1.0 - rs;
        double c_n = 2.0 * (rs * cos_theta + rd / 3.0);

        Vec3 fp = f_coeff * (c_s * sun_hat + c_n * n);
        force_out += fp;

        Vec3 lv = levs.row(i);
        torque_out += lv.cross(fp);
    }
}

} // namespace gnc
