#include <gtest/gtest.h>
#include <cmath>

#include "math/types.hpp"
#include "math/constants.hpp"
#include "environment/gravity.hpp"
#include "environment/density.hpp"
#include "environment/magnetic.hpp"
#include "environment/wind.hpp"
#include "environment/aerodynamics.hpp"
#include "environment/srp.hpp"

using namespace gnc;

// =====================================================================
//  J2 Acceleration
// =====================================================================

TEST(Gravity, J2AccelDirection) {
    // Position on +x axis → gravity should point mostly -x
    Vec3 r(6678137.0, 0.0, 0.0);  // R_E + 300 km
    Vec3 a = j2_accel(r);
    EXPECT_LT(a.x(), 0.0);  // pointing inward
    EXPECT_NEAR(a.y(), 0.0, 1e-4);
    EXPECT_NEAR(a.z(), 0.0, 1e-4);
}

TEST(Gravity, J2AccelMagnitude) {
    Vec3 r(6678137.0, 0.0, 0.0);
    Vec3 a = j2_accel(r);
    double g = a.norm();
    // Should be ~8.9 m/s^2 at 300km altitude
    EXPECT_GT(g, 8.0);
    EXPECT_LT(g, 10.0);
}

// =====================================================================
//  Gravity Gradient Torque
// =====================================================================

TEST(Gravity, GravGradZero) {
    Vec3 r(6678137.0, 0.0, 0.0);
    Vec4 q(0, 0, 0, 1);  // identity
    Vec3 t = gravity_gradient_torque(r, q);
    // Along x-axis with identity rotation → nonzero (off-diagonal inertia)
    // but should be small ~10^-6 Nm
    EXPECT_LT(t.norm(), 1e-3);
}

// =====================================================================
//  Sun / Moon Position
// =====================================================================

TEST(Gravity, SunPosition) {
    // MJD for 2024-03-20 (vernal equinox) → Sun near +x ECI
    double mjd = 60388.5;
    Vec3 r = sun_pos(mjd);
    EXPECT_GT(r.norm(), 1.4e11);   // ~1 AU
    EXPECT_LT(r.norm(), 1.6e11);
}

TEST(Gravity, MoonPosition) {
    double mjd = 60388.5;
    Vec3 r = moon_pos(mjd);
    double d_km = r.norm() / 1000.0;
    EXPECT_GT(d_km, 350000.0);
    EXPECT_LT(d_km, 410000.0);
}

// =====================================================================
//  Third Body Acceleration
// =====================================================================

TEST(Gravity, ThirdBodySmall) {
    Vec3 r(6678137.0, 0.0, 0.0);
    Vec3 rm, a;
    third_body_accel(r, 60388.5, rm, a);
    // Third body perturbation ~10^-6 m/s^2
    EXPECT_LT(a.norm(), 1e-3);
    EXPECT_GT(a.norm(), 1e-10);
}

// =====================================================================
//  Exponential Density
// =====================================================================

TEST(Density, ExpDensityAtRef) {
    Vec3 r(constants::R_EARTH + constants::H_0, 0.0, 0.0);
    double rho = exp_density(r);
    EXPECT_NEAR(rho, constants::RHO_0, constants::RHO_0 * 0.01);
}

TEST(Density, ExpDensityDecay) {
    Vec3 r_low(constants::R_EARTH + 200e3, 0.0, 0.0);
    Vec3 r_high(constants::R_EARTH + 500e3, 0.0, 0.0);
    EXPECT_GT(exp_density(r_low), exp_density(r_high));
}

// =====================================================================
//  Dipole Magnetic Field
// =====================================================================

TEST(Magnetic, DipoleMagnitude) {
    Vec3 r(6678137.0, 0.0, 0.0);
    Vec3 B = dipole_B(r, 0.0);
    double B_mag = B.norm();
    // ~20-60 μT at 300km
    EXPECT_GT(B_mag, 1e-6);
    EXPECT_LT(B_mag, 1e-3);
}

TEST(Magnetic, BTorqueCross) {
    Vec3 m(1, 0, 0);
    Vec3 B(0, 1, 0);
    Vec3 t = B_torque(m, B);
    EXPECT_NEAR(t.x(), 0.0, 1e-15);
    EXPECT_NEAR(t.y(), 0.0, 1e-15);
    EXPECT_NEAR(t.z(), 1.0, 1e-15);
}

// =====================================================================
//  Wind
// =====================================================================

TEST(Wind, SimpleCorotation) {
    Vec3 r(6678137.0, 0.0, 0.0);
    Vec3 v(0.0, 7725.0, 0.0);
    Vec3 v_rel = wind_simple(r, v);
    // v_atm = [-w_e*y, w_e*x, 0] = [0, w_e*r_x, 0]
    double v_atm_y = constants::W_EARTH * r.x();
    EXPECT_NEAR(v_rel.x(), 0.0, 1e-10);
    EXPECT_NEAR(v_rel.y(), v.y() - v_atm_y, 1e-6);
    EXPECT_NEAR(v_rel.z(), 0.0, 1e-10);
}

// =====================================================================
//  SRP
// =====================================================================

TEST(SRP, ShadowFullSun) {
    Vec3 r(6678137.0, 0.0, 0.0);
    Vec3 r_sun_far(1.5e11, 0.0, 0.0);
    double sf = shadow_factor(r, r_sun_far);
    EXPECT_NEAR(sf, 1.0, 1e-6);  // not in shadow
}

TEST(SRP, ShadowUmbra) {
    // Spacecraft directly behind Earth from Sun
    Vec3 r_sun(1.5e11, 0.0, 0.0);   // Sun at +x
    Vec3 r(-6678137.0, 0.0, 0.0);   // SC at -x (behind Earth)
    double sf = shadow_factor(r, r_sun);
    EXPECT_NEAR(sf, 0.0, 0.1);      // in shadow
}

TEST(SRP, ForceZeroInShadow) {
    Vec3 sun_body(1, 0, 0);
    Vec3 f, t;
    srp_force_torque(sun_body, 0.0, f, t);
    EXPECT_EQ(f.norm(), 0.0);
    EXPECT_EQ(t.norm(), 0.0);
}

TEST(SRP, ForceNonzeroInSun) {
    Vec3 sun_body(1, 0, 0);
    Vec3 f, t;
    srp_force_torque(sun_body, 1.0, f, t);
    EXPECT_GT(f.norm(), 0.0);
}

// =====================================================================
//  Aerodynamics (CLL Interp)
// =====================================================================

TEST(Aero, CllInterpZeroVelocity) {
    AeroCoeffs ac;
    ac.cn = Eigen::MatrixXd::Zero(1, 3);
    ac.ct = Eigen::MatrixXd::Zero(1, 3);
    ac.cn(0, 0) = 2.0;  // constant cn = 2

    Vec3 v_body = Vec3::Zero();
    Vec3 f, t;
    cll_interp(v_body, 1e-11, ac, f, t);
    EXPECT_EQ(f.norm(), 0.0);
}
