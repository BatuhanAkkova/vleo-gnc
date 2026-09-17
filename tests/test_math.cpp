#include <gtest/gtest.h>
#include <cmath>

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "math/math_utils.hpp"
#include "math/frames.hpp"
#include "math/interp3d.hpp"
#include "math/integrator.hpp"

using namespace gnc;

// =====================================================================
//  Constants
// =====================================================================

TEST(Constants, EarthMu) {
    EXPECT_NEAR(constants::MU_EARTH, 3.986004418e14, 1e6);
}

TEST(Constants, InertiaSymmetry) {
    const Mat3& I = constants::inertia();
    EXPECT_NEAR(I(0,1), 0.0, 1e-15);
    EXPECT_NEAR(I(1,0), 0.0, 1e-15);
    EXPECT_GT(I(0,0), 0.0);
    EXPECT_GT(I(1,1), 0.0);
    EXPECT_GT(I(2,2), 0.0);
}

TEST(Constants, InvInertiaConsistency) {
    const Mat3 product = constants::inertia() * constants::inv_inertia();
    EXPECT_NEAR((product - Mat3::Identity()).norm(), 0.0, 1e-12);
}

// =====================================================================
//  Quaternion
// =====================================================================

TEST(Quaternion, IdentityDCM) {
    Vec4 q_id(0, 0, 0, 1);
    Mat3 R = quat_to_dcm(q_id);
    EXPECT_NEAR((R - Mat3::Identity()).norm(), 0.0, 1e-14);
}

TEST(Quaternion, ProductIdentity) {
    Vec4 q(0.1, 0.2, 0.3, 0.9);
    q = quat_normalize(q);

    Vec4 q_id(0, 0, 0, 1);
    Vec4 result = quat_product(q, q_id);
    EXPECT_NEAR((result - q).norm(), 0.0, 1e-14);
}

TEST(Quaternion, ProductInverse) {
    Vec4 q(0.1, 0.2, 0.3, 0.9);
    q = quat_normalize(q);

    Vec4 q_inv = quat_conjugate(q);
    Vec4 result = quat_product(q, q_inv);
    // Should be identity
    Vec4 id(0, 0, 0, 1);
    result = quat_ensure_positive_w(result);
    EXPECT_NEAR((result - id).norm(), 0.0, 1e-13);
}

TEST(Quaternion, DCMRoundTrip) {
    // 90° rotation about Z
    double a = M_PI / 4.0;
    Vec4 q(0, 0, std::sin(a), std::cos(a));  // [x,y,z,w]
    Mat3 R = quat_to_dcm(q);

    // Apply to x-hat: should go to y-hat (body sees inertial x as its y)
    Vec3 v_eci(1, 0, 0);
    Vec3 v_body = R * v_eci;
    EXPECT_NEAR(v_body(0), 0.0, 1e-14);
    EXPECT_NEAR(v_body(1), -1.0, 1e-14);
    EXPECT_NEAR(v_body(2), 0.0, 1e-14);
}

TEST(Quaternion, FromRotvec) {
    Vec3 rv(0, 0, M_PI / 2.0);  // 90° about Z
    Vec4 q = quat_from_rotvec(rv);
    q = quat_normalize(q);

    double a = M_PI / 4.0;
    EXPECT_NEAR(q(0), 0.0, 1e-14);
    EXPECT_NEAR(q(1), 0.0, 1e-14);
    EXPECT_NEAR(q(2), std::sin(a), 1e-14);
    EXPECT_NEAR(q(3), std::cos(a), 1e-14);
}

TEST(Quaternion, QuatRotateConsistency) {
    Vec4 q(0.1, 0.2, 0.3, 0.9);
    q = quat_normalize(q);
    Vec3 v(1, 2, 3);

    Vec3 v1 = quat_to_dcm(q) * v;
    Vec3 v2 = quat_rotate(q, v);
    EXPECT_NEAR((v1 - v2).norm(), 0.0, 1e-13);
}

// =====================================================================
//  Binary Search
// =====================================================================

TEST(MathUtils, BinarySearchMiddle) {
    double arr[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    // 2.5 sits between index 2 and 3 → bisect_right returns 3
    int idx = binary_search(arr, 2.5, 5);
    EXPECT_EQ(idx, 3);
}

TEST(MathUtils, BinarySearchExact) {
    double arr[] = {0.0, 1.0, 2.0, 3.0, 4.0};
    // bisect_right(arr, 2.0) → 3  (goes PAST equal elements)
    int idx = binary_search(arr, 2.0, 5);
    EXPECT_EQ(idx, 3);
}

TEST(MathUtils, BinarySearchBounds) {
    double arr[] = {10.0, 20.0, 30.0};
    EXPECT_EQ(binary_search(arr, 5.0, 3), 0);   // below min
    EXPECT_EQ(binary_search(arr, 35.0, 3), 3);  // above max
}

// =====================================================================
//  Frames
// =====================================================================

TEST(Frames, ECItoECEFZeroGMST) {
    Vec3 r_eci(1000.0, 0.0, 0.0);
    Vec3 r_ecef = eci_to_ecef(r_eci, 0.0);
    EXPECT_NEAR((r_ecef - r_eci).norm(), 0.0, 1e-10);
}

TEST(Frames, RoundTrip) {
    Vec3 r_eci(1000.0, 2000.0, 3000.0);
    double gmst = 0.7854;  // ~45 deg
    Vec3 r_ecef = eci_to_ecef(r_eci, gmst);
    Vec3 r_back = ecef_to_eci(r_ecef, gmst);
    EXPECT_NEAR((r_back - r_eci).norm(), 0.0, 1e-10);
}

TEST(Frames, GeodeticEquator) {
    auto g = ecef_to_geodetic(6678.137, 0.0, 0.0);  // R_E + 300 km
    EXPECT_NEAR(g.lat_deg, 0.0, 1e-10);
    EXPECT_NEAR(g.lon_deg, 0.0, 1e-10);
    EXPECT_NEAR(g.alt_km, 300.0, 1.0);  // approximate
}

TEST(Frames, GeodeticPole) {
    auto g = ecef_to_geodetic(0.0, 0.0, 6678.137);
    EXPECT_NEAR(g.lat_deg, 90.0, 1e-10);
    EXPECT_GT(g.alt_km, 0.0);
}

// =====================================================================
//  Interp3D
// =====================================================================

TEST(Interp3D, AxisLookupClamp) {
    double axis[] = {0.0, 1.0, 2.0, 3.0};
    auto lo = axis_lookup(axis, 4, -1.0);
    EXPECT_EQ(lo.idx, 0);
    EXPECT_DOUBLE_EQ(lo.w, 0.0);

    auto hi = axis_lookup(axis, 4, 10.0);
    EXPECT_EQ(hi.idx, 2);
    EXPECT_DOUBLE_EQ(hi.w, 1.0);
}

TEST(Interp3D, TrilinearConstant) {
    // 2×2×2 grid, all ones
    double grid[8];
    std::fill_n(grid, 8, 1.0);
    double v = trilinear(grid, 2, 2, 2, 0, 0, 0, 0.5, 0.5, 0.5);
    EXPECT_NEAR(v, 1.0, 1e-14);
}

TEST(Interp3D, TrilinearCorners) {
    // 2×2×2 grid: value = i + 2*j + 4*k
    double grid[8];
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                grid[i*4 + j*2 + k] = i + 2.0*j + 4.0*k;

    // At corner (0,0,0) → 0
    EXPECT_NEAR(trilinear(grid, 2, 2, 2, 0, 0, 0, 0.0, 0.0, 0.0), 0.0, 1e-14);
    // At corner (1,1,1) → 1+2+4=7
    EXPECT_NEAR(trilinear(grid, 2, 2, 2, 0, 0, 0, 1.0, 1.0, 1.0), 7.0, 1e-14);
    // Center → average of all corners = (0+1+2+3+4+5+6+7)/8 = 3.5
    EXPECT_NEAR(trilinear(grid, 2, 2, 2, 0, 0, 0, 0.5, 0.5, 0.5), 3.5, 1e-14);
}

// =====================================================================
//  RK4 Integrator (simple harmonic oscillator test)
// =====================================================================

TEST(Integrator, RK4SimpleOscillator) {
    // Test RK4 on a simple 13-state where only indices 0,3 are used:
    // x' = v, v' = -x  (SHO with ω = 1)
    // All other states stay zero.
    Vec13 state = Vec13::Zero();
    state(0) = 1.0;  // x0 = 1
    state(3) = 0.0;  // v0 = 0
    // Set quaternion to identity so normalisation doesn't break
    state(9) = 1.0;  // qw = 1

    Vec13 k1, k2, k3, k4, tmp, next;

    auto rhs = [](const Vec13& s, Vec13& dst) {
        dst.setZero();
        dst(0) = s(3);        // dx/dt = v
        dst(3) = -s(0);       // dv/dt = -x
        // quaternion kinematics: d(qw)/dt = 0 (static)
    };

    double dt = 0.001;
    int steps = static_cast<int>(2.0 * M_PI / dt);  // one full period
    for (int i = 0; i < steps; ++i) {
        rk4_step_13(state, dt, rhs, k1, k2, k3, k4, tmp, next);
        state = next;
    }

    // After one period, x ≈ 1, v ≈ 0
    EXPECT_NEAR(state(0), 1.0, 1e-3);
    EXPECT_NEAR(state(3), 0.0, 1e-3);
}
