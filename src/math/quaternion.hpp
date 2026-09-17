#pragma once
/// @file quaternion.hpp
/// Quaternion helpers.  Convention: scalar-last [x,y,z,w].
/// Wraps Eigen::Quaterniond whose coeffs() is already [x,y,z,w].

#include "math/types.hpp"
#include <cmath>

namespace gnc {

// ── Conversion between Vec4 [x,y,z,w] and Eigen::Quaterniond ────────

/// Build a Quat from a Vec4 stored as [x,y,z,w].
inline Quat quat_from_vec4(const Vec4& v) {
    // Eigen ctor: Quaterniond(w, x, y, z)
    return Quat(v(3), v(0), v(1), v(2));
}

/// Extract [x,y,z,w] Vec4 from an Eigen Quaterniond.
inline Vec4 quat_to_vec4(const Quat& q) {
    return q.coeffs();   // already [x,y,z,w]
}

// ── Core operations ──────────────────────────────────────────────────

/// Quaternion product p cross q  (Hamilton convention, scalar-last).
inline Vec4 quat_product(const Vec4& p, const Vec4& q) {
    const Quat qp = quat_from_vec4(p) * quat_from_vec4(q);
    return qp.coeffs();
}

/// Conjugate (inverse for unit quaternion).
inline Vec4 quat_conjugate(const Vec4& q) {
    return Vec4(-q(0), -q(1), -q(2), q(3));
}

/// Normalise in-place and return.
inline Vec4 quat_normalize(const Vec4& q) {
    double n = q.norm();
    // If norm is too small or NaN, return identity but check for NaNs first
    if (std::isnan(n)) return q; // Propagate NaN
    if (n < 1e-12) return Vec4(0, 0, 0, 1);
    return q / n;
}

/// Quaternion -> DCM (rotation matrix: body-from-inertial).
///   v_body = R * v_eci
inline Mat3 quat_to_dcm(const Vec4& q) {
    const double qx = q(0), qy = q(1), qz = q(2), qw = q(3);
    Mat3 R;
    R(0,0) = 1 - 2*(qy*qy + qz*qz);
    R(0,1) = 2*(qx*qy + qw*qz);
    R(0,2) = 2*(qx*qz - qw*qy);
    R(1,0) = 2*(qx*qy - qw*qz);
    R(1,1) = 1 - 2*(qx*qx + qz*qz);
    R(1,2) = 2*(qy*qz + qw*qx);
    R(2,0) = 2*(qx*qz + qw*qy);
    R(2,1) = 2*(qy*qz - qw*qx);
    R(2,2) = 1 - 2*(qx*qx + qy*qy);
    return R;
}

/// Rotate vector v by quaternion q: v_body = R(q) * v.
inline Vec3 quat_rotate(const Vec4& q, const Vec3& v) {
    return quat_to_dcm(q) * v;
}

/// Build quaternion from rotation vector (axis-angle, axis = rv/|rv|, angle = |rv|).
inline Vec4 quat_from_rotvec(const Vec3& rv) {
    double angle = rv.norm();
    if (angle < 1e-14) return Vec4(0, 0, 0, 1);
    Vec3 axis = rv / angle;
    double ha = angle * 0.5;
    double sha = std::sin(ha);
    return Vec4(axis.x() * sha, axis.y() * sha, axis.z() * sha, std::cos(ha));
}

/// Ensure shortest-path (w >= 0).
inline Vec4 quat_ensure_positive_w(const Vec4& q) {
    return (q(3) < 0) ? Vec4(-q) : q;
}

} // namespace gnc
