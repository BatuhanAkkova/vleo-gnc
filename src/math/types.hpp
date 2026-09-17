#pragma once
/// @file types.hpp
/// Fixed-size Eigen type aliases used across the GNC simulation.

#include <Eigen/Dense>

namespace gnc {

// ── Vector typedefs ──────────────────────────────────────────────────
using Vec3   = Eigen::Vector3d;
using Vec4   = Eigen::Vector4d; // quaternion storage [x,y,z,w]
using Vec6   = Eigen::Matrix<double, 6, 1>;
using Vec7   = Eigen::Matrix<double, 7, 1>;
using Vec13  = Eigen::Matrix<double, 13, 1>;
using Vec15  = Eigen::Matrix<double, 15, 1>; // 2*n+1 UKF weights (n=7)

// ── Matrix typedefs ──────────────────────────────────────────────────
using Mat3   = Eigen::Matrix3d;
using Mat6   = Eigen::Matrix<double, 6, 6>;
using Mat63  = Eigen::Matrix<double, 6, 3>;
using Mat36  = Eigen::Matrix<double, 3, 6>;
using Mat7   = Eigen::Matrix<double, 7, 7>;
using MatX   = Eigen::MatrixXd; // dynamic-size fallback

// ── Quaternion ───────────────────────────────────────────────────────
// Eigen stores coeffs as [x,y,z,w] internally (scalar-last).
// Constructor takes (w,x,y,z).  Use coeffs() for [x,y,z,w].
using Quat   = Eigen::Quaterniond;

} // namespace gnc
