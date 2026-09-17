#pragma once
/// @file constants.hpp
/// Physical and spacecraft constants (all SI).

#include "math/types.hpp"
#include <cmath>
#include <array>

namespace gnc {
namespace constants {

// ── Earth Parameters ─────────────────────────────────────────────────
inline constexpr double MU_EARTH  = 3.986004418e14;   // [m^3/s^2]
inline constexpr double R_EARTH   = 6378137.0;        // [m]
inline constexpr double J2_EARTH  = 1.082626e-3;      // [-]
inline constexpr double W_EARTH   = 7.292115e-5;      // [rad/s]

// ── Solar System ─────────────────────────────────────────────────────
inline constexpr double AU        = 149597870700.0;    // [m]
inline constexpr double MU_SUN    = 1.32712440018e20;  // [m^3/s^2]
inline constexpr double MU_MOON   = 4.9048695e12;      // [m^3/s^2]
inline constexpr double P_SUN     = 4.56e-6;           // [N/m^2] SRP at 1 AU
inline constexpr double R_SUN     = 696340000.0;       // [m]
inline constexpr int    LEAP_SEC  = 37;                // [s]

// ── Thermodynamics / Atomic ──────────────────────────────────────────
inline constexpr double BOLTZMANN = 1.380649e-23;      // [J/K]
inline constexpr double AVOGADRO  = 6.02214076e23;     // [mol^-1]
inline constexpr double GAS_CONST = 8.314462618;       // [J/(mol K)]
inline constexpr double MOLAR_MASS_O = 0.015999;       // [kg/mol]
inline constexpr double AMU       = 1.66053906660e-27; // [kg]

// ── Density (exponential model) ──────────────────────────────────────
inline constexpr double RHO_0     = 2e-11;             // [kg/m^3]
inline constexpr double H_0       = 300000.0;          // [m]
inline constexpr double H_SCALE   = 50000.0;           // [m]

// ── Drag / CLL Parameters ───────────────────────────────────────────
inline constexpr double ALPHA_N   = 0.9;
inline constexpr double SIGMA_T   = 0.85;
inline constexpr int    NUM_PARTICLES = 500;
inline constexpr double T_W       = 300.0;             // [K]
inline constexpr double T_A       = 1000.0;            // [K]

// ── Magnetic Parameters ──────────────────────────────────────────────
inline constexpr double M_MAG     = 7.7e22;            // [A m^2]
inline const     double TILT      = 11.5 * M_PI / 180.0;
inline const     double PHI_POLE  = -70.0 * M_PI / 180.0;

// ── Spacecraft Parameters (3U CubeSat) ──────────────────────────────
inline constexpr double MASS      = 4.0;               // [kg]
inline constexpr double SC_L      = 0.3;               // [m]
inline constexpr double SC_W      = 0.1;
inline constexpr double SC_H      = 0.1;

// Inertia
inline constexpr double Ix = (1.0/12.0) * MASS * (SC_W*SC_W + SC_H*SC_H);
inline constexpr double Iy = (1.0/12.0) * MASS * (SC_L*SC_L + SC_H*SC_H);
inline constexpr double Iz = (1.0/12.0) * MASS * (SC_L*SC_L + SC_W*SC_W);

/// 3x3 inertia matrix (diagonal for box).
inline const Mat3& inertia() {
    static const Mat3 I = (Mat3() << Ix, 0, 0,
                                      0, Iy, 0,
                                      0, 0, Iz).finished();
    return I;
}
/// Inverse inertia.
inline const Mat3& inv_inertia() {
    static const Mat3 Iinv = inertia().inverse();
    return Iinv;
}

// ── Plate Model (6 faces) ────────────────────────────────────────────
// Order: +X, -X, +Y, -Y, +Z, -Z
inline constexpr int N_PLATES = 6;

inline constexpr std::array<double, 6> AREAS = {
    SC_W*SC_H, SC_W*SC_H,
    SC_L*SC_H, SC_L*SC_H,
    SC_L*SC_W, SC_L*SC_W
};

/// Plate outward normals (row-major, 6x3).
inline const Eigen::Matrix<double, 6, 3>& normals() {
    static const Eigen::Matrix<double, 6, 3> N =
        (Eigen::Matrix<double, 6, 3>() <<
          1, 0, 0,  -1, 0, 0,
          0, 1, 0,   0,-1, 0,
          0, 0, 1,   0, 0,-1).finished();
    return N;
}

/// Plate geometric centres (row-major, 6x3).
inline const Eigen::Matrix<double, 6, 3>& centers() {
    static const Eigen::Matrix<double, 6, 3> C =
        (Eigen::Matrix<double, 6, 3>() <<
          SC_L/2, 0, 0,  -SC_L/2, 0, 0,
          0, SC_W/2, 0,   0,-SC_W/2, 0,
          0, 0, SC_H/2,   0, 0,-SC_H/2).finished();
    return C;
}

/// Specular reflectivity per plate.
inline constexpr std::array<double, 6> RS = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
/// Diffuse reflectivity per plate.
inline constexpr std::array<double, 6> RD = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

/// Lever arms from CoM (= centres − CoM; CoM at origin).
inline const Eigen::Matrix<double, 6, 3>& levers() {
    // CoM at origin → levers == centres
    return centers();
}

// ── Residual Magnetic Dipole ─────────────────────────────────────────
inline const Vec3& m_body() {
    static const Vec3 m{0.01, 0.01, 0.01};   // [A m^2]
    return m;
}

} // namespace constants
} // namespace gnc
