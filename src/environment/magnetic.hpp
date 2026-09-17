#pragma once
/// @file magnetic.hpp
/// Magnetic field models: dipole and IGRF LUT.

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/frames.hpp"
#include "math/interp3d.hpp"
#include <cmath>
#include <vector>
#include <string>

namespace gnc {

// ── Magnetic torque ──────────────────────────────────────────────────

inline Vec3 B_torque(const Vec3& m_body, const Vec3& B_body) {
    return m_body.cross(B_body);
}

// ── Dipole Magnetic Field (ECI) ──────────────────────────────────────

inline Vec3 dipole_B(const Vec3& r_eci, double gmst) {
    Vec3 r_ecef = eci_to_ecef(r_eci, gmst);
    double r_norm = r_ecef.norm();
    if (r_norm < constants::R_EARTH) return Vec3::Zero();

    double tilt = constants::TILT;
    double phi_pole = constants::PHI_POLE;
    double m_mag = constants::M_MAG;

    double mx = -m_mag * std::sin(tilt) * std::cos(phi_pole);
    double my = -m_mag * std::sin(tilt) * std::sin(phi_pole);
    double mz = -m_mag * std::cos(tilt);

    Vec3 u = r_ecef / r_norm;
    double m_dot_u = mx*u.x() + my*u.y() + mz*u.z();
    double factor = 1e-7 / (r_norm * r_norm * r_norm);

    Vec3 B_ecef(factor * (3.0*m_dot_u*u.x() - mx),
                factor * (3.0*m_dot_u*u.y() - my),
                factor * (3.0*m_dot_u*u.z() - mz));

    return ecef_to_eci(B_ecef, gmst);
}

// ── IGRF LUT ─────────────────────────────────────────────────────────

struct IgrfLUT {
    std::vector<double> alts_km, lats_deg, lons_deg;
    std::vector<double> be, bn, bu;  // flat [lat x lon x alt]
    int n_lat = 0, n_lon = 0, n_alt = 0;
};

/// Trilinear-interpolated IGRF magnetic field → ECI.
inline Vec3 igrf_B(const Vec3& r_eci, double gmst, const IgrfLUT& lut) {
    Vec3 r_ecef = eci_to_ecef(r_eci, gmst);
    auto geo = ecef_to_geodetic(r_ecef.x()/1000.0, r_ecef.y()/1000.0, r_ecef.z()/1000.0);

    // LUT axes: (lat, lon, alt)
    auto a0 = axis_lookup(lut.lats_deg.data(), lut.n_lat, geo.lat_deg);
    auto a1 = axis_lookup_periodic(lut.lons_deg.data(), lut.n_lon, geo.lon_deg);
    auto a2 = axis_lookup(lut.alts_km.data(),  lut.n_alt, geo.alt_km);

    double Be, Bn, Bu;
    trilinear_vec3(lut.be.data(), lut.bn.data(), lut.bu.data(),
                   lut.n_lat, lut.n_lon, lut.n_alt,
                   a0.idx, a1.idx, a2.idx,
                   a0.w, a1.w, a2.w,
                   Be, Bn, Bu,
                   false, true, false); // periodic_dim1

    // ENU → ECEF
    double lat_rad = geo.lat_deg * M_PI / 180.0;
    double lon_rad = geo.lon_deg * M_PI / 180.0;
    double sP = std::sin(lat_rad), cP = std::cos(lat_rad);
    double sL = std::sin(lon_rad), cL = std::cos(lon_rad);

    Vec3 B_ecef(Bn*(-sP*cL) + Be*(-sL) + Bu*(cP*cL),
                Bn*(-sP*sL) + Be*(cL)   + Bu*(cP*sL),
                Bn*(cP)     + Be*(0.0)  + Bu*(sP));

    return ecef_to_eci(B_ecef, gmst);
}

/// Load IGRF LUT from binary file.
IgrfLUT load_igrf_lut(const std::string& bin_path);

} // namespace gnc
