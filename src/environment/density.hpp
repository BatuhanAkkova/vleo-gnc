#pragma once
/// @file density.hpp
/// Atmospheric density models: exponential and NRLMSISE-00 LUT.

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/frames.hpp"
#include "math/interp3d.hpp"
#include <cmath>
#include <vector>
#include <string>

namespace gnc {

// ── Exponential density model ────────────────────────────────────────

inline double exp_density(const Vec3& r_eci) {
    double r_mag = r_eci.norm();
    double h = r_mag - constants::R_EARTH;
    if (h < 0.0) return constants::RHO_0;
    return constants::RHO_0 * std::exp(-(h - constants::H_0) / constants::H_SCALE);
}

// ── NRLMSISE-00 LUT ─────────────────────────────────────────────────

struct DensityLUT {
    std::vector<double> alts_km;   // sorted altitude axis [km]
    std::vector<double> lats_deg;  // sorted latitude axis [deg]
    std::vector<double> lons_deg;  // sorted longitude axis [deg]
    std::vector<double> rho_grid;  // flat row-major [alt x lat x lon]
    int n_alt = 0, n_lat = 0, n_lon = 0;
};

/// Trilinear-interpolated NRLMSISE density.
inline double nrlmsise_density(const Vec3& r_eci, double gmst,
                               const DensityLUT& lut) {
    Vec3 r_ecef = eci_to_ecef(r_eci, gmst);
    auto geo = ecef_to_geodetic(r_ecef.x()/1000.0, r_ecef.y()/1000.0, r_ecef.z()/1000.0);

    auto a0 = axis_lookup(lut.alts_km.data(), lut.n_alt, geo.alt_km);
    auto a1 = axis_lookup(lut.lats_deg.data(), lut.n_lat, geo.lat_deg);
    auto a2 = axis_lookup_periodic(lut.lons_deg.data(), lut.n_lon, geo.lon_deg);

    return trilinear(lut.rho_grid.data(),
                     lut.n_alt, lut.n_lat, lut.n_lon,
                     a0.idx, a1.idx, a2.idx,
                     a0.w, a1.w, a2.w,
                     true); // periodic_dim2
}

/// Load NRLMSISE LUT from binary file.
DensityLUT load_density_lut(const std::string& bin_path);

} // namespace gnc
