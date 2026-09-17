#pragma once
/// @file wind.hpp
/// Wind models: simple co-rotation and HWM14 LUT.

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/frames.hpp"
#include "math/interp3d.hpp"
#include <cmath>
#include <vector>
#include <string>

namespace gnc {

// ── Simple co-rotation wind ──────────────────────────────────────────

inline Vec3 wind_simple(const Vec3& r, const Vec3& v) {
    constexpr double we = constants::W_EARTH;
    Vec3 v_atm(-we * r.y(), we * r.x(), 0.0);
    return v - v_atm;
}

// ── HWM14 LUT ────────────────────────────────────────────────────────

struct WindLUT {
    std::vector<double> alts_km, lats_deg, lons_deg;
    std::vector<double> u_zonal, u_merid;  // flat [alt x lat x lon]
    int n_alt = 0, n_lat = 0, n_lon = 0;
};

/// HWM14-style wind: co-rotation + LUT horizontal wind.
inline Vec3 wind_hwm14(const Vec3& r_eci, const Vec3& v_eci,
                       double gmst, const WindLUT& lut) {
    constexpr double we = constants::W_EARTH;

    // Co-rotation
    Vec3 v_atm(-we * r_eci.y(), we * r_eci.x(), 0.0);

    // Geodetic lookup
    Vec3 r_ecef = eci_to_ecef(r_eci, gmst);
    double ex = r_ecef.x() / 1000.0;
    double ey = r_ecef.y() / 1000.0;
    double ez = r_ecef.z() / 1000.0;
    if (std::abs(ex) < 1e-6 && std::abs(ey) < 1e-6) ex = 1e-6;
    auto geo = ecef_to_geodetic(ex, ey, ez);

    double lat_deg = geo.lat_deg;
    double lon_deg = std::isfinite(geo.lon_deg) ? geo.lon_deg : 0.0;
    double alt_km  = std::isfinite(geo.alt_km)  ? geo.alt_km  : 300.0;

    // Trilinear interp both components
    auto a0 = axis_lookup(lut.alts_km.data(),  lut.n_alt, alt_km);
    auto a1 = axis_lookup(lut.lats_deg.data(), lut.n_lat, lat_deg);
    auto a2 = axis_lookup_periodic(lut.lons_deg.data(), lut.n_lon, lon_deg);

    double uz = trilinear(lut.u_zonal.data(), lut.n_alt, lut.n_lat, lut.n_lon,
                          a0.idx, a1.idx, a2.idx, a0.w, a1.w, a2.w,
                          false, false, true); // periodic_dim2
    double um = trilinear(lut.u_merid.data(), lut.n_alt, lut.n_lat, lut.n_lon,
                          a0.idx, a1.idx, a2.idx, a0.w, a1.w, a2.w,
                          false, false, true); // periodic_dim2

    // ENU → ECEF
    double lat_rad = lat_deg * M_PI / 180.0;
    double lon_rad = lon_deg * M_PI / 180.0;
    double sP = std::sin(lat_rad), cP = std::cos(lat_rad);
    double sL = std::sin(lon_rad), cL = std::cos(lon_rad);

    Vec3 e_E(-sL, cL, 0.0);
    Vec3 e_N(-sP*cL, -sP*sL, cP);
    Vec3 w_ecef = uz * e_E + um * e_N;

    // ECEF → ECI
    Vec3 w_eci = ecef_to_eci(w_ecef, gmst);
    v_atm += w_eci;

    return v_eci - v_atm;
}

/// Load HWM14 LUT from binary file.
WindLUT load_wind_lut(const std::string& bin_path);

} // namespace gnc
