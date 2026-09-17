#pragma once
/// @file frames.hpp
/// ECI - ECEF - geodetic coordinate transforms.

#include "math/types.hpp"
#include <cmath>

namespace gnc {

/// Rotate an ECI vector to ECEF given Greenwich Mean Sidereal Time [rad].
inline Vec3 eci_to_ecef(const Vec3 &r_eci, double gmst) {
  double c = std::cos(gmst), s = std::sin(gmst);
  return Vec3(c * r_eci.x() + s * r_eci.y(), -s * r_eci.x() + c * r_eci.y(),
              r_eci.z());
}

/// Rotate an ECEF vector to ECI.
inline Vec3 ecef_to_eci(const Vec3 &r_ecef, double gmst) {
  double c = std::cos(gmst), s = std::sin(gmst);
  return Vec3(c * r_ecef.x() - s * r_ecef.y(), s * r_ecef.x() + c * r_ecef.y(),
              r_ecef.z());
}

/// Convert ECEF position [km] to geodetic (lat_deg, lon_deg, alt_km).
struct Geodetic {
  double lat_deg;
  double lon_deg;
  double alt_km;
};

inline Geodetic ecef_to_geodetic(double x_km, double y_km, double z_km) {
  constexpr double R_E_km = 6378.137;
  double r = std::sqrt(x_km * x_km + y_km * y_km + z_km * z_km);
  double lat = std::asin(z_km / std::max(r, 1e-12)) * (180.0 / M_PI);
  double lon = std::atan2(y_km, x_km) * (180.0 / M_PI);
  lon = std::atan2(std::sin(lon), std::cos(lon));
  double alt = r - R_E_km;
  return {lat, lon, alt};
}

} // namespace gnc
