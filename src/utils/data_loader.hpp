#pragma once
/// @file data_loader.hpp
/// Binary LUT and JSON data loaders for environment models.

#include "environment/density.hpp"
#include "environment/magnetic.hpp"
#include "environment/wind.hpp"
#include "environment/aerodynamics.hpp"
#include <string>

namespace gnc {

// These are declared in their respective headers, implemented in data_loader.cpp:
// DensityLUT  load_density_lut(const std::string& bin_path);
// IgrfLUT     load_igrf_lut(const std::string& bin_path);
// WindLUT     load_wind_lut(const std::string& bin_path);
// AeroCoeffs  load_aero_coeffs(const std::string& cn_path, const std::string& ct_path);

} // namespace gnc
