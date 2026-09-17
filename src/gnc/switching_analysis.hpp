#pragma once
/// @file switching_analysis.hpp
/// Conservative dwell-time margin calculations for hysteresis switching.

#include "math/types.hpp"
#include <algorithm>

namespace gnc {

struct DwellTimeInputs {
    double omega_lo_rad_s = 5.0 * 3.14159265358979323846 / 180.0;
    double omega_hi_rad_s = 7.0 * 3.14159265358979323846 / 180.0;
    double max_control_torque_nm = 0.002;
    double disturbance_torque_bound_nm = 2e-6;
    Mat3 inertia = Mat3::Identity();
};

struct DwellTimeResult {
    double hysteresis_gap_rad_s = 0.0;
    double angular_accel_bound_rad_s2 = 0.0;
    double required_dwell_s = 0.0;
    double observed_dwell_s = 0.0;
    double margin_s = 0.0;
};

inline DwellTimeResult compute_dwell_time_margin(const DwellTimeInputs& in,
                                                 double observed_dwell_s) {
    DwellTimeResult out;
    out.hysteresis_gap_rad_s = std::max(0.0, in.omega_hi_rad_s - in.omega_lo_rad_s);

    const double min_inertia = std::max(1e-12, in.inertia.eigenvalues().real().minCoeff());
    const double torque_bound =
        std::max(0.0, in.max_control_torque_nm) +
        std::max(0.0, in.disturbance_torque_bound_nm);
    out.angular_accel_bound_rad_s2 = torque_bound / min_inertia;

    if (out.angular_accel_bound_rad_s2 > 0.0) {
        out.required_dwell_s = out.hysteresis_gap_rad_s / out.angular_accel_bound_rad_s2;
    }
    out.observed_dwell_s = observed_dwell_s;
    out.margin_s = observed_dwell_s - out.required_dwell_s;
    return out;
}

} // namespace gnc
