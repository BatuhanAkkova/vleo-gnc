#pragma once
/// @file actuators.hpp
/// Actuator models: ReactionWheel (per-axis), Magnetorquer.

#include <cmath>
#include <algorithm>
#include <random>

namespace gnc {

// ── Reaction Wheel (single axis) ─────────────────────────────────────

struct ReactionWheel {
    double max_torque;
    double max_speed;   // rad/s
    double inertia;     // kg*m^2
    double jitter;      // torque noise σ [Nm]
    double speed = 0.0; // current wheel speed [rad/s]

    ReactionWheel(double mt = 0.002, double ms = 150.0,
                  double I = 1e-4, double j = 0.0)
        : max_torque(mt), max_speed(ms), inertia(I), jitter(j) {}

    /// Returns actual torque applied TO THE BODY (reaction torque).
    double update(double torque_cmd, double dt, std::mt19937& rng) {
        double t_app = std::clamp(torque_cmd, -max_torque, max_torque);

        // Saturation
        if (speed >= max_speed && t_app > 0.0)  t_app = 0.0;
        if (speed <= -max_speed && t_app < 0.0) t_app = 0.0;

        // Jitter
        double t_actual = t_app;
        if (jitter > 0.0) {
            std::normal_distribution<double> nd(0.0, jitter);
            t_actual += nd(rng);
        }

        // Update wheel state
        speed += (t_actual / inertia) * dt;

        return -t_actual;  // reaction torque on body
    }
};

// ── Magnetorquer ─────────────────────────────────────────────────────

struct Magnetorquer {
    double max_dipole;

    explicit Magnetorquer(double md = 0.3) : max_dipole(md) {}

    /// Clamp each component independently.
    void get_dipole(const double m_cmd[3], double m_out[3]) const {
        for (int i = 0; i < 3; ++i)
            m_out[i] = std::clamp(m_cmd[i], -max_dipole, max_dipole);
    }
};

} // namespace gnc
