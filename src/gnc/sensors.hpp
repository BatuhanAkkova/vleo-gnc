#pragma once
/// @file sensors.hpp
/// Sensor models: Gyroscope, Magnetometer, SunSensor, GPS, StarTracker.

#include "math/types.hpp"
#include "math/quaternion.hpp"
#include <cmath>
#include <random>

namespace gnc {

// ── Gyroscope ────────────────────────────────────────────────────────

struct Gyroscope {
    Vec3 bias = Vec3::Zero();
    double rrw_std;   // rate random walk std [rad/s/sqrt(s)]
    double arw_std;   // angle random walk std [rad/sqrt(s)]
    double scale_err; // scale factor error (fraction)
    double dt;

    Gyroscope(double rrw = 1e-5, double arw = 1e-4, double sc = 0.0, double _dt = 0.1)
        : rrw_std(rrw), arw_std(arw), scale_err(sc), dt(_dt) {}

    Vec3 update(const Vec3& true_omega, std::mt19937& rng) {
        std::normal_distribution<double> nd(0.0, 1.0);
        // Bias random walk
        double bw = rrw_std * std::sqrt(dt);
        bias.x() += nd(rng) * bw;
        bias.y() += nd(rng) * bw;
        bias.z() += nd(rng) * bw;
        // ARW noise
        double rn = arw_std / std::sqrt(dt);
        Vec3 noise(nd(rng)*rn, nd(rng)*rn, nd(rng)*rn);
        // Scale factor
        return (1.0 + scale_err) * true_omega + bias + noise;
    }
};

// ── Magnetometer ─────────────────────────────────────────────────────

struct Magnetometer {
    Vec3 bias;
    double noise_std;

    Magnetometer(double ns = 0.2e-6, double bias_std = 50e-9, std::mt19937* rng = nullptr)
        : noise_std(ns) {
        if (rng) {
            std::normal_distribution<double> nd(0.0, bias_std);
            bias = Vec3(nd(*rng), nd(*rng), nd(*rng));
        } else {
            bias = Vec3::Zero();
        }
    }

    Vec3 update(const Vec3& true_B, std::mt19937& rng) {
        std::normal_distribution<double> nd(0.0, noise_std);
        return true_B + bias + Vec3(nd(rng), nd(rng), nd(rng));
    }
};

// ── Sun Sensor ───────────────────────────────────────────────────────

struct SunSensor {
    double noise_std_rad;

    explicit SunSensor(double noise_deg = 0.1)
        : noise_std_rad(noise_deg * M_PI / 180.0) {}

    /// Returns true + measured unit vector, or false if eclipse.
    bool update(const Vec3& true_s_body, Vec3& out, std::mt19937& rng) {
        double s_norm = true_s_body.norm();
        if (s_norm < 1e-9) return false;

        Vec3 u = true_s_body / s_norm;
        std::normal_distribution<double> nd(0.0, noise_std_rad);
        Vec3 noise(nd(rng), nd(rng), nd(rng));
        out = (u + noise).normalized();
        return true;
    }
};

// ── GPS ──────────────────────────────────────────────────────────────

struct GPS {
    double pos_std;
    double vel_std;

    GPS(double ps = 10.0, double vs = 0.1) : pos_std(ps), vel_std(vs) {}

    void update(const Vec3& true_pos, const Vec3& true_vel,
                Vec3& meas_pos, Vec3& meas_vel, std::mt19937& rng) {
        std::normal_distribution<double> ndp(0.0, pos_std);
        std::normal_distribution<double> ndv(0.0, vel_std);
        meas_pos = true_pos + Vec3(ndp(rng), ndp(rng), ndp(rng));
        meas_vel = true_vel + Vec3(ndv(rng), ndv(rng), ndv(rng));
    }
};

// ── Star Tracker ─────────────────────────────────────────────────────

struct StarTracker {
    double noise_rad;
    double max_rate_rad;
    Mat3 R_matrix;  // measurement covariance

    StarTracker(double noise_deg = 0.005, double max_rate_deg = 1.0)
        : noise_rad(noise_deg * M_PI / 180.0),
          max_rate_rad(max_rate_deg * M_PI / 180.0),
          R_matrix(Mat3::Identity() * (noise_rad * noise_rad)) {}

    /// Returns (is_valid, q_meas). q_meas valid only if is_valid == true.
    bool read(const Vec4& q_truth, const Vec3& omega_truth,
              Vec4& q_meas, std::mt19937& rng) {
        if (omega_truth.norm() > max_rate_rad) return false;

        // Small rotation noise
        std::normal_distribution<double> nd(0.0, noise_rad);
        Vec3 noise_vec(nd(rng), nd(rng), nd(rng));
        Vec4 q_noise = quat_from_rotvec(noise_vec);

        // Apply: q_meas = q_noise * q_truth
        q_meas = quat_product(q_noise, q_truth);
        q_meas = quat_normalize(q_meas);

        // Shortest path
        if (q_meas.w() < 0.0) q_meas = -q_meas;
        return true;
    }
};

} // namespace gnc
