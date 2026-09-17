#pragma once
/// @file mekf.hpp
/// Multiplicative Extended Kalman Filter for attitude estimation.

#include "math/types.hpp"
#include "math/quaternion.hpp"
#include <cmath>

namespace gnc {

// Vec6, Mat6, Mat36, Mat63 already defined in types.hpp

// ── MEKF Predict ─────────────────────────────────────────────────────

inline void mekf_predict(Vec4& q, const Vec3& beta, const Vec3& omega_meas,
                         Mat6& P, const Mat6& Q, double dt) {
    // Corrected angular velocity
    Vec3 w = omega_meas - beta;

    // Propagate quaternion: q_new = q x deltaq(w*dt)
    Vec3 dtheta = w * dt;
    Vec4 dq = quat_from_rotvec(dtheta);
    q = quat_product(q, dq);
    q = quat_normalize(q);

    // State transition Jacobian F (6x6)
    //   F = [ -[wx]  -I ]
    //       [   0     0  ]
    Mat6 F = Mat6::Zero();
    F(0,1) =  w.z(); F(0,2) = -w.y();
    F(1,0) = -w.z(); F(1,2) =  w.x();
    F(2,0) =  w.y(); F(2,1) = -w.x();
    F(0,3) = -1.0; F(1,4) = -1.0; F(2,5) = -1.0;

    // Discrete transition: phi = I + F*dt
    Mat6 Phi = Mat6::Identity() + F * dt;

    // Covariance prediction: P = phi P phi^T + Q*dt
    P = Phi * P * Phi.transpose() + Q * dt;

    // Symmetrise
    P = 0.5 * (P + P.transpose());
}

// ── MEKF Update (vector measurement: magnetometer, sun sensor) ───────

inline void mekf_update(Vec4& q, Vec3& beta, Mat6& P,
                        const Vec3& z_meas, const Vec3& z_ref,
                        const Mat3& R) {
    // Predicted measurement: z_pred = DCM(q) * z_ref (ECI -> body)
    Mat3 Rq = quat_to_dcm(q);
    Vec3 z_pred = Rq * z_ref;

    // Innovation
    Vec3 y = z_meas - z_pred;

    // Observation Jacobian: H = [ skew(z_pred) | phi ]
    Mat36 H = Mat36::Zero();
    H(0,1) = -z_pred.z(); H(0,2) =  z_pred.y();
    H(1,0) =  z_pred.z(); H(1,2) = -z_pred.x();
    H(2,0) = -z_pred.y(); H(2,1) =  z_pred.x();

    // Innovation covariance: S = H P H^T + R
    Mat3 S = H * P * H.transpose() + R;

    // Kalman gain: K = P H^T S^-1
    Mat63 K = P * H.transpose() * S.inverse();

    // State correction
    Vec6 dx = K * y;

    // Quaternion update: q+ = q cross [1/2 detlaphi; 1]
    Vec4 dq;
    dq.x() = 0.5 * dx(0);
    dq.y() = 0.5 * dx(1);
    dq.z() = 0.5 * dx(2);
    dq.w() = 1.0;
    q = quat_normalize(quat_product(q, dq));

    // Bias update
    beta += dx.tail<3>();

    // Covariance update (Joseph form): P = (I-KH) P (I-KH)^T + K R K^T
    Mat6 IKH = Mat6::Identity() - K * H;
    P = IKH * P * IKH.transpose() + K * R * K.transpose();
    P = 0.5 * (P + P.transpose());
}

// ── MEKF Update (star tracker quaternion measurement) ────────────────

inline void mekf_update_star(Vec4& q, Vec3& beta, Mat6& P,
                             const Vec4& q_meas, const Mat3& R_st) {
    // Error quaternion: dq = q^-1 * q_meas
    Vec4 q_inv = quat_conjugate(q);
    Vec4 dq_err = quat_product(q_inv, q_meas);

    // Shortest path
    if (dq_err.w() < 0.0) dq_err = -dq_err;

    // Innovation: y = 2 * dq_vec
    Vec3 y(2.0*dq_err.x(), 2.0*dq_err.y(), 2.0*dq_err.z());

    // H = [I_3 | phi_3x3]
    Mat36 H = Mat36::Zero();
    H(0,0) = 1.0; H(1,1) = 1.0; H(2,2) = 1.0;

    // S = P_att + R_st
    Mat3 S = P.block<3,3>(0,0) + R_st;

    // K = P(:,0:3) * S^-1
    Mat63 PH = P.block<6,3>(0,0);
    Mat63 K = PH * S.inverse();

    // Correction
    Vec6 dx = K * y;

    // Quaternion update (right multiplication)
    Vec4 dq;
    dq.x() = 0.5 * dx(0);
    dq.y() = 0.5 * dx(1);
    dq.z() = 0.5 * dx(2);
    dq.w() = 1.0;
    q = quat_normalize(quat_product(q, dq));

    // Bias update
    beta += dx.tail<3>();

    // Joseph form covariance
    Mat6 IKH = Mat6::Identity() - K * H;
    P = IKH * P * IKH.transpose() + K * R_st * K.transpose();
    P = 0.5 * (P + P.transpose());
}

// ── MEKF class wrapper ───────────────────────────────────────────────

struct MEKF {
    Vec4 q;
    Vec3 beta = Vec3::Zero();
    Mat6 P;
    Mat6 Q;
    Mat3 R;
    Mat3 R_st;

    MEKF(const Vec4& q_init, const Mat6& P_init,
         const Mat6& Q_proc, const Mat3& R_meas,
         const Mat3& R_star = Mat3::Identity() * 1e-9)
        : q(q_init), P(P_init), Q(Q_proc), R(R_meas), R_st(R_star) {}

    void predict(const Vec3& omega_meas, double dt) {
        mekf_predict(q, beta, omega_meas, P, Q, dt);
    }

    void update(const Vec3& z_meas, const Vec3& z_ref) {
        mekf_update(q, beta, P, z_meas, z_ref, R);
    }

    void update_star_tracker(const Vec4& q_meas) {
        mekf_update_star(q, beta, P, q_meas, R_st);
    }
};

} // namespace gnc
