#pragma once
/// @file ukf.hpp
/// Unscented Kalman Filter for orbit estimation.
/// 7-state: [r(3), v(3), ln_rho_scale(1)]. Measurement: GPS [r(3), v(3)].

#include "math/types.hpp"
#include "math/constants.hpp"
#include "math/quaternion.hpp"
#include "environment/gravity.hpp"
#include "environment/density.hpp"
#include "environment/wind.hpp"
#include "environment/aerodynamics.hpp"
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>

namespace gnc {

// State/measurement dimensions
constexpr int UKF_N = 7;    // [r, v, ln_rho_scale]
constexpr int UKF_M = 6;    // [r, v]  (GPS)
constexpr int UKF_NS = 2 * UKF_N + 1;  // 15 sigma points

using UVec7  = Eigen::Matrix<double, UKF_N, 1>;
using UMat7  = Eigen::Matrix<double, UKF_N, UKF_N>;
using UVec6  = Eigen::Matrix<double, UKF_M, 1>;
using UMat6  = Eigen::Matrix<double, UKF_M, UKF_M>;
using UMat76 = Eigen::Matrix<double, UKF_N, UKF_M>;

// ── Sigma weights ────────────────────────────────────────────────────

struct UKFWeights {
    double wm[UKF_NS];
    double wc[UKF_NS];
    double lambda;

    void compute(double alpha = 1e-3, double beta = 2.0, double kappa = 0.0) {
        lambda = alpha * alpha * (UKF_N + kappa) - UKF_N;
        double w_i = 0.5 / (UKF_N + lambda);
        for (int i = 0; i < UKF_NS; ++i) {
            wm[i] = w_i;
            wc[i] = w_i;
        }
        wm[0] = lambda / (UKF_N + lambda);
        wc[0] = lambda / (UKF_N + lambda) + (1.0 - alpha * alpha + beta);
    }
};

// ── Sigma point generation ───────────────────────────────────────────

// Row-major so that sigmas.row(i).data() returns contiguous memory
using SigmaMatrix = Eigen::Matrix<double, UKF_NS, UKF_N, Eigen::RowMajor>;
using HSigmaMatrix = Eigen::Matrix<double, UKF_NS, UKF_M, Eigen::RowMajor>;

inline void sigma_points(const UVec7& x, const UMat7& P, double lambda,
                         SigmaMatrix& sigmas) {
    UMat7 A = (UKF_N + lambda) * P;
    // Small perturbation for numerical stability
    for (int i = 0; i < UKF_N; ++i) A(i, i) += 1e-12;

    // Cholesky: A = L L^T
    Eigen::LLT<UMat7> llt(A);
    UMat7 L = llt.matrixL();

    sigmas.row(0) = x.transpose();
    for (int i = 0; i < UKF_N; ++i) {
        sigmas.row(i + 1) = (x + L.col(i)).transpose();
        sigmas.row(i + 1 + UKF_N) = (x - L.col(i)).transpose();
    }
}

// ── UKF Process Model ───────────────────────────────────────────────

/// Process model: propagate orbit sigma point by dt.
/// state = [r(3), v(3), ln_rho_scale]
inline Vec3 onboard_gravity_accel(const Vec3& r, double gmst,
                                  int gravity_degree,
                                  const Eigen::MatrixXd& C_nm,
                                  const Eigen::MatrixXd& S_nm,
                                  SHBuf& sh_buf) {
    if (gravity_degree <= 0) {
        const double r_norm = r.norm();
        return -constants::MU_EARTH * r / (r_norm * r_norm * r_norm);
    }
    if (C_nm.rows() == 0 || S_nm.rows() == 0) {
        return j2_accel(r);
    }

    const int coeff_degree = static_cast<int>(std::min(C_nm.rows(), S_nm.rows())) - 1;
    const int max_degree = std::min(gravity_degree, coeff_degree);
    if (max_degree < 2) return j2_accel(r);

    if (sh_buf.n_max < max_degree) sh_buf = SHBuf(max_degree);

    Vec3 a_grav = Vec3::Zero();
    spherical_harmonics_accel(r, max_degree, C_nm, S_nm, gmst, sh_buf, a_grav);
    return a_grav;
}

inline void ukf_process_model(const double* x_in, double dt, double gmst,
                              const Vec4& q,
                              const AeroCoeffs& ac, double onboard_density_scale,
                              int onboard_gravity_degree,
                              const Eigen::MatrixXd& C_nm,
                              const Eigen::MatrixXd& S_nm,
                              SHBuf& sh_buf, double* x_out) {
    Vec3 r(x_in[0], x_in[1], x_in[2]);
    Vec3 v(x_in[3], x_in[4], x_in[5]);
    double ln_rho_scale = x_in[6];

    Vec3 a_grav = onboard_gravity_accel(
        r, gmst, onboard_gravity_degree, C_nm, S_nm, sh_buf);

    // Density (exponential, scaled)
    double rho = onboard_density_scale * std::exp(ln_rho_scale) * exp_density(r);

    // Wind + aero
    Vec3 v_rel_eci = wind_simple(r, v);
    Mat3 R = quat_to_dcm(q);
    Vec3 v_rel_body = R * v_rel_eci;

    Vec3 f_aero, t_aero;
    cll_interp(v_rel_body, rho, ac, f_aero, t_aero);

    // dv = a_grav + R^T f_aero / MASS
    Vec3 f_eci = R.transpose() * f_aero;
    Vec3 a_total = a_grav + f_eci / constants::MASS;

    // Euler forward
    x_out[0] = r.x() + v.x() * dt;
    x_out[1] = r.y() + v.y() * dt;
    x_out[2] = r.z() + v.z() * dt;
    x_out[3] = v.x() + a_total.x() * dt;
    x_out[4] = v.y() + a_total.y() * dt;
    x_out[5] = v.z() + a_total.z() * dt;
    x_out[6] = ln_rho_scale;  // constant
}

// ── UKF measurement model ────────────────────────────────────────────

/// GPS: z = [r, v]
inline UVec6 ukf_measurement(const double* x) {
    UVec6 z;
    for (int i = 0; i < 6; ++i) z(i) = x[i];
    return z;
}

// ── UKF class ────────────────────────────────────────────────────────

struct UKF {
    UVec7 x;
    UMat7 P;
    UMat7 Q;
    UMat6 R;
    UKFWeights w;
    double alpha, beta, kappa;
    double onboard_density_scale = 1.0;
    int onboard_gravity_degree = 2;
    Eigen::MatrixXd onboard_C_nm;
    Eigen::MatrixXd onboard_S_nm;
    SHBuf onboard_sh_buf{0};

    // Preallocated buffers (RowMajor for contiguous row access)
    SigmaMatrix sigmas;
    HSigmaMatrix h_sigmas;

    UKF(const UVec7& x_init, const UMat7& P_init,
        const UMat7& Q_init, const UMat6& R_init,
        double a = 1e-3, double b = 2.0, double k = 0.0)
        : x(x_init), P(P_init), Q(Q_init), R(R_init),
          alpha(a), beta(b), kappa(k) {
        w.compute(alpha, beta, kappa);
    }

    void predict(double dt, double gmst, const Vec4& q, const AeroCoeffs& ac) {
        w.compute(alpha, beta, kappa);
        sigma_points(x, P, w.lambda, sigmas);

        // Propagate each sigma point
        for (int i = 0; i < UKF_NS; ++i) {
            double buf[UKF_N];
            ukf_process_model(sigmas.row(i).data(), dt, gmst, q, ac,
                              onboard_density_scale,
                              onboard_gravity_degree,
                              onboard_C_nm, onboard_S_nm,
                              onboard_sh_buf, buf);
            for (int j = 0; j < UKF_N; ++j) sigmas(i, j) = buf[j];
        }

        // Predicted mean
        UVec7 x_pred = UVec7::Zero();
        for (int i = 0; i < UKF_NS; ++i)
            x_pred += w.wm[i] * sigmas.row(i).transpose();

        // Predicted covariance
        UMat7 P_pred = Q;
        for (int i = 0; i < UKF_NS; ++i) {
            UVec7 dx = sigmas.row(i).transpose() - x_pred;
            P_pred += w.wc[i] * dx * dx.transpose();
        }
        P_pred = 0.5 * (P_pred + P_pred.transpose());

        x = x_pred;
        P = P_pred;
    }

    void predict(double dt, const Vec4& q, const AeroCoeffs& ac) {
        predict(dt, 0.0, q, ac);
    }

    void update(const UVec6& z) {
        // Measurement sigma points
        for (int i = 0; i < UKF_NS; ++i)
            h_sigmas.row(i) = ukf_measurement(sigmas.row(i).data()).transpose();

        // Predicted measurement mean
        UVec6 z_pred = UVec6::Zero();
        for (int i = 0; i < UKF_NS; ++i)
            z_pred += w.wm[i] * h_sigmas.row(i).transpose();

        // Innovation covariance S and cross-covariance Pxz
        UMat6 S = R;
        UMat76 Pxz = UMat76::Zero();
        for (int i = 0; i < UKF_NS; ++i) {
            UVec6 dz = h_sigmas.row(i).transpose() - z_pred;
            S += w.wc[i] * dz * dz.transpose();

            UVec7 dx = sigmas.row(i).transpose() - x;
            Pxz += w.wc[i] * dx * dz.transpose();
        }

        // Kalman gain K = Pxz * S^{-1}
        UMat76 K = Pxz * S.inverse();

        // State + covariance update
        UVec6 innov = z - z_pred;
        x += K * innov;
        UMat7 KS = K * S * K.transpose();
        P -= KS;
        P = 0.5 * (P + P.transpose());
    }
};

} // namespace gnc
