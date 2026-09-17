#include "gnc/controller.hpp"
#include "math/constants.hpp"
#include <Eigen/Dense>
// Eigen's unsupported CARE solver
#include <unsupported/Eigen/MatrixFunctions>
#include <stdexcept>

namespace gnc {

// Manual CARE solver via Schur decomposition of the Hamiltonian.
// Solves A^T P + P A - P B R^{-1} B^T P + Q = 0.
static Eigen::MatrixXd solve_care(const Eigen::MatrixXd& A,
                                   const Eigen::MatrixXd& B,
                                   const Eigen::MatrixXd& Q,
                                   const Eigen::MatrixXd& R) {
    int n = static_cast<int>(A.rows());
    Eigen::MatrixXd R_inv = R.inverse();
    Eigen::MatrixXd BR = B * R_inv * B.transpose();

    // Build 2n x 2n Hamiltonian matrix
    Eigen::MatrixXd H(2*n, 2*n);
    H.topLeftCorner(n, n) = A;
    H.topRightCorner(n, n) = -BR;
    H.bottomLeftCorner(n, n) = -Q;
    H.bottomRightCorner(n, n) = -A.transpose();

    // Eigendecomposition of the Hamiltonian
    Eigen::EigenSolver<Eigen::MatrixXd> es(H);
    auto eigenvalues = es.eigenvalues();
    auto eigenvectors = es.eigenvectors();

    // Select stable eigenvectors (Re(lambda) < 0)
    Eigen::MatrixXcd U(2*n, n);
    int col = 0;
    for (int i = 0; i < 2*n && col < n; ++i) {
        if (eigenvalues(i).real() < 0.0) {
            U.col(col) = eigenvectors.col(i);
            ++col;
        }
    }
    if (col < n)
        throw std::runtime_error("CARE: insufficient stable eigenvalues");

    // P = U21 * U11^{-1}
    Eigen::MatrixXcd U11 = U.topRows(n);
    Eigen::MatrixXcd U21 = U.bottomRows(n);
    Eigen::MatrixXcd P_complex = U21 * U11.inverse();

    return P_complex.real();
}

LQRController::LQRController(const Mat3& inertia, double q_w, double w_w, double u_w)
    : I_sc(inertia) {
    Mat3 inv_I = inertia.inverse();

    // Linearized A, B
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6, 6);
    A.block<3,3>(0,3) = Mat3::Identity();

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 3);
    B.block<3,3>(3,0) = inv_I;

    // Q, R weights
    Eigen::MatrixXd Q_mat = Eigen::MatrixXd::Zero(6, 6);
    Q_mat(0,0) = q_w; Q_mat(1,1) = q_w; Q_mat(2,2) = q_w;
    Q_mat(3,3) = w_w; Q_mat(4,4) = w_w; Q_mat(5,5) = w_w;

    Eigen::MatrixXd R_mat = Eigen::MatrixXd::Identity(3, 3) * u_w;

    // Solve CARE
    Eigen::MatrixXd P = solve_care(A, B, Q_mat, R_mat);
    K = (R_mat.inverse() * B.transpose() * P).block<3,6>(0,0);
}

std::pair<Vec3, Vec3> LQRController::compute(
    const Vec4& q_est, const Vec3& w_est,
    const Vec4& q_target, const Vec3& w_target,
    const Vec3* B_body, const Vec3* h_sys,
    const Vec3* t_ff_in,
    double max_torque, double k_desat, double max_dip) const {

    const Vec3& ff = t_ff_in ? *t_ff_in : t_ff;
    Vec3 u_rw = lqr_control(q_est, w_est, q_target, w_target, K, ff, max_torque);

    Vec3 m_mtq = Vec3::Zero();
    if (B_body && h_sys) {
        Vec3 m_cmd = momentum_desaturation(*h_sys, *B_body, k_desat);
        m_mtq = limit_dipole(m_cmd, max_dip);
    }
    return {u_rw, m_mtq};
}

} // namespace gnc
