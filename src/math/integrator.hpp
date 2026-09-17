#pragma once
/// @file integrator.hpp
/// Fixed-size RK4 integrator for 13-state dynamics (r, v, q, w).

#include "math/types.hpp"
#include <cmath>

namespace gnc {

/// Generic RK4 step for a 13-element state vector.
///
/// @tparam RhsFunc  callable (const Vec13& state, Vec13& dst) writing the
///                  derivative into dst.  Any captures are allowed.
/// @param state     current state
/// @param dt        time step [s]
/// @param rhs       right-hand-side functor
/// @param k1..k4    preallocated scratch buffers (13-element each)
/// @param tmp       preallocated scratch for intermediate state
/// @param next      output: state at t + dt
template <typename RhsFunc>
inline void rk4_step_13(const Vec13& state, double dt,
                        RhsFunc&& rhs,
                        Vec13& k1, Vec13& k2, Vec13& k3, Vec13& k4,
                        Vec13& tmp, Vec13& next) {
    rhs(state, k1);

    tmp = state + 0.5 * dt * k1;
    rhs(tmp, k2);

    tmp = state + 0.5 * dt * k2;
    rhs(tmp, k3);

    tmp = state + dt * k3;
    rhs(tmp, k4);

    next = state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    // Renormalise quaternion (indices 6..9)
    double qn = 0.0;
    for (int i = 6; i < 10; ++i) qn += next(i) * next(i);
    qn = std::sqrt(qn);
    if (qn > 1e-12) {
        for (int i = 6; i < 10; ++i) next(i) /= qn;
    }
}

} // namespace gnc
