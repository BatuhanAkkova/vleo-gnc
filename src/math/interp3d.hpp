#pragma once
/// @file interp3d.hpp
/// Trilinear interpolation helpers for 3-D lookup tables.

#include "math/math_utils.hpp"
#include <cmath>
#include <algorithm>

namespace gnc {

/// Result of axis lookup: index + fractional weight.
struct AxisLookup { int idx; double w; };

/// Find the lower index and fractional weight for one axis.
/// Clamps to the table bounds.
template <typename ArrayLike>
inline AxisLookup axis_lookup(const ArrayLike& axis, int n, double val) {
    if (val <= axis[0])    return {0, 0.0};
    if (val >= axis[n - 1]) return {n - 2, 1.0};
    int i = binary_search(axis, val, n) - 1;
    i = std::max(0, std::min(i, n - 2));
    double w = (val - axis[i]) / (axis[i + 1] - axis[i]);
    return {i, w};
}

/// Find the lower index and fractional weight for a periodic axis (e.g. longitude).
/// Assumes the axis covers a 360-degree period but might not include the wrap-around point.
template <typename ArrayLike>
inline AxisLookup axis_lookup_periodic(const ArrayLike& axis, int n, double val, double period = 360.0) {
    // Wrap val to [axis[0], axis[0] + period)
    double v = std::fmod(val - axis[0], period);
    if (v < 0) v += period;
    v += axis[0];

    if (v >= axis[n - 1]) {
        double gap = (axis[0] + period) - axis[n - 1];
        return {n - 1, (v - axis[n - 1]) / gap};
    }
    int i = binary_search(axis, v, n) - 1;
    i = std::max(0, std::min(i, n - 2));
    double w = (v - axis[i]) / (axis[i + 1] - axis[i]);
    return {i, w};
}

/// Scalar trilinear interpolation.
/// Grid shape: [dim0 x dim1 x dim2], stored row-major as a flat array.
/// idx order: (i, j, k) for the three axes.
/// Flags p0, p1, p2 enable periodic wrapping for the respective dimensions.
inline double trilinear(const double* grid,
                        int n0, int n1, int n2,
                        int i, int j, int k,
                        double w0, double w1, double w2,
                        bool p0 = false, bool p1 = false, bool p2 = false) {
    auto G = [&](int a, int b, int c) -> double {
        if (p0) a = (a + n0) % n0;
        if (p1) b = (b + n1) % n1;
        if (p2) c = (c + n2) % n2;
        return grid[a * n1 * n2 + b * n2 + c];
    };

    double v000 = G(i,   j,   k  );
    double v100 = G(i+1, j,   k  );
    double v010 = G(i,   j+1, k  );
    double v001 = G(i,   j,   k+1);
    double v110 = G(i+1, j+1, k  );
    double v101 = G(i+1, j,   k+1);
    double v011 = G(i,   j+1, k+1);
    double v111 = G(i+1, j+1, k+1);

    // Along axis-0
    double v00 = v000 * (1 - w0) + v100 * w0;
    double v01 = v001 * (1 - w0) + v101 * w0;
    double v10 = v010 * (1 - w0) + v110 * w0;
    double v11 = v011 * (1 - w0) + v111 * w0;

    // Along axis-1
    double v0 = v00 * (1 - w1) + v10 * w1;
    double v1 = v01 * (1 - w1) + v11 * w1;

    // Along axis-2
    return v0 * (1 - w2) + v1 * w2;
}

/// Three-component vector trilinear interpolation (e.g. magnetic field B(E,N,U)).
/// comp0/comp1/comp2 are separate flat grids of the same shape.
inline void trilinear_vec3(const double* comp0,
                           const double* comp1,
                           const double* comp2,
                           int n0, int n1, int n2,
                           int i, int j, int k,
                           double w0, double w1, double w2,
                           double& out0, double& out1, double& out2,
                           bool p0 = false, bool p1 = false, bool p2 = false) {
    out0 = trilinear(comp0, n0, n1, n2, i, j, k, w0, w1, w2, p0, p1, p2);
    out1 = trilinear(comp1, n0, n1, n2, i, j, k, w0, w1, w2, p0, p1, p2);
    out2 = trilinear(comp2, n0, n1, n2, i, j, k, w0, w1, w2, p0, p1, p2);
}

} // namespace gnc
