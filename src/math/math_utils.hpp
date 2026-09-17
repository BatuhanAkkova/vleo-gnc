#pragma once
/// @file math_utils.hpp
/// Small numeric helpers used across the codebase.

#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

namespace gnc {

/// Binary search: returns the index such that arr[idx-1] <= val < arr[idx].
/// Equivalent to Python's bisect_right.
/// Assumes arr is sorted ascending.
template <typename ArrayLike>
inline int binary_search(const ArrayLike& arr, double val, int n) {
    int lo = 0, hi = n;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (arr[mid] <= val)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

/// Overload for Eigen vectors (uses .size()).
template <typename Derived>
inline int binary_search(const Eigen::DenseBase<Derived>& arr, double val) {
    return binary_search(arr, val, static_cast<int>(arr.size()));
}

/// Clamp a scalar.
inline double clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

} // namespace gnc
