#include <gtest/gtest.h>
#include "math/interp3d.hpp"
#include <vector>

using namespace gnc;

TEST(Interp3D, AxisLookupPeriodic) {
    std::vector<double> axis = {0, 10, 20, 30};
    int n = 4;
    double period = 40.0;

    // Normal lookup
    auto res = axis_lookup_periodic(axis.data(), n, 15.0, period);
    EXPECT_EQ(res.idx, 1);
    EXPECT_NEAR(res.w, 0.5, 1e-10);

    // Wrap around from negative
    res = axis_lookup_periodic(axis.data(), n, -5.0, period);
    EXPECT_EQ(res.idx, 3);
    EXPECT_NEAR(res.w, 0.5, 1e-10);

    // Wrap around from positive
    res = axis_lookup_periodic(axis.data(), n, 45.0, period);
    EXPECT_EQ(res.idx, 0);
    EXPECT_NEAR(res.w, 0.5, 1e-10);

    // Exactly at last element
    res = axis_lookup_periodic(axis.data(), n, 30.0, period);
    EXPECT_EQ(res.idx, 3);
    EXPECT_NEAR(res.w, 0.0, 1e-10);

    // Between last and first (periodic wrap)
    res = axis_lookup_periodic(axis.data(), n, 35.0, period);
    EXPECT_EQ(res.idx, 3);
    EXPECT_NEAR(res.w, 0.5, 1e-10);
}

TEST(Interp3D, TrilinearPeriodic) {
    int n0 = 2, n1 = 2, n2 = 2;
    // Grid: [2, 2, 2]
    // Axis2 is periodic: 0, 10 (period 20)
    // grid[0,0,0] = 1, grid[0,0,1] = 2
    double grid[8] = {
        1.0, 2.0, // (0,0,0), (0,0,1)
        3.0, 4.0, // (0,1,0), (0,1,1)
        5.0, 6.0, // (1,0,0), (1,0,1)
        7.0, 8.0  // (1,1,0), (1,1,1)
    };

    // Interpolate at (0.5, 0.5, 15)
    // k = 1 (10 deg), k+1 = 0 (0 deg) because periodic
    // w2 = (15 - 10) / (20 - 10) = 0.5
    double val = trilinear(grid, n0, n1, n2, 0, 0, 1, 0.5, 0.5, 0.5, false, false, true);
    
    // Manual calculation:
    // v000 = G(0,0,1) = 2.0
    // v100 = G(1,0,1) = 6.0
    // v010 = G(0,1,1) = 4.0
    // v001 = G(0,0,0) = 1.0 (wrapped)
    // v110 = G(1,1,1) = 8.0
    // v101 = G(1,0,0) = 5.0 (wrapped)
    // v011 = G(0,1,0) = 3.0 (wrapped)
    // v111 = G(1,1,0) = 7.0 (wrapped)
    
    // Along axis-0:
    // v00 = 2*(0.5) + 6*(0.5) = 4.0
    // v01 = 1*(0.5) + 5*(0.5) = 3.0
    // v10 = 4*(0.5) + 8*(0.5) = 6.0
    // v11 = 3*(0.5) + 7*(0.5) = 5.0
    
    // Along axis-1:
    // v0 = 4*(0.5) + 6*(0.5) = 5.0
    // v1 = 3*(0.5) + 5*(0.5) = 4.0
    
    // Along axis-2:
    // res = 5*(0.5) + 4*(0.5) = 4.5
    
    EXPECT_NEAR(val, 4.5, 1e-10);
}
