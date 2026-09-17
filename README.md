# GNC C++ Simulation Suite

A C++17 Guidance, Navigation, and Control simulation framework for studying autonomous spacecraft deployment recovery, detumbling, and transition into precision pointing in Very-Low Earth Orbit (VLEO). 

The codebase is built entirely around the three-mode **Baseline MEKF+LQR** architecture validated in our IAC 2026 paper:

> **Autonomous Multi-Mode GNC for VLEO CubeSats: Monte Carlo Validation of High-Rate Deployment Recovery and Precision Pointing**

## System Architecture

The simulated spacecraft starts in a high-rate deployment condition, uses magnetorquer-only detumbling until it is safe to enter a protected transition phase, and then transitions into precision fine pointing. The architecture enforces this progression through three modes:

1. **Detumble (Mode 1)**: Magnetic-only B-dot rate reduction. Reaction wheels are inactive.
2. **Acquisition (Mode 2)**: Capped-geodesic LQR control to reduce attitude and rate errors into a defended local linear region using reaction wheels.
3. **Nominal (Mode 3)**: Precision pointing using geometric PD+I control with continuous magnetic momentum desaturation.

### Mode Manager & Switching Stability
The ModeManager acts as the sole switching authority. Transition from Detumble to Acquisition requires the estimated angular rate to fall below a threshold. Transition from Acquisition to Nominal requires both rate and geodesic attitude errors to settle, alongside a valid star-tracker measurement that persists for a strict **dwell time** to prevent mode chatter under aerodynamic disturbances.

## Estimators and Control Laws

### Multiplicative EKF (MEKF)
Attitude is estimated using a 6-state MEKF tracking small-angle attitude errors and gyro biases. It propagates attitude using bias-corrected gyro measurements and applies quaternion-based innovation updates from star tracker data.

### Orbit UKF
A 7-state Unscented Kalman Filter estimates ECI position, velocity, and a logarithmic atmospheric density scale. It propagates 15 sigma points through spherical harmonic gravity and aerodynamic drag, updated by GPS measurements. The UKF operates independently of the MEKF to prevent cascade failures.

### Control Formulations
- **B-dot**: Standard geomagnetic rate-of-change feedback m = -k * B_dot.
- **LQR**: Solves a continuous-time algebraic Riccati equation for the finite-horizon state [theta, omega]. The perceived geodesic error is capped to prevent early saturation.
- **PD+I**: Geometric feedback with gyroscopic decoupling and integral action for steady-state disturbance rejection.
- **Saturation**: All torque commands are passed through a strict boundary limiter (saturate_torque_preserve_direction) that preserves the intended control vector direction even if single-axis limits are exceeded.

## Truth Environment & Dynamics

Propagated via RK4 integration, the 6-DOF truth state [r, v, q, omega] incorporates:
- **Gravity**: EGM2008 spherical harmonics (truncated to degree/order 8).
- **Atmosphere**: NRLMSISE-00 density model with scalable multipliers.
- **Aerodynamics & SRP**: Interpolated aerodynamic coefficients and solar radiation pressure.
- **Magnetic Field**: IGRF dipole modeling.
- **Sensors**: Models for Gyroscope (bias/walk), Magnetometer, GPS, and Star Tracker (with dynamic rate-based validity dropouts).
- **Actuators**: Hardware-limited reaction wheels and magnetorquers.

## Building and Testing

### Prerequisites
- CMake 3.16 or newer
- C++17 compiler (MSVC, GCC, Clang)
- Python 3 (
umpy, matplotlib, scipy)

Dependencies (Eigen, nlohmann/json, GoogleTest) are automatically fetched by CMake.

### Build Instructions
`powershell
cmake -S . -B build -DBUILD_TESTS=ON
cmake --build build --config Release
`

### Unit Tests
`powershell
.\build\Release\gnc_tests.exe
`

## Running Monte Carlo Campaigns

The architecture is validated through large-scale Monte Carlo campaign grids. The simulation sweeps across environments using JSON configuration files.

`powershell
.\build\Release\gnc_sim.exe data configs\campaign_monte_carlo.json
`

**Scenario Generation:** The Cartesian product of deployment rates, density scales, actuator scales, and run indices is automatically expanded. To ensure apples-to-apples comparisons, the stochastic seed is defined as ase_seed + run_idx.

**Binary Log Format:** Outputs are written to a compact column-major binary format alongside a .json metadata sidecar. Logs contain over 80 columns including truth states, estimated states, actuator outputs, and dwell metrics.
`	ext
int32 n_cols
int32 n_rows
double data[n_cols * n_rows]
`

## Figure Generation

Single-run plotting is deprecated. Use the bundled Python script to parse the Monte Carlo binary logs, evaluate the success metrics (sustained pointing durations, RMS errors), and generate the final publication-grade figures.

`powershell
python generate_figures.py
`
Generated plots are output directly into the paper/figures/ directory.
