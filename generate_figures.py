import csv
import math
import os
import struct
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SCENARIO_SUMMARY = os.path.join(ROOT, "results", "campaign_monte_carlo", "figures", "scenario_summary.csv")
RUN_SUMMARY = os.path.join(ROOT, "results", "campaign_monte_carlo", "figures", "campaign_run_summary.csv")
ENVELOPE_DIR = os.path.join(ROOT, "results", "campaign_monte_carlo")
OUT_DIR = os.path.join(ROOT, "paper", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# Publication styling settings
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif", "Liberation Serif"],
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "figure.titlesize": 12,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})

def read_csv_rows(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def read_binary_log(path):
    with open(path, "rb") as f:
        header = f.read(8)
        if len(header) < 8:
            raise ValueError(f"bad header in {path}")
        n_cols, n_rows = struct.unpack("ii", header)
        data = np.fromfile(f, dtype=np.float64)
        if data.size != n_cols * n_rows:
            raise ValueError(f"bad data size in {path}")
        return data.reshape((n_cols, n_rows)).T

def parse_float(val, default=np.nan):
    try:
        return float(val) if val is not None and str(val).strip() != "" else default
    except ValueError:
        return default

def quat_error_deg(q_true, q_est):
    dot = np.sum(q_true * q_est, axis=1)
    dot = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(2.0 * np.arccos(np.abs(dot)))

def get_text_color(val, vmin, vmax, cmap_name="viridis"):
    """Dynamically determine text color (black or white) based on background luminance."""
    if not np.isfinite(val):
        return "black"
    norm_val = (val - vmin) / (vmax - vmin) if vmax > vmin else 0.5
    cmap = matplotlib.colormaps.get_cmap(cmap_name)
    rgba = cmap(norm_val)
    luminance = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "white" if luminance < 0.5 else "black"

# -------------------------------------------------------------
# Figure 1: Recovery & Sustained Hold Heatmaps (High-Contrast)
# -------------------------------------------------------------
def build_fig1():
    rows = read_csv_rows(SCENARIO_SUMMARY)
    rates = sorted({parse_float(r["rate_deg_s"]) for r in rows})
    x_keys = sorted({(parse_float(r["rho_scale"]), parse_float(r["actuator_scale"])) for r in rows}, key=lambda k: (k[0], -k[1]))
    
    z_entry = np.full((len(rates), len(x_keys)), np.nan)
    z_hold = np.full((len(rates), len(x_keys)), np.nan)
    
    for r in rows:
        i = rates.index(parse_float(r["rate_deg_s"]))
        j = x_keys.index((parse_float(r["rho_scale"]), parse_float(r["actuator_scale"])))
        z_entry[i, j] = parse_float(r["nominal_entry_rate"]) * 100.0
        z_hold[i, j] = parse_float(r["longest_nominal_hold_median_s"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)

    # Left: Entry Rate (Viridis)
    im1 = ax1.imshow(z_entry, aspect="auto", cmap="RdYlGn", origin="lower", vmin=0, vmax=100)
    cbar1 = fig.colorbar(im1, ax=ax1, pad=0.02)
    cbar1.set_label("Nominal Acquisition Rate [%]")
    ax1.set_title("(a) Nominal Acquisition Entry Rate")
    ax1.set_xlabel("Density Scale / Actuator Scale")
    ax1.set_ylabel("Initial Body Rate [deg/s]")
    ax1.set_xticks(range(len(x_keys)))
    ax1.set_xticklabels([f"ρ={rho:g}\nact={act:g}" for rho, act in x_keys], fontsize=8.5)
    ax1.set_yticks(range(len(rates)))
    ax1.set_yticklabels([f"{rate:g}" for rate in rates])

    # Draw cell borders
    ax1.set_xticks(np.arange(-.5, len(x_keys), 1), minor=True)
    ax1.set_yticks(np.arange(-.5, len(rates), 1), minor=True)
    ax1.grid(which="minor", color="white", linestyle='-', linewidth=1)

    for i in range(len(rates)):
        for j in range(len(x_keys)):
            v = z_entry[i, j]
            if np.isfinite(v):
                color = get_text_color(v, 0, 100, "viridis")
                ax1.text(j, i, f"{v:.0f}%", ha="center", va="center", color=color, fontsize=8.5, fontweight="bold")

    # Right: Hold Duration (Plasma)
    vmax_hold = np.nanmax(z_hold)
    im2 = ax2.imshow(z_hold, aspect="auto", cmap="RdYlGn", origin="lower", vmin=0, vmax=vmax_hold)
    cbar2 = fig.colorbar(im2, ax=ax2, pad=0.02)
    cbar2.set_label("Median Continuous Hold [s]")
    ax2.set_title("(b) Median Longest Continuous Nominal Hold")
    ax2.set_xlabel("Density Scale / Actuator Scale")
    ax2.set_xticks(range(len(x_keys)))
    ax2.set_xticklabels([f"ρ={rho:g}\nact={act:g}" for rho, act in x_keys], fontsize=8.5)

    ax2.set_xticks(np.arange(-.5, len(x_keys), 1), minor=True)
    ax2.set_yticks(np.arange(-.5, len(rates), 1), minor=True)
    ax2.grid(which="minor", color="white", linestyle='-', linewidth=1)

    for i in range(len(rates)):
        for j in range(len(x_keys)):
            v = z_hold[i, j]
            if np.isfinite(v):
                color = get_text_color(v, 0, vmax_hold, "plasma")
                ax2.text(j, i, f"{v:.0f}s", ha="center", va="center", color=color, fontsize=8.5, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig1_recovery_heatmaps.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 1 (High-Contrast)")

# -------------------------------------------------------------
# Figure 2: Statistical Uncertainty & Hold Variability (High-Contrast)
# -------------------------------------------------------------
def build_fig2():
    rows = read_csv_rows(SCENARIO_SUMMARY)
    rates = sorted({parse_float(r["rate_deg_s"]) for r in rows})
    x_keys = sorted({(parse_float(r["rho_scale"]), parse_float(r["actuator_scale"])) for r in rows}, key=lambda k: (k[0], -k[1]))
    
    z_ci = np.full((len(rates), len(x_keys)), np.nan)
    z_iqr = np.full((len(rates), len(x_keys)), np.nan)
    
    run_rows = read_csv_rows(RUN_SUMMARY)
    scenario_holds = {}
    for r in run_rows:
        sc = r["scenario"]
        h = parse_float(r["longest_nominal_hold_s"], 0.0)
        scenario_holds.setdefault(sc, []).append(h)

    z = 1.96
    n = 20.0
    for r in rows:
        sc = r["scenario"]
        i = rates.index(parse_float(r["rate_deg_s"]))
        j = x_keys.index((parse_float(r["rho_scale"]), parse_float(r["actuator_scale"])))
        p_hat = parse_float(r["nominal_entry_rate"])
        half_w = (z * math.sqrt((p_hat * (1.0 - p_hat) / n) + (z**2 / (4.0 * n**2)))) / (1.0 + (z**2 / n))
        z_ci[i, j] = half_w * 100.0
        
        holds = scenario_holds.get(sc, [0.0])
        z_iqr[i, j] = np.percentile(holds, 75) - np.percentile(holds, 25)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2), sharey=True)

    # Left: Wilson CI Half-Width (Magma_r)
    vmax_ci = np.nanmax(z_ci)
    im1 = ax1.imshow(z_ci, aspect="auto", cmap="RdYlGn_r", origin="lower", vmin=0, vmax=vmax_ci)
    cbar1 = fig.colorbar(im1, ax=ax1, pad=0.02)
    cbar1.set_label("95% Wilson CI Half-Width [%]")
    ax1.set_title("(a) Nominal Entry Rate Uncertainty (95% CI)")
    ax1.set_xlabel("Density Scale / Actuator Scale")
    ax1.set_ylabel("Initial Body Rate [deg/s]")
    ax1.set_xticks(range(len(x_keys)))
    ax1.set_xticklabels([f"ρ={rho:g}\nact={act:g}" for rho, act in x_keys], fontsize=8.5)
    ax1.set_yticks(range(len(rates)))
    ax1.set_yticklabels([f"{rate:g}" for rate in rates])

    ax1.set_xticks(np.arange(-.5, len(x_keys), 1), minor=True)
    ax1.set_yticks(np.arange(-.5, len(rates), 1), minor=True)
    ax1.grid(which="minor", color="white", linestyle='-', linewidth=1)

    for i in range(len(rates)):
        for j in range(len(x_keys)):
            v = z_ci[i, j]
            if np.isfinite(v):
                color = get_text_color(v, 0, vmax_ci, "magma_r")
                ax1.text(j, i, f"±{v:.1f}%", ha="center", va="center", color=color, fontsize=8, fontweight="bold")

    # Right: IQR Variability (Cividis)
    vmax_iqr = np.nanmax(z_iqr)
    im2 = ax2.imshow(z_iqr, aspect="auto", cmap="RdYlGn_r", origin="lower", vmin=0, vmax=vmax_iqr)
    cbar2 = fig.colorbar(im2, ax=ax2, pad=0.02)
    cbar2.set_label("IQR of Continuous Hold [s]")
    ax2.set_title("(b) Continuous Hold Variability (Interquartile Range)")
    ax2.set_xlabel("Density Scale / Actuator Scale")
    ax2.set_xticks(range(len(x_keys)))
    ax2.set_xticklabels([f"ρ={rho:g}\nact={act:g}" for rho, act in x_keys], fontsize=8.5)

    ax2.set_xticks(np.arange(-.5, len(x_keys), 1), minor=True)
    ax2.set_yticks(np.arange(-.5, len(rates), 1), minor=True)
    ax2.grid(which="minor", color="white", linestyle='-', linewidth=1)

    for i in range(len(rates)):
        for j in range(len(x_keys)):
            v = z_iqr[i, j]
            if np.isfinite(v):
                color = get_text_color(v, 0, vmax_iqr, "cividis")
                ax2.text(j, i, f"{v:.0f}s", ha="center", va="center", color=color, fontsize=8, fontweight="bold")

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig2_statistical_uncertainty.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 2 (High-Contrast)")

# -------------------------------------------------------------
# Figure 3a: Trajectory Recovery (3 Subplots)
# -------------------------------------------------------------
def build_fig3a():
    bin_path = os.path.join(ENVELOPE_DIR, "BaselineMEKF_LQR_rate50_rho6_act0.85", "run_0.bin")
    data = read_binary_log(bin_path)[::50]

    t = data[:, 0]
    q_true = data[:, 7:11]
    w_true = np.rad2deg(data[:, 11:14])
    w_est = np.rad2deg(data[:, 18:21])
    mode = data[:, 55]

    fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharex=True)

    # Subplot (a): GNC Mode
    axs[0].plot(t, mode, 'navy', lw=1.5)
    axs[0].set_yticks([1, 2, 4])
    axs[0].set_yticklabels(['Detumble', 'Acquisition', 'Nominal'])
    axs[0].set_ylabel('GNC Mode')
    axs[0].set_title('(a) GNC Mode Sequence & State Machine Handover (50 deg/s, ρ=6, act=0.85)')

    # Subplot (b): Body Rates
    axs[1].plot(t, np.linalg.norm(w_true, axis=1), 'k-', label='True Total Rate', lw=1.2)
    axs[1].plot(t, w_est[:, 0], 'r--', label='Est ωx', alpha=0.7, lw=0.8)
    axs[1].plot(t, w_est[:, 1], 'g--', label='Est ωy', alpha=0.7, lw=0.8)
    axs[1].plot(t, w_est[:, 2], 'b--', label='Est ωz', alpha=0.7, lw=0.8)
    axs[1].set_ylabel('Body Rate [deg/s]')
    axs[1].legend(loc='upper right', ncol=2)

    # Subplot (c): Pointing Error
    target_q = np.array([0.0, 0.0, 0.0, 1.0])
    truth_pointing = quat_error_deg(q_true, np.tile(target_q, (len(q_true), 1)))
    axs[2].plot(t, truth_pointing, 'crimson', lw=1.2)
    axs[2].set_ylabel('Pointing Error [deg]')
    axs[2].set_xlabel('Time [s]')
    axs[2].set_yscale('log')
    axs[2].set_ylim([1e-1, 200])

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig3a_operational_trajectory.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 3a (Split 1)")

# -------------------------------------------------------------
# Figure 3b: Subsystem Telemetry & Envelopes (3 Subplots)
# -------------------------------------------------------------
def build_fig3b():
    bin_path = os.path.join(ENVELOPE_DIR, "BaselineMEKF_LQR_rate50_rho6_act0.85", "run_0.bin")
    data = read_binary_log(bin_path)[::50]

    t = data[:, 0]
    q_true = data[:, 7:11]
    q_est = data[:, 14:18]
    mekf_p = data[:, 30:36]
    rw_speed = data[:, 46:49]
    mtq_dipole = data[:, 52:55]

    att_err = quat_error_deg(q_true, q_est)
    att_sigma = 3.0 * np.rad2deg(np.sqrt(np.clip(np.sum(mekf_p[:, 0:3], axis=1), 0, None)))

    fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharex=True)

    # Subplot (a): MEKF Error & 3-Sigma
    axs[0].plot(t, att_err, 'b-', label='Attitude Error', lw=1.0)
    axs[0].plot(t, att_sigma, 'r--', label='3σ Bounds', lw=1.0)
    axs[0].set_ylabel('MEKF Error [deg]')
    axs[0].set_yscale('log')
    axs[0].set_ylim([1e-2, 100])
    axs[0].set_title('(b) Subsystem Telemetry & Estimator Envelopes')
    axs[0].legend(loc='upper right')

    # Subplot (b): RW Speeds
    axs[1].plot(t, rw_speed[:, 0], 'r-', label='RW 1', lw=0.9)
    axs[1].plot(t, rw_speed[:, 1], 'g-', label='RW 2', lw=0.9)
    axs[1].plot(t, rw_speed[:, 2], 'b-', label='RW 3', lw=0.9)
    axs[1].axhline(150, color='k', ls=':', label='Limits')
    axs[1].axhline(-150, color='k', ls=':')
    axs[1].set_ylabel('RW Speed [rad/s]')
    axs[1].legend(loc='upper right', ncol=2)

    # Subplot (c): MTQ Dipoles
    axs[2].plot(t, mtq_dipole[:, 0], 'r-', label='MTQ X', lw=0.9)
    axs[2].plot(t, mtq_dipole[:, 1], 'g-', label='MTQ Y', lw=0.9)
    axs[2].plot(t, mtq_dipole[:, 2], 'b-', label='MTQ Z', lw=0.9)
    axs[2].axhline(0.3, color='k', ls=':', label='Limit (±0.3)')
    axs[2].axhline(-0.3, color='k', ls=':')
    axs[2].set_ylabel('MTQ Dipole [A m²]')
    axs[2].set_xlabel('Time [s]')
    axs[2].legend(loc='upper right', ncol=2)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig3b_subsystem_telemetry.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 3b (Split 2)")

# Maintain legacy fig3 image for backwards compatibility
def build_fig3_legacy():
    build_fig3a()
    build_fig3b()

# -------------------------------------------------------------
# Figure 4: MEKF Attitude Estimation Performance
# -------------------------------------------------------------
def build_fig4():
    bin_path = os.path.join(ENVELOPE_DIR, "BaselineMEKF_LQR_rate50_rho6_act0.85", "run_0.bin")
    data = read_binary_log(bin_path)[::20]

    t = data[:, 0]
    q_true = data[:, 7:11]
    q_est = data[:, 14:18]
    bias_est = np.rad2deg(data[:, 27:30])
    mekf_p = data[:, 30:36]

    att_err = quat_error_deg(q_true, q_est)
    att_sigma = 3.0 * np.rad2deg(np.sqrt(np.clip(np.sum(mekf_p[:, 0:3], axis=1), 0, None)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.2))

    # Left: Attitude Error & 3-Sigma
    ax1.plot(t, att_err, 'darkblue', lw=1.2, label='Attitude Error |q_err|')
    ax1.plot(t, att_sigma, 'crimson', lw=1.2, ls='--', label='3σ Covariance Envelope')
    ax1.set_ylabel('Attitude Error [deg]')
    ax1.set_xlabel('Time [s]')
    ax1.set_yscale('log')
    ax1.set_ylim([1e-2, 200])
    ax1.set_title('(a) MEKF Attitude Error & 3σ Bounds')
    ax1.legend(loc='upper right')

    # Right: Gyro Bias Convergence
    ax2.plot(t, bias_est[:, 0], 'r-', label='Bias β_x', lw=1.0)
    ax2.plot(t, bias_est[:, 1], 'g-', label='Bias β_y', lw=1.0)
    ax2.plot(t, bias_est[:, 2], 'b-', label='Bias β_z', lw=1.0)
    ax2.set_ylabel('Estimated Gyro Bias [deg/s]')
    ax2.set_xlabel('Time [s]')
    ax2.set_title('(b) Onboard Gyro Bias State Convergence')
    ax2.legend(loc='upper right')

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig4_mekf_estimation_performance.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 4")

# -------------------------------------------------------------
# Figure 5: Reaction Wheel Actuator Dynamics
# -------------------------------------------------------------
def build_fig5():
    bin_path = os.path.join(ENVELOPE_DIR, "BaselineMEKF_LQR_rate50_rho6_act0.85", "run_0.bin")
    data = read_binary_log(bin_path)[::20]

    t = data[:, 0]
    rw_speed = data[:, 46:49]
    rw_torque = data[:, 49:52] * 1000.0
    mtq_dipole = data[:, 52:55]
    rw_momentum = rw_speed * 1e-4 * 1000.0

    fig, axs = plt.subplots(2, 2, figsize=(12.0, 8.5))

    axs[0,0].plot(t, rw_speed[:, 0], 'r-', label='Wheel 1', lw=1.0)
    axs[0,0].plot(t, rw_speed[:, 1], 'g-', label='Wheel 2', lw=1.0)
    axs[0,0].plot(t, rw_speed[:, 2], 'b-', label='Wheel 3', lw=1.0)
    axs[0,0].axhline(150, color='k', ls='--', label='Speed Capacity (+/-150 rad/s)')
    axs[0,0].axhline(-150, color='k', ls='--')
    axs[0,0].set_ylabel('Wheel Angular Speed [rad/s]')
    axs[0,0].set_xlabel('Time [s]')
    axs[0,0].set_title('(a) Reaction Wheel Speeds')
    axs[0,0].legend(loc='upper right', fontsize=8.5)

    axs[0,1].plot(t, rw_torque[:, 0], 'r-', label='RW Torque X', lw=0.9, alpha=0.8)
    axs[0,1].plot(t, rw_torque[:, 1], 'g-', label='RW Torque Y', lw=0.9, alpha=0.8)
    axs[0,1].plot(t, rw_torque[:, 2], 'b-', label='RW Torque Z', lw=0.9, alpha=0.8)
    axs[0,1].axhline(2.0, color='k', ls='--', label='Max Torque Limit (+/-2.0 mN m)')
    axs[0,1].axhline(-2.0, color='k', ls='--')
    axs[0,1].set_ylabel('Wheel Reaction Torque [mN m]')
    axs[0,1].set_xlabel('Time [s]')
    axs[0,1].set_title('(b) Reaction Wheel Applied Control Torques')
    axs[0,1].legend(loc='upper right', fontsize=8.5)

    axs[1,0].plot(t, mtq_dipole[:, 0], 'r-', label='MTQ Dipole X', lw=0.9)
    axs[1,0].plot(t, mtq_dipole[:, 1], 'g-', label='MTQ Dipole Y', lw=0.9)
    axs[1,0].plot(t, mtq_dipole[:, 2], 'b-', label='MTQ Dipole Z', lw=0.9)
    axs[1,0].axhline(0.3, color='k', ls='--', label='Dipole Limit (+/-0.3 A m^2)')
    axs[1,0].axhline(-0.3, color='k', ls='--')
    axs[1,0].set_ylabel('MTQ Magnetic Dipole [A m^2]')
    axs[1,0].set_xlabel('Time [s]')
    axs[1,0].set_title('(c) Magnetorquer Dipole Commands')
    axs[1,0].legend(loc='upper right', fontsize=8.5)

    axs[1,1].plot(t, rw_momentum[:, 0], 'r-', label='RW Momentum X', lw=1.0)
    axs[1,1].plot(t, rw_momentum[:, 1], 'g-', label='RW Momentum Y', lw=1.0)
    axs[1,1].plot(t, rw_momentum[:, 2], 'b-', label='RW Momentum Z', lw=1.0)
    axs[1,1].axhline(15.0, color='k', ls='--', label='Max Capacity (15 mN m s)')
    axs[1,1].axhline(-15.0, color='k', ls='--')
    axs[1,1].set_ylabel('Stored Wheel Momentum [mN m s]')
    axs[1,1].set_xlabel('Time [s]')
    axs[1,1].set_title('(d) Magnetic Momentum Unloading Profile')
    axs[1,1].legend(loc='upper right', fontsize=8.5)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig5_actuator_dynamics.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 5 (4-panel actuator dynamics)")

# Figure 7: UKF Navigation
def build_fig7():
    bin_path = os.path.join(ENVELOPE_DIR, "BaselineMEKF_LQR_rate50_rho6_act0.85", "run_0.bin")
    data = read_binary_log(bin_path)[::20]

    t = data[:, 0]
    r_true = data[:, 1:4]
    r_est = data[:, 21:24]
    ukf_p = data[:, 36:43]
    rho_true = data[:, 43]
    rho_scale_est = np.exp(data[:, 44])

    pos_err = np.linalg.norm(r_est - r_true, axis=1)
    pos_sigma = 3.0 * np.sqrt(np.clip(np.sum(ukf_p[:, 0:3], axis=1), 0, None))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.2))

    # Left: Orbit Position Error
    ax1.plot(t, pos_err, 'darkblue', lw=1.2, label='Orbit Position Error')
    ax1.plot(t, pos_sigma, 'crimson', lw=1.2, ls='--', label='3σ Position Bounds')
    ax1.set_ylabel('ECI Position Error [m]')
    ax1.set_xlabel('Time [s]')
    ax1.set_title('(a) UKF Orbit Navigation Error & 3σ Bounds')
    ax1.legend(loc='upper right')

    # Right: Density Tracking (Showing unobservability and structural bias without dedicated drag accelerometer)
    ax2_scale = ax2.twinx()
    l1 = ax2.plot(t, rho_true, 'k-', lw=1.2, label='True NRLMSISE Density ρ', alpha=0.7)
    l2 = ax2_scale.plot(t, rho_scale_est, 'crimson', lw=1.2, ls='--', label='Onboard Density Scale Estimate')
    
    ax2.set_ylabel('Atmospheric Mass Density [kg/m³]', color='k')
    ax2_scale.set_ylabel('Onboard Density Scale Factor [-]', color='crimson')
    ax2.set_xlabel('Time [s]')
    ax2.set_title('(b) Onboard Density Scale Structural Bias & Unobservability')

    lines = l1 + l2
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc='upper right')

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig7_orbit_density_ukf.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 7")

# -------------------------------------------------------------
# Figure 8: Aerodynamic Torque Drift Mechanism at 70 deg/s
# -------------------------------------------------------------
def build_fig8():
    rho_scenarios = [
        ("ρ = 1.0 (Low Density)", "BaselineMEKF_LQR_rate70_rho1_act1", "forestgreen"),
        ("ρ = 3.0 (Moderate Density)", "BaselineMEKF_LQR_rate70_rho3_act1", "darkorange"),
        ("ρ = 6.0 (High VLEO Density)", "BaselineMEKF_LQR_rate70_rho6_act1", "crimson"),
    ]

    target_q = np.array([0.0, 0.0, 0.0, 1.0])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    for label, sc, color in rho_scenarios:
        bin_path = os.path.join(ENVELOPE_DIR, sc, "run_0.bin")
        if os.path.exists(bin_path):
            data = read_binary_log(bin_path)[::20]
            t = data[:, 0]
            q_true = data[:, 7:11]
            rw_speed = data[:, 46:49]
            mode = data[:, 55]

            pointing = quat_error_deg(q_true, np.tile(target_q, (len(q_true), 1)))
            rw_h_mag = np.linalg.norm(rw_speed * 1e-4 * 1000.0, axis=1) # mN.m.s
            
            # Find first nominal entry index
            nom_idx = np.where(mode == 4.0)[0]
            if len(nom_idx) > 0:
                t_first = t[nom_idx[0]]
                t_rel = t - t_first
                mask = (t_rel >= -50) & (t_rel <= 1500)
                
                p_plot = pointing[mask].copy()
                rw_plot = rw_h_mag[mask].copy()
                m_plot = mode[mask]
                
                # Insert gaps when not in nominal mode (after the initial capture)
                p_plot[m_plot != 4.0] = np.nan
                # For the first 50 seconds (pre-entry), we can leave it or NaN it.
                # Let's just plot it as is to show the entry, but break the line on exit.
                # Actually, a better way is to only show the FIRST continuous hold.
                first_exit_idx = np.where((m_plot != 4.0) & (t_rel[mask] > 5.0))[0]
                if len(first_exit_idx) > 0:
                    cutoff = first_exit_idx[0]
                    p_plot[cutoff:] = np.nan
                    rw_plot[cutoff:] = np.nan
                
                ax1.plot(t_rel[mask], p_plot, color=color, lw=1.4, label=f"{label}")
                ax2.plot(t_rel[mask], rw_plot, color=color, lw=1.4, label=f"{label}")

    # Panel (a): Pointing Drift
    ax1.axhline(20.0, color='k', ls='--', lw=1.2, label='Nominal Exit Guard (20°)')
    ax1.axhline(15.0, color='gray', ls=':', lw=1.2, label='Nominal Entry Guard (15°)')
    ax1.set_ylabel('True Pointing Error [deg]')
    ax1.set_xlabel('Time Elapsed Since Nominal Entry [s]')
    ax1.set_title('(a) Post-Entry Pointing Error Drift')
    ax1.set_ylim([0, 30])
    ax1.legend(loc='upper left', fontsize=8)

    # Panel (b): RW Momentum Saturation
    ax2.axhline(15.0, color='k', ls='--', lw=1.2, label='RW Momentum Capacity (15 mN m s)')
    ax2.set_ylabel('Stored RW Momentum Magnitude [mN m s]')
    ax2.set_xlabel('Time Elapsed Since Nominal Entry [s]')
    ax2.set_title('(b) Reaction Wheel Momentum Accumulation')
    ax2.set_ylim([0, 20])
    ax2.legend(loc='lower right', fontsize=8)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig8_aerodynamic_drift_mechanism.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 8 (2-Panel Physics Drift Mechanism)")

# -------------------------------------------------------------
# Figure 9: Average Dwell-Time & Mode Switching Stability
# -------------------------------------------------------------
def build_fig9():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.5))

    # Left: Gate Accumulator & Handover Window
    bin_path = os.path.join(ENVELOPE_DIR, "BaselineMEKF_LQR_rate50_rho6_act0.85", "run_0.bin")
    data = read_binary_log(bin_path)[::10]
    t = data[:, 0]
    w_true = np.rad2deg(data[:, 11:14])
    q_true = data[:, 7:11]
    mode = data[:, 55]
    target_q = np.array([0.0, 0.0, 0.0, 1.0])
    pointing = quat_error_deg(q_true, np.tile(target_q, (len(q_true), 1)))
    rate_mag = np.linalg.norm(w_true, axis=1)

    # Find the exact handover index (Acquisition -> Nominal)
    nom_idx = np.where(mode == 4.0)[0]
    if len(nom_idx) > 0:
        idx_entry = nom_idx[0]
        t_entry = t[idx_entry]
        window_mask = (t >= t_entry - 40.0) & (t <= t_entry + 60.0)
        
        t_win = t[window_mask]
        pointing_win = pointing[window_mask]
        rate_win = rate_mag[window_mask]
        mode_win = mode[window_mask]

        # Gate dwell accumulation during window (clamped at 20s upon handover)
        in_gate = (pointing_win < 15.0) & (rate_win < 0.8)
        dwell_acc = np.zeros(len(t_win))
        curr = 0.0
        for i in range(1, len(t_win)):
            if mode_win[i] == 4.0:
                curr = 20.0  # Clamped upon handover to nominal
            elif in_gate[i]:
                curr += (t_win[i] - t_win[i-1])
            else:
                curr = 0.0
            dwell_acc[i] = curr

        ax1.plot(t_win - t_entry, dwell_acc, 'navy', lw=1.8, label='Dwell Accumulator [s]')
        ax1_twin = ax1.twinx()
        ax1_twin.plot(t_win - t_entry, mode_win, 'crimson', lw=1.2, ls='--', label='GNC Mode')
        ax1_twin.set_yticks([2, 4])
        ax1_twin.set_yticklabels(['Acq (2)', 'Nom (4)'])
        ax1_twin.set_ylabel('GNC Mode', color='crimson')

        ax1.axhline(20.0, color='darkgreen', ls=':', lw=1.2, label='Required Dwell Bound (20s)')
        ax1.axvline(0.0, color='black', ls='-', lw=1.0, label='Mode Handover Instant')
        ax1.set_ylabel('Continuous Dwell Time [s]', color='navy')
        ax1.set_xlabel('Time Relative to Mode Handover [s]')
        ax1.set_title('(a) Dwell Gate Accumulator & Mode Handover')
        
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax1_twin.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8)

    # Right: Mode Switching Chatter Histogram
    run_rows = read_csv_rows(RUN_SUMMARY)
    switches = []
    for r in run_rows:
        if 'switch_count' in r and r['switch_count'] not in ['', 'nan']:
            switches.append(int(float(r['switch_count'])))

    bins = np.arange(0, max(switches)+3) - 0.5 if switches else 10
    ax2.hist(switches, bins=bins, color='teal', edgecolor='black', alpha=0.7)
    ax2.set_ylabel('Number of Monte Carlo Runs')
    ax2.set_xlabel('Total Mode Switch Count per Run')
    ax2.set_title('(b) Switching Distribution (Ideal = 2)')
    ax2.xaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=10))

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig9_switching_dwell_stability.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 9 (Handover Window & ADT Margin)")


# -------------------------------------------------------------
# Figure 10: Fixed-Gain Disturbance Rejection Tradeoff
# -------------------------------------------------------------
def build_fig10():
    scenarios_list = [
        ('Benign\n(40$^\circ$/s, $\\rho$=1)', '40.0', '1.0'),
        ('Moderate\n(50$^\circ$/s, $\\rho$=3)', '50.0', '3.0'),
        ('Severe\n(60$^\circ$/s, $\\rho$=6)', '60.0', '6.0'),
        ('Harsh Edge\n(70$^\circ$/s, $\\rho$=6)', '70.0', '6.0')
    ]
    
    # We will read from the three campaigns
    campaigns = [
        ('g1', os.path.join(os.path.join(ROOT, "results"), 'campaign_monte_carlo', 'figures', 'scenario_summary.csv'),
               os.path.join(os.path.join(ROOT, "results"), 'campaign_monte_carlo', 'figures', 'campaign_run_summary.csv')),
        ('g2', os.path.join(os.path.join(ROOT, "results"), 'campaign_tradeoff_g2', 'figures', 'scenario_summary.csv'),
               os.path.join(os.path.join(ROOT, "results"), 'campaign_tradeoff_g2', 'figures', 'campaign_run_summary.csv')),
        ('g3', os.path.join(os.path.join(ROOT, "results"), 'campaign_tradeoff_g3', 'figures', 'scenario_summary.csv'),
               os.path.join(os.path.join(ROOT, "results"), 'campaign_tradeoff_g3', 'figures', 'campaign_run_summary.csv'))
    ]
    
    entry_rates = {'g1': [], 'g2': [], 'g3': []}
    median_holds = {'g1': [], 'g2': [], 'g3': []}
    
    for c_id, sc_path, run_path in campaigns:
        sc_data = {}
        if os.path.exists(sc_path):
            rows = read_csv_rows(sc_path)
            for r in rows:
                if r['actuator_scale'] == '1.0':
                    sc_data[(r['rate_deg_s'], r['rho_scale'])] = float(r['nominal_entry_rate']) * 100.0
                    
        run_data = {}
        if os.path.exists(run_path):
            rows = read_csv_rows(run_path)
            for r in rows:
                if 'act0.85' in r['scenario']: continue
                import re
                m = re.search(r'rate([0-9\.]+)_rho([0-9\.]+)_act1', r['scenario'])
                if m:
                    rate, rho = m.group(1), m.group(2)
                    rate = f"{float(rate):.1f}"
                    rho = f"{float(rho):.1f}"
                    run_data.setdefault((rate, rho), []).append(float(r['longest_nominal_hold_s']) if r['longest_nominal_hold_s'] else 0.0)
                    
        for label, r_target, rho_target in scenarios_list:
            entry_rates[c_id].append(sc_data.get((r_target, rho_target), 0.0))
            holds = run_data.get((r_target, rho_target), [0.0])
            median_holds[c_id].append(np.median(holds) if holds else 0.0)

    x = np.arange(len(scenarios_list))
    width = 0.25

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.2))

    # Left: Nominal Entry Rate
    ax1.bar(x - width, entry_rates['g1'], width, label='Baseline Gain (g1)', color='skyblue', edgecolor='black')
    ax1.bar(x, entry_rates['g2'], width, label='High Gain (g2 - Paper)', color='navy', edgecolor='black')
    ax1.bar(x + width, entry_rates['g3'], width, label='Fixed Integral (g3)', color='salmon', edgecolor='black')
    ax1.set_ylabel('Nominal Entry Rate [%]')
    ax1.set_xticks(x)
    ax1.set_xticklabels([s[0] for s in scenarios_list])
    ax1.set_title('(a) Nominal Acquisition Rate Across Control Settings')
    ax1.set_ylim([0, 115])
    ax1.legend(loc='upper right')

    # Right: Median Continuous Nominal Hold
    ax2.bar(x - width, median_holds['g1'], width, label='Baseline Gain (g1)', color='skyblue', edgecolor='black')
    ax2.bar(x, median_holds['g2'], width, label='High Gain (g2 - Paper)', color='navy', edgecolor='black')
    ax2.bar(x + width, median_holds['g3'], width, label='Fixed Integral (g3)', color='salmon', edgecolor='black')
    ax2.set_ylabel('Median Continuous Nominal Hold [s]')
    ax2.set_xticks(x)
    ax2.set_xticklabels([s[0] for s in scenarios_list])
    ax2.set_title('(b) True Nominal Pointing Endurance')
    ax2.axhline(7200, color='gray', ls='--', lw=1.2, label='Simulation Duration (7200s)')
    ax2.set_ylim([0, 7800])
    ax2.legend(loc='upper right')

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "fig10_fixed_gain_tradeoff.png"), dpi=400, bbox_inches='tight')
    plt.close(fig)
    print("Generated Fig 10 (Real Tradeoff Data)")

def main():
    print("Building publication figures...")
    os.makedirs(OUT_DIR, exist_ok=True)
    
    build_fig1()
    build_fig2()
    build_fig3a()
    build_fig3b()
    build_fig4()
    build_fig5()
    build_fig7()
    build_fig8()
    build_fig9()
    build_fig10()
    print("All publication figures successfully built in paper/figures/")

if __name__ == "__main__":
    main()

