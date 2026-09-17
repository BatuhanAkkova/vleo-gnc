#include "environment/density.hpp"
#include "environment/magnetic.hpp"
#include "environment/wind.hpp"
#include "environment/aerodynamics.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <cmath>

namespace gnc {

// ── Binary LUT format ────────────────────────────────────────────────
// Header: 3 x int32 (n0, n1, n2) = axis sizes
// Then: n0 doubles (axis0), n1 doubles (axis1), n2 doubles (axis2)
// Then: n0*n1*n2 doubles (flat grid)
// For multi-component LUTs: additional n0*n1*n2 blocks follow.

static std::vector<double> read_doubles(std::ifstream& f, int n) {
    std::vector<double> v(n);
    f.read(reinterpret_cast<char*>(v.data()), n * sizeof(double));
    return v;
}

static void read_dims(std::ifstream& f, int& n0, int& n1, int& n2) {
    int32_t d[3];
    f.read(reinterpret_cast<char*>(d), 3 * sizeof(int32_t));
    n0 = d[0]; n1 = d[1]; n2 = d[2];
}

// ── Density LUT ──────────────────────────────────────────────────────

DensityLUT load_density_lut(const std::string& bin_path) {
    std::ifstream f(bin_path, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("Cannot open density LUT: " + bin_path);

    DensityLUT lut;
    read_dims(f, lut.n_alt, lut.n_lat, lut.n_lon);
    lut.alts_km  = read_doubles(f, lut.n_alt);
    lut.lats_deg = read_doubles(f, lut.n_lat);
    lut.lons_deg = read_doubles(f, lut.n_lon);

    // log_rho stored in file -> exponentiate
    int total = lut.n_alt * lut.n_lat * lut.n_lon;
    lut.rho_grid = read_doubles(f, total);
    for (auto& v : lut.rho_grid) v = std::exp(v);

    return lut;
}

// ── IGRF LUT ─────────────────────────────────────────────────────────

IgrfLUT load_igrf_lut(const std::string& bin_path) {
    std::ifstream f(bin_path, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("Cannot open IGRF LUT: " + bin_path);

    IgrfLUT lut;
    read_dims(f, lut.n_lat, lut.n_lon, lut.n_alt);  // stored as (lat, lon, alt)
    lut.lats_deg = read_doubles(f, lut.n_lat);
    lut.lons_deg = read_doubles(f, lut.n_lon);
    lut.alts_km  = read_doubles(f, lut.n_alt);

    int total = lut.n_lat * lut.n_lon * lut.n_alt;
    lut.be = read_doubles(f, total);
    lut.bn = read_doubles(f, total);
    lut.bu = read_doubles(f, total);

    return lut;
}

// ── Wind LUT ─────────────────────────────────────────────────────────

WindLUT load_wind_lut(const std::string& bin_path) {
    std::ifstream f(bin_path, std::ios::binary);
    if (!f.is_open()) throw std::runtime_error("Cannot open Wind LUT: " + bin_path);

    WindLUT lut;
    read_dims(f, lut.n_alt, lut.n_lat, lut.n_lon);
    lut.alts_km  = read_doubles(f, lut.n_alt);
    lut.lats_deg = read_doubles(f, lut.n_lat);
    lut.lons_deg = read_doubles(f, lut.n_lon);

    int total = lut.n_alt * lut.n_lat * lut.n_lon;
    lut.u_zonal = read_doubles(f, total);
    lut.u_merid = read_doubles(f, total);

    return lut;
}

// ── Aero Coefficients ────────────────────────────────────────────────

static Eigen::MatrixXd load_gsi_json(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("Cannot open GSI coeffs: " + path);

    auto data = nlohmann::json::parse(f);
    std::vector<std::array<double,3>> rows;
    for (auto& item : data) {
        if (item.value("type", "") != "coeff") continue;
        rows.push_back({std::stod(item.at("val").get<std::string>()),
                        std::stod(item.at("power_s").get<std::string>()),
                        std::stod(item.at("power_t").get<std::string>())});
    }
    Eigen::MatrixXd mat(static_cast<int>(rows.size()), 3);
    for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
        mat(i, 0) = rows[i][0];
        mat(i, 1) = rows[i][1];
        mat(i, 2) = rows[i][2];
    }
    return mat;
}

AeroCoeffs load_aero_coeffs(const std::string& cn_path, const std::string& ct_path) {
    AeroCoeffs ac;
    ac.cn = load_gsi_json(cn_path);
    ac.ct = load_gsi_json(ct_path);
    return ac;
}

} // namespace gnc
