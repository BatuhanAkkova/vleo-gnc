#include "environment/gravity.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>

namespace gnc {

void load_gravity_coefficients(const std::string& egm_path, int n_max,
                               Eigen::MatrixXd& C_nm, Eigen::MatrixXd& S_nm) {
    std::ifstream f(egm_path);
    if (!f.is_open())
        throw std::runtime_error("Gravity coefficients not found: " + egm_path);

    auto data = nlohmann::json::parse(f);

    C_nm = Eigen::MatrixXd::Zero(n_max + 1, n_max + 1);
    S_nm = Eigen::MatrixXd::Zero(n_max + 1, n_max + 1);

    for (auto& item : data) {
        if (item.value("type", "") != "coeff") continue;
        int n = item.value("n", 0);
        int m = item.value("m", 0);
        if (n <= n_max && m <= n_max) {
            C_nm(n, m) = item.value("C", 0.0);
            S_nm(n, m) = item.value("S", 0.0);
        }
    }
}

} // namespace gnc
