#pragma once
/// @file logger.hpp
/// Binary data logger for simulation results.
/// Format: [n_cols:int32][n_rows:int32][data: col-major doubles]

#include <vector>
#include <string>
#include <fstream>
#include <cstdint>

namespace gnc {

struct Logger {
    int n_cols;
    std::vector<std::vector<double>> columns;  // column-major storage

    explicit Logger(int nc) : n_cols(nc), columns(nc) {}

    /// Append one row of data.
    void log(const double* row, int len) {
        int n = std::min(len, n_cols);
        for (int i = 0; i < n; ++i)
            columns[i].push_back(row[i]);
    }

    int n_rows() const { return columns.empty() ? 0 : static_cast<int>(columns[0].size()); }

    /// Write binary file.
    void save(const std::string& path) const {
        std::ofstream f(path, std::ios::binary);
        int32_t nc = n_cols;
        int32_t nr = n_rows();
        f.write(reinterpret_cast<const char*>(&nc), sizeof(int32_t));
        f.write(reinterpret_cast<const char*>(&nr), sizeof(int32_t));
        for (int c = 0; c < n_cols; ++c) {
            f.write(reinterpret_cast<const char*>(columns[c].data()),
                    nr * sizeof(double));
        }
    }

    void clear() { for (auto& c : columns) c.clear(); }
};

} // namespace gnc
