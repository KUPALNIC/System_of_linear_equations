#include "csr.hpp"
#include <algorithm>
#include <stdexcept>

CSR::CSR(const std::map<std::tuple<int, int>, double>& dok, int rows, int cols) 
    : rows_count(rows), cols_count(cols) {
    this->rows.resize(rows_count + 1, 0);
    std::vector<std::tuple<int, int, double>> sorted_elements;

    for (const auto& element : dok) {
        int i = std::get<0>(element.first);
        int j = std::get<1>(element.first);
        double value = element.second;        
        if (value != 0.0) {
            sorted_elements.push_back(std::make_tuple(i, j, value));
        }
    }
    std::sort(sorted_elements.begin(), sorted_elements.end());

    values.reserve(sorted_elements.size());
    cols.reserve(sorted_elements.size());

    int current_row = 0;
    rows[0] = 0;

    for (const auto& element : sorted_elements) {
        int i = std::get<0>(element);
        int j = std::get<1>(element);
        double value = std::get<2>(element);

        while (current_row < i) {
            current_row++;
            rows[current_row] = values.size();
        }

        values.push_back(value);
        cols.push_back(j);
    }

    while (current_row < rows_count) {
        current_row++;
        rows[current_row] = values.size();
    }
}

double CSR::operator()(int i, int j) const {

    int row_start = rows[i];
    int row_end = rows[i + 1];

    for (int idx = row_start; idx < row_end; idx++) {
        if (cols[idx] == j) {
            return values[idx];
        }
    }

    return 0.0;
}
// Умножение матрицы на вектор
std::vector<double> CSR::operator*(const std::vector<double>& vec) const {
    std::vector<double> result(rows_count, 0.0);
    
    for (int i = 0; i < rows_count; ++i) {
        double sum = 0.0;
        for (int idx = rows[i]; idx < rows[i + 1]; ++idx) {
            sum += values[idx] * vec[cols[idx]];
        }
        result[i] = sum;
    }

    return result;
}

// Сумма CSR матриц
CSR CSR::operator+(const CSR& other) const {

    std::map<std::tuple<int, int>, double> result_dok;

    for (int i = 0; i < rows_count; ++i) {
        for (int idx = rows[i]; idx < rows[i + 1]; ++idx) {
            int j = cols[idx];
            result_dok[std::make_tuple(i, j)] = values[idx];
        }
    }

    for (int i = 0; i < other.rows_count; ++i) {
        for (int idx = other.rows[i]; idx < other.rows[i + 1]; ++idx) {
            int j = other.cols[idx];
            auto key = std::make_tuple(i, j);
            if (result_dok.find(key) != result_dok.end()) {
                result_dok[key] += other.values[idx];
            } else {
                result_dok[key] = other.values[idx];
            }
        }
    }

    return CSR(result_dok, rows_count, cols_count);
}
// Умножение матрицы на скаляр
CSR CSR::operator*(double scalar) const {
    std::map<std::tuple<int, int>, double> result_dok;

    for (int i = 0; i < rows_count; ++i) {
        for (int idx = rows[i]; idx < rows[i + 1]; ++idx) {
            int j = cols[idx];
            result_dok[std::make_tuple(i, j)] = values[idx] * scalar;
        }
    }

    return CSR(result_dok, rows_count, cols_count);
}