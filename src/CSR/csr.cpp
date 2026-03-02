#include "csr.hpp"
#include <algorithm>
#include <stdexcept>

CSR::CSR(const std::map<std::tuple<int, int>, double>& dok, int rows, int cols): rows_count(rows), cols_count(cols) {

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