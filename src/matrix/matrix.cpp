#include "matrix.hpp"

Matrix::Matrix(const std::vector<std::vector<double>>& data_) {    
    rows_count = data_.size();
    cols_count = data_[0].size();
    data = data_;
}

double& Matrix::operator()(int i, int j) {
    return data[i][j];
}

double Matrix::operator()(int i, int j) const {
    return data[i][j];
}

// Сумма матриц
Matrix Matrix::operator+(const Matrix& other) const {
    std::vector<std::vector<double>> result(rows_count, std::vector<double>(cols_count));
    for (int i = 0; i < rows_count; ++i) {
        for (int j = 0; j < cols_count; ++j) {
            result[i][j] = data[i][j] + other(i, j);
        }
    }
    return Matrix(result);
}

// Умножение матрицы на вектор
std::vector<double> Matrix::operator*(const std::vector<double>& vec) const {
    std::vector<double> result(rows_count, 0.0);
    for (int i = 0; i < rows_count; ++i) {
        for (int j = 0; j < cols_count; ++j) {
            result[i] += data[i][j] * vec[j];
        }
    }
    return result;
}

// Умножение матрицы на скаляр
Matrix Matrix::operator*(double scalar) const {
    std::vector<std::vector<double>> result(rows_count, std::vector<double>(cols_count));
    for (int i = 0; i < rows_count; ++i) {
        for (int j = 0; j < cols_count; ++j) {
            result[i][j] = data[i][j] * scalar;
        }
    }
    return Matrix(result);
}