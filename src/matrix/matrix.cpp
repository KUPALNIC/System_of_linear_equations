#include "matrix.hpp"

Matrix::Matrix(const std::vector<double>& data_, int r, int c) {    
    rows_count = r;
    cols_count = c;
    data = data_;
}

Matrix::Matrix(size_t r, size_t c) : rows_count(r), cols_count(c), data(r * c, 0.0) {}

double& Matrix::operator()(int i, int j) {
    return data[i*cols_count+j];
}

double Matrix::operator()(int i, int j) const {
    return data[i*cols_count+j];
}

// Сумма матриц
Matrix Matrix::operator+(const Matrix& other) const {
    std::vector<double> result(rows_count*cols_count);
    for (int i = 0; i < rows_count; ++i) {
        for (int j = 0; j < cols_count; ++j) {
            result[i*cols_count+j] = data[i*cols_count+j] + other(i, j);
        }
    }
    return Matrix(result, rows_count, cols_count);
}

// Умножение матрицы на вектор
std::vector<double> Matrix::operator*(const std::vector<double>& vec) const {
    std::vector<double> result(rows_count, 0.0);
    for (int i = 0; i < rows_count; ++i) {
        for (int j = 0; j < cols_count; ++j) {
            result[i] += data[i*cols_count+j] * vec[j];
        }
    }
    return result;
}

// Умножение матрицы на скаляр
Matrix Matrix::operator*(double scalar) const {
    std::vector<double> result(rows_count*cols_count);
    for (int i = 0; i < rows_count; ++i) {
        for (int j = 0; j < cols_count; ++j) {
            result[i*cols_count+j] = data[i*cols_count+j] * scalar;
        }
    }
    return Matrix(result, rows_count, cols_count);
}

//Транспонирование
Matrix Matrix::transpose() const {
    Matrix result(cols_count, rows_count);
    for (size_t i = 0; i < rows_count; ++i) {
        for (size_t j = 0; j < cols_count; ++j) {
            result(j, i) = (*this)(i, j);
        }
    }

    return result;
}
