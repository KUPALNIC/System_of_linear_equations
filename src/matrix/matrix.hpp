#pragma once
#include <vector>
#include <stdexcept>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    int rows_count;
    int cols_count;

public:
    Matrix(const std::vector<std::vector<double>>& data);

    double operator()(int i, int j) const;
    double& operator()(int i, int j);

    int getRows() const { return rows_count; }
    int getCols() const { return cols_count; }

    // Сумма матриц
    Matrix operator+(const Matrix& other) const;
    
    // Умножение матрицы на вектор
    std::vector<double> operator*(const std::vector<double>& vec) const;
    
    // Умножение матрицы на скаляр
    Matrix operator*(double scalar) const;
};