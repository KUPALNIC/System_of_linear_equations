#pragma once 
#include <vector>
#include <map>
#include <tuple>

class CSR {
private:
    std::vector<double> values;
    std::vector<int> cols;
    std::vector<int> rows;

    int rows_count;
    int cols_count;

public:
    CSR(const std::map<std::tuple<int, int>, double>& dok, int rows, int cols);

    double operator()(int i, int j) const;

    int getRows() const { return rows_count; }
    int getCols() const { return cols_count; }
    std::vector<int> Rows() const {return rows;}
    std::vector<int> Cols() const {return cols;}
    std::vector<double> Values() const {return values;}
    int getNonZeroCount() const { return values.size(); }

    // Умножение матрицы на вектор
    std::vector<double> operator*(const std::vector<double>& vec) const;

    // Сумма матриц
    CSR operator+(const CSR& other) const;

    // Умножение на скаляр
    CSR operator*(double scalar) const;
};