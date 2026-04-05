#include "SymmetricGaussSeidel.hpp"
#include <cmath>

std::vector<double> SymmetricGaussSeidelStep(
    const Matrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x,
    double omega
) {
    const int n = A.getRows();
    std::vector<double> x_half = x;
    
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) sum += A(i, j) * x_half[j];
        }
        x_half[i] = (1.0 - omega) * x[i] + omega * (b[i] - sum) / A(i, i);
    }
    
    std::vector<double> x_new = x_half;
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) sum += A(i, j) * x_new[j];
        }
        x_new[i] = (1.0 - omega) * x_half[i] + omega * (b[i] - sum) / A(i, i);
    }
    
    return x_new;
}

std::vector<double> SymmetricGaussSeidelStep(
    const CSR& A,
    const std::vector<double>& b,
    const std::vector<double>& x,
    double omega
) {
    const int n = A.getRows();
    std::vector<double> x_half = x;
    
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        double diag = 1.0;
        for (int idx = A.Rows()[i]; idx < A.Rows()[i + 1]; ++idx) {
            int j = A.Cols()[idx];
            double val = A.Values()[idx];
            if (j == i) diag = val;
            else sum += val * x_half[j];
        }
        x_half[i] = (1.0 - omega) * x[i] + omega * (b[i] - sum) / diag;
    }
    
    std::vector<double> x_new = x_half;
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        double diag = 1.0;
        for (int idx = A.Rows()[i]; idx < A.Rows()[i + 1]; ++idx) {
            int j = A.Cols()[idx];
            double val = A.Values()[idx];
            if (j == i) diag = val;
            else sum += val * x_new[j];
        }
        x_new[i] = (1.0 - omega) * x_half[i] + omega * (b[i] - sum) / diag;
    }
    
    return x_new;
}