#pragma once

#include "../matrix/matrix.hpp"
#include "../CSR/csr.hpp"
#include <vector>

std::vector<double> SymmetricGaussSeidelStep(
    const Matrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x,
    double omega = 1.0
);

std::vector<double> SymmetricGaussSeidelStep(
    const CSR& A,
    const std::vector<double>& b,
    const std::vector<double>& x,
    double omega = 1.0
);

template<typename MatrixType>
std::vector<double> SymmetricGaussSeidelSolver(
    const MatrixType& A,
    const std::vector<double>& b,
    double eps,
    int max_iter,
    double omega = 1.0
) {
    const int n = A.getRows();
    std::vector<double> x(n, 0.0);
    const double eps_sq = eps * eps;
    
    auto step = [&](const std::vector<double>& x_curr) {
        return SymmetricGaussSeidelStep(A, b, x_curr, omega);
    };
    
    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> x_new = step(x);
        
        double norm_sq = 0.0;
        for (int i = 0; i < n; ++i) {
            double r = x_new[i] - x[i];
            norm_sq += r * r;
        }
        
        x = std::move(x_new);
        if (norm_sq < eps_sq) break;
    }
    
    return x;
}