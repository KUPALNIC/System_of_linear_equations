#include "SimpleIterationSolver.hpp"
#include <cmath>
#include <numeric>

// Степенной метод
template <typename T>
double findMaxEigenvalue(const T& A) {
    int n = A.getRows();
    std::vector<double> v(n, 1.0);
    double lambda = 0.0;
    
    for (int i = 0; i < 50; ++i) {
        std::vector<double> Av = A * v;
        
        double norm = 0.0;
        for (double x : Av) norm += x * x;
        norm = std::sqrt(norm);
        
        lambda = norm; 
        for (int j = 0; j < n; ++j) v[j] = Av[j] / norm;
    }
    return lambda;
}

// Версия для Matrix с авто-tau
std::vector<double> SimpleIterationSolver(const Matrix& A, const std::vector<double>& b, double eps, int max_iter) {
    double tau = 1.0 / findMaxEigenvalue(A);
    return SimpleIterationSolver(A, b, tau, eps, max_iter);
}

// Версия для CSR с авто-tau
std::vector<double> SimpleIterationSolver(const CSR& A, const std::vector<double>& b, double eps, int max_iter) {
    double tau = 1.0 / findMaxEigenvalue(A);
    return SimpleIterationSolver(A, b, tau, eps, max_iter);
}

std::vector<double> SimpleIterationSolver(const Matrix& A, const std::vector<double>& b, double tau, double eps, int max_iter) {
    int n = A.getRows();
    std::vector<double> x(n, 0.0);
    double eps_sq = eps * eps;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> Ax = A * x;
        std::vector<double> r(n);
        double norm_sq = 0.0;

        for (int i = 0; i < n; ++i) {
            r[i] = Ax[i] - b[i];
            norm_sq += r[i] * r[i];
        }

        if (norm_sq < eps_sq) break;

        for (int i = 0; i < n; ++i) {
            x[i] -= tau * r[i];
        }
    }
    return x;
}

std::vector<double> SimpleIterationSolver(const CSR& A, const std::vector<double>& b, double tau, double eps, int max_iter) {
    int n = A.getRows();
    std::vector<double> x(n, 0.0);
    double eps_sq = eps * eps;

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<double> Ax = A * x; 
        std::vector<double> r(n);
        double norm_sq = 0.0;

        for (int i = 0; i < n; ++i) {
            r[i] = Ax[i] - b[i];
            norm_sq += r[i] * r[i];
        }

        if (norm_sq < eps_sq) break;

        for (int i = 0; i < n; ++i) {
            x[i] -= tau * r[i];
        }
    }
    return x;
}
