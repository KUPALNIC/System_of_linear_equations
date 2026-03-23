#include "SimpleIterationSolver.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

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

void fillIndices(std::vector<int>& indices, int n, int step) {
    if (step == n) return;
    int current_size = step;
    for (int i = 0; i < current_size; ++i) {
        indices[2 * current_size - 1 - i] = 2 * n / (2 * current_size) - 1 - indices[i];
    }
    fillIndices(indices, n, 2 * step);
}


std::vector<double> getChebyshevTaus(double lambda_min, double lambda_max, int n) {
    std::vector<double> taus(n);
    std::vector<int> indices(n, 0);
    indices[1] = 1; 
    if (n > 2) fillIndices(indices, n, 2);

    double half_sum = (lambda_max + lambda_min) / 2.0;
    double half_diff = (lambda_max - lambda_min) / 2.0;

    for (int k = 0; k < n; ++k) {
        double cos_val = std::cos(M_PI * (2.0 * indices[k] + 1.0) / (2.0 * n));
        double lambda_k = half_sum + half_diff * cos_val;
        taus[k] = 1.0 / lambda_k;
    }
    return taus;
}

template <typename T>
std::vector<double> ChebyshevIterationSolver(const T& A, const std::vector<double>& b, 
                                             double lambda_min, double lambda_max, 
                                             double eps, int max_iter, int n_roots = 16) {
    int n = A.getRows();
    std::vector<double> x(n, 0.0);
    double eps_sq = eps * eps;

    std::vector<double> taus = getChebyshevTaus(lambda_min, lambda_max, n_roots);

    for (int iter = 0; iter < max_iter; ) {
        for (int k = 0; k < n_roots && iter < max_iter; ++k, ++iter) {
            std::vector<double> Ax = A * x;
            double norm_sq = 0.0;
            
            for (int i = 0; i < n; ++i) {
                double r_i = Ax[i] - b[i];
                norm_sq += r_i * r_i;
                x[i] -= taus[k] * r_i;
            }
            if (norm_sq < eps_sq) return x;
        }
    }
    return x;
}

std::vector<double> ChebyshevIterationSolverAuto(const Matrix& A, const std::vector<double>& b, double eps, int max_iter) {
    double l_max = findMaxEigenvalue(A);
    double l_min = l_max * 1e-3;
    return ChebyshevIterationSolver(A, b, l_min, l_max, eps, max_iter);
}
