#include "GaussSeidelSolver.hpp"
#include <cmath>


double norm(const std::vector<double>& r) {
    double sum = 0.0;
    for (double val : r) sum += val * val;
    return std::sqrt(sum);
}

std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> res = a;
    for (size_t i = 0; i < a.size(); ++i) res[i] -= b[i];
    return res;
}


std::vector<double> GaussSeidelSolver(
    const Matrix& A,
    const std::vector<double>& b,
    double eps,
    int max_iter)
{
    int n = A.getRows();
    std::vector<double> x(n, 0.0);

    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) sum += A(i, j) * x[j];
            }
            x[i] = (b[i] - sum) / A(i, i);
        }

        if (norm(A * x - b) < eps) {
            break;
        }
    }
    return x;
}


std::vector<double> GaussSeidelSolver(
    const CSR& A,
    const std::vector<double>& b,
    double eps,
    int max_iter)
{
    int n = A.getRows();
    std::vector<double> x(n, 0.0);

    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            double diag = 1.0;
            for (int idx = A.Rows()[i]; idx < A.Rows()[i+1]; ++idx) {
                int j = A.Cols()[idx];
                double val = A.Values()[idx];
                if (j == i) diag = val;
                else sum += val * x[j];
            }
            x[i] = (b[i] - sum) / diag;
        }

        if (norm(A * x - b) < eps) {
            break;
        }
    }
    return x;
}