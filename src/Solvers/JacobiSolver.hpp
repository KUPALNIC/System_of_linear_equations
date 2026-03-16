#pragma once

#include "../matrix/matrix.hpp"
#include "../CSR/csr.hpp"
#include <vector>

std::vector<double> JacobiSolver(
    const Matrix& A,
    const std::vector<double>& b,
    double eps,
    int max_iter
);

std::vector<double> JacobiSolver(
    const CSR& A,
    const std::vector<double>& b,
    double eps,
    int max_iter
);