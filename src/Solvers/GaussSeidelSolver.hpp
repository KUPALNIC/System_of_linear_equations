#pragma once

#include "../matrix/matrix.hpp"
#include "../CSR/csr.hpp"
#include <vector>

std::vector<double> GaussSeidelSolver(
    const Matrix& A,
    const std::vector<double>& b,
    double eps,
    int max_iter
);

std::vector<double> GaussSeidelSolver(
    const CSR& A,
    const std::vector<double>& b,
    double eps,
    int max_iter
);