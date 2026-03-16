#pragma once

#include "../matrix/matrix.hpp"
#include "../CSR/csr.hpp"
#include <vector>

std::vector<double> SimpleIterationSolver(
    const Matrix& A,
     const std::vector<double>& b,
      double tau,
       double eps,
        int max_iter
);

std::vector<double> SimpleIterationSolver(
    const CSR& A,
     const std::vector<double>& b,
      double tau,
       double eps,
        int max_iter
);


std::vector<double> SimpleIterationSolver(
    const Matrix& A,
     const std::vector<double>& b,
      double eps,
       int max_iter
);

std::vector<double> SimpleIterationSolver(
    const CSR& A,
     const std::vector<double>& b,
      double eps,
       int max_iter
);

