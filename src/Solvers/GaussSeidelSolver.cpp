#include "GaussSeidelSolver.hpp"
#include <cmath>

std::vector<double> GaussSeidelSolver(
    const Matrix& A,
    const std::vector<double>& b,
    double eps,
    int max_iter)
{
    int n = A.getRows();

    std::vector<double> x(n, 0.0);

    for (int iter = 0; iter < max_iter; iter++)
    {
        double diff = 0.0;

        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;

            for (int j = 0; j < n; j++)
            {
                if (j != i)
                    sum += A(i,j) * x[j];
            }

            double old = x[i];

            x[i] = (b[i] - sum) / A(i,i);

            double d = x[i] - old;
            diff += d*d;
        }

        if (diff < eps*eps)
            break;
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

    for (int iter = 0; iter < max_iter; ++iter)
    {
        double diff = 0.0;

        for (int i = 0; i < n; ++i)
        {
            double sum = 0.0;
            double diag = 0.0;

            int row_start = A.Rows()[i];
            int row_end   = A.Rows()[i+1];

            for (int idx = row_start; idx < row_end; ++idx)
            {
                int j = A.Cols()[idx];
                double val = A.Values()[idx];

                if (j == i)
                    diag = val;
                else
                    sum += val * x[j];
            }

            double old = x[i];

            x[i] = (b[i] - sum) / diag;

            double d = x[i] - old;
            diff += d*d;
        }

        if (diff < eps*eps)
            break;
    }

    return x;
}