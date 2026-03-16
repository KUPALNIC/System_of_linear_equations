#include "JacobiSolver.hpp"
#include <cmath>

std::vector<double> JacobiSolver(
    const Matrix& A,
    const std::vector<double>& b,
    double eps,
    int max_iter)
{
    int n = A.getRows();

    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n, 0.0);

    for (int iter = 0; iter < max_iter; iter++)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;

            for (int j = 0; j < n; j++)
            {
                if (j != i)
                    sum += A(i,j) * x[j];
            }

            double diag = A(i,i);

            x_new[i] = (b[i] - sum) / diag;
        }

        std::vector<double> r = A * x;

        double norm = 0.0;
        
        for (int i = 0; i < n; i++)
        {
            double ri = r[i] - b[i];
            norm += ri * ri;
        }
        
        if (norm < eps*eps)
            break;

        x = x_new;
    }
    return x;
}




std::vector<double> JacobiSolver(
    const CSR& A,
    const std::vector<double>& b,
    double eps,
    int max_iter)
{
    int n = A.getRows();

    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n, 0.0);

    for (int iter = 0; iter < max_iter; ++iter)
    {
        for (int i = 0; i < n; ++i)
        {
            double sum = 0.0;
            double diag = 0.0;

            int row_start = A.Rows()[i];
            int row_end = A.Rows()[i + 1];

            for (int idx = row_start; idx < row_end; ++idx)
            {
                int j = A.Cols()[idx];
                double value = A.Values()[idx];

                if (j == i)
                    diag = value;
                else
                    sum += value * x[j];
            }

            x_new[i] = (b[i] - sum) / diag;
        }

        std::vector<double> r = A * x;

        double norm = 0.0;
        
        for (int i = 0; i < n; i++)
        {
            double ri = r[i] - b[i];
            norm += ri * ri;
        }

        if (norm < eps*eps)
            break;
        
        x = x_new;
    }

    return x;
}