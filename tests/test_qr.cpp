#include <gtest/gtest.h>

#include "../src/matrix/matrix.hpp"
#include "../src/matrix/QR/qr_solve.hpp"

#include <vector>
#include <cmath>

double eps = 1e-6;

bool vec_eq(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size())
        return false;

    for (size_t i = 0; i < a.size(); ++i)
        if (std::abs(a[i] - b[i]) > eps)
            return false;

    return true;
}


TEST(QR, simple_sys)
{
    Matrix A({
        2, 1,
        1, 3
    }, 2, 2);

    std::vector<double> b = {1, 2};

    QRSolve s(A);

    std::vector<double> x = s.solve(b);

    std::vector<double> ax = A * x;

    EXPECT_TRUE(vec_eq(ax, b));
}


TEST(QR, diag_sys)
{
    Matrix A({
        2, 0, 0,
        0, 3, 0,
        0, 0, 4
    }, 3, 3);

    std::vector<double> b = {2, 6, 8};

    QRSolve s(A);

    std::vector<double> x = s.solve(b);

    std::vector<double> exp = {1, 2, 2};

    EXPECT_TRUE(vec_eq(x, exp));
}


TEST(QR, sym_sys)
{
    Matrix A({
        2,-1,0,
       -1,2,-1,
        0,-1,2
    }, 3, 3);

    std::vector<double> b = {1,0,1};

    QRSolve s(A);

    std::vector<double> x = s.solve(b);

    std::vector<double> ax = A * x;

    EXPECT_TRUE(vec_eq(ax, b));
}


TEST(QR, identity)
{
    Matrix A({
        1,0,0,
        0,1,0,
        0,0,1
    },3,3);

    std::vector<double> b = {5,7,9};

    QRSolve s(A);

    std::vector<double> x = s.solve(b);

    EXPECT_TRUE(vec_eq(x,b));
}