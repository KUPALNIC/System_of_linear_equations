#include <gtest/gtest.h>
#include "../src/Solvers/JacobiSolver.hpp"
#include "../src/CSR/csr.hpp"
#include "../src/matrix/matrix.hpp"
#include <map>
#include <tuple>
#include <vector>
#include <cmath>

TEST(JacobiSolverDense, Solve2x2)
{
    Matrix A(2,2);

    A(0,0) = 4;
    A(0,1) = 1;
    A(1,0) = 2;
    A(1,1) = 3;

    std::vector<double> b = {1,2};

    auto x = JacobiSolver(A, b, 1e-8, 1000);

    EXPECT_NEAR(x[0], 0.1, 1e-3);
    EXPECT_NEAR(x[1], 0.6, 1e-3);
}

TEST(JacobiSolverDense, Solve3x3)
{
    Matrix A(3,3);

    A(0,0)=10; A(0,1)=1; A(0,2)=1;
    A(1,0)=2;  A(1,1)=10; A(1,2)=1;
    A(2,0)=2;  A(2,1)=2;  A(2,2)=10;

    std::vector<double> b = {12,13,14};

    auto x = JacobiSolver(A,b,1e-8,1000);

    EXPECT_NEAR(x[0],1.0,1e-3);
    EXPECT_NEAR(x[1],1.0,1e-3);
    EXPECT_NEAR(x[2],1.0,1e-3);
}

TEST(JacobiSolverDense, ResidualSmall)
{
    Matrix A(2,2);

    A(0,0)=5;
    A(0,1)=1;
    A(1,0)=2;
    A(1,1)=6;

    std::vector<double> b = {6,8};

    auto x = JacobiSolver(A,b,1e-8,1000);

    auto Ax = A * x;

    double residual = 0.0;

    for(int i=0;i<2;i++)
    {
        double r = Ax[i] - b[i];
        residual += r*r;
    }

    residual = std::sqrt(residual);

    EXPECT_LT(residual,1e-6);
}

TEST(JacobiSolverDense, DiagonalMatrix)
{
    Matrix A(3,3);

    A(0,0)=5;
    A(1,1)=7;
    A(2,2)=9;

    std::vector<double> b = {10,14,18};

    auto x = JacobiSolver(A,b,1e-10,100);

    EXPECT_NEAR(x[0],2.0,1e-8);
    EXPECT_NEAR(x[1],2.0,1e-8);
    EXPECT_NEAR(x[2],2.0,1e-8);
}



//test CSR matrix solving by Jacobi

TEST(JacobiSolverCSR, Solve2x2)
{
    std::map<std::tuple<int,int>, double> dok;

    dok[{0,0}] = 4;
    dok[{0,1}] = 1;
    dok[{1,0}] = 2;
    dok[{1,1}] = 3;

    CSR A(dok, 2, 2);

    std::vector<double> b = {1,2};

    auto x = JacobiSolver(A, b, 1e-8, 1000);

    EXPECT_NEAR(x[0], 0.1, 1e-3);
    EXPECT_NEAR(x[1], 0.6, 1e-3);
}

TEST(JacobiSolverCSR, Solve3x3)
{
    std::map<std::tuple<int,int>, double> dok;

    dok[{0,0}] = 10;
    dok[{0,1}] = 1;
    dok[{0,2}] = 1;

    dok[{1,0}] = 2;
    dok[{1,1}] = 10;
    dok[{1,2}] = 1;

    dok[{2,0}] = 2;
    dok[{2,1}] = 2;
    dok[{2,2}] = 10;

    CSR A(dok,3,3);

    std::vector<double> b = {12,13,14};

    auto x = JacobiSolver(A,b,1e-8,1000);

    EXPECT_NEAR(x[0],1.0,1e-3);
    EXPECT_NEAR(x[1],1.0,1e-3);
    EXPECT_NEAR(x[2],1.0,1e-3);
}

TEST(JacobiSolverCSR, ResidualSmall)
{
    std::map<std::tuple<int,int>, double> dok;

    dok[{0,0}] = 5;
    dok[{0,1}] = 1;

    dok[{1,0}] = 2;
    dok[{1,1}] = 6;

    CSR A(dok,2,2);

    std::vector<double> b = {6,8};

    auto x = JacobiSolver(A,b,1e-8,1000);

    auto Ax = A * x;

    double residual = 0.0;

    for(int i=0;i<2;i++)
    {
        double r = Ax[i] - b[i];
        residual += r*r;
    }

    residual = std::sqrt(residual);

    EXPECT_LT(residual,1e-6);
}
