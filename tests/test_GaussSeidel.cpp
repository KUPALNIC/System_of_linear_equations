#include <gtest/gtest.h>
#include "../src/Solvers/GaussSeidelSolver.hpp"
#include "../src/matrix/matrix.hpp"
#include "../src/CSR/csr.hpp"

#include <map>
#include <tuple>
#include <vector>
#include <cmath>

TEST(GaussSeidelDense, Solve2x2)
{
    Matrix A(2,2);

    A(0,0)=4;
    A(0,1)=1;
    A(1,0)=2;
    A(1,1)=3;

    std::vector<double> b = {1,2};

    auto x = GaussSeidelSolver(A,b,1e-8,1000);

    EXPECT_NEAR(x[0],0.1,1e-3);
    EXPECT_NEAR(x[1],0.6,1e-3);
}

TEST(GaussSeidelDense, Solve3x3)
{
    Matrix A(3,3);

    A(0,0)=10; A(0,1)=1; A(0,2)=1;
    A(1,0)=2;  A(1,1)=10; A(1,2)=1;
    A(2,0)=2;  A(2,1)=2;  A(2,2)=10;

    std::vector<double> b = {12,13,14};

    auto x = GaussSeidelSolver(A,b,1e-8,1000);

    EXPECT_NEAR(x[0],1.0,1e-3);
    EXPECT_NEAR(x[1],1.0,1e-3);
    EXPECT_NEAR(x[2],1.0,1e-3);
}

TEST(GaussSeidelDense, ResidualSmall)
{
    Matrix A(2,2);

    A(0,0)=5;
    A(0,1)=1;
    A(1,0)=2;
    A(1,1)=6;

    std::vector<double> b = {6,8};

    auto x = GaussSeidelSolver(A,b,1e-8,1000);

    auto Ax = A * x;

    double r = 0;

    for(int i=0;i<2;i++)
    {
        double d = Ax[i]-b[i];
        r += d*d;
    }

    r = std::sqrt(r);

    EXPECT_LT(r,1e-6);
}


// CSR

TEST(GaussSeidelCSR, Solve2x2)
{
    std::map<std::tuple<int,int>,double> dok;

    dok[{0,0}]=4;
    dok[{0,1}]=1;
    dok[{1,0}]=2;
    dok[{1,1}]=3;

    CSR A(dok,2,2);

    std::vector<double> b={1,2};

    auto x = GaussSeidelSolver(A,b,1e-8,1000);

    EXPECT_NEAR(x[0],0.1,1e-3);
    EXPECT_NEAR(x[1],0.6,1e-3);
}

TEST(GaussSeidelCSR, ResidualSmall)
{
    std::map<std::tuple<int,int>,double> dok;

    dok[{0,0}]=5;
    dok[{0,1}]=1;
    dok[{1,0}]=2;
    dok[{1,1}]=6;

    CSR A(dok,2,2);

    std::vector<double> b={6,8};

    auto x = GaussSeidelSolver(A,b,1e-8,1000);

    auto Ax = A*x;

    double r=0;

    for(int i=0;i<2;i++)
    {
        double d = Ax[i]-b[i];
        r += d*d;
    }

    r = std::sqrt(r);

    EXPECT_LT(r,1e-6);
}
