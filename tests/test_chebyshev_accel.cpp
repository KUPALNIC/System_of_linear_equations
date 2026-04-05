#include <gtest/gtest.h>
#include "../src/Solvers/ChebyshevAcceleration.hpp"
#include "../src/Solvers/SymmetricGaussSeidel.hpp"
#include "../src/matrix/matrix.hpp"
#include "../src/CSR/csr.hpp"

#include <map>
#include <tuple>
#include <vector>
#include <cmath>

// Вспомогательная: норма невязки ||Ax - b||
double residual_norm(const Matrix& A, const std::vector<double>& x, const std::vector<double>& b) {
    auto Ax = A * x;
    double r = 0;
    for (size_t i = 0; i < b.size(); ++i) {
        double d = Ax[i] - b[i];
        r += d * d;
    }
    return std::sqrt(r);
}

double residual_norm(const CSR& A, const std::vector<double>& x, const std::vector<double>& b) {
    auto Ax = A * x;
    double r = 0;
    for (size_t i = 0; i < b.size(); ++i) {
        double d = Ax[i] - b[i];
        r += d * d;
    }
    return std::sqrt(r);
}

// ==================== Matrix: Symmetric Gauss-Seidel + Chebyshev ====================

TEST(ChebyshevMatrix, SGS_Converges_3x3) {
    Matrix A(3, 3);
    A(0,0) = 4; A(0,1) = -1; A(0,2) = 0;
    A(1,0) = -1; A(1,1) = 4; A(1,2) = -1;
    A(2,0) = 0; A(2,1) = -1; A(2,2) = 4;

    std::vector<double> b = {1, 2, 3};
    std::vector<double> x0(3, 0.0);

    // Лямбда для шага SGS
    auto sgs_step = [&](const std::vector<double>& x) {
        return SymmetricGaussSeidelStep(A, b, x, 1.0);
    };

    double rho = 0.8; // оценка спектрального радиуса
    auto x = ChebyshevAccelerate(sgs_step, x0, rho, 1e-8, 1000);

    double r = residual_norm(A, x, b);
    EXPECT_LT(r, 1e-6);
}

TEST(ChebyshevMatrix, SGS_Converges_4x4) {
    Matrix A(4, 4);
    A(0,0)=10; A(0,1)=-1; A(0,2)=0;  A(0,3)=0;
    A(1,0)=-1; A(1,1)=10; A(1,2)=-1; A(1,3)=0;
    A(2,0)=0;  A(2,1)=-1; A(2,2)=10; A(2,3)=-1;
    A(3,0)=0;  A(3,1)=0;  A(3,2)=-1; A(3,3)=10;

    std::vector<double> b = {9, 8, 8, 9}; // решение ~ [1,1,1,1]
    std::vector<double> x0(4, 0.0);

    auto sgs_step = [&](const std::vector<double>& x) {
        return SymmetricGaussSeidelStep(A, b, x, 1.0);
    };

    double rho = 0.85;
    auto x = ChebyshevAccelerate(sgs_step, x0, rho, 1e-8, 1000);

    EXPECT_NEAR(x[0], 1.0, 1e-3);
    EXPECT_NEAR(x[1], 1.0, 1e-3);
    EXPECT_NEAR(x[2], 1.0, 1e-3);
    EXPECT_NEAR(x[3], 1.0, 1e-3);
}

TEST(ChebyshevMatrix, AutoRho_Converges) {
    Matrix A(3, 3);
    A(0,0) = 4; A(0,1) = -1; A(0,2) = 0;
    A(1,0) = -1; A(1,1) = 4; A(1,2) = -1;
    A(2,0) = 0; A(2,1) = -1; A(2,2) = 4;

    std::vector<double> b = {1, 2, 3};
    std::vector<double> x0(3, 0.0);

    auto sgs_step = [&](const std::vector<double>& x) {
        return SymmetricGaussSeidelStep(A, b, x, 1.0);
    };

    // Авто-версия сама оценит спектральный радиус
    auto x = ChebyshevAccelerateAuto(sgs_step, x0, 1e-8, 1000);

    double r = residual_norm(A, x, b);
    EXPECT_LT(r, 1e-6);
}

// ==================== CSR: Symmetric Gauss-Seidel + Chebyshev ====================

TEST(ChebyshevCSR, SGS_Converges_3x3) {
    std::map<std::tuple<int,int>, double> dok;
    dok[{0,0}] = 4; dok[{0,1}] = -1;
    dok[{1,0}] = -1; dok[{1,1}] = 4; dok[{1,2}] = -1;
    dok[{2,1}] = -1; dok[{2,2}] = 4;

    CSR A(dok, 3, 3);
    std::vector<double> b = {1, 2, 3};
    std::vector<double> x0(3, 0.0);

    auto sgs_step = [&](const std::vector<double>& x) {
        return SymmetricGaussSeidelStep(A, b, x, 1.0);
    };

    double rho = 0.8;
    auto x = ChebyshevAccelerate(sgs_step, x0, rho, 1e-8, 1000);

    double r = residual_norm(A, x, b);
    EXPECT_LT(r, 1e-6);
}

TEST(ChebyshevCSR, SGS_Converges_4x4) {
    std::map<std::tuple<int,int>, double> dok;
    dok[{0,0}] = 10; dok[{0,1}] = -1;
    dok[{1,0}] = -1; dok[{1,1}] = 10; dok[{1,2}] = -1;
    dok[{2,1}] = -1; dok[{2,2}] = 10; dok[{2,3}] = -1;
    dok[{3,2}] = -1; dok[{3,3}] = 10;

    CSR A(dok, 4, 4);
    std::vector<double> b = {9, 8, 8, 9};
    std::vector<double> x0(4, 0.0);

    auto sgs_step = [&](const std::vector<double>& x) {
        return SymmetricGaussSeidelStep(A, b, x, 1.0);
    };

    double rho = 0.85;
    auto x = ChebyshevAccelerate(sgs_step, x0, rho, 1e-8, 1000);

    EXPECT_NEAR(x[0], 1.0, 1e-3);
    EXPECT_NEAR(x[1], 1.0, 1e-3);
    EXPECT_NEAR(x[2], 1.0, 1e-3);
    EXPECT_NEAR(x[3], 1.0, 1e-3);
}

TEST(ChebyshevCSR, AutoRho_Converges) {
    std::map<std::tuple<int,int>, double> dok;
    dok[{0,0}] = 4; dok[{0,1}] = -1;
    dok[{1,0}] = -1; dok[{1,1}] = 4; dok[{1,2}] = -1;
    dok[{2,1}] = -1; dok[{2,2}] = 4;

    CSR A(dok, 3, 3);
    std::vector<double> b = {1, 2, 3};
    std::vector<double> x0(3, 0.0);

    auto sgs_step = [&](const std::vector<double>& x) {
        return SymmetricGaussSeidelStep(A, b, x, 1.0);
    };

    auto x = ChebyshevAccelerateAuto(sgs_step, x0, 1e-8, 1000);

    double r = residual_norm(A, x, b);
    EXPECT_LT(r, 1e-6);
}