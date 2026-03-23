#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <map>
#include <tuple>
#include "../src/Solvers/ChebishevIterationSolver.cpp"
#include "../src/matrix/matrix.hpp"
#include "../src/CSR/csr.hpp"

// Вспомогательная функция для вычисления нормы разности векторов
double vector_norm_diff(const std::vector<double>& v1, const std::vector<double>& v2) {
    double sum = 0;
    for (size_t i = 0; i < v1.size(); ++i) {
        double d = v1[i] - v2[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}

class ChebyshevTest : public ::testing::Test {
protected:
    int n = 4;
    std::vector<double> x_true = {1.0, -1.0, 2.0, 0.5};
    double eps = 1e-7;
    
    std::map<std::tuple<int, int>, double> get_test_map() {
        std::map<std::tuple<int, int>, double> dok;
        for (int i = 0; i < n; ++i) {
            dok[{i, i}] = 10.0; 
            if (i > 0) dok[{i, i-1}] = 1.0;
            if (i < n-1) dok[{i, i+1}] = 1.0;
        }
        return dok;
    }
};

TEST_F(ChebyshevTest, MatrixConvergence) {
    auto dok = get_test_map();
    Matrix A(n, n);
    for (auto const& [pos, val] : dok) {
        A(std::get<0>(pos), std::get<1>(pos)) = val;
    }

    std::vector<double> b = A * x_true;

    double l_min = 8.0;
    double l_max = 12.0;

    std::vector<double> x_res = ChebyshevIterationSolver(A, b, l_min, l_max, eps, 1000, 8);

    EXPECT_NEAR(vector_norm_diff(x_res, x_true), 0.0, eps * 100);
}

TEST_F(ChebyshevTest, CSRConvergence) {
    CSR A(get_test_map(), n, n);
    std::vector<double> b = A * x_true;
    
    double l_min = 8.0;
    double l_max = 12.0;

    std::vector<double> x_res = ChebyshevIterationSolver(A, b, l_min, l_max, eps, 1000, 16);

    EXPECT_NEAR(vector_norm_diff(x_res, x_true), 0.0, eps * 100);
}

TEST_F(ChebyshevTest, StabilityTestLargeN) {
    CSR A(get_test_map(), n, n);
    std::vector<double> b = A * x_true;
    
    std::vector<double> x_res = ChebyshevIterationSolver(A, b, 8.0, 12.0, eps, 1000, 64);

    for (double val : x_res) {
        EXPECT_FALSE(std::isnan(val));
        EXPECT_FALSE(std::isinf(val));
    }
    EXPECT_NEAR(vector_norm_diff(x_res, x_true), 0.0, eps * 100);
}