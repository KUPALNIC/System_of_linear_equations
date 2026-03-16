#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <map>
#include <tuple>
#include "../src/Solvers/SimpleIterationSolver.hpp"

void ExpectVectorNear(const std::vector<double>& actual, const std::vector<double>& expected, double tolerance) {
    ASSERT_EQ(actual.size(), expected.size());
    for (size_t i = 0; i < actual.size(); ++i) {
        EXPECT_NEAR(actual[i], expected[i], tolerance);
    }
}

// Тест для плотной матрицы
TEST(SimpleIterationTest, DenseMatrixAutoTau) {
    Matrix A(2, 2);
    A(0, 0) = 4.0; A(0, 1) = 1.0;
    A(1, 0) = 1.0; A(1, 1) = 3.0;

    std::vector<double> b = {5.0, 4.0};
    std::vector<double> expected = {1.0, 1.0};

    double eps = 1e-7;
    int max_iter = 1000;

    std::vector<double> result = SimpleIterationSolver(A, b, eps, max_iter);
    ExpectVectorNear(result, expected, 1e-5);
}

// Тест для CSR матрицы (через std::map/DOK)
TEST(SimpleIterationTest, CSRMatrixAutoTau) {
    // Используем ваш конструктор: std::map<std::tuple<int, int>, double>
    std::map<std::tuple<int, int>, double> dok;
    dok[{0, 0}] = 10.0;
    dok[{0, 1}] = 2.0;
    dok[{1, 0}] = 2.0;
    dok[{1, 1}] = 10.0;

    CSR A(dok, 2, 2);

    std::vector<double> b = {12.0, 12.0};
    std::vector<double> expected = {1.0, 1.0};

    double eps = 1e-8;
    int max_iter = 1000;

    std::vector<double> result = SimpleIterationSolver(A, b, eps, max_iter);
    ExpectVectorNear(result, expected, 1e-6);
}

// Проверка корректности работы
TEST(SimpleIterationTest, EigenvalueCalculation) {
    Matrix A(2, 2);
    A(0, 0) = 2.0; A(0, 1) = 0.0;
    A(1, 0) = 0.0; A(1, 1) = 1.0;

    std::vector<double> b = {2.0, 1.0};
    std::vector<double> result = SimpleIterationSolver(A, b, 1e-9, 500);

    EXPECT_NEAR(result[0], 1.0, 1e-7);
    EXPECT_NEAR(result[1], 1.0, 1e-7);
}

// Тест на разреженную систему большой размерности
TEST(SimpleIterationTest, LargeSparseSystem) {
    int N = 100;
    std::map<std::tuple<int, int>, double> dok;

    for(int i = 0; i < N; ++i) {
        dok[{i, i}] = 5.0; // Диагональная матрица
    }
    
    CSR A(dok, N, N);
    std::vector<double> b(N, 10.0);
    std::vector<double> expected(N, 2.0);

    std::vector<double> result = SimpleIterationSolver(A, b, 1e-7, 1000);

    for(int i = 0; i < N; ++i) {
        EXPECT_NEAR(result[i], 2.0, 1e-5);
    }
}
