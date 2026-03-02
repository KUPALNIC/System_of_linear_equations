#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include "../src/tomas.hpp"

// Тест 1: Простая система с известным решением
TEST(ThreeMatrixTest, SimpleSystem) {
    std::vector<double> a = {0, -1, -1};
    std::vector<double> b = {2, 2, 2};
    std::vector<double> c = {-1, -1, 0};
    std::vector<double> d = {1, 0, 1};
    
    three_matrix matrix(a, b, c);
    std::vector<double> solution = matrix.solve(d);
    
    EXPECT_EQ(solution.size(), 3);
    EXPECT_NEAR(solution[0], 1.0, 1e-6);
    EXPECT_NEAR(solution[1], 1.0, 1e-6);
    EXPECT_NEAR(solution[2], 1.0, 1e-6);
}

// Тест 2: Единичная матрица (диагональная)
TEST(ThreeMatrixTest, IdentitySystem) {
    std::vector<double> a = {0, 0, 0};
    std::vector<double> b = {1, 1, 1};
    std::vector<double> c = {0, 0, 0};
    std::vector<double> d = {2, 3, 4};
    
    three_matrix matrix(a, b, c);
    std::vector<double> solution = matrix.solve(d);
    
    EXPECT_EQ(solution.size(), 3);
    EXPECT_NEAR(solution[0], 2.0, 1e-6);
    EXPECT_NEAR(solution[1], 3.0, 1e-6);
    EXPECT_NEAR(solution[2], 4.0, 1e-6);
}

// Тест 3: Размер вектора не совпадает с матрицей
TEST(ThreeMatrixTest, MismatchedDiagonalSizes) {
    std::vector<double> a = {0, -1};
    std::vector<double> b = {2, 2, 2};
    std::vector<double> c = {-1, -1, 0};
    
    EXPECT_THROW(three_matrix matrix(a, b, c), std::invalid_argument);
}

// Тест 4: Неправильный размер вектора решения
TEST(ThreeMatrixTest, WrongRightHandSideSize) {
    std::vector<double> a = {0, -1, -1};
    std::vector<double> b = {2, 2, 2};
    std::vector<double> c = {-1, -1, 0};
    std::vector<double> d = {1, 0};  // Неверный размер
    
    three_matrix matrix(a, b, c);
    EXPECT_THROW(matrix.solve(d), std::invalid_argument);
}

// Тест 5: Система размером 5x5
TEST(ThreeMatrixTest, LargerSystem) {
    std::vector<double> a = {0, -1, -1, -1, -1};
    std::vector<double> b = {3, 3, 3, 3, 3};
    std::vector<double> c = {-1, -1, -1, -1, 0};
    std::vector<double> d = {2, 1, 1, 1, 2};
    
    three_matrix matrix(a, b, c);
    std::vector<double> solution = matrix.solve(d);
    
    EXPECT_EQ(solution.size(), 5);
    for (size_t i = 0; i < solution.size(); ++i) {
        EXPECT_FALSE(std::isnan(solution[i]));
        EXPECT_FALSE(std::isinf(solution[i]));
    }
}

// Тест 6: Система размером 2x2
TEST(ThreeMatrixTest, TwoByTwoSystem) {
    std::vector<double> a = {0, -2};
    std::vector<double> b = {4, 4};
    std::vector<double> c = {-2, 0};
    std::vector<double> d = {2, 2};
    
    three_matrix matrix(a, b, c);
    std::vector<double> solution = matrix.solve(d);
    
    EXPECT_EQ(solution.size(), 2);
    EXPECT_NEAR(solution[0], 0.5, 1e-6);
    EXPECT_NEAR(solution[1], 0.5, 1e-6);
}

// Тест 7: Система с большими числами
TEST(ThreeMatrixTest, LargeNumbers) {
    std::vector<double> a = {0, -100, -100};
    std::vector<double> b = {200, 200, 200};
    std::vector<double> c = {-100, -100, 0};
    std::vector<double> d = {100, 0, 100};
    
    three_matrix matrix(a, b, c);
    std::vector<double> solution = matrix.solve(d);
    
    EXPECT_EQ(solution.size(), 3);
    EXPECT_NEAR(solution[0], 1.0, 1e-5);
    EXPECT_NEAR(solution[1], 1.0, 1e-5);
    EXPECT_NEAR(solution[2], 1.0, 1e-5);
}

// Тест 8: Система с отрицательными значениями
TEST(ThreeMatrixTest, NegativeValues) {
    std::vector<double> a = {0, -1, -1};
    std::vector<double> b = {2, 2, 2};
    std::vector<double> c = {-1, -1, 0};
    std::vector<double> d = {-4, -2, -4};
    
    three_matrix matrix(a, b, c);
    std::vector<double> solution = matrix.solve(d);
    
    EXPECT_EQ(solution.size(), 3);
    EXPECT_NEAR(solution[0], -2.0, 1e-6);
    EXPECT_NEAR(solution[1], -2.0, 1e-6);
    EXPECT_NEAR(solution[2], -2.0, 1e-6);
}

// Тест 9: Проверка диагонального преобладания - матрица с преобладанием
TEST(ThreeMatrixTest, DiagonalDominanceCheck_Valid) {
    std::vector<double> a = {0, -1, -1};
    std::vector<double> b = {2, 2, 2};
    std::vector<double> c = {-1, -1, 0};
    
    three_matrix matrix(a, b, c);
    EXPECT_TRUE(matrix.is_diagonally_dominant());
}

// Тест 10: Проверка диагонального преобладания - нарушение условия
TEST(ThreeMatrixTest, DiagonalDominanceCheck_Invalid) {
    std::vector<double> a = {0, -5, -1};
    std::vector<double> b = {2, 2, 2};
    std::vector<double> c = {-5, -1, 0};
    
    // Должно выбросить исключение при создании матрицы
    EXPECT_THROW(three_matrix matrix(a, b, c), std::invalid_argument);
}

// Тест 11: Граничный случай - диагональное преобладание на границе
TEST(ThreeMatrixTest, DiagonalDominanceCheck_Boundary) {
    std::vector<double> a = {0, -1, -1};
    std::vector<double> b = {1, 2, 2};
    std::vector<double> c = {-1, -1, 0};
    
    three_matrix matrix(a, b, c);
    EXPECT_TRUE(matrix.is_diagonally_dominant());
}

// Тест 12: Строгое диагональное преобладание
TEST(ThreeMatrixTest, StrictDiagonalDominance) {
    std::vector<double> a = {0, -1, -1};
    std::vector<double> b = {5, 5, 5};
    std::vector<double> c = {-1, -1, 0};
    
    three_matrix matrix(a, b, c);
    EXPECT_TRUE(matrix.is_diagonally_dominant());
}