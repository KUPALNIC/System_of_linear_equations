#include <gtest/gtest.h>
#include "../CSR.h"
#include "../vector_oper.h"
#include <map>
#include <tuple>
#include <vector>

// Тест конструктора и получения элементов
TEST(CSRTest, ConstructorAndAccess) {
    std::map<std::tuple<int, int>, double> dok = {
        {{0, 0}, 1.0},
        {{0, 2}, 2.0},
        {{1, 1}, 3.0},
        {{2, 0}, 4.0},
        {{2, 2}, 5.0}
    };
    
    CSR matrix(dok, 3, 3);

    // Проверка существующих элементов
    EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix(0, 2), 2.0);
    EXPECT_DOUBLE_EQ(matrix(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(matrix(2, 0), 4.0);
    EXPECT_DOUBLE_EQ(matrix(2, 2), 5.0);

    // Проверка нулевых элементов
    EXPECT_DOUBLE_EQ(matrix(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(matrix(2, 1), 0.0);
}

// Тест операций из vector_oper (умножение CSR на вектор)
TEST(VectorOperTest, CSRVectorMultiplication) {
    std::map<std::tuple<int, int>, double> dok = {
        {{0, 0}, 1.0}, {{0, 1}, 2.0},
        {{1, 0}, 3.0}, {{1, 1}, 4.0}
    };
    CSR A(dok, 2, 2);
    std::vector<double> x = {1.0, 1.0};
    
    // Ожидаемый результат: [1*1 + 2*1, 3*1 + 4*1] = [3, 7]
    std::vector<double> result = A * x;
    
    ASSERT_EQ(result.size(), 2);
    EXPECT_DOUBLE_EQ(result[0], 3.0);
    EXPECT_DOUBLE_EQ(result[1], 7.0);
}

// Тест на пустую матрицу
TEST(CSRTest, EmptyMatrix) {
    std::map<std::tuple<int, int>, double> empty_dok;
    CSR matrix(empty_dok, 3, 3);
    
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            EXPECT_DOUBLE_EQ(matrix(i, j), 0.0);
        }
    }
}
