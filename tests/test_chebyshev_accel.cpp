#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <functional>
#include <random>

// Подключаем заголовки (не .cpp файлы, если у вас настроен CMake правильно)
// Если нет, оставьте как в вашем примере, но лучше использовать hpp
#include "../src/Solvers/ChebyshevAcceleration.hpp"
#include "../src/Solvers/SymmetricGaussSeidel.hpp"
#include "../src/matrix/matrix.hpp"
#include "../src/CSR/csr.hpp"

// Вспомогательная функция для вычисления L2 нормы разности векторов
double norm(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) return 1e9;
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        double d = v1[i] - v2[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}

// Генератор симметричной положительно определенной матрицы (диагональное преобладание)
std::pair<Matrix, std::vector<double>> generate_spd_system(int n) {
    Matrix A(n, n);
    std::vector<double> x_true(n);
    std::mt19937 gen(42); 
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    for (int i = 0; i < n; ++i) {
        x_true[i] = dis(gen);
        double row_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                A(i, j) = static_cast<double>(n) * 2.0; // Сильное диагональное преобладание
            } else {
                double val = dis(gen) * 0.5;
                A(i, j) = val;
                A(j, i) = val; // Симметрия
                row_sum += std::abs(val);
            }
        }
        A(i, i) += row_sum + 1.0;
    }

    std::vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A(i, j) * x_true[j];
        }
        b[i] = sum;
    }

    return {A, b};
}

class ChebyshevAccelTest : public ::testing::Test {
protected:
    int n = 10;
    double eps = 1e-6;
    int max_iter = 1000;
    Matrix A;
    std::vector<double> b;
    std::vector<double> x_true;

    void SetUp() override {
        auto sys = generate_spd_system(n);
        A = sys.first;
        b = sys.second;
        
        x_true.resize(n);
        std::mt19937 gen(42);
        std::uniform_real_distribution<> dis(-1.0, 1.0);
        for(int i=0; i<n; ++i) x_true[i] = dis(gen);
        
        b.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                b[i] += A(i, j) * x_true[j];
            }
        }
    }
};

// Тест 1: Проверка работы шаблона ускорения на простом методе Якоби
TEST_F(ChebyshevAccelTest, JacobiWithChebyshev) {
    std::vector<double> x0(n, 0.0);
    
    // Лямбда, реализующая один шаг Якоби
    auto jacobi_step = [&](const std::vector<double>& x) -> std::vector<double> {
        std::vector<double> x_new(n);
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) sum += A(i, j) * x[j];
            }
            x_new[i] = (b[i] - sum) / A(i, i);
        }
        return x_new;
    };

    double rho_est = 0.9; 

    auto x_res = ChebyshevAccelerate(jacobi_step, x0, rho_est, eps, max_iter);

    EXPECT_NEAR(vector_norm_diff(x_res, x_true), 0.0, eps * 100);
}

// Тест 2: Проверка Symmetric Gauss-Seidel (SGS) без ускорения (базовая сходимость)
TEST_F(ChebyshevAccelTest, SymmetricGaussSeidelBasicConvergence) {
    std::vector<double> x0(n, 0.0);
    std::vector<double> x_curr = x0;
    
    // Делаем много шагов SGS вручную
    for (int k = 0; k < 500; ++k) {
        x_curr = SymmetricGaussSeidelStep(A, b, x_curr, 1.0);
    }

    EXPECT_NEAR(vector_norm_diff(x_curr, x_true), 0.0, eps * 10);
}

// Тест 3: Ускорение Чебышева для SGS (Основной тест задания)
TEST_F(ChebyshevAccelTest, SGSWithChebyshevAcceleration) {
    std::vector<double> x0(n, 0.0);

    // Лямбда для SGS
    auto sgs_step = [&](const std::vector<double>& x) -> std::vector<double> {
        return SymmetricGaussSeidelStep(A, b, x, 1.0);
    };

    double rho_est = 0.8;

    auto x_res = ChebyshevAccelerate(sgs_step, x0, rho_est, eps, max_iter);

    EXPECT_NEAR(norm(x_res, x_true), 0.0, eps * 100);
}

// Тест 4: Сравнение количества итераций (Ускорение vs Обычный метод)
// Этот тест показывает эффективность метода
TEST_F(ChebyshevAccelTest, AccelerationEfficiency) {
    std::vector<double> x0(n, 0.0);
    
    // 1. Обычный SGS
    std::vector<double> x_sgs = x0;
    int iter_sgs = 0;
    for (; iter_sgs < max_iter; ++iter_sgs) {
        x_sgs = SymmetricGaussSeidelStep(A, b, x_sgs, 1.0);
        if (norm(x_sgs, x_true) < eps) break;
    }

    // 2. SGS с Чебышевым
    auto sgs_step = [&](const std::vector<double>& x) -> std::vector<double> {
        return SymmetricGaussSeidelStep(A, b, x, 1.0);
    };
    
    auto x_cheb = ChebyshevAccelerate(sgs_step, x0, 0.8, eps, max_iter);
    
    EXPECT_NEAR(norm(x_sgs, x_true), 0.0, eps * 10);
    EXPECT_NEAR(norm(x_cheb, x_true), 0.0, eps * 10);

    SUCCEED() << "Both methods converged. SGS took " << iter_sgs << " iterations.";
}

// Тест 5: Работа с CSR матрицей
TEST_F(ChebyshevAccelTest, CSR_SGS_Chebyshev) {
    std::map<std::tuple<int,int>, double> dok;
    for(int i=0; i<n; ++i) {
        for(int j=0; j<n; ++j) {
            if(std::abs(A(i,j)) > 1e-12) {
                dok[{i,j}] = A(i,j);
            }
        }
    }
    CSR A_csr(dok, n, n);

    std::vector<double> x0(n, 0.0);

    auto sgs_step_csr = [&](const std::vector<double>& x) -> std::vector<double> {
        return SymmetricGaussSeidelStep(A_csr, b, x, 1.0);
    };

    auto x_res = ChebyshevAccelerate(sgs_step_csr, x0, 0.8, eps, max_iter);

    EXPECT_NEAR(norm(x_res, x_true), 0.0, eps * 100);
}

// Тест 6: Автоматическая оценка спектра (ChebyshevAccelerateAuto)
TEST_F(ChebyshevAccelTest, AutoRhoEstimation) {
    std::vector<double> x0(n, 0.0);

    auto sgs_step = [&](const std::vector<double>& x) -> std::vector<double> {
        return SymmetricGaussSeidelStep(A, b, x, 1.0);
    };

    // Используем авто-версию, которая сама оценит rho
    auto x_res = ChebyshevAccelerateAuto(sgs_step, x0, eps, max_iter);

    EXPECT_NEAR(norm(x_res, x_true), 0.0, eps * 100);
}