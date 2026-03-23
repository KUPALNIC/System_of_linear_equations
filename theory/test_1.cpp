#include <vector>
#include "../src//matrix/matrix.hpp"
#include "../src/Solvers/JacobiSolver.hpp"
#include <cmath>
#include <vector>
#include <iostream>

double calculate_residual(const Matrix& A, const std::vector<double>& x, const std::vector<double>& b) {
    std::vector<double> Ax = A * x;
    double norm_sq = 0.0;
    for (size_t i = 0; i < b.size(); ++i) {
        double r = Ax[i] - b[i];
        norm_sq += r * r;
    }
    return std::sqrt(norm_sq);
}


int main(){
    std::vector<double> v{77.0, 0.0, 0.0, 0.0, 0.0, 58.0, 0.0, 0.0, 0.0, 0.0, 54.0, 0.0, 0.0, 0.0, 0.0, 84.0};
    std::vector<double> b{7.0, 8.0, 5.0, 9.0};
    Matrix A(v, 4, 4);
    std::vector<double> x{0.0, 0.0, 0.0, 0.0};

    int total_iter = 1;
    double eps = calculate_residual(A, x, b);
    double required_eps = eps/1000;
    while (calculate_residual(A, x, b) - required_eps > 10e-15){
        std::vector<double> x_new = JacobiSolver(A, b, eps, total_iter);
        total_iter++;
        x = x_new;
        std::cout << "количество итераций ... " << total_iter << std::endl;
    }
    std::cout << "Итоговое количество итераций: " << total_iter << std::endl;  
}