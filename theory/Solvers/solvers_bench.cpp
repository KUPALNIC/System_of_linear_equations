#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include <map>

#include "../src/matrix/matrix.hpp"
#include "../src/CSR/csr.hpp"
#include "../src/Solvers/JacobiSolver.hpp"
#include "../src/Solvers/SimpleIterationSolver.hpp"
#include "../src/Solvers/GaussSeidelSolver.hpp"
#include "../src/Solvers/ChebishevIterationSolver.cpp"

double calculate_residual(const CSR& A, const std::vector<double>& x, const std::vector<double>& b) {
    std::vector<double> Ax = A * x;
    double norm_sq = 0.0;
    for (size_t i = 0; i < b.size(); ++i) {
        double r = Ax[i] - b[i];
        norm_sq += r * r;
    }
    return std::sqrt(norm_sq);
}

int main() {
    int N = 400; 
    std::map<std::tuple<int, int>, double> dok;
    std::vector<double> b(N, 1.0);
    double eps = 1e-15; 
    int total_iters = 300; 

    // Генерация матрицы с диагональным преобладанием
    for (int i = 0; i < N; ++i) {
        dok[{i, i}] = 2.5;
        if (i > 0) dok[{i, i - 1}] = -1.0;
        if (i < N - 1) dok[{i, i + 1}] = -1.0;
    }
    CSR A(dok, N, N);

    double l_max = findMaxEigenvalue(A); 
    double l_min = 2.5 - 2.0 * std::cos(M_PI / (N + 1)); 

    std::ofstream out("convergence.csv");
    out << "Iteration,Method,Time_ms,Residual\n";

    for (int k = 1; k <= total_iters; ++k) {
        
        // 1. Gauss-Seidel
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> x_gs = GaussSeidelSolver(A, b, eps, k); 
        auto end = std::chrono::high_resolution_clock::now();
        double t_gs = std::chrono::duration<double, std::milli>(end - start).count();
        out << k << ",GaussSeidel," << t_gs << "," << calculate_residual(A, x_gs, b) << "\n";

        // 2. Jacobi
        start = std::chrono::high_resolution_clock::now();
        std::vector<double> x_jac = JacobiSolver(A, b, eps, k);
        end = std::chrono::high_resolution_clock::now();
        double t_jac = std::chrono::duration<double, std::milli>(end - start).count();
        out << k << ",Jacobi," << t_jac << "," << calculate_residual(A, x_jac, b) << "\n";

        // 3. Simple Iteration
        start = std::chrono::high_resolution_clock::now();
        std::vector<double> x_si = SimpleIterationSolver(A, b, eps, k);
        end = std::chrono::high_resolution_clock::now();
        double t_si = std::chrono::duration<double, std::milli>(end - start).count();
        out << k << ",SimpleIteration," << t_si << "," << calculate_residual(A, x_si, b) << "\n";

        // 4. Chebyshev
        start = std::chrono::high_resolution_clock::now();
        std::vector<double> x_ch = ChebyshevIterationSolver(A, b, l_min, l_max, eps, k, 8);
        end = std::chrono::high_resolution_clock::now();
        double t_ch = std::chrono::duration<double, std::milli>(end - start).count();
        out << k << ",Chebyshev," << t_ch << "," << calculate_residual(A, x_ch, b) << "\n";

        if (k % 50 == 0) std::cout << "Processed " << k << " iterations..." << std::endl;
    }

    out.close();
    std::cout << "Results saved to convergence.csv" << std::endl;

    return 0;
}
