#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include <map>
#include <tuple>
#include <functional>
#include <string>  

#include "../../src/matrix/matrix.hpp"
#include "../../src/CSR/csr.hpp"
#include "../../src/Solvers/JacobiSolver.hpp"
#include "../../src/Solvers/GaussSeidelSolver.hpp"
#include "../../src/Solvers/SymmetricGaussSeidel.hpp"
#include "../../src/Solvers/ChebyshevAcceleration.hpp"

std::string get_source_directory(const char* filepath) {
    std::string path(filepath);
    size_t last_sep = path.find_last_of("/\\");
    return (last_sep == std::string::npos) ? "." : path.substr(0, last_sep);
}

double calculate_residual(const CSR& A, const std::vector<double>& x, const std::vector<double>& b) {
    std::vector<double> Ax = A * x;
    double norm_sq = 0.0;
    for (size_t i = 0; i < b.size(); ++i) {
        double r = Ax[i] - b[i];
        norm_sq += r * r;
    }
    return std::sqrt(norm_sq);
}

std::vector<double> JacobiStep(const CSR& A, const std::vector<double>& b, const std::vector<double>& x) {
    int n = A.getRows();
    std::vector<double> x_new(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int idx = A.Rows()[i]; idx < A.Rows()[i + 1]; ++idx) {
            int j = A.Cols()[idx];
            if (i != j) sum += A.Values()[idx] * x[j];
        }
        // Диагональный элемент: ищем в строке i элемент с колонкой i
        double diag = 1.0;
        for (int idx = A.Rows()[i]; idx < A.Rows()[i + 1]; ++idx) {
            if (A.Cols()[idx] == i) {
                diag = A.Values()[idx];
                break;
            }
        }
        x_new[i] = (b[i] - sum) / diag;
    }
    return x_new;
}

int main() {
    const int N = 500;
    const double eps = 1e-10;
    const int max_iters = 250;
    const double rho_est = 0.95;

    std::map<std::tuple<int, int>, double> dok;
    std::vector<double> b(N, 1.0);
    for (int i = 0; i < N; ++i) {
        dok[{i, i}] = 4.0;
        if (i > 0) dok[{i, i - 1}] = -1.0;
        if (i < N - 1) dok[{i, i + 1}] = -1.0;
    }
    CSR A(dok, N, N);

    // Формируем путь к выходному файлу рядом с исходником
    std::string out_path = get_source_directory(__FILE__) + "/convergence.csv";
    std::ofstream out(out_path);
    
    if (!out.is_open()) {
        std::cerr << "Error: cannot open " << out_path << " for writing\n";
        return 1;
    }
    
    out << "Iteration,Method,Time_ms,Residual\n";

    std::vector<double> x0(N, 0.0);

    for (int k = 1; k <= max_iters; ++k) {
        
        // 1. Jacobi
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> x_jac = x0;
        for (int iter = 0; iter < k; ++iter) {
            x_jac = JacobiStep(A, b, x_jac);
        }
        auto end = std::chrono::high_resolution_clock::now();
        double t_jac = std::chrono::duration<double, std::milli>(end - start).count();
        out << k << ",Jacobi," << t_jac << "," << calculate_residual(A, x_jac, b) << "\n";

        // 2. Gauss-Seidel
        start = std::chrono::high_resolution_clock::now();
        std::vector<double> x_gs = GaussSeidelSolver(A, b, eps, k);
        end = std::chrono::high_resolution_clock::now();
        double t_gs = std::chrono::duration<double, std::milli>(end - start).count();
        out << k << ",GaussSeidel," << t_gs << "," << calculate_residual(A, x_gs, b) << "\n";

        // 3. Symmetric GS
        start = std::chrono::high_resolution_clock::now();
        std::vector<double> x_sgs = x0;
        for (int iter = 0; iter < k; ++iter) {
            x_sgs = SymmetricGaussSeidelStep(A, b, x_sgs, 1.0);
        }
        end = std::chrono::high_resolution_clock::now();
        double t_sgs = std::chrono::duration<double, std::milli>(end - start).count();
        out << k << ",SymmetricGS," << t_sgs << "," << calculate_residual(A, x_sgs, b) << "\n";

        // 4. Symmetric GS + Chebyshev
        start = std::chrono::high_resolution_clock::now();
        auto sgs_step = [&](const std::vector<double>& x) {
            return SymmetricGaussSeidelStep(A, b, x, 1.0);
        };
        std::vector<double> x_cheb = ChebyshevAccelerateAuto(sgs_step, x0, eps, k);
        end = std::chrono::high_resolution_clock::now();
        double t_cheb = std::chrono::duration<double, std::milli>(end - start).count();
        out << k << ",SymmetricGS_Chebyshev," << t_cheb << "," << calculate_residual(A, x_cheb, b) << "\n";

        if (k % 50 == 0) {
            std::cout << "Processed " << k << " iterations..." << std::endl;
        }
    }

    out.close();
    std::cout << "Results saved to " << out_path << std::endl;
    return 0;
}