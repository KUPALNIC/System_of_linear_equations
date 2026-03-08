#include <fstream> 
#include <chrono>
#include <vector>
#include <map>
#include <random>
#include "../src/matrix/matrix.hpp" 
#include "../src/CSR/csr.hpp"

using namespace std::chrono;

void run_benchmarks(const std::string& filename) {
    std::ofstream file(filename);

    file << "Size,Sparsity,Dense_ms,CSR_ms\n";

    std::vector<int> sizes = {100, 500, 1000, 3000, 5000, 7500, 10000};
    std::vector<double> sparsities = {0.5, 0.9, 0.95, 0.99, 0.999};

    std::default_random_engine gen(42); 
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int n : sizes) {        
        for (double s : sparsities) {
            std::map<std::tuple<int, int>, double> dok;
            std::vector<double> matrix_data(n*n, 0.0);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (dist(gen) > s) {
                        double v = 1.0; 
                        dok[{i, j}] = v;
                        matrix_data[i*n + j] = v;
                    }
                }
            }

            Matrix matrix_mat(matrix_data, n, n);
            CSR csr_mat(dok, n, n);
            std::vector<double> x(n, 1.0);

            // Замер matrix
            auto s1 = high_resolution_clock::now();
            auto res1 = matrix_mat * x; 
            auto e1 = high_resolution_clock::now();
            double matrix_time = duration<double, std::milli>(e1 - s1).count();

            // Замер CSR
            auto s2 = high_resolution_clock::now();
            auto res2 = csr_mat * x;
            auto e2 = high_resolution_clock::now();
            double csr_time = duration<double, std::milli>(e2 - s2).count();

            // Запись в файл
            file << n << "," << s << "," << matrix_time << "," << csr_time << "\n";
        }
    }

    file.close();
}

int main() {
    run_benchmarks("benchmark_results.csv");
    return 0;
}
