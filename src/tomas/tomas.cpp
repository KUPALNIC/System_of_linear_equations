#include "tomas.hpp"
#include <stdexcept>
#include <cmath>

three_matrix::three_matrix(std::vector<double> a_, std::vector<double> b_, std::vector<double> c_) {
    if(a_.size() != b_.size() || b_.size() != c_.size()) {
        throw std::invalid_argument("Размеры диагоналей не совпадают!");
    }
    a = a_;
    b = b_;
    c = c_;
    
    check_diagonal_dominance();
}

void three_matrix::check_diagonal_dominance() {
    int n = b.size();
    
    for (int i = 0; i < n; i++) {
        double diagonal = std::abs(b[i]);
        double sum_off_diagonal = 0.0;
        if (i > 0) {
            sum_off_diagonal += std::abs(a[i]);
        }
        if (i < n - 1) {
            sum_off_diagonal += std::abs(c[i]);
        }
        if (diagonal < sum_off_diagonal) {
            throw std::invalid_argument(
                "Нарушено условие диагонального преобладания на строке "
            );
        }
    }
}

bool three_matrix::is_diagonally_dominant() const {
    int n = b.size();
    for (int i = 0; i < n; i++) {
        double diagonal = std::abs(b[i]);
        double sum_off_diagonal = 0.0;
        if (i > 0) {
            sum_off_diagonal += std::abs(a[i]);
        }
        if (i < n - 1) {
            sum_off_diagonal += std::abs(c[i]);
        }
        if (diagonal < sum_off_diagonal) {
            return false;
        }
    }

    return true;
}

std::vector<double> three_matrix::solve(std::vector<double> d) {
    int n = b.size();
    if(d.size() != n) {
        throw std::invalid_argument("Размер вектора d не совпадает с размером матрицы!");
    }
    
    std::vector<double> x(n);
    std::vector<double> p(n), q(n);
    p[0] = -c[0]/b[0];
    q[0] = d[0]/b[0];
    
    // Прямой метод
    for (int i=0; i<n-1; i++) {
        p[i+1] = -c[i]/(b[i]+a[i]*p[i]);
        q[i+1] = (d[i] - a[i]*q[i])/(b[i]+a[i]*p[i]);
    }
    
    // Обратный метод
    x[n-1] = (d[n-1] - a[n-1] * q[n-1])/(a[n-1]*p[n-1] + b[n-1]);

    for(int i = n-2; i >= 0; i--) {
        x[i] = p[i+1] * x[i+1] + q[i+1];
    }
    return x;
}