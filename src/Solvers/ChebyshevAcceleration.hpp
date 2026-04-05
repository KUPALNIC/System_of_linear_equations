#pragma once

#include <vector>
#include <cmath>
#include <functional>

template<typename StepFunc>
std::vector<double> ChebyshevAccelerate(
    StepFunc step,
    const std::vector<double>& x0,
    double rho,
    double eps,
    int max_iter
) {
    const int n = x0.size();
    const double eps_sq = eps * eps;
    
    std::vector<double> y_prev = x0;
    std::vector<double> y_curr = step(y_prev);
    
    double w = 2.0 / (2.0 - rho * rho);
    
    for (int k = 1; k < max_iter; ++k) {
        std::vector<double> y_next(n);
        std::vector<double> stepped = step(y_curr);
        
        for (int i = 0; i < n; ++i) {
            y_next[i] = w * (stepped[i] - y_prev[i]) + y_prev[i];
        }
        
        double norm_sq = 0.0;
        for (int i = 0; i < n; ++i) {
            double r = stepped[i] - y_curr[i];
            norm_sq += r * r;
        }
        
        if (norm_sq < eps_sq) return y_next;
        
        y_prev = std::move(y_curr);
        y_curr = std::move(y_next);
        
        w = 1.0 / (1.0 - rho * rho * w / 4.0);
    }
    
    return y_curr;
}

template<typename StepFunc>
std::vector<double> ChebyshevAccelerateAuto(
    StepFunc step,
    const std::vector<double>& x0,
    double eps,
    int max_iter,
    int power_iter = 30
) {
    const int n = x0.size();
    std::vector<double> v(n, 1.0);
    
    for (int i = 0; i < power_iter; ++i) {
        std::vector<double> Av = step(v);
        double norm = 0.0;
        for (double x : Av) norm += x * x;
        norm = std::sqrt(norm);
        if (norm < 1e-15) break;
        for (int j = 0; j < n; ++j) v[j] = Av[j] / norm;
    }
    
    std::vector<double> Av = step(v);
    double rho_est = 0.0;
    for (int i = 0; i < n; ++i) rho_est += v[i] * Av[i];
    rho_est = std::abs(rho_est);
    if (rho_est >= 1.0) rho_est = 0.99;
    
    return ChebyshevAccelerate(step, x0, rho_est, eps, max_iter);
}