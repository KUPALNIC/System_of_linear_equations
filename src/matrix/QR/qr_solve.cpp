#include "qr_solve.hpp"

QRSolve::QRSolve(const Matrix& A) : qr(A) {}

std::vector<double> QRSolve::back(std::vector<double>& b) const {
    const Matrix& R = qr.getR();
    int n = R.getRows();
    std::vector<double> x(n);

    for (int i = n-1; i >= 0; --i) {
        double s = b[i];
        for (int j = i+1; j < n; ++j)
            s -= R(i,j) * x[j];
        x[i] = s / R(i,i);
    }

    return x;
}

std::vector<double> QRSolve::solve(std::vector<double> b) {
    qr.apply(b);
    return back(b);
}