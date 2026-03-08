#include "qr.hpp"
#include <cmath>

double QR::norm(const std::vector<double>& x) const {
    double s = 0.0;
    for (double a : x)
        s += a * a;
    return std::sqrt(s);
}

QR::QR(const Matrix& A) : R(A) {
    n = R.getRows();
    v.resize(n);

    for (int k = 0; k < n - 1; ++k) {
        int m = n - k;
        std::vector<double> x(m);

        for (int i = 0; i < m; ++i)
            x[i] = R(k + i, k);

        double nx = norm(x);
        if (nx == 0.0)
            continue;
        
        v[k] = x;

        if (x[0] >= 0)
            v[k][0] += nx;
        else
            v[k][0] -= nx;

        double nv = norm(v[k]);

        for (double& a : v[k])
            a /= nv;

        for (int j = k; j < n; ++j) {

            double dot = 0.0;

            for (int i = 0; i < m; ++i)
                dot += v[k][i] * R(k + i, j);

            for (int i = 0; i < m; ++i)
                R(k + i, j) -= 2 * v[k][i] * dot;
        }
    }
}

const Matrix& QR::getR() const {
    return R;
}

void QR::apply(std::vector<double>& b) const {

    for (int k = 0; k < n - 1; ++k) {

        int m = n - k;

        if (v[k].empty())
            continue;

        double dot = 0.0;

        for (int i = 0; i < m; ++i)
            dot += v[k][i] * b[k + i];

        for (int i = 0; i < m; ++i)
            b[k + i] -= 2 * v[k][i] * dot;
    }
}