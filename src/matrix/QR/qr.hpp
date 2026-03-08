#pragma once
#include "../matrix.hpp"
#include <vector>

class QR {
private:
    Matrix R;
    std::vector<std::vector<double>> v;
    int n;

    double norm(const std::vector<double>& x) const;

public:
    QR(const Matrix& A);

    const Matrix& getR() const;
    void apply(std::vector<double>& b) const;
};