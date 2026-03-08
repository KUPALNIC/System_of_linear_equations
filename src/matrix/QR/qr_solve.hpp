#pragma once
#include "qr.hpp"

class QRSolve {
private:
    QR qr;
    std::vector<double> back(std::vector<double>& b) const;

public:
    QRSolve(const Matrix& A);
    std::vector<double> solve(std::vector<double> b);
};