#pragma once
#include <vector>

class Matrix{
    private:
        int  nro_;
        int ny_;
        std::vector<double> m_;
        
        public:
            Matrix(int rows, int cols, const std::vector<double>& matrix): nx_(rows), ny_(cols){}

            double operator() ( int i, int j) const {
                return m_[i*ny_ + j]
            } 
};