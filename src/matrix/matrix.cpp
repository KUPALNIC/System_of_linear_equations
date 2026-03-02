#include <stdexcept>
#include <vector>
#include "matrix.hpp"

class Matrix{
    private:
        int  nx_;
        int ny_;
        std::vector<double> m_;
        
        public:
            Matrix(int rows, int cols, const std::vector<double>& matrix): nx_(rows), ny_(cols){
                if (matrix.size() != rows*cols){
                    throw std::invalid_argument("Некорректный размер массива")
                }

                m_ = matrix
            }

            double operator() ( int i, int j) const {
                return m_[i*ny_ + j];
            } 
};