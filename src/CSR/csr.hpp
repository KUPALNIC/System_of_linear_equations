#pragma once 
#include <vector>
#include <map>
#include <tuple>

class CSR {
private:
    std::vector<double> values;     
    std::vector<int> cols;          
    std::vector<int> rows;       

    int rows_count;
    int cols_count;

public:
    CSR(const std::map<std::tuple<int, int>, double>& dok, int rows, int cols);

    double operator()(int i, int j) const;

    
};