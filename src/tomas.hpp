#ifndef TOMAS_HPP
#define TOMAS_HPP

#include <vector>

class three_matrix {
private:
    std::vector<double> c;
    std::vector<double> a;
    std::vector<double> b;
    
    void check_diagonal_dominance();

public:
    three_matrix(std::vector<double> a_, std::vector<double> b_, std::vector<double> c_);
    std::vector<double> solve(std::vector<double> d);
    
    bool is_diagonally_dominant() const;
};

#endif