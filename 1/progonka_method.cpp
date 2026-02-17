#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
using namespace std;

class three_matrix(){
    private:
    vector<double> c;
    vector<double> a;
    vector<double> b;
    
    public:
    three_matrix(vector<double> a_, vector<double> b_, vector<double> c_){
        if(a_.size() != b_.size() || b_.size() != c_.size()){
            throw invalid_argument("Размеры диагоналей не совпадают!");
        }
        a=a_;
        b=b_;
        c=c_;
    }

    vector<double> solve(vector<double> d){
        if(d.size() != n){
            throw invalid_argument("Размер вектора d не совпадает с размером матрицы!");
        }
        int n = d.size();
        vector<double> x(n);
        vector<double> p(n), q(n);
        p[0] = -c[0]/b[0];
        q[0] = d[0]/b[0];
        //Прямой метод
        for (int i=1; i<n; i++){
            p[i+1] = -c[i]/(b[i]+a[i]*p[i]);
            q[i+1] = (d[i] - a[i]*q[i])/(b[i]+a[i]*p[i]);
        }
        //Обратный метод
        x[n-1] = (d[n-1] - a[n-1] * q[n-1]);

        for(int i = n-2; i >= 0; i--){
            x[i] = p[i+1] * x[i+1] + q[i+1];
        }
        cout << x.size() << endl;
        return x;

    }

    void test_simple_system(){
        cout << "\n=== Тест 1: Простая система ===" << endl;
        vector<double> a = {0, -1, -1};
        vector<double> b = {2, 2, 2};
        vector<double> c = {-1, -1, 0};
        vector<double> d = {1, 0, 1};
        
        three_matrix matrix(a, b, c);
        vector<double> solution = matrix.solve(d);
        if(abs(solution[0] - 1.0) < 1e-6 && 
           abs(solution[1] - 1.0) < 1e-6 && 
           abs(solution[2] - 1.0) < 1e-6){
            cout << "Тест пройден" << endl;
        } else {
            cout << "Тест не пройден" << endl;
        }
    }
    
    void test_identity_system(){
        cout << "\n=== Тест 2: Единичная система ===" << endl;
        vector<double> a = {0, 0, 0};
        vector<double> b = {1, 1, 1};
        vector<double> c = {0, 0, 0};
        vector<double> d = {2, 3, 4};
        
        three_matrix matrix(a, b, c);
        vector<double> solution = matrix.solve(d);
        if(abs(solution[0] - 2.0) < 1e-6 && 
           abs(solution[1] - 3.0) < 1e-6 && 
           abs(solution[2] - 4.0) < 1e-6){
            cout << "Тест пройден" << endl;
        } else {
            cout << "Тест не пройден" << endl;
        }
    }
};


int main(){
    try {
        test_simple_system();
        test_identity_system();
    } catch(const exception& e){
        cerr << "Неожиданная ошибка: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}