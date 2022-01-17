#include <stdio.h>
#include <iostream>
#include <tuple>
#include <complex>
#include <math.h>
#include <chrono>

using std::tuple;
using std::complex;
using std::cout;
using namespace std::chrono;
template <typename T> std::string type_name();

template <typename R, typename C>
class System_2D{

    private:
        R* real_data;
        C* complex_data;
        int j_len,i_len;

    public:
    
        System_2D(int dx, int dy):
            real_data(new R[dx*dy]), complex_data(new C[(dx/2+1)*dy]), i_len(dx), j_len(dy) {}

        ~System_2D(){
            delete[] real_data;
            delete[] complex_data;
        }

        R& operator()(int  i, int j){
            return real_data[i_len*j+i];
        }

        R* real() {
            return real_data;
        }

        C* complex() {
            return complex_data;
        }

        tuple<int,int> get_dimensions() {
            return tuple<int,int>{i_len,j_len};
        }

        void print(double rescale){
            for (auto j = 0; j < j_len; j++) {
                for (auto i = 0; i < i_len; i++) {
                    cout << (*this)(i,j)/rescale << " ";
                }        
                cout << "\n"; 
            }   
            cout << "\n";
        }
};