#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <complex>
#include <math.h>
#include <chrono>

using std::tuple;
using std::complex;
using std::cout;
using std::setw;
using std::setprecision;
using namespace std::chrono;
template <typename T> std::string type_name();

void init_arguments(int argc, char const* argv[], int* dx, int* dy, int* iteration_count);

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

        void print_R(){

            //printing limits to avoid horrible formatting
            if (j_len > 20 || i_len > 20) return;
            cout << "Real space: \n";
            cout << "    ";
            for (auto j = 0; j < i_len; j++) {
                cout << setw(3) << j << ":";
            }
            cout << '\n';
            for (auto j = 0; j < j_len; j++) {
                cout << setw(3) << j <<": ";
                for (auto i = 0; i < i_len; i++) {
                    cout << setw(3) << (*this)(i,j) << " ";
                }        
                cout << "\n"; 
            }   
            cout << "\n";
        }

        void print_C(){
            if (j_len > 20 || i_len > 20) return;
            cout << "Complex space: \n";
            double real;
            double imag;
            cout << "    ";
            for (auto j = 0; j < j_len; j++) {
                cout << setw(5) << j << ":" << setw(4) << " ";
            }
            cout << '\n';
            for (auto i=0; i < i_len/2+1; i++) {
                cout << setw(3) << i <<": ";
                for (auto j = 0; j < j_len; j++) {
                    real = complex_data[j_len*i+j].real();
                    if (real < 1.0e-14) real = 0;
                    imag = complex_data[j_len*i+j].imag();
                    if (imag < 1.0e-14) imag = 0;
                    cout << "(" << setw(3) << setprecision(3) << real << "," << setw(3) << setprecision(3) << imag << ") ";
                }
                cout << "\n";
            }
            cout << "\n";
        }
};