#include <tuple>

void launch_kernel();

template <typename T>
class system_2D{

    private:
        T* data;
        int j_len,i_len;

    public:
        system_2D(int w, int h):
            data(new T[w*h]), j_len(w), i_len(h) {}

        ~system_2D(){
            delete[] data;
        }

        T& operator()(int i, int j){
            return data[i_len*j+i];
        }
};