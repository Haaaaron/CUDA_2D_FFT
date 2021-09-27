#ifndef __TwoDMPIPFC_H__
#define __TwoDMPIPFC_H__
typedef struct {
    int lx;
    int ly;
    int lxhp; 
    int lyhp;
    int lhx;
    int lhy;
    double pi;
    double sq3;
    double at;
    double qqt;
    double dt; 
    double VV;
    double dx;
    double dy;
    double r;
    double pm;
} Parameters;

typedef struct {
    int myid;
    int numprocs;
    int lyr,lxr;
    int rcounts, scounts;
    int sumscounts, sumrcounts;
    int* rdispls;
    int* sdispls;
    int* imin,* imax;
    int* jmin,* jmax;
    int iminmyid, imaxmyid;
    int jminmyid, jmaxmyid;
    int j_len;
    int i_len;
} Init_rank_parameters;

//Structs for FFTW arrays and plans
typedef struct {
    fftw_plan* forward_plan;
    fftw_plan* backward_plan;
    fftw_complex* complex_arr;
    double* real_arr;
} fft_real_complex;

typedef struct {
    fftw_plan* forward_plan;
    fftw_plan* backward_plan;
    fftw_complex* f_arr;
    fftw_complex* b_arr;
} fft_complex_complex;

typedef struct {
    double start_time;
    double* total_time;
} Timer;

typedef struct {
    fft_real_complex psi;
    fft_real_complex zero;
    fft_real_complex lxhp;
    fft_complex_complex psic;
    double complex* send;
    double complex* recv;
    struct Time {
        Timer fft_plan;
        Timer before_MPI;
        Timer during_MPI;
        Timer after_MPI;
        Timer run_time;
        Timer without_init;
    } time;

} fft_builds;

// Structs fo Physical system
typedef struct {
    double* lxhp,* zero,* rest;
} Phys_coeff;

typedef struct {
    // kx_zero nad kx_final are a waste of memory if rank isn't 0
    fftw_complex* kx_zero,* kx_final,* rest;
} Complex_field;

typedef struct {
    Phys_coeff energy;
    Phys_coeff ntc;
    Phys_coeff psi;
} Phys_all_coeff;

typedef struct {
    Complex_field complex_buffer;
    Complex_field complex_nt;
    Complex_field complex_psi;
    double* psi;
    double* nt;
    double* real_buffer;
    double* final_energy;
    double total_energy;
    Phys_all_coeff coeff;
} Physical_system;



extern void coeff_init(Parameters params, Init_rank_parameters info, Physical_system* system);
extern void fill_2d_phys_arr(Parameters params, Init_rank_parameters info, Phys_all_coeff coefficients, int start, int end);
extern double energy(Parameters params, Init_rank_parameters info, Physical_system system, fft_builds builds);
extern double* update(Parameters params, Init_rank_parameters info, Physical_system system, fft_builds builds);

extern void fft_init(Parameters params, Init_rank_parameters* info, fft_builds* plans);
extern void fft_forward(Parameters params, Init_rank_parameters info, fft_builds plans, double* input_arr, Complex_field output_arr);
extern void fft_backward(Parameters params, Init_rank_parameters info, fft_builds builds, Complex_field input_arr, double* output_arr);

extern void init_params(int argc, char* argv[], int* size_x, int* size_y, int* iteration_count);
extern void output_time(struct Time time, Init_rank_parameters info);

//Macros for debugging
#define deb(arr, end, type)                 \
  for (int l = 0; l < (end); l++) {         \
    type(arr[l],l);                         \
  }                                         \
printf("\n");//exit(0);
#define rdeb(x,i) printf("%.16e \n",x,i)
#define cdeb(x,i) printf("%.16e %.16e %d\n",creal(x),cimag(x),i)
#define fcdeb(x,i) printf("%.16e %.16e %d\n",x[0],x[1],i)

#endif