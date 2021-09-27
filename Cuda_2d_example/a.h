#ifndef A_H
#define A_H

//#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <time.h>
//#include "hdf5.h"

typedef struct {
   cufftDoubleReal *h;
   cufftDoubleReal *d;
} cufftrealvector;

typedef struct {
   cufftDoubleComplex *h;
   cufftDoubleComplex *d;
} cufftcomplexvector;

typedef struct {
   double *h;
   double *d;
} doublevector;


typedef struct {
   double2 *h;
   double2 *d;
} double2vector;

typedef struct {
   clock_t *ttttime;
   cufftrealvector enemat;
   cufftrealvector redone;
   cufftrealvector redtwo;
   cufftrealvector redthr;
   cufftrealvector totene;
   cufftrealvector psione;
   cufftrealvector dxxpsione;
   cufftrealvector dyypsione;
   cufftrealvector nntone;
   cufftcomplexvector psionek;
   cufftcomplexvector dxxpsionek;
   cufftcomplexvector dyypsionek;
   cufftcomplexvector nntonek;
   cufftHandle D2Z;
   cufftHandle Z2D;

} field;

typedef struct {
   double dt;
   double dx;
   double dy;
   double dz;
   double dkx;
   double dky;
   double dkz;
   double Lx;
   double Ly;
   /*double kappaone;
   double kappatwo;
   double Cmax;*/
   double kmax;
   double alpha;
   double pmone;
   double qtone;
   double atone;
   int sssone;
   int ssstwo;
   int sssthr;
   int nx;
   int ny;
   int totsize;
   int totsize_invspa;
   double epsilon_notone;
   double vvvone;
   double tttone;
   double oneqone;
   double flength;
   double ELiqOneMode(){return (epsilon_notone+pow(oneqone,4))*0.5*(pmone*pmone)+(1.0/3.0)*tttone*(pmone*pmone)*pmone+0.25*vvvone*(pmone*pmone)*(pmone*pmone);};
//   double ELiqTwoModes(){return (-epsilon+(1+Rnot)*((Qone*Qone)*(Qone*Qone)+Rone))*0.5*(rho*rho)+0.25*(rho*rho)*(rho*rho);};
//   double ELiqThreeModes(){return (-epsilon+(1+Rnot)*((Qone*Qone)*(Qone*Qone)+Rone)*((Qtwo*Qtwo)*(Qone*Qtwo)+Rtwo))*0.5*(rho*rho)+0.25*(rho*rho)*(rho*rho);};
} modelparams;

#define tpb 256


#endif
