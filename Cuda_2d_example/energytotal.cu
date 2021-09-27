#include "a.h"

__global__ void
//__launch_bounds__( tpb, 2 )

Lop(cufftDoubleComplex *psik, cufftDoubleComplex *dxxpsik, cufftDoubleComplex *dyypsik, int nx, int ny,int totsize_invspa,int totsize, double epsilon_not, double knot, double lx,  double ly)
{
    int ibl=blockIdx.y+blockIdx.x*gridDim.y;
    int ind=threadIdx.x+blockDim.x*ibl;
    double pi=acos(-1.0);
    double invvol=1.0/(double)totsize;
    double dkx=(2.0*pi/lx);
    double dky=(2.0*pi/ly);
    if(ind<totsize_invspa)
    {
        int ix=ind/(ny/2+1);
        int iy=ind%(ny/2+1);
        cufftDoubleComplex pppk=psik[ind],dxxpppk,dyypppk;
        pppk.x=pppk.x*invvol;
        pppk.y=pppk.y*invvol;
        double Ky=iy*dky;
        double Kx=dkx*(ix<=nx/2?ix:ix-nx);
        dxxpppk.x=-(Kx*Kx)*pppk.x;
        dxxpppk.y=-(Kx*Kx)*pppk.y;
        dyypppk.x=-(Ky*Ky)*pppk.x;
        dyypppk.y=-(Ky*Ky)*pppk.y;


        dxxpsik[ind]=dxxpppk;
        dyypsik[ind]=dyypppk;
    }
}


__global__ void energy_mat_double(cufftDoubleReal *psione, cufftDoubleReal *dxxpsione, cufftDoubleReal *dyypsione, cufftDoubleReal *ene, double epsilon_one, double ktone, double tttone, double vvvone, double dx, double dy, int nx, int ny,int totsize)
{
    int ibl=blockIdx.y+blockIdx.x*gridDim.y;
    int ind=threadIdx.x+blockDim.x*ibl;
    if(ind<totsize)
    {
        cufftDoubleReal pppone=psione[ind];
        cufftDoubleReal dxxpppone=dxxpsione[ind];
        cufftDoubleReal dyypppone=dyypsione[ind];
        ene[ind]=0.5*epsilon_one*pppone*pppone+0.5*(ktone*ktone*pppone+dxxpppone+dyypppone)*(ktone*ktone*pppone+dxxpppone+dyypppone)+(1.0/3.0)*tttone*(pppone*pppone)*pppone+0.25*vvvone*(pppone*pppone)*(pppone*pppone);
    }
}

__global__ void reduction_one(cufftDoubleReal *ene, cufftDoubleReal *redcone, int totsize)
{
    int ibl=blockIdx.y+blockIdx.x*gridDim.y;
    int ind=threadIdx.x+blockDim.x*ibl;
    __shared__ cufftDoubleReal shene[512];
    shene[threadIdx.x]=0;
    if(ind<totsize)
    {
        shene[threadIdx.x]=ene[ind];
    }
    __syncthreads();
    for(int s=256;s>0;s>>=1)
    {
        if(threadIdx.x<s)
        {shene[threadIdx.x]+=shene[threadIdx.x+s];}
        __syncthreads();
    }
    if(threadIdx.x==0)
    {
        redcone[ibl]=shene[0];
    }
}

__global__ void reduction_two(cufftDoubleReal *redcone, cufftDoubleReal *redctwo, int tsize)
{
    int ibl=blockIdx.y+blockIdx.x*gridDim.y;
    int ind=threadIdx.x+blockDim.x*ibl;
    __shared__ cufftDoubleReal shene[512];
    shene[threadIdx.x]=0;
    if(ind<tsize)
    {
        shene[threadIdx.x]=redcone[ind];
    }
    __syncthreads();
    for(int s=256;s>0;s>>=1)
    {
        if(threadIdx.x<s)
        {shene[threadIdx.x]+=shene[threadIdx.x+s];}
        __syncthreads();
    }
    if(threadIdx.x==0)
    {
        redctwo[ibl]=shene[0];
    }
}
__global__ void reduction_thr(cufftDoubleReal *redctwo, cufftDoubleReal *redcthr, int tsize)
{
    int ibl=blockIdx.y+blockIdx.x*gridDim.y;
    int ind=threadIdx.x+blockDim.x*ibl;
    __shared__ cufftDoubleReal shene[512];
    shene[threadIdx.x]=0;
    if(ind<tsize)
    {
        shene[threadIdx.x]=redctwo[ind];
    }
    __syncthreads();
    for(int s=256;s>0;s>>=1)
    {
        if(threadIdx.x<s)
        {shene[threadIdx.x]+=shene[threadIdx.x+s];}
        __syncthreads();
    }
    if(threadIdx.x==0)
    {
        redcthr[ibl]=shene[0];
    }
}


__global__ void final_ene(cufftDoubleReal *redcthr, cufftDoubleReal *totene, int totsize)
{

    __shared__ cufftDoubleReal shene[512];

    shene[threadIdx.x]=redcthr[threadIdx.x];

    __syncthreads();
    for(int s=256;s>0;s>>=1)
    {
        if(threadIdx.x<s)
        {shene[threadIdx.x]+=shene[threadIdx.x+s];}
        __syncthreads();
    }
    if(threadIdx.x==0)
    {
        totene[0]=shene[0]/totsize;
        //printf("%lf \n", shene[0]/totsize);
    }
}

extern "C" __host__ void energyall(field gassed,modelparams parameters)
{
        //printf("\n Start of energy computing. \n \n");

        cufftExecD2Z(gassed.D2Z,gassed.psione.d,gassed.psionek.d);

        Lop<<<(parameters.totsize_invspa+tpb-1)/tpb,tpb>>>(gassed.psionek.d, gassed.dxxpsionek.d, gassed.dyypsionek.d, parameters.nx, parameters.ny, parameters.totsize_invspa, parameters.totsize,  parameters.epsilon_notone, parameters.oneqone, parameters.Lx,  parameters.Ly);


        cufftExecZ2D(gassed.Z2D,gassed.dxxpsionek.d,gassed.dxxpsione.d);
        cufftExecZ2D(gassed.Z2D,gassed.dyypsionek.d,gassed.dyypsione.d);

        cudaMemset(gassed.enemat.d, 0,  sizeof(cufftDoubleReal)*parameters.totsize);

        energy_mat_double<<<(parameters.totsize+tpb-1)/tpb,tpb>>>(gassed.psione.d, gassed.dxxpsione.d, gassed.dyypsione.d, gassed.enemat.d, parameters.epsilon_notone, parameters.oneqone, parameters.tttone, parameters.vvvone, parameters.dx, parameters.dy, parameters.nx, parameters.ny, parameters.totsize);

        cudaMemset(gassed.redone.d, 0, sizeof(cufftDoubleReal)*parameters.sssone);
        cudaMemset(gassed.redtwo.d, 0, sizeof(cufftDoubleReal)*parameters.ssstwo);
        cudaMemset(gassed.redthr.d, 0, sizeof(cufftDoubleReal)*512);
        cudaMemset(gassed.totene.d, 0, sizeof(cufftDoubleReal));
        reduction_one<<<parameters.sssone,512>>>(gassed.enemat.d,gassed.redone.d, parameters.totsize);
//        printf("%d \n", (parameters.totsize+512-1)/512);
        reduction_two<<<parameters.ssstwo,512>>>(gassed.redone.d,gassed.redtwo.d, parameters.sssone);
        reduction_thr<<<512,512>>>(gassed.redtwo.d,gassed.redthr.d, parameters.ssstwo);
        final_ene<<<1,512>>>(gassed.redthr.d,gassed.totene.d, parameters.totsize);
        cudaMemcpy(gassed.totene.h, gassed.totene.d, sizeof(double),cudaMemcpyDeviceToHost);
        cudaMemcpy(gassed.enemat.h, gassed.enemat.d, sizeof(cufftDoubleReal)*parameters.totsize,cudaMemcpyDeviceToHost);
        /*int ccc=0; long double ene_host=0.0;
        for(int i=0;i<parameters.nx; i++)
        {
          long double localenei=0.0;
          for(int j=0;j<parameters.ny; j++)
          {
            localenei+=gassed.enemat.h[ccc];
            ccc++;
          }
          ene_host+=localenei/parameters.ny;
        }
        gassed.totene.h[0]=ene_host/parameters.nx;*/
        cudaDeviceSynchronize();
        //printf("\n End of energy computing. \n \n");

}
