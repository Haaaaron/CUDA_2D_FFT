#include "a.h"

__global__ void nonlinear(cufftDoubleReal *psiiii, cufftDoubleReal *nnt, double ttt, double vvv, double dx, double dy, int nx, int ny,int totsize)
{
    int ibl=blockIdx.y+blockIdx.x*gridDim.y;
    int ind=threadIdx.x+blockDim.x*ibl;
    if(ind<totsize)
    {
        cufftDoubleReal pppi=psiiii[ind];
        nnt[ind]=ttt*(pppi*pppi)+vvv*(pppi*pppi)*pppi;
    }
}

__global__ void k_update_psi_unc(cufftDoubleComplex *psik, cufftDoubleComplex *nntk, int nx, int ny,int totsize_invspa,int totsize, double epsilon_not, double knot, double dt, double lx,  double ly, double pm)
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
        cufftDoubleComplex pppk=psik[ind];
        cufftDoubleComplex nnnk=nntk[ind];
        double Ky=iy*dky;
        double Kx=dkx*(ix<=nx/2?ix:ix-nx);
        double klapl=-(Kx*Kx+Ky*Ky);
        double Lfact=(knot*knot)-(Kx*Kx+Ky*Ky);
        double factor=1.0/(1-klapl*dt*(epsilon_not+Lfact*Lfact));
        if(ind==0)
        {
                pppk.x=pm;
                pppk.y=0;
        }
        if(ind!=0)
        {
                pppk.x=factor*(pppk.x*invvol+klapl*dt*(nnnk.x)*invvol);
                pppk.y=factor*(pppk.y*invvol+klapl*dt*(nnnk.y)*invvol);
        }
        psik[ind]=pppk;
    }
}

extern "C" __host__ void unc_constr_update(field gassed,modelparams parameters, int istep, int nsteps)
{

        nonlinear<<<(parameters.totsize+tpb-1)/tpb,tpb>>>(gassed.psione.d, gassed.nntone.d, parameters.tttone, parameters.vvvone, parameters.dx, parameters.dy, parameters.nx, parameters.ny, parameters.totsize);

        cufftExecD2Z(gassed.D2Z,gassed.psione.d,gassed.psionek.d);

        cufftExecD2Z(gassed.D2Z,gassed.nntone.d,gassed.nntonek.d);

        k_update_psi_unc<<<(parameters.totsize_invspa+tpb-1)/tpb,tpb>>>(gassed.psionek.d, gassed.nntonek.d, parameters.nx, parameters.ny, parameters.totsize_invspa, parameters.totsize,  parameters.epsilon_notone, parameters.oneqone, parameters.dt, parameters.Lx,  parameters.Ly,  parameters.pmone);

        cufftExecZ2D(gassed.Z2D,gassed.psionek.d,gassed.psione.d);
}

extern "C" __host__ void advance_in_time(field gassed,modelparams parameters, int istep, int nsteps)
{
        for(int it=0;it<nsteps;it++)
        {
               unc_constr_update(gassed, parameters, istep*nsteps+it, nsteps);
        }
}
