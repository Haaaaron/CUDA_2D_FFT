#include "a.h"

__global__ void hexfield(double *psi, double qt, double qteq, double At, double pm, double dx, double dy, int nx, int ny,int totsize, double alpha, double ynotshift)
{
    int ibl=blockIdx.y+blockIdx.x*gridDim.y;
    int ind=threadIdx.x+blockDim.x*ibl;
    if(ind<totsize)
    {
        int ix=ind/(ny);
        int iy=ind%(ny);
        double xx=(ix*dx)*cos(alpha)-(iy*dy+ynotshift)*sin(alpha);
        double yy=(ix*dx)*sin(alpha)+(iy*dy+ynotshift)*cos(alpha);
        //double xx=(ix*dx);
        //double yy=(iy*dy);
        double ppp=At*(cos(qt*xx)*cos(qteq*yy/sqrt(3.0))-0.5*cos(2*qteq*yy/sqrt(3.0)))+pm;
        psi[ind]=ppp;

    }
}

extern "C" __host__ void initfield(field gassed,modelparams parameters)
{
        printf("\n Create Initial configurations. \n \n");
        //FILE *fphhhone, *fphhhtwo, *fppsione, *fppsitwo;
        /*fphhhone=fopen("hhhone.in","r");
        fphhhtwo=fopen("hhhtwo.in","r");
        fppsione=fopen("psione.in","r");
        fppsitwo=fopen("psitwo.in","r");*/

        printf("\nCreating Initial Configuration with OMA for psi.\n");
                double At=0.179482078213165*4;

                /*srand(time(NULL));
                double rr=1.0*(rand() & 100000)/100000.;*/

                double ynotshiftone=0.0* (parameters.atone/parameters.flength)/2.0;

                hexfield<<<(parameters.totsize+tpb-1)/tpb,tpb>>>(gassed.psione.d, parameters.qtone*parameters.flength, parameters.qtone*parameters.flength, At, parameters.pmone, parameters.dx, parameters.dy, parameters.nx, parameters.ny, parameters.totsize,  parameters.alpha, ynotshiftone);

/*
                hid_t file, space, dset; // handles
                herr_t      status;
                hsize_t dims[1]=   {parameters.totsize};
                int ndims;


                file = H5Fopen ("confpsione.in", H5F_ACC_RDONLY, H5P_DEFAULT);
                dset = H5Dopen (file, "none", H5P_DEFAULT);


                space = H5Dget_space (dset);
                ndims = H5Sget_simple_extent_dims (space, dims, NULL);

                status = H5Dread (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            gassed.psione.h);

                status = H5Dclose (dset);
                status = H5Sclose (space);
                status = H5Fclose (file);

          cudaMemcpy(gassed.psione.d, gassed.psione.h, sizeof(double)*parameters.totsize,cudaMemcpyHostToDevice); */
           //exit(0);
}
