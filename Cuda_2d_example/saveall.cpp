#include "a.h"
extern "C" __host__ void energyall(field gassed,modelparams parameters);

__host__ void savedata(field gassed, modelparams parameters, int n)
{

   FILE *fp, *Fpt;
   char filename[18];
   int ccc;
   energyall(gassed, parameters);

   /*int ret,ccc;
   hid_t file, space, dset; // handles
   herr_t      status;
   hsize_t dims[1]=   {parameters.totsize};

    // Create a new file using the default properties
   sprintf(filename,"psione%08d.h5", (n%2));
   file= H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   // Create dataspace. Setting maximum size to NULL sets the maximum size to the current size.
   space=H5Screate_simple(1, dims,NULL);

   cudaMemcpy(gassed.psione.h, gassed.psione.d, sizeof(cufftDoubleReal)*parameters.totsize,cudaMemcpyDeviceToHost);
   dset = H5Dcreate (file, "none", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, gassed.psione.h);

   status = H5Dclose (dset);
   status = H5Sclose (space);
   status = H5Fclose (file);


   sprintf(filename,"eneloc%08d.h5", (n%2));
   file= H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   // Create dataspace. Setting maximum size to NULL sets the maximum size to the current size.
   space=H5Screate_simple(1, dims,NULL);

   dset = H5Dcreate (file, "ene", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, gassed.enemat.h);
   status = H5Dclose (dset);
   status = H5Sclose (space);
   status = H5Fclose (file);
*/

   ccc=0;
   long double pm_one_host=0.0;
   for(int i=0;i<parameters.nx; i++)
   {
     long double localpmonei=0.0;
     for(int j=0;j<parameters.ny; j++)
     {
       localpmonei+=gassed.psione.h[ccc];
       ccc++;
     }
     pm_one_host+=localpmonei/parameters.ny;
   }
   pm_one_host=pm_one_host/parameters.nx;
   cudaDeviceSynchronize();

   fp=fopen("stats.txt","a");
   fprintf(fp,"%06d %.16e %.16e %.16e %.16e %.16e\n", n, parameters.ELiqOneMode(),  gassed.totene.h[0], gassed.totene.h[0]-parameters.ELiqOneMode(),
   (double) pm_one_host);
   fclose(fp);
   printf("%06d %.16e %.16e %.16e %.16e %.16e\n",n, parameters.ELiqOneMode(), gassed.totene.h[0], gassed.totene.h[0]-parameters.ELiqOneMode(),
   (double) pm_one_host);
   Fpt=fopen("inint.in", "w");
   fprintf(Fpt," %d\n", n%2);
   fclose(Fpt);

//   delete(fp);
//   delete(Fpt);

}
