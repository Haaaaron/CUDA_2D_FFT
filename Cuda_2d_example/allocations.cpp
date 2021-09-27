#include "a.h"

__host__ void parameterinit(modelparams *parameters)
{
   parameters[0].oneqone=1.0;
   parameters[0].qtone=parameters[0].oneqone*sqrt(3.0)/2.0;
   parameters[0].atone=2*acos(-1.0)/parameters[0].qtone;


  parameters[0].nx=424; //3...
   parameters[0].ny=742;
   int kang=26;
   int lang=0;

  /* parameters[0].nx=2011; //3...
   parameters[0].ny=3484;
   int kang=124;
   int lang=0;*/
/*     parameters[0].nx=840; //3...
   parameters[0].ny=1440;
   int kang=75;
   int lang=1;
*/
/*   parameters[0].nx=1680; //3...
   parameters[0].ny=2880;
   int kang=144;
   int lang=1;*/
/*   parameters[0].nx=2304;
   parameters[0].ny=4096;
   int kang=204;
   int lang=1;
*/
   /*parameters[0].nx=216; //7.
   parameters[0].ny=384;
   int kang=20;
   int lang=1;*/
   /*parameters[0].nx=540; //3...
   parameters[0].ny=960;
   int kang=45;
   int lang=1;*/
/*   parameters[0].nx=4320;
   parameters[0].ny=7680;
   int kang=358;
   int lang=1;*/
   /*parameters[0].nx=2304;
   parameters[0].ny=4096;
   int kang=204;
   int lang=1;*/
   /*parameters[0].nx=1152; //1.4 ...
   parameters[0].ny=2048;
   int kang=104;
   int lang=1;*/

   parameters[0].alpha=atan((sqrt(3.0)*(lang+0.5))/(kang+0.5)); //25*acos(-1.0)/180.0;
   parameters[0].flength=0.99767*(1.0+0.001*0); //between -50 and 50
   parameters[0].Lx=(parameters[0].atone/parameters[0].flength)*sqrt(3.0*(lang+0.5)*(lang+0.5)+(kang+0.5)*(kang+0.5)); //parameters[0].at*24.0*cos(parameters[0].alpha); //parameters[0].dx*parameters[0].nx;
   parameters[0].Ly=(parameters[0].atone/parameters[0].flength)*sqrt((1.5)*(1.5)+3.0*(kang+0.5)*(kang+0.5)); //parameters[0].at*14.0*sqrt(3.0)*cos(parameters[0].alpha);//parameters[0].dy*parameters[0].ny;
   parameters[0].dx=parameters[0].Lx/parameters[0].nx;
   parameters[0].dy=parameters[0].Ly/parameters[0].ny;
   parameters[0].dt=0.25;
   parameters[0].pmone=0;
   parameters[0].epsilon_notone=-0.15;
   parameters[0].vvvone=1.0;
   parameters[0].tttone=0.5/sqrt(0.98/3.0);
   parameters[0].totsize=parameters[0].nx*parameters[0].ny;
   //parameters[0].totsize_pad=parameters[0].totsize; //parameters[0].nx*2*(parameters[0].ny/2+1);
   parameters[0].totsize_invspa=parameters[0].nx*(parameters[0].ny/2+1);
   parameters[0].sssone=(parameters[0].totsize+512-1)/512;
   parameters[0].sssone=(parameters[0].sssone>512?parameters[0].sssone:512);
   parameters[0].ssstwo=(parameters[0].sssone+512-1)/512;
   parameters[0].ssstwo=(parameters[0].ssstwo>512?parameters[0].ssstwo:512);
   parameters[0].sssthr=(parameters[0].ssstwo+512-1)/512;
   printf("%d %d %d \n", parameters[0].sssone, parameters[0].ssstwo, parameters[0].sssthr);
   if(parameters[0].totsize>pow(512,4) || parameters[0].sssthr>= 512)
   {printf("\n Alert!!!!!!! One more layer needed for the reduction of the energy matrix\n ");}

}

__host__ void pointerinit(field *gassed,modelparams *parameters)
{
   parameterinit(parameters);

   cudaMalloc((void**)&gassed[0].enemat.d, sizeof(cufftDoubleReal)*parameters[0].totsize);

   cudaMalloc((void**)&gassed[0].redone.d, sizeof(cufftDoubleReal)*parameters[0].sssone);

   cudaMalloc((void**)&gassed[0].redtwo.d, sizeof(cufftDoubleReal)*parameters[0].ssstwo);

   cudaMalloc((void**)&gassed[0].redthr.d, sizeof(cufftDoubleReal)*512);

   gassed[0].totene.h=new cufftDoubleReal[1];
   cudaMalloc((void**)&gassed[0].totene.d, sizeof(cufftDoubleReal));

   gassed[0].enemat.h=new cufftDoubleReal[parameters[0].totsize];
   cudaMalloc((void**)&gassed[0].enemat.d, sizeof(cufftDoubleReal)*parameters[0].totsize);

   gassed[0].psione.h=new cufftDoubleReal[parameters[0].totsize];
   cudaMalloc((void**)&gassed[0].psione.d, sizeof(cufftDoubleReal)*parameters[0].totsize);
   cudaMalloc((void**)&gassed[0].dxxpsione.d, sizeof(cufftDoubleReal)*parameters[0].totsize);
   cudaMalloc((void**)&gassed[0].dyypsione.d, sizeof(cufftDoubleReal)*parameters[0].totsize);
   cudaMalloc((void**)&gassed[0].psionek.d, sizeof(cufftDoubleComplex)*parameters[0].totsize_invspa);
   cudaMalloc((void**)&gassed[0].dxxpsionek.d, sizeof(cufftDoubleComplex)*parameters[0].totsize_invspa);
   cudaMalloc((void**)&gassed[0].dyypsionek.d, sizeof(cufftDoubleComplex)*parameters[0].totsize_invspa);

   cudaMalloc((void**)&gassed[0].nntone.d, sizeof(cufftDoubleReal)*parameters[0].totsize);
   cudaMalloc((void**)&gassed[0].nntonek.d, sizeof(cufftDoubleComplex)*parameters[0].totsize_invspa);


   cufftPlan2d(&gassed[0].D2Z,parameters[0].nx,parameters[0].ny,CUFFT_D2Z);
   cufftPlan2d(&gassed[0].Z2D,parameters[0].nx,parameters[0].ny,CUFFT_Z2D);

   printf(" Allocations done \n Lx= %lf  \n Ly= %lf \n", parameters[0].Lx, parameters[0].Ly);
   printf(" totalsize= %d  \n totalsiz inverse= %d \n", parameters[0].totsize, parameters[0].totsize_invspa);
   printf(" Allocations done. angle= %lf \n", parameters[0].alpha*180/acos(-1.0));
   printf(" dx= %22.16lf;  \n dy= %22.16lf; \n", parameters[0].dx, parameters[0].dy);
   FILE *fp1;
   fp1=fopen("parameters", "w");
   fprintf(fp1,"Angle= %lf \n", parameters[0].alpha*180/acos(-1.0));
   fprintf(fp1,"totalsize= %d  \n totalsiz inverse= %d \n", parameters[0].totsize, parameters[0].totsize_invspa);
   fprintf(fp1,"nx= %d;   ny= %d; \n", parameters[0].nx, parameters[0].ny);
   fprintf(fp1,"Lx= %lf;  \n Ly= %lf; \n", parameters[0].Lx, parameters[0].Ly);
   fprintf(fp1,"dx= %22.16lf;\ndy= %22.16lf; \n", parameters[0].dx, parameters[0].dy);
   fprintf(fp1,"atone= %22.16lf ; qtone= %22.16lf;\n", parameters[0].atone, parameters[0].qtone);
   fprintf(fp1,"deltaB_1= %22.16lf;  \n", parameters[0].epsilon_notone);
   fprintf(fp1,"Bx_1= %22.16lf;   \n", 1.0);
   fprintf(fp1,"Bx_2= %22.16lf;  \n", 1.0);
   fprintf(fp1,"tau_1= %22.16lf;   \n",  parameters[0].tttone);
   fprintf(fp1,"v_1= %22.16lf;   \n", parameters[0].vvvone);
   /*fprintf(fp1,"kappa_1= %22.16lf;   \n", parameters[0].kappaone);
   fprintf(fp1,"kappa_2= %22.16lf;   \n", parameters[0].kappatwo);*/

   fclose(fp1);
   //exit(0);
}
