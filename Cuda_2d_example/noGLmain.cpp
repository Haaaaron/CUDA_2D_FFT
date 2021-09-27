#include "a.h"

field globalgassed[1];
modelparams globalparameters[1];

__host__ void savedata(field gassed, modelparams parameters, int n);
__host__ void pointerinit(field *gassed,modelparams *parameters);
//extern "C" __host__ void fft_tests(field gassed, modelparams parameters);
extern "C" __host__ void initfield(field gassed, modelparams parameters);
extern "C" __host__ void advance_in_time(field gassed,modelparams parameters, int i, int nsteps);

// Main program
int main(int argc, char** argv)
{
  cudaSetDevice(0);

  double max_time=5*24*3600;

  field gassed[1];
  modelparams parameters[1];

  cudaDeviceSynchronize();
  printf("\n Test \n");

  time_t *ttttime=new time_t[3];
  time(&ttttime[0]);
  pointerinit(&gassed[0],&parameters[0]);
  initfield(gassed[0], parameters[0]);
  savedata(gassed[0], parameters[0], 0);
  time(&ttttime[1]);
  double time_to_init=difftime(ttttime[1], ttttime[0]);
  printf("Time spent to initialize  %lf. \n", time_to_init);
  for(int in=1;in<=1000;in++)
  {
        advance_in_time(gassed[0], parameters[0],in-1, 1000);
        savedata(gassed[0], parameters[0], in);
        time(&ttttime[2]);
        double time_per_step=difftime(ttttime[2],ttttime[1])/in; // time per step in seconds
        printf("\n So far, on average  %lf seconds spent per step \n", time_per_step);
        double time_elapsed=difftime(ttttime[2],ttttime[0]);
        double time_left=max_time-time_elapsed;
        printf("Step %04d. Time from start %lf, time left %lf. time per step %lf\n",in, time_elapsed, time_left, time_per_step);
        if(time_left < 3.0*time_per_step)
        {
             printf("\n Time to exit. On average  %lf  hours spent per step \n", time_per_step/3600);
             exit(0);
        }
}

  cudaDeviceSynchronize();
  printf("\n The End??? \n \n");
  /*cudaThreadExit();
  cudaDeviceReset();*/
}
