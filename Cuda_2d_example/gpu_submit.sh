#!/bin/bash

#SBATCH --time=5-00:00:00          ## wallclock time hh:mm:ss
#SBATCH --gres=gpu:1

module load cuda
module load hdf5/1.10.2

## run my GPU accelerated executable, note the --gres
srun --gres=gpu:1 ./graphene.exe
