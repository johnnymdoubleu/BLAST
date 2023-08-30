#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N NimbleRun              
#$ -cwd                  
#$ -l h_rt=23:59:59 
#$ -l h_vmem=128G
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem

# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load R modules
module load igmm/apps/R/4.2.2
module load phys/compilers/gcc/9.3.0
# module load igmm/apps/jags/4.3.0
module load igmm/libs/lapack/3.5.0
module load igmm/libs/eigen/3.4.0
module load igmm/apps/hdf5/1.8.13
module load igmm/libs/gsl/2.6
module load igmm/libs/gdal/3.2.0
module load igmm/apps/sqlite3/3.33.0

# Run the program
Rscript ./Simulation/BayesianPsplines/1shot_listofmatrices.R

