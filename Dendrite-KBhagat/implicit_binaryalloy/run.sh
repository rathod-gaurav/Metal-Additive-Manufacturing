#!/bin/sh
#SBATCH --partition=compphys2026 # default "shared", if not specified
#SBATCH --time=0-23:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1               # require 1 nodes
#SBATCH --ntasks-per-node=1    # cpus per node (by default, "ntasks"="cpus")
#SBATCH --mem=16000             # RAM per node in megabytes
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
# Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output

# Load modules and dealii directory
module load gcc/13.2.0
module load cmake/3.27.9-gcc-13.2.0
module load openmpi/5.0.3-gcc-13.2.0
module load openblas/0.3.26-gcc-13.2.0-nsrzk5g

export DEAL_II_DIR=/software/groups/rudraraju_group/libs/dealii-9.7.1

# User Script
cmake -DCMAKE_BUILD_TYPE=Release .
cmake -B build
cmake --build build -j4
time ./build/main