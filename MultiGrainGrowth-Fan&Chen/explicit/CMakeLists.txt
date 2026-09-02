cmake_minimum_required(VERSION 3.16)
project(Fan&Chen CXX)

# Force your preferred build type here (e.g., "Debug" or "Release")
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type" FORCE)

# Point CMake to your deal.II installation path
find_package(deal.II 9.7.1 REQUIRED
  HINTS /home/grathod/Downloads/dealii/dealii/installation_files
#   HINTS $ENV{DEAL_II_DIR}
)

include_directories("include")

# Initialize cached variables required by deal.II
DEAL_II_INITIALIZE_CACHED_VARIABLES()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -march=native")

# Define your executable and source file
add_executable(main main.cpp)

# Link deal.II and all dependencies (MPI, PETSc, Trilinos, p4est)
DEAL_II_SETUP_TARGET(main)
