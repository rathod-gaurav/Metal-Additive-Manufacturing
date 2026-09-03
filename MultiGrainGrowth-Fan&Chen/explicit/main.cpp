//Explicit Fan&Chen - DealII

#include <deal.II/grid/tria.h> //triangulation
#include <deal.II/dofs/dof_handler.h> //enumeration of degrees of freedom
#include <deal.II/grid/grid_generator.h> //grid generation

#include <deal.II/fe/fe_q.h> //Lagrange finite elements
#include <deal.II/dofs/dof_tools.h> //DoFHandler tools

#include <deal.II/fe/fe_values.h> //used to assemble matrix using quadrature on each cell
#include <deal.II/base/quadrature_lib.h> //quadrature rules

//need thiese three for treatment of boundary values
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

//linear algebra
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/base/table.h>
// #include <deal.II/lac/solver_cg.h>
// #include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/precondition.h>

//output
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <filesystem>

//random number generator
#include <random>

using namespace dealii;

#include <FanChen.hpp>
// #include <OutputWriter.hpp>

int main(){
    constexpr unsigned int Nsd = 2; //2 for 2D, 3 for 3D
    constexpr unsigned int BfOrder = 1; //1 for linear, 2 for quadratic

    unsigned int quadOrder = 2; //quadrature order for numerical integration
    
    //Problem parameters
    const double x_ll = 0.0; //lower left corner of the domain
    const double x_ul = 2.0; //upper right corner of the domain
    const unsigned int n_refinements = 9; //number of global refinements in the mesh grid

    unsigned int p = 36; //number of distinct grain orientations in the domain

    const unsigned int NT = 5; //number of time steps
    const double dt = 0.25; //time step size

    const double L = 1.0; //mobility
    const double kappa = 2.0; //gradient energy coefficient
    const double alpha = 1.0; //kinetic equation coefficient alpha
    const double beta = 1.0; //kinetic equation coefficient beta
    const double gamma = 1.0; //kinetic equation coefficient gamma

    OutputWriter<Nsd,BfOrder> output_writer("output");
    FanChen<Nsd, BfOrder> problem(
        x_ll, x_ul, n_refinements, //for creating mesh grid
        p, //number of distinct grain orientations in the domain
        dt, NT, //time step size, number of time steps
        L, kappa, //mobility, gradient energy coefficient //both assumed constant here
        alpha, beta, gamma, //kinetic equation coefficients
        quadOrder, //quadrature order for numerical integration
        output_writer //output writer object
    );

    problem.run();
}