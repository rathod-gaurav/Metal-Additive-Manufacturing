// Pure melt dendrite evolution (Implicit]) - KBhagat 2023

#include <deal.II/grid/tria.h> //triangulation
#include <deal.II/dofs/dof_handler.h> //enumeration of degrees of freedom
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h> //grid generation

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/differentiation/ad.h>

//linear algebra
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

using namespace dealii;
#include "PureMelt.hpp"

int main(){
    constexpr unsigned int Nsd = 2;
    constexpr unsigned int BfOrder = 1;

    //Mesh
    const double x_ll = -250.0;
    const double x_ul = 250.0;
    const unsigned int n_elements = 250; //number of elements per side

    //NR parameters
    const unsigned int NT = 10; //total number of timesteps
    const double dt = 0.01; //timestep size
    const unsigned int N_NR = 10; //total number of newton iterates per timestep
    const double epsilon_NR = 1e-6; //NR convergence tolerance


    PureMelt<Nsd,BfOrder> problem(
        x_ll, x_ul, n_elements,
        NT, dt, N_NR, epsilon_NR
    );

    problem.run();

}