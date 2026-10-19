// Pure melt dendrite evolution (Implicit]) - KBhagat 2023

#include <deal.II/grid/tria.h> //triangulation
#include <deal.II/dofs/dof_handler.h> //enumeration of degrees of freedom
#include <deal.II/grid/grid_generator.h> //grid generation

using namespace dealii;
#include "PureMelt.hpp"

int main(){
    constexpr unsigned int Nsd = 1;
    constexpr unsigned int BfOrder = 2;

    //Mesh
    const double x_ll = -250.0;
    const double x_ul = 250.0;
    const unsigned int n_elements = 250; //number of elements per side

    PureMelt<Nsd,BfOrder> problem(
        x_ll, x_ul, n_elements
    );

    problem.run();

}