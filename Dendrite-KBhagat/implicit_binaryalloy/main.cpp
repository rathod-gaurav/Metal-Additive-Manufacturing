// Binary Alloy dendrite evolution (Implicit]) - KBhagat 2023

#include <random>

#include <deal.II/grid/tria.h> //triangulation
#include <deal.II/dofs/dof_handler.h> //enumeration of degrees of freedom
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h> //grid generation
#include <deal.II/base/point.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/fe/fe_values.h> //used to assemble matrix using quadrature on each cell
#include <deal.II/base/quadrature_lib.h> //quadrature rules

#include <deal.II/differentiation/ad.h>

//linear algebra
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

//output
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <cmath>

//Adaptive meshing
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

using namespace dealii;
#include "InitialConditions.hpp"
#include "BinaryAlloy.hpp"
#include "ConstitutiveModel.hpp"
#include "OutputWriter.hpp"

int main(){
    constexpr unsigned int Nsd = 2;
    constexpr unsigned int BfOrder = 1;
    const unsigned int quadOrder = BfOrder + 1;

    //Mesh
    const double x_ll = -300.0;
    const double x_ul = 300.0;
    const unsigned int n_elements = 50; //number of elements per side

    //NR parameters
    const unsigned int NT = 100000; //total number of timesteps
    double dt = 0.0005; //timestep size
    const unsigned int N_NR = 20; //total number of newton iterates per timestep
    const double epsilon_NR = 1e-6; //NR convergence tolerance

    //Field variables of the problem (component index inside FESystem)
    static constexpr unsigned int UC = 0; //uc : composition
    static constexpr unsigned int PHI = 1; //phi : phase field order parameter
    static constexpr unsigned int N_FIELDS = 2; //total two fields in this problem

    //Initial conditions - randomly distributed nucleation seeds in the domain
    const double supersaturation = -1.0;
    const unsigned int Nseeds = 5; //number of nucleation seeds
    std::vector<NucleationSeed<Nsd>> seeds(Nseeds);
    //for multiple equiaxed case
    seeds = generate_random_seeds<Nsd>(Nseeds, x_ll, x_ul);
    //for single columnar case
    // seeds[0].center = Point<Nsd>(0.0, x_ll+10.0); //seed at the bottom of the domain
    // seeds[0].radius = 5.0; //seed radius
    //for multiple columnar case
    // for(unsigned int i = 0 ; i < Nseeds ; i++){
    //     double x_pos = i*((x_ul - x_ll)/(Nseeds-1)) + x_ll; //evenly spaced seeds along the bottom of the domain
    //     seeds[i].center = Point<Nsd>(x_pos, x_ll + 5.0); //seed at the bottom of the domain
    //     seeds[i].radius = 5.0; //seed radius
    // }

    //Adaptive Meshing parameters
    const unsigned int max_level = 5; //maximum refinement level
    const unsigned int min_level = 0; //minimum refinement level
    const double delta_phi = 1e-3; //threshold for refinement (|phi| < 1 - delta_phi)
    const double delta_coarsen = 1e-4; //threshold for coarsening (|phi| > 1 - delta_coarsen)

    //Adaptive time stepping parameters
    const double dt_max = 5*dt; //maximum timestep size
    const double dt_min = 1e-6; //minimum timestep size
    const double dt_reduction_factor = 0.5; //factor to reduce dt when NR fails to converge
    const double dt_increase_factor = 2.0; //factor to increase dt when NR converges quickly
    const unsigned int residual_divergence_factor = 2; //factor to check for divergence of residual

    //code here - for constitutive model parameters, if any
    const double D_l = 2.4*1e-9; //diffusivity of solute in liquid 
    const double D_s = 1.15*1e-12; //diffusivity of solute in solid
    const double k = 0.14; //partition coefficient
    const double eps4 = 0.01; //anisotropy strength
    const double Gamma_g = 2.36*1e-7; //gibbs-thomson coefficient
    const double m_l = -3.5; //liquidus slope
    const double c0 = 4.0; //initial solute concentration
    const double G = 0.0; //thermal gradient
    const double Vp = 0.0; //pulling velocity
    const double lambda0 = 1.058*1e-6; //interface thickness
    const double a1 = 0.8839; //model constant 1
    const double a2 = 0.6267; //model constant 2

    ConstitutiveModel model(//params here
        D_l, D_s, k,
        eps4, Gamma_g, m_l, c0, 
        G, Vp, 
        lambda0, a1, a2
    );

    OutputWriter<Nsd,BfOrder> output_writer("output");

    BinaryAlloy<Nsd,BfOrder> problem(
        x_ll, x_ul, n_elements, quadOrder,
        UC, PHI, N_FIELDS, supersaturation, seeds,
        NT, dt, N_NR, epsilon_NR,
        model, 
        output_writer,
        max_level, min_level, delta_phi, delta_coarsen,
        dt_max, dt_min, dt_reduction_factor, dt_increase_factor, residual_divergence_factor
    );

    problem.run();

}