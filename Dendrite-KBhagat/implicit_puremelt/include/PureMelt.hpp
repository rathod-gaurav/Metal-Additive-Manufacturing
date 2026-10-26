#pragma once

#include "OutputWriter.hpp"
#include "QPFields.hpp"

template<unsigned int Nsd, unsigned int BfOrder>
class PureMelt{
    public:
        PureMelt(
            const double x_ll, const double x_ul, const unsigned int n_elements, const unsigned int quadOrder, //mesh
            const unsigned int N_FIELDS, const double undercooling, Point<Nsd> seed_center, const double seed_radius, //initial conditions
            const unsigned int NT, const double dt, const unsigned int N_NR, const double epsilon_NR, //newton solver
            OutputWriter<Nsd,BfOrder>& output_writer
        );
        void run();
    
    private:
        const double x_ll_, x_ul_;
        const unsigned int n_elements_;
        const unsigned int quadOrder_;
        const unsigned int N_FIELDS_;
        const double undercooling_;
        Point<Nsd> seed_center_;
        const double seed_radius_;
        const unsigned int NT_;
        const double dt_;
        const unsigned int N_NR_;
        const double epsilon_NR_;
        OutputWriter<Nsd,BfOrder>& output_writer_;

        void make_grid();
        void setup_system();
        void assemble_system();
        void compute_element(
            const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, 
            FullMatrix<double>& Jlocal, Vector<double>& Rlocal,
            const std::vector<types::global_dof_index>& local_dof_indices
        );
        template<typename Number>
        void evaluate_fields(
            const FEValues<Nsd>& fe_values, const std::vector<Number>& local_dof_values, QPFields<Nsd,Number>& fields
        ) const;
        void solve();

        Triangulation<Nsd> triangulation;
        FESystem<Nsd> fe; //two variables per node
        DoFHandler<Nsd> dof_handler;
        unsigned int Nt;

        SparsityPattern sparsity_pattern;
        SparseMatrix<double> NR_jacobian;
        Vector<double> sol_k, sol_n; //solution at current iterate, solution at previous timestep
        Vector<double> NR_update;
        Vector<double> NR_residual;

};

#include "PureMelt.tpp"
#include "MakeGrid.tpp"
#include "SetupSystem.tpp"
#include "Solver.tpp"

