#pragma once

#include "OutputWriter.hpp"

template<unsigned int Nsd, unsigned int BfOrder>
class FanChen{
    public:
        FanChen(
            const double x_ll, const double x_ul, const unsigned int n_refinements, //for creating mesh grid
            const unsigned int p, //number of distinct grain orientations in the domain
            const double dt, const unsigned int NT, //time step size, number of time steps
            const double L, const double kappa, //mobility, gradient energy coefficient //both assumed constant here
            const double alpha, const double beta, const double gamma, //kinetic equation
            const unsigned int quadOrder, //quadrature order for numerical integration
            OutputWriter<Nsd,BfOrder>& output_writer //pointer to output writer object
        );
        void run();
    
    private:
        const double x_ll_, x_ul_;
        const unsigned int n_refinements_;
        const unsigned int p_;
        const double dt_;
        const unsigned int NT_;
        const double L_, kappa_;
        const double alpha_, beta_, gamma_;
        const unsigned int quadOrder_;
        OutputWriter<Nsd,BfOrder>& output_writer_;

        void make_grid();
        void setup_system();
        void compute_element(const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, FullMatrix<double>& Mlocal, FullMatrix<double>& Klocal, std::vector<types::global_dof_index>& local_dof_indices);
        void compute_element_F(const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, Vector<double>& Flocal, std::vector<types::global_dof_index>& local_dof_indices);
        void assemble_system();
        void assemble_system_F();
        void solve();
        void post_process();

        Triangulation<Nsd> triangulation;
        const FE_Q<Nsd> fe;
        DoFHandler<Nsd> dof_handler;
        unsigned int Nt; //total number of degrees of freedom

        SparsityPattern sparsity_pattern;
        SparseMatrix<double> Mglobal, Kglobal;
        Vector<double> Fglobal;
        Table<2, double> eta_n, eta_np1;
        Vector<double> eta_ni, eta_np1i, eta_n2, RHS;
        Vector<double> phi;
};

#include <FanChen.tpp>

#include <MakeGrid.tpp>
#include <SetupSystem.tpp>
#include <Solver.tpp>