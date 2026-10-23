#pragma once

template<unsigned int Nsd, unsigned int Nne>
class PureMelt{
    public:
        PureMelt(
            const double x_ll, const double x_ul, const unsigned int n_elements, //mesh
            const double undercooling, const Point<Nsd> seed_center, const double seed_radius, //initial conditions
            const unsigned int NT, const double dt, const unsigned int N_NR, const double epsilon_NR //newton solver
        );
        void run();
    
    private:
        const double x_ll_, x_ul_;
        const unsigned int n_elements_;
        const double undercooling_;
        const Point<Nsd> seed_center_;
        const double seed_radius_;
        const unsigned int NT_;
        const double dt_;
        const unsigned int N_NR_;
        const double epsilon_NR_;

        void make_grid();
        void setup_system();
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

