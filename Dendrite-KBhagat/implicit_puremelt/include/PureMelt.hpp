#pragma once

template<unsigned int Nsd, unsigned int Nne>
class PureMelt{
    public:
        PureMelt(
            const double x_ll, const double x_ul, const unsigned int n_elements
        );
        void run();
    
    private:
        const double x_ll_, x_ul_;
        const unsigned int n_elements_;

        void make_grid();
        void setup_system();

        Triangulation<Nsd> triangulation;
        FESystem<Nsd> fe; //two variables per node
        DoFHandler<Nsd> dof_handler;
        unsigned int Nt;

        SparsityPattern sparsity_pattern;
        SparseMatrix<double> system_jacobian;
        Vector<double> solution; //solution at current timestep
        Vector<double> solution_old; //solution at previous timestep
        Vector<double> newton_update;
        Vector<double> system_residual;


};

#include "PureMelt.tpp"
#include "MakeGrid.tpp"
#include "SetupSystem.tpp"
