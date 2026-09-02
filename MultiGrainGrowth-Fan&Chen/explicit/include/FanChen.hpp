#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
class FanChen{
    public:
        FanChen(
            const double x_ll, const double x_ul, const unsigned int n_refinements, //for creating mesh grid
            const unsigned int p //number of distinct grain orientations in the domain
        );
        void run();
    
    private:
        const double x_ll_, x_ul_;
        const unsigned int n_refinements_;
        const unsigned int p_;

        void make_grid();
        void setup_system();
        // void compute_element();
        // void compute_element_F();
        // void assemble_system();
        // void assemble_system_F();
        // void solve();

        Triangulation<Nsd> triangulation;
        const FE_Q<Nsd> fe;
        DoFHandler<Nsd> dof_handler;
        unsigned int Nt; //total number of degrees of freedom

        SparsityPattern sparsity_pattern;
        SparseMatrix<double> Mglobal, Kglobal;
        Vector<double> Fglobal;
        Table<2, double> eta_n, eta_np1;
        Vector<double> eta_ni, eta_np1i, eta_n2, RHS;
};

#include <FanChen.tpp>

#include <MakeGrid.tpp>
#include <SetupSystem.tpp>