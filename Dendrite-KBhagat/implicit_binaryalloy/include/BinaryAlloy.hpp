#pragma once

#include "OutputWriter.hpp"
#include "QPFields.hpp"
#include "InitialConditions.hpp"
#include "ConstitutiveModel.hpp"

template<unsigned int Nsd, unsigned int BfOrder>
class BinaryAlloy{
    public:
        BinaryAlloy(
            const double x_ll, const double x_ul, const unsigned int n_elements, const unsigned int quadOrder, //mesh
            const unsigned int UC, const unsigned int PHI, const unsigned int N_FIELDS, 
            const double c0, const std::vector<NucleationSeed<Nsd>>& seeds, //initial conditions
            const unsigned int NT, const double dt, const unsigned int N_NR, const double epsilon_NR, //newton solver
            const ConstitutiveModel& model, 
            OutputWriter<Nsd,BfOrder>& output_writer
        );
        void run();
    
    private:
        const double x_ll_, x_ul_;
        const unsigned int n_elements_;
        const unsigned int quadOrder_;
        const unsigned int UC_;
        const unsigned int PHI_;
        const unsigned int N_FIELDS_;
        const double c0_;
        const std::vector<NucleationSeed<Nsd>> seeds_;
        const unsigned int NT_;
        const double dt_;
        const unsigned int N_NR_;
        const double epsilon_NR_;
        const ConstitutiveModel& model_;
        OutputWriter<Nsd,BfOrder>& output_writer_;

        void make_grid();
        void setup_system();
        void assemble_system(double t);
        void compute_element(
            double t,
            const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, 
            FullMatrix<double>& Jlocal, Vector<double>& Rlocal,
            const std::vector<types::global_dof_index>& local_dof_indices
        );
        UClate<typename Number>
        void evaluate_fields(
            const FEValues<Nsd>& fe_values, const std::vector<Number>& local_dof_values, QPFields<Nsd,Number>& fields
        ) const;
        UClate<typename Number>
        void element_residual(
            double t, const FEValues<Nsd>& fe_values, const QPFields<Nsd,Number>& fields, std::vector<Number>& Rlocal
        ) const;
        void solve();

        //--- debugging: AD Jacobian vs finite differences ---
        void verify_jacobian(const unsigned int n_elems_to_check);

        Triangulation<Nsd> triangulation;
        FESystem<Nsd> fe; //two variables per node
        DoFHandler<Nsd> dof_handler;
        unsigned int Nt;
        std::vector<FEValuesExtractors::Scalar> field_ext; //field_ext[f] picks fielf f

        SparsityPattern sparsity_pattern;
        SparseMatrix<double> NR_jacobian;
        Vector<double> sol_k, sol_n; //solution at current iterate, solution at previous timestep
        Vector<double> NR_update;
        Vector<double> NR_residual;

};

#include "BinaryAlloy.tpp"
#include "MakeGrid.tpp"
#include "SetupSystem.tpp"
#include "VerifyJacobian.tpp"
#include "Solver.tpp"

