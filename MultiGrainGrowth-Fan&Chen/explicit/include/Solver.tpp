#pragma once

#include "AssembleSystem.tpp"
#include "AssembleSystemF.tpp"
#include "PostProcess.tpp"

template<unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::solve(){
    std::cout << "Initiating solver..." << std::endl;

    std::cout << "-----------------------------" << std::endl;
    Mglobal = 0.0;
    Kglobal = 0.0;
    assemble_system();

    PreconditionJacobi<SparseMatrix<double>> prec; //PreconditionSSOR in case of SSOR
    prec.initialize(Mglobal, 1.0);

    // SparseDirectUMFPACK directsolver;
    // directsolver.initialize(Mglobal);

    double t = dt_;
    for(unsigned int timestep = 1 ; timestep < NT_ ; timestep++){
        //Assemble eta_n2
        for (unsigned int k = 0; k < p_; ++k){
            const double *const row = &eta_n[k][0];
            for (unsigned int j = 0; j < Nt; ++j){
                eta_n2(j) += row[j] * row[j];
            }
        }

        for(unsigned int i = 0 ; i < p_ ; i++){
            //copy ith row from eta_n to eta_ni
            std::copy(&eta_n[i][0], &eta_n[i][0] + Nt, eta_ni.begin());

            Fglobal = 0.0;
            assemble_system_F();

            //RHS = (Mglobal - dt_*L_*kappa_*Kglobal)*eta_ni - dt_*L_*Mglobal*Fglobal;
            Kglobal.vmult(RHS, eta_ni);
            RHS *= -1*dt_*L_*kappa_;
            Mglobal.vmult_add(RHS, eta_ni);
            Fglobal *= -1*dt_*L_;
            Mglobal.vmult_add(RHS, Fglobal);

            // directsolver.vmult(eta_np1i, RHS);

            eta_np1i = eta_ni;
            SolverControl control(200, 1e-10 * RHS.l2_norm());
            SolverCG<Vector<double>> cgsolver(control);
            cgsolver.solve(Mglobal, eta_np1i, RHS, prec);

            // std::cout << control.last_step() << std::endl;

            std::copy(eta_np1i.begin(), eta_np1i.end(), &eta_np1[i][0]);
        }
        std::cout << "Solve completed for timestep " << timestep << std::endl;
        std::cout << "-----------------------------" << std::endl;
        
        post_process();
        output_writer_.write_vtu(dof_handler, phi, timestep);

        std::swap(eta_n, eta_np1);
        t += dt_;

    }

    std::cout << "Solve completed." << std::endl;

    output_writer_.write_pvd();
}