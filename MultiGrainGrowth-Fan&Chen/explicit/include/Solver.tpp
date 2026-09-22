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

    post_process(); //to store the initial conditions in result file
    output_writer_.write_vtu(dof_handler, phi, 0);

    double t = dt_;
    for(unsigned int timestep = 1 ; timestep < NT_ ; timestep++){
        //Assemble eta_n2
        // for (unsigned int k = 0; k < p_; k++){
        //     const double *const row = &eta_n[k][0];
        //     for (unsigned int j = 0; j < Nt; j++){
        //         eta_n2(j) += row[j] * row[j];
        //     }
        // }

        //debug
        double max_eta = 0.0;
        double Fglobal_norm = 0.0;

        for(unsigned int i = 0 ; i < p_ ; i++){
            // std::cout << "solver initiated for timestep " << timestep << " | Solve: " << i << std::endl;
            
            //copy ith row from eta_n to eta_ni
            std::copy(&eta_n[i][0], &eta_n[i][0] + Nt, eta_ni.begin());

            Fglobal = 0.0;
            assemble_system_F();

            // std::cout << "Assembled Fglobal for timestep " << timestep << " | Solve: " << i << std::endl;

            // constraints.distribute(eta_ni);

            //RHS = (Mglobal - dt_*L_*kappa_*Kglobal)*eta_ni - dt_*L_*Fglobal;
            Kglobal.vmult(RHS, eta_ni);
            RHS *= -1*dt_*L_*kappa_;
            Mglobal.vmult_add(RHS, eta_ni);
            Fglobal *= -1*dt_*L_;
            RHS += Fglobal;

            // directsolver.vmult(eta_np1i, RHS);

            eta_np1i = eta_ni;
            SolverControl control(1000, 1e-10 * RHS.l2_norm());
            SolverCG<Vector<double>> cgsolver(control);
            cgsolver.solve(Mglobal, eta_np1i, RHS, prec);

            // constraints.distribute(eta_np1i);

            //debug
            for(unsigned int i = 0 ; i < p_ ; i++){
                max_eta = std::max(max_eta, eta_np1i.linfty_norm());
                Fglobal_norm = std::max(Fglobal_norm, Fglobal.l2_norm());
            }

            std::cout << "Timestep: " << timestep << " | Solve: " << i << " | Iterations: " << control.last_step() << " | Residual: " << control.last_value() << std::endl;

            std::copy(eta_np1i.begin(), eta_np1i.end(), &eta_np1[i][0]);
        }
        std::cout << "Solve completed for timestep " << timestep << " | Max Eta: " << max_eta << " | Max Fglobal: " << Fglobal_norm << std::endl;
        std::cout << "-----------------------------" << std::endl;
        std::swap(eta_n, eta_np1);

        post_process();
        if(timestep%10 == 0){
            output_writer_.write_vtu(dof_handler, phi, timestep);
        }
        
        t += dt_;

    }

    std::cout << "Solve completed. Congratulations!!" << std::endl;

    output_writer_.write_pvd();
}