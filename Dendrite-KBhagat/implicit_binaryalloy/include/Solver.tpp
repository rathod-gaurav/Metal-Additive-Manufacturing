#pragma once

#include "AssembleSystem.tpp"
// #include "VerifyJacobian.tpp"

template <unsigned int Nsd, unsigned int BfOrder>
void BinaryAlloy<Nsd,BfOrder>::solve(){

    std::cout << "Initiating solver..." << std::endl;

    double t = 0.0;

    std::ofstream pvd_output(output_writer_.get_output_dir() + "/final_solution.pvd");
    output_writer_.initiate_pvd(pvd_output); //initiated pvd file with headers
    output_writer_.write_vtu_and_pvd(dof_handler, sol_n, {"uc", "phi"}, 0, pvd_output);

    for(unsigned int timestep = 1 ; timestep <= NT_ ; timestep++){
        t += dt_;
        sol_k = sol_n; //initial guess for the new time level

        //NR loops
        bool converged = false;
        for(unsigned int k = 0 ; k < N_NR_ ; k++){

            assemble_system(t); //uses sol_k (and sol_n) -> NR_jacobian = dR/dd, NR_residual = R

            double error = NR_residual.l2_norm();
            std::cout << "timestep: " << timestep << " | iteration: " << k << " | residual: " << error << std::endl;

            if(error < epsilon_NR_){
                converged = true;
                std::cout << "convergence achieved for timestep " << timestep << " in " << k << " iterations. Total time:" << t << std::endl;
                if (k <= 4) {// easy step: grow dt
                    dt_ = std::min(dt_increase_factor_*dt_, dt_max_);
                    std::cout << "increasing dt to " << dt_ << std::endl;
                }
                break;
            }

            //Solve J dd = -R
            NR_residual *= -1.0;
            NR_update = 0.0;
            SparseDirectUMFPACK directsolver;
            directsolver.initialize(NR_jacobian);
            directsolver.vmult(NR_update, NR_residual);

            constraints.distribute(NR_update); //distributes the update to the constrained dofs

            sol_k += NR_update;
        }

        if(!converged){
            std::cout << "WARNING: Newton did not converge in " << N_NR_ << " iterations at timestep " << timestep << ". Executing timestep refinement." << std::endl;
            sol_k = sol_n;
            dt_ *= dt_reduction_factor_; //halve the timestep
            if (dt_ < dt_min_){ std::cerr << "dt below dt_min, abort!\n"; break; }
            std::cout << "reducing dt to " << dt_ << std::endl;
            timestep--;
            continue; //retry the same timestep with a smaller dt
        }

        sol_n = sol_k; //accept the new time level

        if(timestep % 50 == 0){ refine_mesh(); } //refine mesh every 50 timesteps

        if(timestep % 10 == 0){
            output_writer_.write_vtu_and_pvd(dof_handler, sol_n, {"u", "phi"}, timestep, pvd_output);
        }
    }

    output_writer_.finish_pvd(pvd_output); //finishes pvd file with footer

    std::cout << "Congratulations! Solve Completed." << std::endl;
}
