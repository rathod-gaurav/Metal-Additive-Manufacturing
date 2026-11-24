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
        unsigned int n_iter = 0;
        double error0 = 0.0; //initial residual for divergence check
        double allowed_divergence = 0.0; //threshold for divergence

        for(unsigned int k = 0 ; k < N_NR_ ; k++){

            assemble_system(t); //uses sol_k (and sol_n) -> NR_jacobian = dR/dd, NR_residual = R

            double error = NR_residual.l2_norm();
            std::cout << "timestep: " << timestep << " | iteration: " << k << " | residual: " << error << std::endl;
            if(k == 0){ //store initial residual for divergence check
                error0 = error; 
                allowed_divergence = std::pow(error0, residual_divergence_factor_); 
            } 

            if(!std::isfinite(error) || error > allowed_divergence){ //if residual diverges
                std::cerr << "ERROR: Non-finite residual encountered at timestep " << timestep << ", iteration " << k << std::endl;
                break;
            }

            if(error < epsilon_NR_){ //if residual is converged!
                converged = true;
                n_iter = k;
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
            if (dt_ < dt_min_){ std::cerr << "dt below dt_min, Abort!\n"; break; }
            std::cout << "reducing dt to " << dt_ << std::endl;
            timestep--;
            continue; //retry the same timestep with a smaller dt
        }

        sol_n = sol_k; //accept the new time level
        std::cout << "convergence achieved for timestep " << timestep << " in " << n_iter << " iterations. Total time:" << t << std::endl;
        

        //debug - check max phi and its location
        std::vector<types::global_dof_index> dofs(fe.n_dofs_per_cell());
        double phi_max = -1e300;  Point<Nsd> x_max;  double h_max = 0;  int lev_max = 0;
        for (const auto& cell : dof_handler.active_cell_iterators()){
            cell->get_dof_indices(dofs);
            for (unsigned int i = 0; i < dofs.size(); ++i){
                if (fe.system_to_component_index(i).first == PHI_ && sol_n(dofs[i]) > phi_max){
                    phi_max = sol_n(dofs[i]);
                    x_max = cell->center();  h_max = cell->diameter()/std::sqrt(2.0);  lev_max = cell->level();
                }
            }
        }
        std::cout << "max phi = " << phi_max << " at " << x_max << " | dx = " << h_max << " | level = " << lev_max << std::endl;
        //debug - end        
        
        if (n_iter <= 3) {// easy step: grow dt
            dt_ = std::min(dt_increase_factor_*dt_, dt_max_);
            std::cout << "increasing dt to " << dt_ << std::endl;
        }

        if(timestep % 10 == 0){ refine_mesh(); } //refine mesh every 20 timesteps

        if(timestep % 10 == 0){
            output_writer_.write_vtu_and_pvd(dof_handler, sol_n, {"uc", "phi"}, timestep, pvd_output);
        }
    }

    output_writer_.finish_pvd(pvd_output); //finishes pvd file with footer

    std::cout << "Congratulations! Solve Completed." << std::endl;
}
