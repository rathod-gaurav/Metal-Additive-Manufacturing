#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
void BinaryAlloy<Nsd,BfOrder>::flag_cells_for_refinement() {
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for(const typename DoFHandler<Nsd>::active_cell_iterator &elem : dof_handler.active_cell_iterators()){
        elem->clear_refine_flag();
        elem->clear_coarsen_flag();
        elem->get_dof_indices(local_dof_indices);

        double min_abs_phi = 1.0e10; //smallest |phi| on this cell
        for (unsigned int i = 0; i < dofs_per_cell; ++i){
            if (fe.system_to_component_index(i).first == PHI_){
                min_abs_phi = std::min(min_abs_phi, std::abs(sol_n(local_dof_indices[i])));
            }
        }

        const bool in_band = (min_abs_phi < 1.0 - delta_phi_); //Step 1
        const bool in_bulk = (min_abs_phi > 1.0 - delta_coarsen_); //Step 2

        if (in_band && elem->level() < static_cast<int>(max_level_)){
            elem->set_refine_flag();
        }
        else if (in_bulk && elem->level() > static_cast<int>(min_level_)){
            elem->set_coarsen_flag();
        }
    }
}

template<unsigned int Nsd, unsigned int BfOrder>
void BinaryAlloy<Nsd,BfOrder>::refine_mesh() {
    flag_cells_for_refinement(); //steps 1 and 2
    
    triangulation.prepare_coarsening_and_refinement(); //3.1
    SolutionTransfer<Nsd, Vector<double>> transfer(dof_handler);
    transfer.prepare_for_coarsening_and_refinement(sol_n); //3.2
    
    const Vector<double> sol_old = sol_n; //keep the old data: sol_n is resized below
    
    triangulation.execute_coarsening_and_refinement(); //3.3
    setup_system_dofs(); //3.4 (NO initial conditions here)

    transfer.interpolate(sol_old, sol_n); //3.5 (deal.II 9.5 signature)
    constraints.distribute(sol_n); //3.6
    sol_k = sol_n; //3.7
    
    std::cout << "remesh: cells = " << triangulation.n_active_cells() << " | DoFs = " << dof_handler.n_dofs() << std::endl;
}