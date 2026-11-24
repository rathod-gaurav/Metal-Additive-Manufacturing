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
    
    //debug - check coarsening
    auto count = [&](){ 
        unsigned int c = 0;
        for (const auto& cell : triangulation.active_cell_iterators()){c += cell->coarsen_flag_set();}
        return c;     
    };
        
    std::cout << "coarsen flags: " << count();

    triangulation.prepare_coarsening_and_refinement(); //3.1

    std::cout << " -> " << count() << std::endl; //debug - end

    SolutionTransfer<Nsd, Vector<double>> transfer(dof_handler);
    transfer.prepare_for_coarsening_and_refinement(sol_n); //3.2

    triangulation.execute_coarsening_and_refinement(); //3.3
    setup_system_dofs(); //3.4 (NO initial conditions here)

    transfer.interpolate(sol_n); //3.5 (deal.II 9.5 signature)
    constraints.distribute(sol_n); //3.6
    sol_k = sol_n; //3.7
    
    std::cout << "remesh: cells = " << triangulation.n_active_cells() << " | DoFs = " << dof_handler.n_dofs() << std::endl;
}

//check total solute upon mesh refinement (should be conserved)
template <unsigned int Nsd, unsigned int BfOrder>
double BinaryAlloy<Nsd,BfOrder>::compute_total_solute() const {
    const QGauss<Nsd> quad(quadOrder_);
    FEValues<Nsd> fe_values(fe, quad, update_values | update_JxW_values);
    const unsigned int n_q = quad.size();
    std::vector<double> uc(n_q), phi(n_q);

    const double k   = model_.partition_coefficient();
    const double cl0 = model_.c0()/k;                          // reference liquid composition
    double C = 0.0;

    for (const auto& cell : dof_handler.active_cell_iterators()){
        fe_values.reinit(cell);
        fe_values[field_ext[UC_]].get_function_values(sol_n, uc);
        fe_values[field_ext[PHI_]].get_function_values(sol_n, phi);
        for (unsigned int q = 0; q < n_q; ++q){
            const double c = 0.5*cl0*(1.0 - phi[q] + k*(1.0 + phi[q]))*(1.0 + (1.0 - k)*uc[q]);
            C += c*fe_values.JxW(q);
        }
    }
    return C;
}