#pragma once

#include "ComputeElement.tpp"

template<unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::assemble_system(){
    const QGauss<Nsd> quadrature_formula(quadOrder_);

    // the following class handles three things at once - finite element(basis functions), quadrature, and mapping from parent to real domains
    FEValues<Nsd> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values); //the list of what kind of information we need on each cell is given as a collection of flags as the third argument to the constructor of FEValues class 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    FullMatrix<double> Mlocal(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> Klocal(dofs_per_cell, dofs_per_cell);

    for(const typename DoFHandler<Nsd>::active_cell_iterator &elem : dof_handler.active_cell_iterators()){
        Mlocal = 0.0;
        Klocal = 0.0;

        elem->get_dof_indices(local_dof_indices);
        compute_element(elem, fe_values, Mlocal, Klocal, local_dof_indices);

        for(const unsigned int i : fe_values.dof_indices()){
            for(const unsigned int j : fe_values.dof_indices()){
                Mglobal.add(local_dof_indices[i], local_dof_indices[j], Mlocal(i,j));
                Kglobal.add(local_dof_indices[i], local_dof_indices[j], Klocal(i,j));
            }
        }
    }

    // std::cout << "System assembly for Mglobal and Kglobal completed" << std::endl;
}