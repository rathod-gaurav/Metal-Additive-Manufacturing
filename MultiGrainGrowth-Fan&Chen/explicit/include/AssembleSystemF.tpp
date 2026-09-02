#pragma once

#include "ComputeElementF.tpp"

template <unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::assemble_system_F(){
    const QGauss<Nsd> quadrature_formula(quadOrder_);

    // the following class handles three things at once - finite element(basis functions), quadrature, and mapping from parent to real domains
    FEValues<Nsd> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values); //the list of what kind of information we need on each cell is given as a collection of flags as the third argument to the constructor of FEValues class

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    Vector<double> Flocal(dofs_per_cell);

    for(const typename DoFHandler<Nsd>::active_cell_iterator &elem : dof_handler.active_cell_iterators()){
        Flocal = 0.0;

        elem->get_dof_indices(local_dof_indices);
        compute_element_F(elem, fe_values, Flocal, local_dof_indices);

        for(const unsigned int i : fe_values.dof_indices()){
            Fglobal(local_dof_indices[i]) += Flocal(i);
        }
    }

    // std::cout << "System assembly for Fglobal completed" << std::endl;
}