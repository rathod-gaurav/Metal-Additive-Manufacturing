#pragma once

#include "ComputeElement.tpp"

template<unsigned int Nsd, unsigned int BfOrder>
void PureMelt<Nsd,BfOrder>::assemble_system(){
    NR_jacobian = 0.0;
    NR_residual = 0.0;

    const QGauss<Nsd> quadrature_formula(quadOrder_);

    // the following class handles three things at once - finite element(basis functions), quadrature, and mapping from parent to real domains
    FEValues<Nsd> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    FullMatrix<double> Jlocal(dofs_per_cell, dofs_per_cell);
    Vector<double> Rlocal(dofs_per_cell);

    for(const typename DoFHandler<Nsd>::active_cell_iterator &elem : dof_handler.active_cell_iterators()){
        Jlocal = 0.0;
        Rlocal = 0.0;

        elem->get_dof_indices(local_dof_indices);
        compute_element(elem, fe_values, Jlocal, Rlocal, local_dof_indices);

        for(const unsigned int i : fe_values.dof_indices()){
            NR_residual(local_dof_indices[i]) += Rlocal(i);
            for(const unsigned int j : fe_values.dof_indices()){
                NR_jacobian.add(local_dof_indices[i], local_dof_indices[j], Jlocal(i,j));
            }
        }
    }
}