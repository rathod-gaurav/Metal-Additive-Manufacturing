#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::compute_element(const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, FullMatrix<double>& Mlocal, FullMatrix<double>& Klocal, std::vector<types::global_dof_index>& local_dof_indices){
    fe_values.reinit(elem);

    for(const unsigned int q_index : fe_values.quadrature_point_indices()){
        for(const unsigned int i : fe_values.dof_indices()){
            for(const unsigned int j : fe_values.dof_indices()){
                Klocal(i,j) += (fe_values.shape_grad(i,q_index)*fe_values.shape_grad(j,q_index)*fe_values.JxW(q_index));
                Mlocal(i,j) +=  (fe_values.shape_value(i,q_index)*fe_values.shape_value(j,q_index)*fe_values.JxW(q_index));
            }
        }
    }
}