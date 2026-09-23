#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::compute_element_F(const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, Vector<double>& Flocal, std::vector<types::global_dof_index>& local_dof_indices){
    fe_values.reinit(elem);

    for(const unsigned int A : fe_values.dof_indices()){
        eta_ni_local(A) = eta_ni(local_dof_indices[A]);
        for(unsigned int j = 0 ; j < p_ ; j++){
            eta_nj_local(j*dofs_per_cell + A) = eta_n[j][local_dof_indices[A]];
        }
    }

    for(const unsigned int q_index : fe_values.quadrature_point_indices()){
        double eta = 0.0;

        for(const unsigned int A : fe_values.dof_indices()){
            eta += fe_values.shape_value(A,q_index)*eta_ni_local(A);
        }

        double eta2_sum = 0.0;
        for(unsigned int j = 0 ; j < p_ ; j++){
            double eta_j = 0.0;
            for(const unsigned int A : fe_values.dof_indices()){
                eta_j += fe_values.shape_value(A,q_index)*eta_nj_local(j*dofs_per_cell + A);
            }
            eta2_sum += eta_j*eta_j;
        }

        double prefactor = (-1*alpha_*eta + beta_*eta*eta*eta + 2*gamma_*eta*(eta2_sum - eta*eta))*fe_values.JxW(q_index);
        for(const unsigned int A : fe_values.dof_indices()){
            Flocal(A) += prefactor*fe_values.shape_value(A,q_index);
        }
        
    }
}