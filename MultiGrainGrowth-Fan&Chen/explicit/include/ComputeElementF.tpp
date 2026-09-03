#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::compute_element_F(const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, Vector<double>& Flocal, std::vector<types::global_dof_index>& local_dof_indices){
    fe_values.reinit(elem);

    for(const unsigned int q_index : fe_values.quadrature_point_indices()){
        double etaA = 0.0;
        double etaA_2 = 0.0;
        for(const unsigned int i : fe_values.dof_indices()){
            etaA += fe_values.shape_value(i,q_index)*eta_ni(local_dof_indices[i]);
            etaA_2 += fe_values.shape_value(i,q_index)*eta_n2(local_dof_indices[i]);
        }

        for(const unsigned int i : fe_values.dof_indices()){
            Flocal(i) += -1*alpha_*etaA + beta_*pow(etaA,3) + 2*gamma_*etaA*(etaA_2 - etaA*etaA);
        }
    }
}