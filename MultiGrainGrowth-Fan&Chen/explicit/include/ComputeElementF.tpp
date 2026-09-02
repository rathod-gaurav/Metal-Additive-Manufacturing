#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::compute_element_F(const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, Vector<double>& Flocal, std::vector<types::global_dof_index>& local_dof_indices){
    fe_values.reinit(elem);

    double etaA = 0.0;
    double etaA_2 = 0.0;

    for(const unsigned int i : fe_values.dof_indices()){
        Flocal(i) += -1*alpha_*etaA + beta_*pow(etaA,3) + 2*gamma_*etaA*(etaA_2 - etaA*etaA);
    }
}