#pragma once

template <unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::compute_element_F(const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, Vector<double>& Flocal, std::vector<types::global_dof_index>& local_dof_indices){
    fe_values.reinit(elem);

    for(const unsigned int q_index : fe_values.quadrature_point_indices()){
        double eta = 0.0;
        Vector<double> eta_nj(Nt);
        eta_nj = 0.0;

        std::cout << "initiating first for loop inside computeElementF" << std::endl;

        for(const unsigned int A : fe_values.dof_indices()){
            eta += fe_values.shape_value(A,q_index)*eta_ni(local_dof_indices[A]);
        }

        std::cout << "initiating second for loop inside computeElementF" << std::endl;
        double eta2_sum = 0.0;
        for(unsigned int j = 0 ; j < p_ ; j++){
            double eta_j = 0.0;
            std::copy(&eta_n[j][0], &eta_n[j][0] + Nt, eta_nj.begin());
            for(const unsigned int A : fe_values.dof_indices()){
                eta_j += fe_values.shape_value(A,q_index)*eta_nj(local_dof_indices[A]);
            }
            eta2_sum += eta_j*eta_j;
        }

        std::cout << "initiating third for loop inside computeElementF" << std::endl;

        for(const unsigned int A : fe_values.dof_indices()){
            Flocal(A) += (-1*alpha_*eta + beta_*pow(eta,3) + 2*gamma_*eta*(eta2_sum - eta*eta))*fe_values.shape_value(A,q_index)*fe_values.JxW(q_index);
        }
        
        
        
        
        // double etaA = 0.0;
        // double etaA_2 = 0.0;
        // for(const unsigned int i : fe_values.dof_indices()){//i = sum over all nodes
        //     etaA += fe_values.shape_value(i,q_index)*eta_ni(local_dof_indices[i]);

        //     for()//j = sum over all p field variables
        //     etaA_2 += fe_values.shape_value(i,q_index)*eta_n2(local_dof_indices[i]);
        // }

        // for(const unsigned int i : fe_values.dof_indices()){
        //     // std::cout << "i : " << i << std::endl;
        //     Flocal(i) += (-1*alpha_*etaA + beta_*pow(etaA,3) + 2*gamma_*etaA*(etaA_2 - etaA*etaA))*fe_values.shape_value(i,q_index)*fe_values.JxW(q_index);
        // }
    }
}