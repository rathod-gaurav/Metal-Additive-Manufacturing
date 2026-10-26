#pragma once

template<unsigend int Nsd, unsigned int BfOrder>
template<typename Number>
void PureMelt<Nsd,BfOrder>::evaluate_fields(
    const FEValues<Nsd>& fe_values, const std::vector<Number>& local_dof_values, QPFields<Nsd,Number>& fields
) const {
    for(unsigned int f = 0; f < N_FIELDS_; f++){
        //current iterate : carries derivatives when Number = ADNumber
        fe_values[field_ext[f]].get_function_values_from_local_dof_values(local_dof_values, fields.val[f]);
        fe_values[field_ext[f]].get_function_gradients_from_local_dof_values(local_dof_values, fields.grad[f]);

        //previous iterate : plain doubles
        fe_values[field_ext[f]].get_function_values(sol_n, fields.val_old[f]);
        fe_values[field_ext[f]].get_function_gradients(sol_n, fields.grad_old[f]);
    }
}