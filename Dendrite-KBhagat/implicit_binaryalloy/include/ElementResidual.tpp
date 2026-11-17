#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
template<typename Number>
void PureMelt<Nsd,BfOrder>::element_residual(
    const FEValues<Nsd>& fe_values, const QPFields<Nsd,Number>& fields, std::vector<Number>& Rlocal
) const {
    const double D = model_.diffusivity();

    for(const unsigned int q : fe_values.quadrature_point_indices()){
        //fields at this quadrature point
        const Number&                u        = fields.val[TEMP_][q];
        const Number&                phi      = fields.val[PHI_][q];
        const Tensor<1,Nsd,Number>&  grad_u   = fields.grad[TEMP_][q];
        const Tensor<1,Nsd,Number>&  grad_phi = fields.grad[PHI_][q];
        const double                 u_n      = fields.val_old[TEMP_][q];
        const double                 phi_n    = fields.val_old[PHI_][q];

        //code here

        const double JxW = fe_values.JxW(q);

        for(const unsigned int i : fe_values.dof_indices()){
            //shape functions of each field (zero if dof i belongs to another field)
            const double          N_u       = fe_values[field_ext[TEMP_]].value(i,q);
            const Tensor<1,Nsd>   gradN_u   = fe_values[field_ext[TEMP_]].gradient(i,q);
            const double          N_phi     = fe_values[field_ext[PHI_]].value(i,q);
            const Tensor<1,Nsd>   gradN_phi = fe_values[field_ext[PHI_]].gradient(i,q);

            //Eq. 6a : temperature
            const Number R_u   = ;
            //Eq. 6b : phase field
            const Number R_phi = ;

            Rlocal[i] += (R_u + R_phi)*JxW;
        }
    }
}