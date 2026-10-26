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

        //nonlinear constitutive quantities
        Number lambda, tau;
        Tensor<1,Nsd,Number> dlambda_dq;
        model_.interface_coefficients(grad_phi, lambda, dlambda_dq, tau);

        const Number dfdphi   = model_.bulk_driving_force(phi, u);
        const Number hprime   = model_.latent_heat_slope(phi);

        //time derivatives (backward Euler)
        const Number dudt     = (u - u_n)/dt_;
        const Number dphidt   = (phi - phi_n)/dt_;

        //gradient + anisotropy flux: lambda^2 grad phi + |grad phi|^2 lambda dlambda/dq
        const Number grad_phi2 = grad_phi*grad_phi;
        const Number lambda2   = lambda*lambda;
        const Number coef      = grad_phi2*lambda;
        const Tensor<1,Nsd,Number> phi_flux = lambda2*grad_phi + coef*dlambda_dq;

        const double JxW = fe_values.JxW(q);

        for(const unsigned int i : fe_values.dof_indices()){
            //shape functions of each field (zero if dof i belongs to another field)
            const double          N_u       = fe_values[field_ext[TEMP_]].value(i,q);
            const Tensor<1,Nsd>   gradN_u   = fe_values[field_ext[TEMP_]].gradient(i,q);
            const double          N_phi     = fe_values[field_ext[PHI_]].value(i,q);
            const Tensor<1,Nsd>   gradN_phi = fe_values[field_ext[PHI_]].gradient(i,q);

            //Eq. 6a : temperature
            const Number R_u   = N_u*(dudt - 0.5*hprime*dphidt) + D*(gradN_u*grad_u);
            //Eq. 6b : phase field
            const Number R_phi = N_phi*(tau*dphidt + dfdphi) + gradN_phi*phi_flux;

            Rlocal[i] += (R_u + R_phi)*JxW;
        }
    }
}