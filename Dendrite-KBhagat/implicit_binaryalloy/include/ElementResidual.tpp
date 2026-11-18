#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
template<typename Number>
void BinaryAlloy<Nsd,BfOrder>::element_residual(
    double t, const FEValues<Nsd>& fe_values, const QPFields<Nsd,Number>& fields, std::vector<Number>& Rlocal
) const {
    const double D_ = model_.diffusivity_non_dim();
    const double k = model_.partition_coefficient();
    const double lambda0 = model_.lambda0();
    const double v_ = model_.pulling_speed_non_dim();
    const double lT_ = model_.thermal_length_non_dim();
    const double xi_ = model_.xi();

    for(const unsigned int q : fe_values.quadrature_point_indices()){
        //fields at this quadrature point
        const Number&                uc        = fields.val[TEMP_][q];
        const Number&                phi      = fields.val[PHI_][q];
        const Tensor<1,Nsd,Number>&  grad_uc   = fields.grad[TEMP_][q];
        const Tensor<1,Nsd,Number>&  grad_phi = fields.grad[PHI_][q];
        const double                 uc_n      = fields.val_old[TEMP_][q];
        const double                 phi_n    = fields.val_old[PHI_][q];

        //nonlinear constitutive quantities
        Number as;
        Tensor<1,Nsd,Number> dasdq;
        model_.anisotropy<Nsd,Number>(grad_phi, as, dasdq);

        //time derivatives (backward euler)
        const Number dudt = (uc - uc_n)/dt_;
        const Number dphidt = (phi - phi_n)/dt_;

        //some required coefficients that multiply with dudt and dphidt
        const Number coef1 = ((1.0 + k)/2.0) + ((1.0 - k)/2.0)*phi;
        const Number coef2 = 1.0 + (1.0 - k)*uc;

        //uc flux : the grad_uc term in the weak form of the temperature equation
        const Tensor<1,Nsd,Number> flux_uc = grad_uc + (grad_phi/grad_phi.norm())*0.25*sqrt(2.0)*lambda0*coef2*dphidt;
        //phi flux : all grad_phi terms in the weak form of the phase field equation
        const Number as2 = as*as;
        const Number grad_phi2 = grad_phi*grad_phi;
        const Number as_dasdq = as*dasdq;
        const Tensor<1,Nsd,Number> flux_phi = as2*grad_phi + grad_phi2*as_dasdq;

        //first integral in phi residual
        const Point<Nsd> p = fe_values.quadrature_point(q);
        const double y = p[1];
        const Number coef3 = (y - v_*t)/lT_;
        const Number coef4 = 1.0 - (1.0 - k)*coef3;
        const Number coef5 = uc + coef3;
        const Number dfdphi = phi - phi*phi*phi - xi_*(1.0 - phi*phi)*(1.0 - phi*phi)*coef5;


        const double JxW = fe_values.JxW(q);

        for(const unsigned int i : fe_values.dof_indices()){
            //shape functions of each field (zero if dof i belongs to another field)
            const double          N_uc       = fe_values[field_ext[TEMP_]].value(i,q);
            const Tensor<1,Nsd>   gradN_uc   = fe_values[field_ext[TEMP_]].gradient(i,q);
            const double          N_phi     = fe_values[field_ext[PHI_]].value(i,q);
            const Tensor<1,Nsd>   gradN_phi = fe_values[field_ext[PHI_]].gradient(i,q);

            //Eq. 6a : temperature
            const Number R_uc = N_uc*(coef1*dudt - coef2*0.5*dphidt) + gradN_uc*flux_uc*D_*(0.5*(1.0-phi));
            //Eq. 6b : phase field
            const Number R_phi = N_phi*(coef4*as2*dphidt - dfdphi) + gradN_phi*flux_phi;

            Rlocal[i] += (R_uc + R_phi)*JxW;
        }
    }
}