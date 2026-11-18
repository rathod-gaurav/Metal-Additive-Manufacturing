#pragma once

//---------------------------------------------------------------------------
// DEBUG  -  check the AD Jacobian against finite differences
//
// For the first n elements:  J_AD(i,j)  vs  [R(d + h e_j) - R(d)]_i / h
// using the double version of evaluate_fields() + element_residual().
// A relative error ~1e-6 or smaller means the AD Jacobian is consistent
// with the residual (i.e. Newton will converge quadratically).
// Run once after changing ElementResidual.tpp or ConstitutiveModel.hpp.
//---------------------------------------------------------------------------

template <unsigned int Nsd, unsigned int BfOrder>
void BinaryAlloy<Nsd,BfOrder>::verify_jacobian(const unsigned int n_elems_to_check){
    const QGauss<Nsd> quadrature_formula(quadOrder_);
    FEValues<Nsd> fe_values(fe, quadrature_formula, update_values | update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q = quadrature_formula.size();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    FullMatrix<double> Jlocal(dofs_per_cell, dofs_per_cell);
    Vector<double> Rlocal(dofs_per_cell);

    //residual of one element with plain doubles
    auto residual_double = [&](const std::vector<double>& d, std::vector<double>& R){
        QPFields<Nsd,double> fields(N_FIELDS, n_q);
        evaluate_fields<double>(fe_values, d, fields);
        std::fill(R.begin(), R.end(), 0.0);
        element_residual<double>(fe_values, fields, R);
    };

    //check elements that cut the interface (where all terms are active)
    unsigned int checked = 0;
    double worst = 0.0;
    for(const typename DoFHandler<Nsd>::active_cell_iterator &elem : dof_handler.active_cell_iterators()){
        if(checked == n_elems_to_check){ break; }
        if(std::abs(elem->center().norm() - seed_radius_) > 1.0){ continue; }

        elem->get_dof_indices(local_dof_indices);
        Jlocal = 0.0; Rlocal = 0.0;
        compute_element_ad(elem, fe_values, Jlocal, Rlocal, local_dof_indices); //also reinits fe_values

        std::vector<double> d(dofs_per_cell), R0(dofs_per_cell), R1(dofs_per_cell);
        for(unsigned int i = 0 ; i < dofs_per_cell ; i++){ d[i] = sol_k(local_dof_indices[i]); }
        residual_double(d, R0);

        const double h = 1e-7;
        double max_err = 0.0, max_J = 0.0;
        for(unsigned int j = 0 ; j < dofs_per_cell ; j++){
            std::vector<double> dp = d;
            dp[j] += h;
            residual_double(dp, R1);
            for(unsigned int i = 0 ; i < dofs_per_cell ; i++){
                const double J_fd = (R1[i] - R0[i])/h;
                max_err = std::max(max_err, std::abs(J_fd - Jlocal(i,j)));
                max_J   = std::max(max_J, std::abs(Jlocal(i,j)));
            }
        }
        const double rel = max_err/max_J;
        worst = std::max(worst, rel);
        std::cout << "verify_jacobian: element " << elem->active_cell_index()
                  << " | max |J_AD - J_FD| / max |J_AD| = " << rel << std::endl;
        checked++;
    }
    std::cout << "verify_jacobian: worst relative difference = " << worst
              << (worst < 1e-5 ? "  -> OK" : "  -> CHECK THE RESIDUAL/MODEL") << std::endl;
}
