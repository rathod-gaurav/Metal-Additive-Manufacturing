#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
double BinaryAlloy<Nsd,BfOrder>::max_interface_speed() const {
    const QGauss<Nsd> quad(quadOrder_);
    FEValues<Nsd> fe_values(fe, quad, update_values | update_gradients);
    const unsigned int n_q = quad.size();
    std::vector<double> phi(n_q), phi_old(n_q);
    std::vector<Tensor<1,Nsd>> grad_phi(n_q);

    double v_max = 0.0;
    for (const auto& cell : dof_handler.active_cell_iterators()){
        fe_values.reinit(cell);
        fe_values[field_ext[PHI_]].get_function_values(sol_k, phi);
        fe_values[field_ext[PHI_]].get_function_values(sol_n, phi_old);
        fe_values[field_ext[PHI_]].get_function_gradients(sol_k, grad_phi);
        for (unsigned int q = 0; q < n_q; ++q)
            if (std::abs(phi[q]) < 0.5 && grad_phi[q].norm() > 1e-6)
                v_max = std::max(v_max, std::abs(phi[q] - phi_old[q])/dt_ / grad_phi[q].norm());
    }
    return v_max;
}