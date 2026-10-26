#pragma once

#include "ADConfig.hpp"
#include "EvaluateFields.tpp"
#include "ElementResidual.tpp"

template<unsigned int Nsd, unsigned int BfOrder>
void PureMelt<Nsd,BfOrder>::compute_element(
    const typename DoFHandler<Nsd>::active_cell_iterator& elem, FEValues<Nsd>& fe_values, 
    FullMatrix<double>& Jlocal, Vector<double>& Rlocal,
    const std::vector<types::global_dof_index>& local_dof_indices
){
    using namespace AD;

    fe_values.reinit(elem);
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q = fe_values.n_quadrature_points;

    //Step 1: fresh helper for this element
    ADHelper ad_helper(dofs_per_cell, dofs_per_cell);

    //Step 2: register the current iterate as the independent variables
    ad_helper.register_dof_values(sol_k, local_dof_indices);
    const std::vector<ADNumber>& dof_values_ad = ad_helper.get_sensitive_dof_values();

    //Step 3: fields at the quadrature points, as AD Numbers
    QPFields<Nsd,ADNumber> fields(N_FIELDS_, n_q);
    evaluate_fields<ADNumber>(fe_values, dof_values_ad, fields);

    //Step 4: compute the element residual as AD numbers
    std::vector<ADNumber> residual_ad(dofs_per_cell, ADNumber(0.0)); //Creates the output: a vector with one entry per local DoF (8), each an AD number set to zero.
    element_residual<ADNumber>(fe_values, fields, residual_ad);
}