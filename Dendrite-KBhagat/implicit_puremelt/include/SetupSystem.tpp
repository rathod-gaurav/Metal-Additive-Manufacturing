#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
void PureMelt<Nsd,BfOrder>::setup_system(){
    dof_handler.distribute_dofs(fe);
    Nt = dof_handler.n_dofs();
    std::cout << "DoFs per vertex: " << fe.n_dofs_per_vertex() << std::endl;
    std::cout << "DoFs per cell: " << fe.n_dofs_per_cell() << std::endl;
    std::cout << "Total number of degrees of freedom: " << Nt << std::endl;

    // automatic differentiation setup
    constexpr Differentiation::AD::NumberTypes ADTypeCode = Differentiation::AD::NumberTypes::sacado_dfad;
    using ADHelper = Differentiation::AD::ResidualLinearization<ADTypeCode, double>;
    using ADNumber = typename ADHelper::ad_type;

    DynamicSparsityPattern dsp(Nt);
    DoFTools::make_sparsity_pattern(dof_handler,dsp); //coupled both u and phi
    sparsity_pattern.copy_from(dsp);

    system_jacobian.reinit(sparsity_pattern);
    system_residual.reinit(Nt);
    solution.reinit(Nt);
    solution_old.reinit(Nt);
    newton_update.reinit(Nt);
}