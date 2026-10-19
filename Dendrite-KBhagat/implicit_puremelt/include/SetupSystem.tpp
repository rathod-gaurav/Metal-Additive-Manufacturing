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
}