#pragma once

#include "InitialConditions.hpp"

template<unsigned int Nsd, unsigned int BfOrder>
void PureMelt<Nsd,BfOrder>::setup_system(){
    dof_handler.distribute_dofs(fe);
    Nt = dof_handler.n_dofs();
    std::cout << "DoFs per vertex: " << fe.n_dofs_per_vertex() << std::endl;
    std::cout << "DoFs per cell: " << fe.n_dofs_per_cell() << std::endl;
    std::cout << "Total number of degrees of freedom: " << Nt << std::endl;

    DynamicSparsityPattern dsp(Nt);
    DoFTools::make_sparsity_pattern(dof_handler,dsp); //coupled both u and phi
    sparsity_pattern.copy_from(dsp);

    NR_jacobian.reinit(sparsity_pattern);
    NR_residual.reinit(Nt);
    sol_k.reinit(Nt);
    sol_n.reinit(Nt);
    NR_update.reinit(Nt);

    VectorTools::interpolate(dof_handler, InitialConditions<Nsd>(undercooling_, seed_center_, seed_radius_), sol_n);
    sol_k = sol_n; //guess for Newton's first iteration
}