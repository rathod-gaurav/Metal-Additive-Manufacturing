#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::setup_system(){
    dof_handler.distribute_dofs(fe);
    Nt = dof_handler.n_dofs();

    std::cout << "Number of degrees of freedom: " << Nt << std::endl;

    DynamicSparsityPattern dsp(Nt);
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    Mglobal.reinit(sparsity_pattern);
    Kglobal.reinit(sparsity_pattern);
    Fglobal.reinit(Nt);

    eta_ni.reinit(Nt);
    eta_np1i.reinit(Nt);
    eta_n2.reinit(Nt);
    RHS.reinit(Nt);

    eta_n.reinit(p_, Nt, true);
    eta_np1.reinit(p_, Nt, true);

    phi.reinit(Nt);

    //initial conditions
    // code to fill eta_n here
    const std::map<types::global_dof_index, Point<Nsd>> dof_locations_map = DoFTools::map_dofs_to_support_points(MappingQ1<Nsd>(), dof_handler);
    
    std::default_random_engine gen(123);
    std::uniform_real_distribution<double> dist(-0.001, 0.001);
    for(unsigned int i = 0 ; i < p_ ; i++){
        double *const row = &eta_n[i][0];
        for(unsigned int j = 0 ; j < Nt ; j++){
            row[j] = dist(gen);
        }
    }

    std::cout << "System setup successful. All global system matrices and vectors initialized" << std::endl;
}