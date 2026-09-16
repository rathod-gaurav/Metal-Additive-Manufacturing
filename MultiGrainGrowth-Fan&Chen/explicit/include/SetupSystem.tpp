#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::setup_system(){

    DynamicSparsityPattern dsp(Nt);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); //false means don't keep constrained DoFs in the sparsity pattern
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
    dof_locations_map = DoFTools::map_dofs_to_support_points(MappingQ1<Nsd>(), dof_handler);
    //random noise between -0.001 and 0.001
    // std::default_random_engine gen(123);
    // std::uniform_real_distribution<double> dist(-0.001, 0.001);
    // for(unsigned int i = 0 ; i < p_ ; i++){
    //     double *const row = &eta_n[i][0];
    //     for(unsigned int j = 0 ; j < Nt ; j++){
    //         row[j] = dist(gen);
    //     }
    // }

    //vertical line
    for(unsigned int i = 0 ; i < p_ ; i++){
        double *const row = &eta_n[i][0];
        for(unsigned int j = 0 ; j < Nt ; j++){
            if(dof_locations_map[j][0] > 192.0){
                row[j] = 1.0;
            }
            else{
                row[j] = -1.0;
            }
        }
    }

    // for(unsigned int i = 0 ; i < p_ ; i++){
    //     std::copy(&eta_n[i][0], &eta_n[i][0] + Nt, eta_ni.begin());
    //     constraints.distribute(eta_ni);
    //     std::copy(eta_ni.begin(), eta_ni.end(), &eta_n[i][0]);
    // }

    std::cout << "System setup successful. All global system matrices and vectors initialized" << std::endl;
}