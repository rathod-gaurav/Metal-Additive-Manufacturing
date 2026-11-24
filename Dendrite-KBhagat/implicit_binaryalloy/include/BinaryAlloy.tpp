#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
BinaryAlloy<Nsd,BfOrder>::BinaryAlloy(
    const double x_ll, const double x_ul, const unsigned int n_elements, const unsigned int quadOrder,
    const unsigned int UC, const unsigned int PHI, const unsigned int N_FIELDS, 
    const double supersaturation, const std::vector<NucleationSeed<Nsd>>& seeds,
    const unsigned int NT, const double dt, const unsigned int N_NR, const double epsilon_NR,
    const ConstitutiveModel& model, 
    OutputWriter<Nsd,BfOrder>& output_writer,
    const unsigned int max_level, const unsigned int min_level, const double delta_phi, const double delta_coarsen //adaptive mesh refinement
) :
    x_ll_(x_ll), x_ul_(x_ul), n_elements_(n_elements), quadOrder_(quadOrder),
    UC_(UC), PHI_(PHI), N_FIELDS_(N_FIELDS), 
    supersaturation_(supersaturation), seeds_(seeds),
    NT_(NT), dt_(dt), N_NR_(N_NR), epsilon_NR_(epsilon_NR),
    model_(model),
    output_writer_(output_writer),
    max_level_(max_level), min_level_(min_level), delta_phi_(delta_phi), delta_coarsen_(delta_coarsen),
    triangulation(Triangulation<Nsd>::limit_level_difference_at_vertices),
    fe(FE_Q<Nsd>(BfOrder), N_FIELDS_), //2 for 2 fields per node
    dof_handler(triangulation)
{
    for(unsigned int f = 0 ; f < N_FIELDS_ ; f++){
        field_ext.emplace_back(f);
    }
}

template<unsigned int Nsd, unsigned int BfOrder>
void BinaryAlloy<Nsd,BfOrder>::run(){
    make_grid();
    setup_system_dofs();
    apply_initial_conditions();
    for (unsigned int cycle = 0; cycle < max_level_ - min_level_; ++cycle){
        flag_cells_for_refinement();
        triangulation.prepare_coarsening_and_refinement();
        triangulation.execute_coarsening_and_refinement();
        setup_system_dofs();
        apply_initial_conditions(); //exact IC on the finer mesh
    }

    //check total solute upon mesh refinement (should be conserved)
    const double C0 = compute_total_solute(sol_n);
    for (unsigned int i = 1; i <= 3; ++i){
        refine_mesh();
        const double C = compute_total_solute(sol_n);
        std::cout << "remesh " << i << " | C = " << C << " | rel. change = " << std::abs(C - C0)/C0 << std::endl;
        std::ofstream pvd_output_test(output_writer_.get_output_dir() + "/total_solute_test.pvd");
        output_writer_.initiate_pvd(pvd_output_test); //initiated pvd file with headers
        output_writer_.write_vtu_and_pvd(dof_handler, sol_n, {"uc","phi"}, i, pvd_output_test);
    }
    return;

    verify_jacobian(5); //check 5 elements for AD vs FD Jacobian
    solve();
}
