#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
BinaryAlloy<Nsd,BfOrder>::BinaryAlloy(
    const double x_ll, const double x_ul, const unsigned int n_elements, const unsigned int quadOrder,
    const unsigned int UC, const unsigned int PHI, const unsigned int N_FIELDS, 
    const double c0, const std::vector<NucleationSeed<Nsd>>& seeds,
    const unsigned int NT, const double dt, const unsigned int N_NR, const double epsilon_NR,
    const ConstitutiveModel& model, 
    OutputWriter<Nsd,BfOrder>& output_writer
) :
    x_ll_(x_ll), x_ul_(x_ul), n_elements_(n_elements), quadOrder_(quadOrder),
    UC_(UC), PHI_(PHI), N_FIELDS_(N_FIELDS), 
    c0_(c0), seeds_(seeds),
    NT_(NT), dt_(dt), N_NR_(N_NR), epsilon_NR_(epsilon_NR),
    model_(model),
    output_writer_(output_writer),
    fe(FE_Q<Nsd>(BfOrder), N_FIELDS_), //2 for 2 fields per node
    dof_handler(triangulation)
{
    for(unsigned int f = 0 ; f < N_FIELDS_ ; f++){
        field_ext.emplace_back(f);
    }
}

UClate<unsigned int Nsd, unsigned int BfOrder>
void BinaryAlloy<Nsd,BfOrder>::run(){
    make_grid();
    setup_system();
    verify_jacobian(5); //check 5 elements for AD vs FD Jacobian
    solve();
}
