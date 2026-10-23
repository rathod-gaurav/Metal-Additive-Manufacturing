#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
PureMelt<Nsd,BfOrder>::PureMelt(
    const double x_ll, const double x_ul, const unsigned int n_elements,
    const double undercooling, Point<Nsd> seed_center, const double seed_radius,
    const unsigned int NT, const double dt, const unsigned int N_NR, const double epsilon_NR,
    OutputWriter<Nsd,BfOrder>& output_writer
) :
    x_ll_(x_ll), x_ul_(x_ul), n_elements_(n_elements),
    undercooling_(undercooling), seed_center_(seed_center), seed_radius_(seed_radius),
    NT_(NT), dt_(dt), N_NR_(N_NR), epsilon_NR_(epsilon_NR),
    output_writer_(output_writer),
    fe(FE_Q<Nsd>(BfOrder), 2), //2 for 2 fields per node
    dof_handler(triangulation)
{}

template<unsigned int Nsd, unsigned int BfOrder>
void PureMelt<Nsd,BfOrder>::run(){
    make_grid();
    setup_system();
    solve();
}
