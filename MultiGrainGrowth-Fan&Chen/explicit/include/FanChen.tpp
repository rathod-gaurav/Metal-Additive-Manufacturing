#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
FanChen<Nsd,BfOrder>::FanChen(
    const double x_ll, const double x_ul, const unsigned int n_refinements, //for creating mesh grid
    const unsigned int p, //number of distinct grain orientations in the domain
    const double dt, const unsigned int NT, //time step size, number of time steps
    const double L, const double kappa, //mobility, gradient energy coefficient //both assumed constant here
    const double alpha, const double beta, const double gamma, //kinetic equation coefficients
    const unsigned int quadOrder
):
    x_ll_(x_ll),
    x_ul_(x_ul),
    n_refinements_(n_refinements),
    p_(p),
    dt_(dt),
    NT_(NT),
    L_(L),
    kappa_(kappa),
    alpha_(alpha),
    beta_(beta),
    gamma_(gamma),
    quadOrder_(quadOrder),
    fe(BfOrder),
    dof_handler(triangulation)
{}

template<unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::run(){
    make_grid();
    setup_system();
    solve();
}