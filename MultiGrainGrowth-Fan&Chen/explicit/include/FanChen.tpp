#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
FanChen<Nsd,BfOrder>::FanChen(
    const double x_ll, const double x_ul, const unsigned int n_refinements, //for creating mesh grid
    const unsigned int p //number of distinct grain orientations in the domain
):
    x_ll_(x_ll),
    x_ul_(x_ul),
    n_refinements_(n_refinements),
    p_(p),
    fe(BfOrder),
    dof_handler(triangulation)
{}

template<unsigned int Nsd, unsigned int BfOrder>
void FanChen<Nsd,BfOrder>::run(){
    make_grid();
    setup_system();
}