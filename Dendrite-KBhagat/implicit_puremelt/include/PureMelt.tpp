#pragma once

template<unsigned int Nsd, unsigned int BfOrder>
PureMelt<Nsd,BfOrder>::PureMelt(
    const double x_ll, const double x_ul, const unsigned int n_elements
) :
    x_ll_(x_ll), x_ul_(x_ul), n_elements_(n_elements)
{}

template<unsigned int Nsd, unsigned int BfOrder>
void PureMelt<Nsd,BfOrder>::run(){
    make_grid();
}
